#![feature(iter_next_chunk, iter_map_windows, iter_collect_into)]

use std::{
    collections::HashMap,
    env::args,
    fs::{read_dir, File},
    io::{stdout, BufRead, BufReader, Write},
    path::{Path, PathBuf},
};

use anyhow::{anyhow, Context, Result};
use colorous::{Gradient, TURBO};
use futures::FutureExt;
use macroquad::{
    prelude::*,
    ui::{
        hash, root_ui,
        widgets::{Checkbox, Editbox, Label},
        Ui,
    },
    Window,
};
use miniquad::window::request_quit;
use nucleo::{
    pattern::{CaseMatching, Normalization, Pattern},
    Config, Matcher,
};

// TODO(refact): split into separate files
// TODO(feat): add ability to launch without arguments and open directory through UI
// TODO(feat): add ability to change controls (config file? env vars? arguments?)
// TODO(feat): add ability to change colour scheme
// TODO(feat): draggable/editable bounding box points
// TODO(feat): add support for discrete features
// TODO(feat): add support for using segmentations instead of centroids
// TODO(feat): show top N matches for feature names in selector dialog
// TODO(feat): add edit history / ability to undo
// TODO(perf): add some kind of caching for feature colouring
// TODO(refact): improve error handling

#[allow(non_camel_case_types)]
type float = f64;
#[allow(non_upper_case_globals)]
const tau: float = std::f64::consts::TAU;

const DOT_SIZE: f32 = 2.;
const BB_LINE_ACTIVE_WIDTH: f32 = 3.;
const BB_LINE_FINISHED_WIDTH: f32 = 1.;
const LEGEND_SIZE: usize = 64;
const LEGEND_SIZE_F32: f32 = LEGEND_SIZE as f32;
const BORDER_PADDING: f32 = 3.;
const TEXT_SIZE: f32 = 18.;
const TEXT_SIZE_U16: u16 = TEXT_SIZE as u16;
const CONTINUOUS_PALLETE: Gradient = TURBO;
const MOV_SPEED: float = 384.;
const ROT_SPEED: float = tau / 8.;
const ZOM_SPEED: float = 2.56;
const DIALOG_BOX_W: f32 = 200.;
const DIALOG_BOX_H: f32 = 80.;
const EDIT_BOX_W: f32 = 194.;
const EDIT_BOX_H: f32 = 24.;
const CHECKBOX_MAGIC_H: f32 = 19.;
const CHECKBOX_MAGIC_X_SHIFT: f32 = 35.;

type ExprMap = HashMap<String, float>;
type ExprMapBundle = (ExprMap, (String, TextDimensions));
type CellCoordsBundle = (String, float, float);

fn read_coords<P: AsRef<Path>>(inp_path: P) -> Result<Vec<CellCoordsBundle>> {
    let fpath = inp_path.as_ref().join("cells.tsv");
    let mut lines = BufReader::new(File::open(&fpath)?).lines();

    let header = lines
        .next()
        .ok_or_else(|| anyhow!("coord file {} is empty", inp_path.as_ref().display()))??;
    let header = header
        .split('\t')
        .next_chunk()
        .map_err(|_| anyhow!("less than 3 columns in header for coord file {}", inp_path.as_ref().display()))?;
    if header != ["cell", "x", "y"] {
        return Err(anyhow!(
            "expected coord file {fpath:?} header to be [\"cell\", \"x\", \"y\"], got {header:?}"
        ));
    }

    lines
        .enumerate()
        .map(|(i, l)| {
            let l = l?;
            let [cell, x, y] = l
                .split('\t')
                .next_chunk()
                .map_err(|_| anyhow!("less than 3 columns in coord file {fpath:?} at line {}", i + 1))?;
            Ok((cell.to_string(), x.parse::<float>()?, y.parse::<float>()?))
        })
        .collect()
}

fn read_feat<P: AsRef<Path>, S: AsRef<str>>(inp_path: P, gene: S) -> Result<ExprMap> {
    let fpath = inp_path.as_ref().join("feat").join(format!("{}.tsv", gene.as_ref()));
    let mut lines = BufReader::new(File::open(&fpath)?).lines();

    let header = lines
        .next()
        .ok_or_else(|| anyhow!("feature file file {} is empty", inp_path.as_ref().display()))??;
    let header = header.split('\t').next_chunk().map_err(|_| {
        anyhow!(
            "less than 2 columns in header for feature file {}",
            inp_path.as_ref().display()
        )
    })?;
    if header != ["cell", "expr"] {
        return Err(anyhow!(
            "expected feature file {fpath:?} header to be [\"cell\", \"expr\"], got {header:?}"
        ));
    }

    lines
        .enumerate()
        .map(|(i, l)| {
            let l = l?;
            let [cell, val] = l
                .split('\t')
                .next_chunk()
                .map_err(|_| anyhow!("less than 2 columns in feature file {fpath:?} at line {}", i + 1))?;
            Ok((cell.to_string(), val.parse::<float>()?))
        })
        .collect()
}

fn write_out<W: Write>(out_sink: &mut W, cells: &[CellCoordsBundle], bb_list: &[Vec<(float, float)>]) -> Result<()> {
    let mut out = vec![["cell", "within_bb"].join("\t")];
    cells
        .iter()
        .map(|&(ref cell, x, y)| [cell, bb_list.iter().any(|bb| within_bb(bb, (x, y))).to_string().as_str()].join("\t"))
        .collect_into(&mut out);
    out_sink.write_all(out.join("\n").as_bytes())?;
    Ok(())
}

fn draw_cells<'a, I: Iterator<Item = &'a CellCoordsBundle>>(
    c_iter: I,
    transform: &Transform,
    expr: &Option<ExprMapBundle>,
    m_x: f32,
    m_y: f32,
    w: f32,
    h: f32,
) {
    let c_iter = c_iter
        .map(|&(ref cell, x, y)| {
            let (x, y) = transform.to_pixl(x, y);
            (cell, x + m_x, y + m_y)
        })
        .filter(|&(_cell, x, y)| (x > 0.) && (x < w) && (y > 0.) && (y < h));

    match expr {
        Some((hm, _)) => c_iter.for_each(|(cell, x, y)| {
            draw_rectangle(x, y, DOT_SIZE, DOT_SIZE, {
                let c = CONTINUOUS_PALLETE.eval_continuous(hm.get(cell).cloned().unwrap_or(0.));
                Color::from_rgba(c.r, c.g, c.b, 255)
            });
        }),
        None => c_iter.for_each(|(_cell, x, y)| {
            draw_rectangle(x, y, DOT_SIZE, DOT_SIZE, WHITE);
        }),
    };
}

fn extend_with_first<T>(v: &[T]) -> impl Iterator<Item = &T> {
    v.iter().cycle().take(v.len() + 1)
}

fn within_bb(poly: &[(float, float)], (x, y): (float, float)) -> bool {
    let mut winding_num = 0;
    extend_with_first(poly)
        .map_windows(|e| *e)
        .for_each(|[&(prev_x, prev_y), &(curr_x, curr_y)]| {
            if curr_y <= y {
                if prev_y > y && ((prev_x - curr_x) * (y - curr_y) > (x - curr_x) * (prev_y - curr_y)) {
                    winding_num += 1;
                }
            } else if prev_y <= y && ((prev_x - curr_x) * (y - curr_y) < (x - curr_x) * (prev_y - curr_y)) {
                winding_num -= 1;
            }
        });

    winding_num != 0
}

fn pos_rot(theta: float, x: float, y: float) -> (float, float) {
    let (t_cos, t_sin) = (theta.cos(), theta.sin());
    (x * t_cos - y * t_sin, x * t_sin + y * t_cos)
}

fn neg_rot(theta: float, x: float, y: float) -> (float, float) {
    let (t_cos, t_sin) = (theta.cos(), theta.sin());
    (x * t_cos + y * t_sin, y * t_cos - x * t_sin)
}

struct Transform {
    shift_x: float,
    shift_y: float,
    zoom: float,
    theta: float,
    inv_x: bool,
    inv_y: bool,
    frame_time: float,
}

impl Transform {
    fn to_pixl(&self, x: float, y: float) -> (f32, f32) {
        let (x, y) = pos_rot(self.theta, x, y);
        let (x, y) = (x + self.shift_x, y + self.shift_y);
        let (x, y) = (if self.inv_x { -x } else { x }, if self.inv_y { -y } else { y });
        let (x, y) = (x * self.zoom, y * self.zoom);
        let (x, y) = (x as f32, y as f32);
        (x, y)
    }
    fn to_real(&self, x: f32, y: f32) -> (float, float) {
        let (x, y) = (x as float, y as float);
        let (x, y) = (x / self.zoom, y / self.zoom);
        let (x, y) = (if self.inv_x { -x } else { x }, if self.inv_y { -y } else { y });
        let (x, y) = (x - self.shift_x, y - self.shift_y);
        let (x, y) = neg_rot(self.theta, x, y);
        (x, y)
    }
    fn mov(&mut self, x: float, y: float) {
        let mult = self.zoom.recip() * self.frame_time;
        if self.inv_x {
            self.shift_x -= x * mult;
        } else {
            self.shift_x += x * mult;
        }
        if self.inv_y {
            self.shift_y -= y * mult;
        } else {
            self.shift_y += y * mult;
        }
    }
    fn zoom_inc(&mut self) {
        self.zoom = (self.zoom * ZOM_SPEED.powf(self.frame_time)).min(float::MAX);
    }
    fn zoom_dec(&mut self) {
        self.zoom = (self.zoom / ZOM_SPEED.powf(self.frame_time)).max(0.);
    }
    fn theta_inc(&mut self) {
        let scaled_rot = ROT_SPEED * self.frame_time;
        if self.inv_x ^ self.inv_y {
            self.theta = (self.theta + scaled_rot).rem_euclid(tau);
            (self.shift_x, self.shift_y) = pos_rot(scaled_rot, self.shift_x, self.shift_y);
        } else {
            self.theta = (self.theta - scaled_rot).rem_euclid(tau);
            (self.shift_x, self.shift_y) = neg_rot(scaled_rot, self.shift_x, self.shift_y);
        };
    }
    fn theta_dec(&mut self) {
        let scaled_rot = ROT_SPEED * self.frame_time;
        if self.inv_x ^ self.inv_y {
            self.theta = (self.theta - scaled_rot).rem_euclid(tau);
            (self.shift_x, self.shift_y) = neg_rot(scaled_rot, self.shift_x, self.shift_y);
        } else {
            self.theta = (self.theta + scaled_rot).rem_euclid(tau);
            (self.shift_x, self.shift_y) = pos_rot(scaled_rot, self.shift_x, self.shift_y);
        };
    }
}

fn draw_lines<'a, I: Iterator<Item = &'a (float, float)>>(
    coords: I,
    transform: &Transform,
    mid_x: f32,
    mid_y: f32,
    width: f32,
    color: Color,
) {
    coords
        .map(|&(x, y)| {
            let (x, y) = transform.to_pixl(x, y);
            (x + mid_x, y + mid_y)
        })
        .map_windows(|e| *e)
        .for_each(|[(x1, y1), (x2, y2)]| draw_line(x1, y1, x2, y2, width, color))
}

fn draw_status(inp: &str, height: &mut f32) {
    *height += measure_text(inp, None, TEXT_SIZE_U16, 1.).height + BORDER_PADDING;
    draw_text(inp, BORDER_PADDING, *height, TEXT_SIZE, YELLOW);
}

fn draw_bool_opt(pos: &mut Vec2, label: &str, ui: &mut Ui, opt: &mut bool) {
    let pos_c = *pos;
    pos.y += CHECKBOX_MAGIC_H + BORDER_PADDING;
    Checkbox::new(hash!()).size(vec2(0., 0.)).pos(pos_c).label(label).ui(ui, opt)
}

struct StatusConfig {
    draw_fps: bool,
    draw_rot: bool,
    draw_flip: bool,
}

enum UIState {
    CellSelection,
    FeatureInput,
    HelpMenu,
}

fn main() {
    Window::from_config(
        Conf {
            window_title: "scs".to_string(),
            high_dpi: true,
            ..Default::default()
        },
        ui().map(|f| {
            if let Err(error) = f {
                println!("ERROR: {:?}", error)
            };
        }),
    );
}

async fn ui() -> Result<()> {
    let mut args = args().skip(1);
    let inp_path = args
        .next()
        .ok_or_else(|| anyhow!("input path must be supplied as argument"))?;
    let out_path = args
        .next()
        .map(|p| {
            File::create(&p)
                .with_context(|| format!("trying to open output file {p}"))
                .map(|e| (p, e))
        })
        .transpose()?;

    let mut bb: Option<Vec<(f64, f64)>> = None;
    let mut bb_list = Vec::new();

    let valid_genes = read_dir(PathBuf::from(&inp_path).join("feat"))
        .context("trying to scan feature file directory")?
        .map(|e| {
            let f = e?.path();
            Ok(f.file_stem()
                .ok_or_else(|| anyhow!("could not get file stem for feature file {}", f.display()))?
                .to_str()
                .ok_or_else(|| anyhow!("file stem is not valid UTF-8 for feature file {}", f.display()))?
                .chars()
                .collect::<String>())
        })
        .collect::<Result<Vec<_>>>()
        .context("trying to load all feature names")?;

    let cells = read_coords(&inp_path).context("trying to read cells.tsv file")?;

    // initialize shift to average x/y values
    let mut transform = {
        let mut iter = cells.iter().map(|&(ref _cell, x, y)| (x, y));
        let (init_x, init_y) = iter.next().ok_or_else(|| anyhow!("no cells in cells.tsv file"))?;
        let (mut avg_x, mut avg_y, min_x, min_y, max_x, max_y) = iter.fold(
            (init_x, init_y, init_x, init_y, init_x, init_y),
            |(a_x, a_y, min_x, min_y, max_x, max_y), (x, y)| {
                (a_x + x, a_y + y, min_x.min(x), min_y.min(y), max_x.max(x), max_y.max(y))
            },
        );
        let l = cells.len() as float;
        (avg_x, avg_y) = (-avg_x / l, -avg_y / l);
        Transform {
            shift_x: avg_x,
            shift_y: avg_y,
            zoom: ((screen_width() as f64) / (max_x - min_x)).min((screen_height() as f64) / (max_y - min_y)),
            theta: 0.,
            inv_x: false,
            inv_y: false,
            frame_time: get_frame_time() as float,
        }
    };

    let mut ui_state = UIState::CellSelection;

    let mut input = "".to_string();
    let mut prev_input = input.clone();
    let mut pat = Pattern::parse(&input, CaseMatching::Ignore, Normalization::Smart);
    let mut matcher = Matcher::new({
        let mut c = Config::DEFAULT;
        c.prefer_prefix = true;
        c
    });

    let mut feat_match: (Vec<(&String, _)>, _) = (Vec::new(), 0);
    let mut expr_bundle: Option<ExprMapBundle> = None;

    let mut show_only_filtered = false;
    let mut filtered = Vec::new();

    let mut status_config = StatusConfig {
        draw_fps: false,
        draw_rot: false,
        draw_flip: false,
    };

    let colour_scale = (0..LEGEND_SIZE)
        .map(|i| CONTINUOUS_PALLETE.eval_rational(i, LEGEND_SIZE))
        .map(|c| Color::from_rgba(c.r, c.g, c.b, 255))
        .collect::<Vec<_>>();

    let mut quit_attempt = false;
    let mut quit_success = false;
    prevent_quit();

    loop {
        clear_background(BLACK);

        let w = screen_width();
        let h = screen_height();
        let (m_x, m_y) = (w / 2., h / 2.);

        if is_quit_requested() {
            quit_attempt = true;
        };
        if quit_success {
            break;
        }
        if quit_attempt {
            let quit_size = measure_text("CONFIRM QUIT?YESNO", None, TEXT_SIZE_U16, 1.);
            let exit_dialog_size = vec2(quit_size.width + 9., quit_size.height + 15.);
            root_ui().window(
                hash!(hash!(format!("{}{}", m_x, m_y), hash!())),
                vec2(m_x, m_y) - exit_dialog_size / 2.,
                exit_dialog_size,
                |ui| {
                    ui.label(None, "CONFIRM QUIT?");
                    ui.same_line(0.);
                    if ui.button(None, "YES") {
                        quit_success = true;
                    }
                    ui.same_line(0.);
                    if ui.button(None, "NO") {
                        quit_attempt = false;
                    };
                },
            );
            next_frame().await;
            continue;
        };

        transform.frame_time = get_frame_time() as float;

        match ui_state {
            UIState::CellSelection => {
                get_keys_pressed().into_iter().for_each(|k| match k {
                    KeyCode::Space => {
                        ui_state = UIState::FeatureInput;
                    }
                    KeyCode::H => {
                        ui_state = UIState::HelpMenu;
                    }
                    KeyCode::Escape => {
                        request_quit();
                    }
                    _ => {}
                });

                if matches!(ui_state, UIState::CellSelection) {
                    if is_mouse_button_pressed(MouseButton::Left) {
                        let (mouse_x, mouse_y) = mouse_position();
                        bb.get_or_insert_with(Vec::new)
                            .push(transform.to_real(mouse_x - m_x, mouse_y - m_y));
                    }

                    let (mut x_mov, mut y_mov): (float, float) = (0., 0.);
                    get_keys_pressed().into_iter().for_each(|k| match k {
                        KeyCode::Enter => {
                            if let Some(v) = bb.take() {
                                bb_list.push(v);
                            }
                            cells
                                .iter()
                                .filter(|&&(ref _cell, x, y)| bb_list.iter().any(|bb| within_bb(bb, (x, y))))
                                .cloned()
                                .collect_into(&mut filtered);
                        }
                        KeyCode::F => {
                            show_only_filtered = !show_only_filtered;
                        }
                        KeyCode::C => {
                            transform.inv_x = !transform.inv_x;
                        }
                        KeyCode::V => {
                            transform.inv_y = !transform.inv_y;
                        }
                        KeyCode::R => {
                            feat_match.0.clear();
                            feat_match.1 = 0;
                            expr_bundle = None;
                        }
                        _ => {}
                    });
                    get_keys_down().into_iter().for_each(|k| match k {
                        KeyCode::A => {
                            x_mov += 1.;
                        }
                        KeyCode::D => {
                            x_mov -= 1.;
                        }
                        KeyCode::W => {
                            y_mov += 1.;
                        }
                        KeyCode::S => {
                            y_mov -= 1.;
                        }
                        KeyCode::Z => transform.zoom_inc(),
                        KeyCode::X => transform.zoom_dec(),
                        KeyCode::Q => transform.theta_inc(),
                        KeyCode::E => transform.theta_dec(),
                        _ => {}
                    });
                    let l_mov = x_mov.hypot(y_mov);
                    if l_mov > 0. {
                        let norm = MOV_SPEED / l_mov;
                        transform.mov(x_mov * norm, y_mov * norm);
                    }
                }
            }
            UIState::FeatureInput => {
                for k in get_keys_pressed().into_iter() {
                    match k {
                        KeyCode::Enter => {
                            if !feat_match.0.is_empty() {
                                let sel_feat = feat_match.0[feat_match.1].0;
                                let expr_map = read_feat(&inp_path, sel_feat).context("trying to read feature file")?;
                                let max_recip = expr_map
                                    .values()
                                    .cloned()
                                    .reduce(|l, r| l.max(r))
                                    .ok_or_else(|| anyhow!("empty expression map - should never happen"))?
                                    .recip();
                                expr_bundle = Some(expr_map.into_iter().map(|(k, v)| (k, v * max_recip)).collect())
                                    .map(|e| (e, (sel_feat.clone(), measure_text(sel_feat, None, TEXT_SIZE_U16, 1.))));
                            }
                            input.clear();
                            ui_state = UIState::CellSelection;
                            break;
                        }
                        KeyCode::Escape => {
                            input.clear();
                            ui_state = UIState::CellSelection;
                            break;
                        }
                        _ => {}
                    }
                }

                if matches!(ui_state, UIState::FeatureInput) {
                    get_keys_pressed().into_iter().for_each(|k| {
                        #[allow(clippy::single_match)]
                        match k {
                            KeyCode::Tab => {
                                let l = feat_match.0.len();
                                if l != 0 {
                                    feat_match.1 += 1;
                                    feat_match.1 %= l;
                                }
                            }
                            _ => {}
                        }
                    });

                    // detect input event into editbox using cached previous value
                    if input.ne(&prev_input) {
                        prev_input = input.clone();
                        pat.reparse(&input, CaseMatching::Ignore, Normalization::Smart);
                        // TODO(perf): determine if it's worthwhile for matching to be moved to a background thread
                        let mut matches = pat.match_list(&valid_genes, &mut matcher);
                        matches.retain(|(_m, s)| 0.ne(s));
                        matches.sort_unstable_by(|lhs, rhs| rhs.1.cmp(&lhs.1));
                        feat_match = (matches, 0);
                    }

                    let dialog_size = vec2(DIALOG_BOX_W, DIALOG_BOX_H);
                    let window_pos = vec2(m_x, m_y) - dialog_size / 2.;
                    root_ui().window(hash!(format!("{}{}", m_x, m_y), hash!()), window_pos, dialog_size, |ui| {
                        let edit_box_id = hash!();
                        ui.label(None, "select a gene:");
                        ui.separator();
                        Editbox::new(edit_box_id, vec2(EDIT_BOX_W, EDIT_BOX_H))
                            .multiline(false)
                            .ui(ui, &mut input);
                        ui.set_input_focus(edit_box_id);
                        ui.separator();
                        if !feat_match.0.is_empty() {
                            ui.label(None, &format!("selected: {}", feat_match.0[feat_match.1].0));
                        } else {
                            ui.label(None, "no match found");
                        }
                    });
                }
            }
            UIState::HelpMenu => {
                get_keys_pressed().into_iter().for_each(|k| {
                    #[allow(clippy::single_match)]
                    match k {
                        KeyCode::Escape => {
                            ui_state = UIState::CellSelection;
                        }
                        _ => {}
                    }
                });

                if matches!(ui_state, UIState::HelpMenu) {
                    let help_txt = [
                        "general:",
                        "\t- use mouse left click select create bounding box vertices",
                        "\t- use enter to finish bounding box",
                        "\t- use F to toggle showing only selected cells",
                        "\t- use W/A/S/D to move along the x/y axis",
                        "\t- use Q/E to rotate",
                        "\t- use Z/X to zoom",
                        "\t- use C/V to flip x/y axis",
                        "\t- use space to bring up feature selection dialog",
                        "\t- use R to reset colour",
                        "\t- use H to bring up this menu",
                        "\t- use escape to quit",
                        "feature selection:",
                        "\t- type to select a feature",
                        "\t- use enter to confirm a selected feature",
                        "\t- use tab to cycle through matches",
                        "\t- use escape to exit feature selection dialog",
                        "help dialog:",
                        "\t- use escape to exit help dialog",
                    ];
                    let help_txt_dims = measure_text(
                        help_txt
                            .iter()
                            .reduce(|a, c| if c.len() > a.len() { c } else { a })
                            .ok_or_else(|| anyhow!("help string empty - should never happen"))?,
                        None,
                        TEXT_SIZE_U16,
                        1.,
                    );
                    let base_height = (CHECKBOX_MAGIC_H + BORDER_PADDING) * (std::mem::size_of::<StatusConfig>() as f32);
                    let dialog_size = vec2(
                        help_txt_dims.width + BORDER_PADDING * 2.,
                        base_height + (help_txt_dims.height + BORDER_PADDING) * (help_txt.len() as f32) + 3.,
                    );
                    let window_pos = vec2(m_x, m_y) - dialog_size / 2.;
                    let mut coords = vec2(
                        window_pos.x + CHECKBOX_MAGIC_X_SHIFT + BORDER_PADDING,
                        window_pos.y - 1. + BORDER_PADDING,
                    );
                    root_ui().window(
                        hash!(format!("{}{}", window_pos.x, window_pos.y), hash!()),
                        window_pos,
                        dialog_size,
                        |ui| {
                            draw_bool_opt(&mut coords, "show fps counter", ui, &mut status_config.draw_fps);
                            draw_bool_opt(&mut coords, "show rotation angle", ui, &mut status_config.draw_rot);
                            draw_bool_opt(&mut coords, "show x/y axis flip status", ui, &mut status_config.draw_flip);
                            help_txt.iter().fold(base_height, |y, &line| {
                                Label::new(line).position(vec2(BORDER_PADDING, y)).ui(ui);
                                y + help_txt_dims.height + BORDER_PADDING
                            });
                        },
                    );
                }
            }
        }

        if show_only_filtered {
            draw_cells(filtered.iter(), &transform, &expr_bundle, m_x, m_y, w, h);
        } else {
            draw_cells(cells.iter(), &transform, &expr_bundle, m_x, m_y, w, h);
        }

        if let Some(v) = &bb {
            draw_lines(v.iter(), &transform, m_x, m_y, BB_LINE_ACTIVE_WIDTH, YELLOW);
        }

        bb_list
            .iter()
            .for_each(|v| draw_lines(extend_with_first(v), &transform, m_x, m_y, BB_LINE_FINISHED_WIDTH, YELLOW));

        if let Some((_, (s, size))) = &expr_bundle {
            colour_scale
                .iter()
                .rev()
                .enumerate()
                .for_each(|(i, &c)| draw_rectangle(w - (i as f32) - BORDER_PADDING, BORDER_PADDING, 1., size.height, c));
            draw_text(
                s,
                w - size.width - LEGEND_SIZE_F32 - BORDER_PADDING,
                BORDER_PADDING + size.height,
                TEXT_SIZE,
                YELLOW,
            );
        }

        let mut height = 0.;
        if status_config.draw_fps {
            draw_status(&format!("fps: {}", get_fps()), &mut height);
        }
        if status_config.draw_rot {
            draw_status(&format!("rotation angle: {:.2}", transform.theta.to_degrees()), &mut height);
        }
        if status_config.draw_flip {
            draw_status(
                &format!(
                    "x/y: {}/{}",
                    if transform.inv_x { "flip" } else { "orig" },
                    if transform.inv_y { "flip" } else { "orig" }
                ),
                &mut height,
            );
        }

        next_frame().await
    }

    match out_path {
        Some((p, mut f)) => write_out(&mut f, &cells, &bb_list).with_context(|| format!("trying to write output to {p}")),
        None => write_out(&mut stdout(), &cells, &bb_list).context("trying to write output to standard out"),
    }?;
    Ok(())
}
