# TODO: deduplicate dev deps (currently just cargo-dist)

{
  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixpkgs-unstable";
  };
  outputs = { self, nixpkgs }:
    {
      devShells."aarch64-darwin".default = with nixpkgs.legacyPackages."aarch64-darwin"; mkShell {
        buildInputs = [
          darwin.libiconv
        ] ++ (
          builtins.filter
            (e: (builtins.tryEval e).success)
            (builtins.attrValues pkgs.darwin.apple_sdk.frameworks)
        );
        packages = [
          cargo-dist
        ];
        shellHook = "exec zsh";
      };
      devShells."x86_64-linux".default = with nixpkgs.legacyPackages."x86_64-linux"; mkShell rec {
        buildInputs = [
          libxkbcommon
          libGL
          xorg.libXcursor
          xorg.libXrandr
          xorg.libXi
          xorg.libX11
          zlib
          stdenv.cc.cc.lib
        ];
        packages = [
          cargo-dist
        ];
        LD_LIBRARY_PATH = "${lib.makeLibraryPath buildInputs}";
        shellHook = "exec zsh";
      };
    };
}
