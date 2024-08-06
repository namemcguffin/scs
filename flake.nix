# TODO: deduplicate dev deps (currently just cargo-dist)

{
  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixpkgs-unstable";
  };
  outputs = { self, nixpkgs }:
    {
      devShells."aarch64-darwin".default = with nixpkgs.legacyPackages."aarch64-darwin"; mkShell {
        buildInputs = [
          darwin.apple_sdk.frameworks.Foundation
          darwin.apple_sdk.frameworks.AppKit
          darwin.apple_sdk.frameworks.Vision
          darwin.apple_sdk.frameworks.AVFoundation
          darwin.apple_sdk.frameworks.MetalKit
          darwin.libiconv
        ];
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
        ];
        packages = [
          cargo-dist
        ];
        LD_LIBRARY_PATH = "${lib.makeLibraryPath buildInputs}";
        shellHook = "exec zsh";
      };
    };
}
