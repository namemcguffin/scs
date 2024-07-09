# TODO: add linux output
{
  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixpkgs-unstable";
  };
  outputs = { self, nixpkgs }:
    let
      system = "aarch64-darwin";
      pkgs = nixpkgs.legacyPackages.${system};
    in
    {
      devShells.${system}.default = pkgs.mkShell {
        buildInputs = [
          pkgs.darwin.apple_sdk.frameworks.Foundation
          pkgs.darwin.apple_sdk.frameworks.AppKit
          pkgs.darwin.apple_sdk.frameworks.Vision
          pkgs.darwin.apple_sdk.frameworks.AVFoundation
          pkgs.darwin.apple_sdk.frameworks.MetalKit
          pkgs.darwin.libiconv
        ];
        shellHook = "exec zsh";
      };
    };
}
