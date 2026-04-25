# pic-nix-foot

Thin wrapper around [pic-nix](https://github.com/amanotk/pic-nix) that applies
the ion background B-field subtraction hack to the `foot` example.

The physics setup, source code, and scenario configurations all live in
`pic-nix` under `pic/example/foot/`. This repository only provides:

- `patches/ion-bg-subtraction.patch` — modifies the particle push to subtract
  the uniform background magnetic field from the Lorentz force acting on ions
  (but not electrons). Ions still feel the fluctuating component.
- `setup.sh` — clones pic-nix, checks out a given ref, and applies the patch.

See `pic-nix/pic/example/foot/README.md` for the physical parameters and
scenario descriptions.


# Usage

```
$ git clone git@github.com:amanotk/pic-nix-foot.git
$ cd pic-nix-foot
$ ./setup.sh
```

This clones `pic-nix` into `./pic-nix` and applies the patch.
Build the `foot` target yourself:
```
$ cmake -S pic-nix -B build
$ cmake --build build --target foot
```

## Options

```
$ ./setup.sh -r <ref> -d <directory>
```

- `-r, --ref REF` — pic-nix branch or tag to checkout (default: `develop`)
- `-d, --dir DIR` — pic-nix clone destination (default: `./pic-nix`)
