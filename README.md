# pic-nix-foot

Thin wrapper around [pic-nix](https://github.com/amanotk/pic-nix) that applies
the ion background B-field subtraction hack to the `foot` example.

The physics setup, source code, and scenario configurations all live in
`pic-nix` under `pic/example/foot/`. This repository only provides:

- `patches/ion-bg-subtraction.patch` — modifies the particle push to subtract
  the uniform background magnetic field from the Lorentz force acting on ions
  (but not electrons). Ions still feel the fluctuating component.
- `setup.sh` — clones pic-nix, checks out a given ref, and applies the patch.

The patch is generated against the `develop` branch of pic-nix.
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


# References
- Hoshino, M., & Terasawa, T. (1985). Numerical Study of the Upstream Wave
  Excitation Mechanism, 1. Nonlinear Phase Bunching of Beam Ions.
  *Journal of Geophysical Research*, *90*(A1), 57–64.
  https://doi.org/10.1029/JA090iA01p00057
- Matsukiyo, S., & Scholer, M. (2003). Modified two-stream instability in
  the foot of high Mach number quasi-perpendicular shocks.
  *Journal of Geophysical Research: Space Physics*, *108*(A12), 1459.
  https://doi.org/10.1029/2003JA010080
- Matsukiyo, S., & Scholer, M. (2006). On microinstabilities in the foot of
  high Mach number perpendicular shocks.
  *Journal of Geophysical Research: Space Physics*, *111*(6), 1–10.
  https://doi.org/10.1029/2005JA011409
- Amano, T., & Hoshino, M. (2009). Nonlinear evolution of Buneman instability
  and its implication for electron acceleration in high Mach number collisionless
  perpendicular shocks. *Physics of Plasmas*, *16*(10), 102901.
  https://doi.org/10.1063/1.3240336
- Muschietti, L., & Lembège, B. (2013). Microturbulence in the electron
  cyclotron frequency range at perpendicular supercritical shocks.
  *Journal of Geophysical Research: Space Physics*, *118*(5), 2267–2285.
  https://doi.org/10.1002/jgra.50224
