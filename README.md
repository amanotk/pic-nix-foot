# pic-nix-foot

Thin wrapper around [pic-nix](https://github.com/amanotk/pic-nix) that applies
the ion background B-field subtraction hack and builds the `foot` example.

The physics setup, source code, and scenario configurations all live in
`pic-nix` under `pic/example/foot/`. This repository only provides:

- `patches/ion-bg-subtraction.patch` — modifies the particle push to subtract
  the uniform background magnetic field from the Lorentz force acting on ions
  (but not electrons). Ions still feel the fluctuating component.
- `build.sh` — clones pic-nix, applies the patch, and builds the `foot` target.

See `pic-nix/pic/example/foot/README.md` for the physical parameters and
scenario descriptions.


# Usage

## Build

Clone and build:
```
$ git clone git@github.com:amanotk/pic-nix-foot.git
$ cd pic-nix-foot
$ ./build.sh
```

This will clone `pic-nix` into the current directory (or use an existing one
if `PICNIX_DIR` is set), apply the patch, and build only the `foot` target.

Environment variables:
- `PICNIX_DIR` — path to an existing pic-nix checkout (default: auto-clone)
- `PICNIX_REF` — git ref to checkout in pic-nix (default: `develop`)
- `CMAKE_CXX_COMPILER` — C++ compiler (default: `mpicxx`)
- `CMAKE_CXX_FLAGS` — compiler flags (default: `-O2 -fopenmp`)

## Run

After building, go to a scenario directory and run:
```
$ cd build/pic/example/foot/<scenario>
$ mpiexec -n 8 ./main.out -e 86400 -t <duration> -c config.toml
```

Set `PICNIX_DIR` to the pic-nix root and run the quicklook script:
```
$ python quicklook.py data/profile.msgpack
```

### Scenarios

| Scenario | Description | Suggested `-t` |
|----------|-------------|----------------|
| buneman  | Electrostatic ion-electron beam instability | 200 |
| alfven   | Electromagnetic ion-cyclotron beam instability | 5000 |
| weibel   | Weibel instability with magnetized electrons | 5000 |


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
