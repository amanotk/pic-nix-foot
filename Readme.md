# About
This setup models the dynamics of a collisionless shock transition layer with a three-component homogeneous plasma with the periodic boundary condition in all three directions. The system consists of the reflected ions, the incoming core (upstream) ions and the background electrons, all represented by isotropic Maxwellian distributions in the rest frame of each component. The simulation frame corresponds to the rest frame of electrons.


# Physical Parameters
The following parameters should be defined in the confirugation file:
- `cc` : speed of light $c$
- `wp` : electron plasma frequency $\omega_{pe}$
- `mime` : ion-to-electron mass ratio $m_i/m_e$
- `mach` : Alfven Mach number $M_A$
- `theta` : polar angle of the ambient magnetic field with respect to the x axis $\theta$
- `phi` : azimuthul angle of the ambient magnetic field with respect to the x axis $\phi$
- `sigma` : electron cyclotron-to-plasma frequency squared $\sigma = \Omega_{ce}^2/\omega_{pe}^2$
- `alpha` : density of the reflected ion beam normalized to the total density $\alpha = n_r/n_0$
- `betae` : electron plasma beta $\beta_e = 2 v_{th,e}^2/V_{A,e}^2$
- `betai` : core ion plasma beta $\beta_i = 2 v_{th,i}^2/V_{A,i}^2$
- `betar` : reflected ion plasma beta $\beta_r = 2 v_{th,r}^2/V_{A,i}^2$

The ambient magnetic feld is normalized to unity $B_0 = 1$. The three components of the ambient magnetic field are given by
```math
\begin{aligned}
B_{0,x} &= B_0 \cos \theta \\
B_{0,y} &= B_0 \sin \theta \cos \phi \\
B_{0,z} &= B_0 \sin \theta \sin \phi
\end{aligned}
```

Note that the electron velocity in the shock rest frame is written as $V_e = -(1 - 2\alpha) V_s$ where $V_s$ is the three shock speed calculated from the four shock speed $U_s = M_A V_{A,i}$.
The density ratio $\alpha$ is defined in the shock rest frame, which we need to convert to the simulation frame to yield
```math
\alpha' = \alpha \frac{1 + (1 - 2\alpha) V_s^2}{1 - (1 - 2\alpha)^2 V_s^2}
```
Similarly, the core and reflected ion drift velocities are given by
```math
\begin{aligned}
  V_i' = \frac{-2 \alpha V_s}{1 - (1 - 2\alpha) V_s^2}\\
  V_r' = \frac{+2 (1-\alpha) V_s}{1 + (1 - 2\alpha) V_s^2}
\end{aligned}
```
These drift velocities are always parallel to the x axis.
It is easy to check that the above quantities satisfy the charge and current neutrality conditions.  
Note that the electron and ion Alfven speeds are defined by $V_{A,e} = B_0 / \sqrt{n_0 m_e}$ and $V_{A,i} = B_0 / \sqrt{n_0 m_i}$, respectively.


# Compiling and Executing the Code
The procedure is very similar to the description available [here](https://github.com/amanotk/pic-nix).  
Note that the [pic-nix](https://github.com/amanotk/pic-nix) repository will automatically be cloned in the build directory with the following procedure.

## Clone
Clone the repository via:
```
$ git clone git@github.com:amanotk/pic-nix-foot.git
```

## Compile
Compile the code to obtain an executable `main.out` as follows:
```
$ cd pic-nix-foot
$ cmake -S . -B build \
	-DCMAKE_CXX_COMPILER=mpicxx \
	-DCMAKE_CXX_FLAGS="-O3 -fopenmp"
$ cmake --build build
```
The executable file will be found in the `build` directory.


## Run

### Example 1: `buneman`
This example provides a setup for the Buneman instability, an electrostatic ion-electron beam instability.  

Go to `buneman` directory:
```
$ cd buneman
```
and run the code, for instance, via:
```
$ export OMP_NUM_THREADS=2
$ mpiexec -n 8 ../main.out -e 86400 -t 200 -c config.json
```
The data files are written to `data` directory by default.  

Set the environment variable `PICNIX_DIR` to the path to the [pic-nix](https://github.com/amanotk/pic-nix) directory.
Then, run the following command to examine the simulation results:
```
$ python batch.py data/profile.msgpack
```
You will find image files generated in the same directory.


### Example 2: `aflven`
This example provides a setup for an instability where an ion beam propagating parallel to the ambient magnetic field generates Alfven waves via the cyclotron resonance.  

Run with the same procedure with `buneman`, but this time with the following command to run up to $\omega_{pe} t = 5000$
```
$ mpiexec -n 16 ../main.out -e 86400 -t 5000 -c config.json
```
The rest is the same as `buneman`.
For details of the problem, see Hoshino & Terasawa (1985).


# References
- <div class="csl-entry">Hoshino, M., &#38; Terasawa, T. (1985). Numerical Study of the Upstream Wave Excitation Mechanism, 1. Nonlinear Phase Bunching of Beam Ions. <i>Journal of Geophysical Research</i>, <i>90</i>(A1), 57–64. https://doi.org/10.1029/JA090iA01p00057</div>
- <div class="csl-entry">Matsukiyo, S., &#38; Scholer, M. (2003). Modified two-stream instability in the foot of high Mach number quasi-perpendicular shocks. <i>Journal of Geophysical Research: Space Physics</i>, <i>108</i>(A12), 1459. https://doi.org/10.1029/2003JA010080</div>
- <div class="csl-entry">Matsukiyo, S., &#38; Scholer, M. (2006). On microinstabilities in the foot of high Mach number perpendicular shocks. <i>Journal of Geophysical Research: Space Physics</i>, <i>111</i>(6), 1–10. https://doi.org/10.1029/2005JA011409</div>
- <div class="csl-entry">Amano, T., &#38; Hoshino, M. (2009). Nonlinear evolution of Buneman instability and its implication for electron acceleration in high Mach number collisionless perpendicular shocks. <i>Physics of Plasmas</i>, <i>16</i>(10), 102901. https://doi.org/10.1063/1.3240336</div>
- <div class="csl-entry">Muschietti, L., &#38; Lembège, B. (2013). Microturbulence in the electron cyclotron frequency range at perpendicular supercritical shocks. <i>Journal of Geophysical Research: Space Physics</i>, <i>118</i>(5), 2267–2285. https://doi.org/10.1002/jgra.50224</div>
