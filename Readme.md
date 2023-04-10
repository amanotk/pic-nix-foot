# About

This setup models the dynamics of a collisionless shock transition layer with a three-component homogeneous plasma with the periodic boundary condition in all three directions. The system consists of the reflected ions, the incoming core (upstream) ions and the background electrons, all represented by isotropic Maxwellian distributions in the rest frame of each component. The simulation frame corresponds to the rest frame of electrons.

# Compiling and Executing the Code

## Clone Repositories
Clone repositories to local directory via:
```
$ git clone git@github.com:amanotk/pic-nix.git
$ git clone git@github.com:amanotk/pic-nix-foot.git
$ export PICNIX_DIR=${PWD}/pic-nix
```
The environment variable `PICNIX_DIR` will be used to compile the code and also as a search path for running diagnostic python scripts.

## Compile
A proper C++ compiler and its compiler flags should be specified with `cmake`.  
The following example assumes `mpicxx` as a compiler with OpenMP enabled.
```
$ cd pic-nix-foot
$ mkdir build
$ cd build
$ cmake .. \
	-DPICNIX_DIR=${PICNIX_DIR} \
	-DCMAKE_CXX_COMPILER=mpicxx \
	-DCMAKE_CXX_FLAGS="-O3 -fopenmp"
$ make
```
You will find `main.out` as an executable in the working directory.

## Run
You can now execute `main.out` using `mpiexec` (or `mpirun`). Here is an example:
```
$ export OMP_NUM_THREADS=4
$ mpiexec -n 4 ./main.out -e 86400 -t 100 -c config.json
```
Some options are:
- `-c` or `--config` : configuration file
- `-t` or `--tmax`   : maximum physical time in simulation unit (default 1.79769e+308)
- `-e` or `--emax`   : maximum elapsed time in sec (default 3600)

Note that you may specify the environment variable `OMP_NUM_THREADS` to control the number of threads per process for OpenMP parallelization.


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
With these parameters, the core and reflected ion drift velocities are given by
```math
\begin{aligned}
	V_{d,i}/V_{A,i} &= - M_A \alpha, \\
	V_{d,r}/V_{A,i} &= + M_A (1 - \alpha).
\end{aligned}
```
These drift velocities are always parallel to the x axis.  
Note that the electron and ion Alfven speeds are defined by $V_{A,e} = B_0 / \sqrt{n_0 m_e}$ and $V_{A,i} = B_0 / \sqrt{n_0 m_i}$, respectively.

# Examples
## `aflven`
An ion beam propagating parallel to the ambient magnetic field can generate an Alfven wave via the cyclotron resonance.  
The subdirectory `alfven` provides an example setup for this problem in 1D.  

Go to `alfven` directory:
```
$ cd alfven
```
and run the coe, for instance, via:
```
$ mpiexec -n 16 ../main.out -e 86400 -t 5000 -c config.json
```
Then, run the following command on the same directory to generate plots to examine the simulation results:
```
$ python batch.py profile.msgpack
```
For details of the problem, see, Hoshino & Terasawa (1985).


# References
- <div class="csl-entry">Hoshino, M., &#38; Terasawa, T. (1985). Numerical Study of the Upstream Wave Excitation Mechanism, 1. Nonlinear Phase Bunching of Beam Ions. <i>Journal of Geophysical Research</i>, <i>90</i>(A1), 57–64. https://doi.org/10.1029/JA090iA01p00057</div>
- <div class="csl-entry">Matsukiyo, S., &#38; Scholer, M. (2003). Modified two-stream instability in the foot of high Mach number quasi-perpendicular shocks. <i>Journal of Geophysical Research: Space Physics</i>, <i>108</i>(A12), 1459. https://doi.org/10.1029/2003JA010080</div>
- <div class="csl-entry">Matsukiyo, S., &#38; Scholer, M. (2006). On microinstabilities in the foot of high Mach number perpendicular shocks. <i>Journal of Geophysical Research: Space Physics</i>, <i>111</i>(6), 1–10. https://doi.org/10.1029/2005JA011409</div>
- <div class="csl-entry">Amano, T., &#38; Hoshino, M. (2009). Nonlinear evolution of Buneman instability and its implication for electron acceleration in high Mach number collisionless perpendicular shocks. <i>Physics of Plasmas</i>, <i>16</i>(10), 102901. https://doi.org/10.1063/1.3240336</div>
- <div class="csl-entry">Muschietti, L., &#38; Lembège, B. (2013). Microturbulence in the electron cyclotron frequency range at perpendicular supercritical shocks. <i>Journal of Geophysical Research: Space Physics</i>, <i>118</i>(5), 2267–2285. https://doi.org/10.1002/jgra.50224</div>
