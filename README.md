# UQCreator - a HPC uncertainty quantification framework for hyperbolic equations

The UQCreator is a powerful Uncertainty Quantification framework for any kind of user defined hyperbolic conservation equation. Even though the framework in its core is written in C++ the user can comfortably define the uncertainies, the set of equations and initial parameters from the python interface without loosing any performance.

The core features of the UQCreator are:
- intrusive methods: stochastic Galerkin, (adaptive) IPM
- non-intrusive methods: (sparse-grid) stochastic collocation, MLMC
- MPI and OpenMP support
- support for high dimensional stochastic spaces
- triangular 2-D meshes


## Installation:
The code is available through multiple channels with no difference in the versions between `conda`, `docker` and the current master branch.

Currently only Linux devices are supported.

### Conda-forge
`conda` is an open source package and environment management system commonly used for installing multiple versions of different software packages for seamless interoperability. `conda-forge` is a community driven channel of additionally installable packages.

Having `conda` installed on your system, you can just install the framework using:

     conda install uqcreator -c conda-forge

---
### Docker
`docker` is a virtualization-like environment that bundles a lightweight fully functional operating system into a single shipable container.
If `docker` is available on your system, calling:

     docker run -ti --rm -v $(pwd):/mnt -w /mnt --shm-size=1g uqcreator:latest

will automatically download latest the UQCreator version and start a terminal session with the current working directory as a shared folder with the container.

---
### Singularity
Similarly to `docker`, `singularity` is a virtualization-like environment, but with focus on HPC environments as it does not require a running daemon and no root permissons, which generally are unwanted privileges on a HPC infrastructure.

     Coming soon

--- 
### Installing from source:
Dependencies:
- Compiler with C++17 support
- cmake >= v3.5
- LAPACK
- vtk 
- git
- ninja or make

Initialize all submodules (internet connections required):

     git submodule update --init

Afterwards compile the code with CMake and a suitable generator.
In case of `ninja`:

     mkdir build && cd build
     cmake -G Ninja -DCMAKE_BUILD_TYPE=Release ../src
     ninja && ninja install
---
## Citing
When using/refering to this work in your research, please consider giving proper attribution by citing the following publication:

- Jonas Kusch, Jannick Wolters, Martin Frank, "[Intrusive acceleration strategies for uncertainty quantification for hyperbolic systems of conservation laws](https://doi.org/10.1016/j.jcp.2020.109698)", Journal of Computational Physics, 419, S. 109698, 2020


