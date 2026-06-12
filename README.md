# SF-DDFT
We developed SF-DDFT (Stochastic Fluctuations in Dynamic Density Functional Theory), a Fortran 90 software package that incorporates stochastic noise into DDFT using a pseudospectral method.

The code incorporates stochastic fluctuations into the DDFT equations to capture mesoscale phenomena beyond mean-field descriptions. Time evolution is performed using a pseudospectral method, which enables efficient and accurate evaluation of spatial derivatives in Fourier space while handling nonlinear terms in real space. The package is suitable for studying time-dependent behavior of inhomogeneous polymer and soft-matter systems, including phase separation, interfacial dynamics, and fluctuation-driven processes. Its modular structure allows flexible control of physical parameters, numerical resolution, and noise strength, making it a useful tool for exploring stochastic effects in density-based continuum models.

# Overview of the SF-DDFT package
1. Block copolymers in melt (`Main_Diblock`)   
2. Homopolymer Blend in melt (`Main_blend`)

# Directory Structure
.
├── General_Routines
│   ├── mod_iterate.f90
│   ├── mod_mt.f90
│   ├── mod_propagate.f90
│   └── sub_wrap.f90
├── Geometry_Modules
│   └── Bulk
│       ├── fftw3.f
│       └── mod_fourier.f90
├── Interaction_Modules
│   ├── mod_FloryHuggins.f90
│   ├── mod_FloryHuggins_Solution.f90
│   └── mod_onecomponent_virialexpansion.f90
├── Main_Blend
│   ├── direction.py
│   ├── extract_last_sample.f90
│   ├── extract.sh
│   ├── input_general
│   ├── input_interactions
│   ├── input_mobilities_d
│   ├── main_blend.f90
│   ├── Makefile
│   ├── Mobilityentangledprint
│   ├── mod_global.f90
│   ├── run.sh
│   └── S_q.py
├── Main_Diblock
│   ├── direction.py
│   ├── extract_last_sample.f90
│   ├── input_general
│   ├── input_interactions
│   ├── input_mobilities_d
│   ├── input_mobilities-KG400
│   ├── main_multiblock.f90
│   ├── Makefile
│   ├── Mobility-delta-print.txt
│   ├── mod_global.f90
│   ├── plots.ipynb
│   ├── run.sh
│   └── S_q.py
├── Polymer_Modules
│   ├── mod_diblock.f90
│   ├── mod_homopolymer.f90
│   ├── mod_multiblock.f90
│   ├── mod_solvent.f90
│   └── mod_triblock.f90
└── README.md

The code can be used in one-, two-, and three-dimensional systems.

# Usage

Clone the repository:
```
git clone https://github.com/jafarcheraghalizadeh/SF-DDFT.git
cd DDFT
```

Build and run using the Main_Blend code as an example:
```
cd /Main_Blend
bash run.sh
```
It runs 20 independent simulations using `nohup` in the background for the interaction parameter $\chi N = 0$, starting from random field configurations.



# Important parameters in the code

1. `KT`  in mod_global shows the noise strength.

2. Line 15 in `input_general`: set to 4 to read the initial configuration from a file; set to 3 to generate a random field (uploaded version is set to 3).



# Information for Other files

1. `S_q.py` :

   For calculating the structure factor for given `KT` and different time steps
2. `direction.py` :

    For analysing the structure tensor
3. `plots.ipynb`

   For plotting Energy, The Structure factors $S(q)$, and the average growth rate $\bar{R}(q)$.
4. `extract_last_sample.f90`
   
   It extracts a sample at a specific time step. This can be done by running `extract.sh`
# License
Jafar Cheraghalizadeh (jcheragh@uni-mainz.de)

Friederike Schmid (friederike.schmid@uni-mainz.de)

Johannes Gutenberg-Universität Mainz, Institute of Physics, Schmid group
