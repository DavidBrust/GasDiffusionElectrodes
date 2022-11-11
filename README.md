# GasDiffusionElectrodes
Project package for implementation of a multiphysics simulation of gas
diffusion electrodes (GDEs) for CO2 electrolysis.

In particular, the 1D, isothermal, steady state model of a GDE used in 
CO2 electrolysis is to be reimplemented, presented in:
[Weng, L. C., et al. (2018).](https://pubs.rsc.org/en/content/articlehtml/2018/cp/c8cp01319e)
(Open Access)

The simulation will utilize the finite volume method as implemented in the
[VoronoiFVM.jl](https://github.com/j-fu/VoronoiFVM.jl) package.

The code development takes place in Pluto notebooks, located in the notebooks 
subfolder.