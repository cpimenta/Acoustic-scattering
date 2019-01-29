# Acoustic-scattering

A boundary element method is developed to solve the non-homogeneous Helmholtz equation subjected to boundary conditions related to the plate vibration. This model represent a fluid-structure interaction problem to investigate how porosity and elasticity can affect the noise scattered at far-field. Poroelastic materiais can be studied by the combination of porous and elastic effects.

For an overview of the technical details, please see the following article:

"Numerical solution of acoustic scattering by finite perforated elastic plates", Proc. R. Soc. A 472 (2016)20150767 
https://royalsocietypublishing.org/doi/10.1098/rspa.2015.0767

## SET UP

All integral equations is implemented in fortran and wrapped to python. Use f2py to create a module from fortran and use in python.
The files betas.txt and modes.txt is an example of a modal basis for a clamped plate at leading edge and free at trailing edge. Different geometries and boundary conditions can be used in the code.
