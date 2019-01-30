# Acoustic Scattering by Poroelastic Plates

A boundary element method is developed to solve the non-homogeneous Helmholtz equation subjected to boundary conditions related to the vibration of a poroelastic plate. This model represents a fluid-structure interaction problem where one can investigate how porosity and elasticity affect the noise scattered by the plate trailing edge.

For an overview of the technical details, please see the following article:

"Numerical solution of acoustic scattering by finite perforated elastic plates", Proc. R. Soc. A 472 (2016)20150767 
https://royalsocietypublishing.org/doi/10.1098/rspa.2015.0767

## SET UP

The boundary integral equations are implemented in fortran and wrapped to python. Use f2py to create a module from fortran and use in python. The files 'betas.txt' and 'modes.txt' are examples of a modal basis for a plate clamped at the leading edge and free at the trailing edge. Other boundary conditions can be used in the code.
