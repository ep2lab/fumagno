FUMAGNO
=======

[![DOI](https://zenodo.org/badge/86265405.svg)](https://zenodo.org/badge/latestdoi/86265405)

Matlab code that implements a 3D fully magnetized magnetic nozzle model.
You can find all the details of the model in [TBD].

## Installation

Installation requires simply that you 
[download FUMAGNO](https://github.com/ep2lab/fumagno/archive/master.zip) 
and add the base directory (the one that contains the `+fumagno` directory) to 
your Matlab path.

### Dependencies

A recent version of Matlab is needed to run the code. 
FUMAGNO has been developed in Matlab 2016a Academic version. 

FUMAGNO depends on other Matlab packages that you can download from my GitHub
account:
[magnetic_field](https://github.com/ep2lab/magnetic_field),
[fluid_plasma](https://github.com/ep2lab/fluid_plasma),
[utilities](https://github.com/ep2lab/utilities),
[constants_and_units](https://github.com/ep2lab/constants_and_units)
and [akiles](https://github.com/ep2lab/akiles). These packages must be installed and added to your Matlab path before
running FUMAGNO.


**Important comment**: I frequently update these repositories and some of the
FUMAGNO dependencies in those packages may break. If you detect any such
problem, please contact me through our 
[website](http://ep2lab.uc3m.es/).

## Usage

This package solves to different problems, both in the Fully Magnetized Ions Limit (FMIL):
 1. The fluid equations in the plasma expansion, considering that electrons cool down isotropically (I-FUMAGNO).
 2. The fluid-kinetic model for plasma expansion, using a kinetic approach to solve the cooling of electrons (K-FUMAGNO).

These models are both contained in the same repository as they use mostly the same matlab functions.

### Fluid Model

I-FUMAGNO is an improved version of the prior model FUMAGNO. It solves the problem by interpolation as follow:

1. Create a magnetic_field object with the 3D field of the nozzle 
to be studied. The object must be of a subclass of `magnetic_field.element_3d`
2. Create the arrays with the points where the plasma properties will be
calculated. This can be done  by generating `X0,Y0`, the points
at the initial plane of the magnetic lines of interest. The function `x0y0_direct`is used to compute the remaining arrays.
3. Create the a `fluid_plasma` object and the initial condition functions
for the potential, velocities and densities in the inital plane.
4. Use `flow_solver` to compute the solution of the plasma properties. The solver contains the following: 
    4.1 A function `library` that computes an interpolation library for a random vector of plasma densities 
    4.2 A function `interpolation` that finds the solution for each magnetic line by interpolation 
    4.3 The solution is postprocessed in order to impose the conditions at the throat 
5. Use the output as you see fit (save, plot, etc)

![Example workflow diagram](/docs/figs/fumagno-workflow.png "FUMAGNO example workflow")

See the tests in `tests/fumagno_test.m` for examples of calls to these functions.

## Contributing

If you have any comments for improvement or 
are interested in contributing to the continued 
development of this or any of my other codes, you can contact us
through our [website](http://ep2.uc3m.es/). 

For updates and news, follow us on Twitter: [@ep2lab.](https://twitter.com/ep2lab).

## Acknowledging 

This program is released as open source in the hope that it will be useful to
other people. If you find FUMAGNO useful and/or use it in any of your works, we kindly ask you
to acknowledge it by citing the following article (preferred):

> **Mario Merino and Eduardo Ahedo**, 
"_Contactless steering of a plasma jet with a 3D magnetic nozzle_", 
Plasma Sources Science and Technology 26, 095001 (2017) 
[DOI](http://iopscience.iop.org/article/10.1088/1361-6595/aa8061)

and/or by citing the code directly as:

> Mario Merino (2017). mariomerinomartinez/fumagno: First release DOI:10.5281/zenodo.593787

## License

Copyright (c) 2017 Mario Merino and Eduardo Ahedo. The software is released as open source with the [MIT License](LICENSE.md).

