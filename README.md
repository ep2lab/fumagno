FUMAGNO
=======

[![DOI](https://zenodo.org/badge/86265405.svg)](https://zenodo.org/badge/latestdoi/86265405)

Matlab code that implements a 3D fully magnetized magnetic nozzle model. The code allows to simulate the properties of the expanding plasma in a prescribed magnetic field, under the assumption of full magnetization, quasineutrality, and no collisions.
Plasma potential, density velocity profiles can be simulated from the provided boundary condition upstream.

The code was developed at UC3M as part of the activities of the [EP2 research group](https://ep2.uc3m.es/), by [Mario Merino](https://mariomerino.uc3m.es/), with funding from the Spanish R&D National Plan (grant number ESP2016-75887-P).

You can find all the details of the model in *Mario Merino and Eduardo Ahedo, "Contactless steering of a plasma jet with a 3D magnetic nozzle", Plasma Sources Science and Technology 26, 095001 (2017)*.

## Installation

Installation requires simply that you 
[download FUMAGNO](https://github.com/ep2lab/fumagno/archive/master.zip) 
and add the base directory (the one that contains the `+fumagno` directory) to 
your Matlab path.

### Dependencies

A recent version of Matlab is needed to run the code. 
FUMAGNO has been developed in Matlab 2016a Academic version. 

FUMAGNO depends on other of our Matlab packages that you can download below:
account:
[magnetic_field](https://github.com/ep2lab/magnetic_field),
[fluid_plasma](https://github.com/ep2lab/fluid_plasma),
[utilities](https://github.com/ep2lab/utilities)
and
[constants_and_units](https://github.com/ep2lab/constants_and_units).
These packages must be installed and added to your Matlab path before
running FUMAGNO. 

## Usage

The basic workflow with FUMAGNO is as follows:

1. Create a magnetic_field object with the 3D field of the nozzle 
to be studied. The object must be of a subclass of `magnetic_field.element_3d`
2. Create the arrays with the points where the plasma properties will be
calculated. This can be done in two ways: by generating `X0,Y0`, the points
at the initial plane of the magnetic lines of interest, or by generating
`X,Y,Z`. The functions `x0y0_direct`, `x0y0_to_plane` and `x0y0_inverse` are 
used to compute the remaining arrays
3. Create the a `fluid_plasma` object and the initial condition functions
for the potential, velocities and densities in the inital plane.
4. Use `flow_solver` to compute the solution of the plasma properties
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

