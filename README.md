FUMAGNO
=======

Matlab code that implements the fully magnetized magnetic nozzle model
presented in meri16a -- Fully magnetized plasma flow in a magnetic nozzle.

Workflow examples
-----------------

The basic workflow with FUMAGNO is as follows. 

1. Create a magnetic_field object with the 3D field of the nozzle 
to be studied. The object must be of a subclass of magnetic_field.element_3d
2. Create the arrays with the points where the plasma properties will be
calculated. This can be done in two ways: by generating X0,Y0, the points
at the initial plane of the magnetic lines of interest, or by generating
X,Y,Z. The functions x0y0_direct and x0y0_inverse are used to compute the
remaining arrays
3. Create the a fluid_plasma object and the initial condition functions
for the potential, velocities and densities in the inital plane.
4. Use flow_solver to compute the solution of the plasma properties

![Example workflow diagram](/docs/figs/fumagno-workflow.svg "FUMAGNO example workflow")

See the tests in tests/fumagno_test.m for examples of calls to these functions.

Dependencies
------------

For FUMAGNO to work, the matlab packages 
[magnetic_field](https://github.com/mariomerinomartinez/magnetic_field)
and
[fluid_plasma](https://github.com/mariomerinomartinez/fluid_plasma),
must be in the user's path.

Updates
-----

This readme was last updated: 20170326
