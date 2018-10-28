%{
----------------- K-Fumagno ----------------------
This code explains how to obtain the solution of the fluid-kinetic problem
in the FMIL, by using both, akiles2d and fumagno packages.

Akiles2d: solves the fluid-kinetic problem in the UMIL
Fumagno:  solves the fluid problem (isotropic electrons) in the FMIL

It has been demonstrated that the model behind akiles2d is suitable for
both, the FMIL and the UMIL. Then, akiles2d can be use to solve the FMIL
for an arbitrary magnetic line.

As in the FMIL, the equations can be solved for each line independently,
the solution can be interpolated for a set of magnetic lines (previously
defined by fumagno).

The result is the fluid-kinetic solution in the FMIL, for a certain
magnetic geometry.

Author: Judit Nuez
Date: 05/10/2018

%}

clear all; close all; clc

userdata.ions.model = 'cold'; % the ion model to use
userdata.ions.mu = Inf;
userdata.guess.npoints = 100; 
userdata.guess.h = [linspace(1,400,userdata.guess.npoints-1),Inf].'; % column; independent variable: plume characteristic radius at each test point. The first value must be 1; the final value must be infinity
userdata.guess.r = zeros(1,userdata.guess.npoints).'; % column; corresponding values of the radius for each test point
userdata.guess.phi = linspace(0,-4,userdata.guess.npoints).'; % column; potential at each test point. Must be 0 at origin
userdata.guess.ne00p = 0.51; % density of the (vz > 0) electrons at the origin
userdata.ions.chi= 0.02;

% ------------------------------- AKILES -----------------------------------% 
% Akiles2d is solved for a randorm 'h' vector
 
             [data,solution] = akiles2d.akiles2d(fullfile('fixtures/akiles2d_simrc.m'),userdata); 

% ------------------------------- FUMAGNO ----------------------------------%           
  
% Compute gamma to get the same total potential drop than the kinetic
% solution

Tad0 = (solution.electrons.Tz(1)+solution.electrons.Tr(1)+solution.electrons.Ttheta(1))./3;
gamma = -solution.phi(end-1)*(-solution.phi(end-1)-Tad0)^(-1);

% Definition of the magnetic field
  field = magnetic_field.loop_3d('RL',10,'I',1,'ZL',0);
  electrons = fluid_plasma.species('label','e','m',0,'q',-1,'T0',1,'gamma',gamma);
  ions = fluid_plasma.species('label','ions','m',1,'q',1,'gamma',gamma);
  plasma = fluid_plasma.plasma('electrons',{electrons},'ions',ions);
  odeoptions = odeset;
  odeoptions.AbsTol = 1e-8;

% 2D definition of the initial conditions
  X0 = [0 0.1 0.2
        0.8 0.6 0.7];
  Y0 = [0 0.1 0.2
        0.8 0.6 0.7];   

% Magnetic geometry propagation

        [O,I]   = fumagno.propagate.x0y0_direct('field',field,'x0',X0,'y0',Y0,'ds',1,'n_steps',[100 0],'odeoptions',odeoptions); % array call
        
% Solve Fumagno (isothermal) to compare with Fumagno (Kinetic)        
        
        [Of,If] = fumagno.solve.flow_solver('plasma',plasma,'B',O.B,'B0',O.B0,'X0',O.X0,'Y0',O.Y0,'Nlib',[1000 0]); % a more complete call, but using default initial conditions

% -------------------------------K-FUMAGNO ----------------------------------%
% Get the solution for the Fluid-kinetic problem in the FMIL
% Interpolation to find the solution for each line
        
        [Kinetic_solution] = fumagno.kinetic.interpolation (O,Of,solution);
        
  