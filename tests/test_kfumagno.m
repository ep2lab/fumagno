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


function tests = test_kfumagno
    clc;
    tests = functiontests(localfunctions);     
end

%----------------------------------------------------------------------

function test_geometry_propagation(~)
    
    field = magnetic_field.loop_3d('RL',5,'I',1,'ZL',0,'axis',[0,0,1]);
    odeoptions.AbsTol = 1e-8;
    
    X0 = [0 0.1 0.2; 0.8 0.6 0.7];
    Y0 = [0 0.1 0.2; 0.8 0.6 0.7];   
    
    [O,I] = fumagno.propagate.x0y0_direct('field',field,'x0',X0,'y0',Y0,'ds',0.5,'n_steps',[100 40],'odeoptions',odeoptions);                

    userdata.ions.model = 'cold'; % the ion model to use
    userdata.ions.mu = Inf;
    userdata.guess.npoints = 100; 
    userdata.guess.h = [linspace(1,100,userdata.guess.npoints-1),Inf].'; % column; independent variable: plume characteristic radius at each test point. The first value must be 1; the final value must be infinity
    userdata.guess.r = zeros(1,userdata.guess.npoints).'; % column; corresponding values of the radius for each test point
    userdata.guess.phi = linspace(0,-4,userdata.guess.npoints).'; % column; potential at each test point. Must be 0 at origin
    userdata.guess.ne00p = 0.51; % density of the (vz > 0) electrons at the origin
    userdata.ions.chi= 0.02;
   
    path = 'fixtures/akiles2d_simrc.m'; 
    
    [Kinetic_solution] = fumagno.solve.kinetic_solver (userdata,path,'B',O.B,'B0',O.B0,'X0',O.X0,'Y0',O.Y0);
        
    
end

