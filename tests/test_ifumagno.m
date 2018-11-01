%{
Compendium of all tests for the function/class in the name of this file.

NOTES:
* The functions to be tested must be in the Matlab path. You can run the
  tests by executing 'runtests'. 
* Your working directory must be the
  directory where this test file is contained.  
%}
function tests = test_ifumagno
    clc;
    tests = functiontests(localfunctions);     
end

%----------------------------------------------------------------------

function test_geometry_propagation(~)
    field = magnetic_field.loop_3d('RL',5,'I',1,'ZL',0,'axis',[0,0,1]);
    odeoptions.AbsTol = 1e-8;
    
    [O,I] = fumagno.propagate.x0y0_direct('field',field,'x0',0.2,'y0',0.1,'ds',0.5,'n_steps',[40 100],'odeoptions',odeoptions); 
    [O,I] = fumagno.propagate.x0y0_direct('field',field,'x0',[0.2, 0.3],'y0',[0.1, 0.5],'ds',0.5,'n_steps',[40 100],'odeoptions',odeoptions); 
    [O,I] = fumagno.propagate.x0y0_direct('field',field,'x0',[0.2 0.3; 0.1 0.4],'y0',[0.1 0.3; 0.2 0.4],'ds',0.5,'n_steps',[40 100],'odeoptions',odeoptions); 
end


function test_solver(~)
    field = magnetic_field.loop_3d('RL',5,'I',1,'ZL',0,'axis',[0,0,1]);
    odeoptions.AbsTol = 1e-8;
    electrons = fluid_plasma.species('label','e','m',0,'q',-1,'T0',1,'gamma',1.15);
    ions = fluid_plasma.species('label','ions','m',1,'q',1,'gamma',1.15);
    plasma = fluid_plasma.plasma('electrons',{electrons},'ions',ions);
    
    [O,I] = fumagno.propagate.x0y0_direct('field',field,'x0',[0.2 0.3; 0.1 0.4],'y0',[0.1 0.3; 0.2 0.4],'ds',0.5,'n_steps',[100 40],'odeoptions',odeoptions);                
    [Of,If] = fumagno.solve.flow_solver('plasma',plasma,'B',O.B,'B0',O.B0,'X0',O.X0,'Y0',O.Y0,'Nlib',[100000 1000]); 
end
