%{
Compendium of all tests for the function/class in the name of this file.
You can run the tests by executing runtests. You must add the package to
your path first. 

%----------------------------------------------------------------------
Author: Mario Merino
Date: 20170326
%----------------------------------------------------------------------
%}
function tests = fumagno_test
    tests = functiontests(localfunctions);
end

%----------------------------------------------------------------------
%----------------------------------------------------------------------
%----------------------------------------------------------------------
 
function test_x0y0_direct(t)    
    % set up required fixture data
    field = magnetic_field.loop_3d('axis',[0,0,1]);
    odeoptions = odeset;
    odeoptions.AbsTol = 1e-8;
    
    [O,I] = fumagno.propagate.x0y0_direct; % default call
    [O,I] = fumagno.propagate.x0y0_direct('field',field,'x0',0.1,'y0',0.2); % scalar call
    [O,I] = fumagno.propagate.x0y0_direct('field',field,'x0',[0.1;0.15],'y0',[0.2;0.22]); % vector call
    [O,I] = fumagno.propagate.x0y0_direct('field',field,'x0',rand(3,4),'y0',rand(3,4)); % array call
    [O,I] = fumagno.propagate.x0y0_direct('field',field,'x0',rand(3,4),'y0',rand(3,4),'ds',0.01,'n_steps',53,'odeoptions',odeoptions); % full call
    [O,I] = fumagno.propagate.x0y0_direct('n_steps',1); % pathologic cases
    [O,I] = fumagno.propagate.x0y0_direct('n_steps',2); 
end 

function test_x0y0_inverse(t)    
    % set up required fixture data
    field = magnetic_field.loop_3d('axis',[0,0,1]);
    odeoptions = odeset;
    odeoptions.AbsTol = 1e-8;
    
    [O,I] = fumagno.propagate.x0y0_inverse; % default call
    [O,I] = fumagno.propagate.x0y0_inverse('field',field,'x',0.1,'y',0.2,'z',3); % scalar call
    [O,I] = fumagno.propagate.x0y0_inverse('field',field,'x',[0.1;0.15],'y',[0.2;0.22],'z',[0.9;0.99]); % vector call
    [O,I] = fumagno.propagate.x0y0_inverse('field',field,'x',rand(3,4),'y',rand(3,4),'z',rand(3,4)+2,'odeoptions',odeoptions); % array call
end 
 
function test_flow_solver(t)    
    % set up required fixture data
    field = magnetic_field.loop_3d('axis',[0,0,1]);
    plasma = fluid_plasma.plasma;
    odeoptions = odeset;
    odeoptions.AbsTol = 1e-8;
    [O,I] = fumagno.propagate.x0y0_inverse('field',field,'x',rand(3,4),'y',rand(3,4),'z',rand(3,4)+2,'odeoptions',odeoptions); % array call
    
    [Of,If] = fumagno.solve.flow_solver; % default call
    [Of,If] = fumagno.solve.flow_solver('plasma',plasma,'B',O.B,'B0',O.B0,'X0',O.X0,'Y0',O.Y0); % a more complete call, but using default initial conditions
    
    [O,I] = fumagno.propagate.x0y0_direct('field',field,'x0',rand(3,4)*0.2,'y0',rand(3,4)*0.2,'ds',0.1,'n_steps',13,'odeoptions',odeoptions); % full call
    [Of,If] = fumagno.solve.flow_solver('plasma',plasma,'B',O.B,'B0',O.B0,'X0',O.X0,'Y0',O.Y0); % a more complete call, but using default initial conditions
    [Of,If] = fumagno.solve.flow_solver('plasma',plasma,'B',0.1,'B0',1,'X0',1,'Y0',1); % a point that is outside the default initial condition
end 
