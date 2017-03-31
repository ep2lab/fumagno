%{
Direct solver for the magnetic streamlines in a 3D magnetic nozzle. 
The code propagates magnetic lines from initial x0,y0 positions at z=0
to create arrays of points x,y,z that contain the geometry of the lines.

INPUT: (can be name:value pairs or a structure) 
* field: an object of a subclass of magnetic_field.element_3d, with the
  desired 3D magnetic field to use
* X0,Y0: arrays with the coordinates of the initial points
* ds: the streamwise step size. Defaults to 0.05
* n_steps: the number of steps to calculate, including initial point on
  each streamline.  Defaults to 100
* odeoptions: options for the integrator of magnetic streamlines, as
  those created by odeset 

OUTPUT:
* O: structure containing the following output fields:
    - X,Y,Z,BX,BY,BZ,B: position and magnetic field arrays of computed
      points. These arrays add a new last dimension to those in x0,y0.
    - X0,Y0,BX0,BY0,BZ0,B0: corresponding values of x0,y0 and the
      magnetic field at the initial plane, for each of the computed
      points
* I: structure containing all effective inputs, after adding any
  defaults to missing variables. 

%----------------------------------------------------------------------
Author: Mario Merino
Date: 20170326
%----------------------------------------------------------------------
%}
function [O,I] = x0y0_direct(varargin)

%% Parse input
p = inputParser;
p.addParameter('field',magnetic_field.loop_3d,@(x)isa(x,'magnetic_field.element_3d'));
p.addParameter('X0',0,@isnumeric);
p.addParameter('Y0',0,@isnumeric);
p.addParameter('ds',0.05,@isnumeric);
p.addParameter('n_steps',100,@isnumeric);
p.addParameter('odeoptions',odeset,@isstruct);


% Validate and parse input
p.parse(varargin{:}); % check all, and assign defaults to p.Results as needed.
 
% Place all input to the program in structure 'input' (returned as last output)
I = p.Results;

% Clear temporary variables
clear p 
 
%% Find X,Y,Z by streamline propagation  
[O.X,O.Y,O.Z] = I.field.streamline_3d(I.X0,I.Y0,I.X0.*0,I.ds,I.n_steps,I.odeoptions);
 
%% Compute trivial output components
O.X0 = O.X; % Allocate
O.Y0 = O.X; 
O.BX0 = O.X;
O.BY0 = O.X;
O.BZ0 = O.X;
O.B0 = O.X;
[BX0,BY0,BZ0] = I.field.field_3d(I.X0,I.Y0,I.X0.*0); % temporary values
B0 = sqrt(BX0.^2+BY0.^2+BZ0.^2);
n_X0 = numel(I.X0);
for i =0:I.n_steps-1
    O.X0(1+n_X0*i:n_X0*(i+1)) = I.X0; 
    O.Y0(1+n_X0*i:n_X0*(i+1)) = I.Y0;
    O.BX0(1+n_X0*i:n_X0*(i+1)) = BX0; 
    O.BY0(1+n_X0*i:n_X0*(i+1)) = BY0;
    O.BZ0(1+n_X0*i:n_X0*(i+1)) = BZ0; 
    O.B0(1+n_X0*i:n_X0*(i+1)) = B0;
end 

%% Compute BX, BY, BZ, B at X,Y,Z
[O.BX,O.BY,O.BZ] = I.field.field_3d(O.X,O.Y,O.Z);
O.B = sqrt(O.BX.^2+O.BY.^2+O.BZ.^2);
  
