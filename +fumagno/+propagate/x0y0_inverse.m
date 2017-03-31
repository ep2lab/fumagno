%{
Inverse solver for the magnetic streamlines in a 3D magnetic nozzle. 
The code finds the coordinates x0,y0 at z=0 of the magnetic lines that
pass by the given points x,y,z (arrays)

INPUT: (can be name:value pairs or a structure) 
* field: an object of a subclass of magnetic_field.element_3d, with the
  desired 3D magnetic field to use
* X,Y,Z: arrays with the coordinates of the points where the properties
  of the plasma are queried
* direction: direction to propagate magnetic streamlines until they
  intersect with initial plane. -1 backward (default), 1 forward
* maxlength: maximum distance for propagation of B lines before dropping
  the integration (will raise error if z=0 is not found)
* odeoptions: options for the integrator of magnetic streamlines, as
  those created by odeset 

OUTPUT:
* O: structure containing the following output fields:
    - X0,Y0,BX0,BY0,BZ0,B0: arrays for the corresponding x0,y0
      coordinates for each magnetic streamline at the initial plane, and
      the corresponding magnetic field there
    - X,Y,Z,BX,BY,BZ,B: position and magnetic field arrays
* I: structure containing all effective inputs, after adding any
  defaults to missing variables. 

%----------------------------------------------------------------------
Author: Mario Merino
Date: 20170326
%----------------------------------------------------------------------
%}
function [O,I] = x0y0_inverse(varargin)
 
%% Parse input
p = inputParser;
p.addParameter('field',magnetic_field.loop_3d,@(x)isa(x,'magnetic_field.element_3d'));
p.addParameter('X',0,@isnumeric);
p.addParameter('Y',0,@isnumeric);
p.addParameter('Z',1,@isnumeric);
p.addParameter('direction',-1,@isnumeric);
p.addParameter('maxlength',1e6,@isnumeric);
p.addParameter('odeoptions',odeset,@isstruct);

% Validate and parse input
p.parse(varargin{:}); % check all, and assign defaults to p.Results as needed.
 
% Place all input to the program in structure 'input' (returned as last output)
I = p.Results;

% Clear temporary variables
clear p  
 
%% Compute trivial output components
O.X = I.X;
O.Y = I.Y;
O.Z = I.Z;
[O.BX,O.BY,O.BZ] = I.field.field_3d(O.X,O.Y,O.Z);
O.B = sqrt(O.BX.^2+O.BY.^2+O.BZ.^2);

%% Find X0, Y0 by streamline propagation and intersection with z = 0
I.odeoptions.Events = @odeevents; % Add the events function to odepotions
I.odeoptions.MaxStep = 1; 
O.X0 = I.X.*0; % Allocate
O.Y0 = I.X.*0;
for i = 1:numel(I.X)
    if O.Z(i) == 0 % the point given is already at the initial plane
        O.X0(i) = I.X(i);
        O.Y0(i) = I.Y(i);
    else
        sol = ode45(@(t,PositionVector)propagateB(I.field,I.direction,PositionVector),[0,I.maxlength],[O.X(i);O.Y(i);O.Z(i)],I.odeoptions);
        O.X0(i) = sol.ye(1);
        O.Y0(i) = sol.ye(2);
    end
end    

%% Compute BX0, BY0, BZ0, B0 at X0, Y0
[O.BX0,O.BY0,O.BZ0] = I.field.field_3d(O.X0,O.Y0,O.X0.*0);
O.B0 = sqrt(O.BX0.^2+O.BY0.^2+O.BZ0.^2);
 
end % end main function

%----------------------------------------------------------------------
%----------------------------------------------------------------------
%----------------------------------------------------------------------

function DerivativeVector = propagateB(field,direction,PositionVector)
% This function is used by ode45 to propagate the magnetic lines
% backward (or forward, according to direction) until they intersect
% with the initial plane. The magnetic_field.element_3d streamline_3d
% function is not used here as we require being able to invert field
% direction
    x = PositionVector(1);
    y = PositionVector(2);
    z = PositionVector(3);
    [Bx,By,Bz] = field.field_3d(x,y,z);
    DerivativeVector = direction*[Bx;By;Bz]/sqrt(Bx^2 + By^2 + Bz^2);
end
function [value,isterminal,direction] = odeevents(~,PositionVector)
% Event function used by ode45 to stop the integration when z = 0
    value = PositionVector(3);
    isterminal = 1;
    direction = 0;    
end 
