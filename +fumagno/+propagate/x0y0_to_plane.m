%{
Direct solver for the magnetic streamlines in a 3D magnetic nozzle where
the computed points are at fixez z = const planes
The code propagates magnetic lines from initial x0,y0 positions at z = 0
to create arrays of points x,y,z that contain the geometry of the lines.

Warning: the code assumes that streamlines do intersect with the
selected planes without checking.

INPUT: (can be name:value pairs or a structure) 
* field: an object of a subclass of magnetic_field.element_3d, with the
  desired 3D magnetic field to use
* X0,Y0: arrays with the coordinates of the initial points
* z_planes: vector of values of z where the streamlines will be
  intersected to create new points
* direction: direction to propagate magnetic streamlines until they
  intersect with desired planes. -1 backward, 1 forward (default)
* maxlength: maximum distance for propagation of B lines before dropping
  the integration (will raise error if z=0 is not found)
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
function [O,I] = x0y0_to_plane(varargin)

%% Parse input
p = inputParser;
p.addParameter('field',magnetic_field.loop_3d,@(x)isa(x,'magnetic_field.element_3d'));
p.addParameter('X0',0,@isnumeric);
p.addParameter('Y0',0,@isnumeric);
p.addParameter('z_planes',[0;1;2],@isnumeric);
p.addParameter('direction',1,@isnumeric);
p.addParameter('maxlength',1e6,@isnumeric);
p.addParameter('odeoptions',odeset,@isstruct);

% Validate and parse input
p.parse(varargin{:}); % check all, and assign defaults to p.Results as needed.
 
% Place all input to the program in structure 'input' (returned as last output)
I = p.Results;
I.z_planes = I.z_planes(:); % force column

% Clear temporary variables
clear p 
 
%% Propagate streamlines until intersection with z planes

% Allocate output
n_X0 = numel(I.X0);
n_z_planes = length(I.z_planes);
input_size = size(I.X0); % find non-singleton size of input
for i = length(input_size):-1:1
    if input_size(i) == 1
        input_size = input_size(1:i-1);
    else
        break
    end
end
if isempty(input_size) % allocate
    O.X = zeros([n_z_planes,1]);
else
    O.X = zeros(cat(2,input_size,n_z_planes));
end 
O.Y = O.X;
O.Z = O.X;
I.odeoptions.MaxStep = 1; 
for i = 1:n_z_planes
    I.odeoptions.Events = @(t,PositionVector)odeevents(t,PositionVector,I.z_planes(i)); % Add the events function to odepotions
    O.Z(1+n_X0*(i-1):n_X0*i) = I.z_planes(i);
    if I.z_planes(i) == 0 % the plane coincides with the initial plane
        O.X(1+n_X0*(i-1):n_X0*i) = I.X0;
        O.Y(1+n_X0*(i-1):n_X0*i) = I.Y0;
    else
        for j = 1:n_X0
            sol = ode45(@(t,PositionVector)propagateB(I.field,I.direction,PositionVector),[0,I.maxlength],[I.X0(j);I.Y0(j);I.X0(j).*0],I.odeoptions);
            O.X(j+n_X0*(i-1)) = sol.ye(1);
            O.Y(j+n_X0*(i-1)) = sol.ye(2);
        end
    end
end    

%% Compute BX, BY, BZ, B at X,Y,Z
[O.BX,O.BY,O.BZ] = I.field.field_3d(O.X,O.Y,O.Z);
O.B = sqrt(O.BX.^2+O.BY.^2+O.BZ.^2);

%% Compute trivial outputs
O.X0 = O.X; % allocate
O.Y0 = O.X;
O.BX0 = O.X;
O.BY0 = O.X;
O.BZ0 = O.X;
O.B0 = O.X;
[BX0,BY0,BZ0] = I.field.field_3d(I.X0,I.Y0,I.X0.*0); % temporary values
B0 = sqrt(BX0.^2+BY0.^2+BZ0.^2);
for i =0:n_z_planes-1
    O.X0(1+n_X0*i:n_X0*(i+1)) = I.X0; 
    O.Y0(1+n_X0*i:n_X0*(i+1)) = I.Y0;
    O.BX0(1+n_X0*i:n_X0*(i+1)) = BX0; 
    O.BY0(1+n_X0*i:n_X0*(i+1)) = BY0;
    O.BZ0(1+n_X0*i:n_X0*(i+1)) = BZ0; 
    O.B0(1+n_X0*i:n_X0*(i+1)) = B0;
end 
 
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
function [value,isterminal,direction] = odeevents(~,PositionVector,z_plane)
% Event function used by ode45 to stop the integration when z = 0
    value = PositionVector(3)-z_plane;
    isterminal = 1;
    direction = 0;    
end 
