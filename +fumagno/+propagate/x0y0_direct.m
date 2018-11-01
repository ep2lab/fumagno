%{
Direct solver for the magnetic streamlines in a 3D magnetic nozzle. 
The code propagates magnetic lines from initial x0,y0 positions at z=0
to create arrays of points x,y,z that contain the geometry of the lines.
INPUT: (can be name:value pairs or a structure) 
* field: an object of a subclass of magnetic_field.element_3d, with the
  desired 3D magnetic field to use
* X0,Y0: arrays with the coordinates of the initial points
* ds: the streamwise step size. Defaults to 0.05
* n_steps: 2-dim vector of steps to calculate, including initial point on
  each streamline.  Defaults to [100 0]
        n_steps(1) = steps towards the divergent part.
        n_steps(2) = steps towards the convergent part.
  the number total of points will be n_steps(1)+n_steps(2)-1
* odeoptions: options for the integrator of magnetic streamlines, as
  those created by odeset 
OUTPUT:
* O: structure containing the following output fields:
    - X,Y,Z,BX,BY,BZ,B: position and magnetic field arrays of computed
      points. These arrays add a new last dimension to those in x0,y0.
    - X0,Y0,BX0,BY0,BZ0,B0: corresponding values of x0,y0 and the
      magnetic field at the throat plane, for each of the computed
      points
* I: structure containing all effective inputs, after adding any
  defaults to missing variables. 
%----------------------------------------------------------------------
Author: Mario Merino 
Date: 20170326

Function update
Author: Judit Nuez
Update date: 20181007
%----------------------------------------------------------------------
%}

function [O,I] = x0y0_direct(varargin)

%% Parse input
p = inputParser;
p.addParameter('field',magnetic_field.loop_3d,@(x)isa(x,'magnetic_field.element_3d'));
p.addParameter('X0',0,@isnumeric);
p.addParameter('Y0',0,@isnumeric);
p.addParameter('ds',0.05,@isnumeric);
p.addParameter('n_steps',[100 0],@isnumeric);
p.addParameter('odeoptions',odeset,@isstruct);

% Validate and parse input
p.parse(varargin{:}); % check all, and assign defaults to p.Results as needed.
 
% Place all input to the program in structure 'input' (returned as last output)
I = p.Results;

% Clear temporary variables
clear p 
 
%% Find X,Y,Z by streamline propagation  
if I.n_steps(1) == 0 % Convergent MN    
    
    %Geometry propagation
    [X2,Y2,Z2] = I.field.streamline_3d(I.X0,I.Y0,I.X0*0,-I.ds,I.n_steps(2),I.odeoptions);
    
    %Size adjustments
    sz = size(X2);
    nd = ndims(X2);

    if nd == 2 %X0 and Y0 have 1 or 2 dimensions
      
        if sz(2) == 1 %X0 and Y0 are scalar
            X2 = X2';Y2 = Y2';Z2 = Z2';  
        end 
        
        sz = size(X2);
        O.X2 = zeros(1,sz(1),I.n_steps(2));
        O.Y2 = zeros(1,sz(1),I.n_steps(2));
        O.Z2 = zeros(1,sz(1),I.n_steps(2));
        
        for i3 = 1:I.n_steps(2)                                
            O.X2(1,:,i3) = X2(:,i3);
            O.Y2(1,:,i3) = Y2(:,i3);
            O.Z2(1,:,i3) = Z2(:,i3);                 
        end
         
    elseif nd == 3 %X0 and Y0 are a 2dim array  
        O.X2 = X2; O.Y2 = Y2; O.Z2 = Z2;    
    end                          
    
    %Invert the order of the elements (convergent part)
    O.X = flip(O.X2,3); O.Y = flip(O.Y2,3); O.Z = flip(O.Z2,3);

    
elseif I.n_steps(2) == 0 %Divergent MN

    [X1,Y1,Z1] = I.field.streamline_3d(I.X0,I.Y0,I.X0.*0,I.ds,I.n_steps(1),I.odeoptions);

    %dimensions adjustment
    sz = size(X1);
    nd = ndims(X1);

    if nd == 2 %X0 and Y0 have dimension 1 or 2
    
        if sz(2) == 1 %X0 and Y0 are scalar
           X1 = X1'; Y1 = Y1'; Z1 = Z1';
        end 
        
        sz = size(X1);
        O.X1 = zeros(1,sz(1),I.n_steps(1));
        O.Y1 = zeros(1,sz(1),I.n_steps(1));
        O.Z1 = zeros(1,sz(1),I.n_steps(1));
        
        for i3 = 1:I.n_steps(1)   
            O.X1(1,:,i3) = X1(:,i3);
            O.Y1(1,:,i3) = Y1(:,i3);
            O.Z1(1,:,i3) = Z1(:,i3);                 
        end
         
    elseif nd == 3 %X0 and Y0 are a 2dim array        
        O.X1 = X1 ; O.Y1 = Y1; O.Z1 = Z1;    
    end                          
 
    O.X = O.X1; O.Y = O.Y1;O.Z = O.Z1;
    
else %Convergent-divergent MN   
    [X1,Y1,Z1] = I.field.streamline_3d(I.X0,I.Y0,I.X0.*0,I.ds,I.n_steps(1),I.odeoptions); %divergent part
    [X2,Y2,Z2] = I.field.streamline_3d(I.X0,I.Y0,I.X0.*0,-I.ds,I.n_steps(2),I.odeoptions);%convergent part

    %Dimensions adjustment
    sz = size(X1);
    nd = ndims(X1);

    if nd == 2 %X0 and Y0 have dimension 1 or 2
    
        if sz(2) == 1 
           X1 = X1';Y1 = Y1';Z1 = Z1';         
           X2 = X2';Y2 = Y2';Z2 = Z2';  
        end 
        
        sz = size(X2);
        O.X1 = zeros(1,sz(1),I.n_steps(1));
        O.Y1 = zeros(1,sz(1),I.n_steps(1));
        O.Z1 = zeros(1,sz(1),I.n_steps(1));
        O.X2 = zeros(1,sz(1),I.n_steps(1));
        O.Y2 = zeros(1,sz(1),I.n_steps(1));
        O.Z2 = zeros(1,sz(1),I.n_steps(1));
        
        for i3 = 1:I.n_steps(1) 
            O.X1(1,:,i3) = X1(:,i3);
            O.Y1(1,:,i3) = Y1(:,i3);
            O.Z1(1,:,i3) = Z1(:,i3);
        end
        
        for i3 = 1:I.n_steps(2) 
            O.X2(1,:,i3) = X2(:,i3);
            O.Y2(1,:,i3) = Y2(:,i3);
            O.Z2(1,:,i3) = Z2(:,i3);
        end
             
    elseif nd == 3  %X0 and Y0 are a 2dim array    
         O.X1 = X1; O.Y1 = Y1; O.Z1 = Z1;
         O.X2 = X2; O.Y2 = Y2; O.Z2 = Z2;
    end      
         
    %concatenate
    O.X = cat(3,flip(O.X2,3),O.X1(:,:,2:end)); 
    O.Y = cat(3,flip(O.Y2,3),O.Y1(:,:,2:end));
    O.Z = cat(3,flip(O.Z2,3),O.Z1(:,:,2:end));
end

%% Compute trivial output components

O.X0 = O.X; % Allocate
O.Y0 = O.X; 
O.BX0 = O.X;
O.BY0 = O.X;
O.BZ0 = O.X;
O.B0 = O.X;

[O.BX,O.BY,O.BZ] = I.field.field_3d(O.X,O.Y,O.Z);
O.B = sqrt(O.BX.^2+O.BY.^2+O.BZ.^2);

[Bth, nth]=max(O.B(1,1,:));

BX0 = O.BX (:,:,nth);
BY0 = O.BY (:,:,nth);
BZ0 = O.BZ (:,:,nth);
B0  = O.B(:,:,nth);

n_X0 = numel(I.X0);

for i =0:I.n_steps(1)+I.n_steps(2)-2
     O.X0(1+n_X0*i:n_X0*(i+1)) = O.X(:,:,nth); 
     O.Y0(1+n_X0*i:n_X0*(i+1)) = O.Y(:,:,nth);
     O.BX0(1+n_X0*i:n_X0*(i+1)) = BX0; 
     O.BY0(1+n_X0*i:n_X0*(i+1)) = BY0;
     O.BZ0(1+n_X0*i:n_X0*(i+1)) = BZ0; 
     O.B0(1+n_X0*i:n_X0*(i+1)) = B0;
end 


end 
