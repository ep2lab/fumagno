%{
This is the solver for the fully magnetized quasineutral plasma flow in
a magnetic nozzle. The functions x0y0_direct and x0y0_inverse can be
used first to generate some of the inputs.
INPUT: (can be name:value pairs or a structure) 
* B,B0: arrays with magnetic field at each input point (B), and the
  corresponding initial condition point (B0)
* X0,Y0: arrays with corresponding initial conditions. This and the
  previous input can be obtained one from the other with x0y0_direct or
  x0y0_inverse 
* plasma: a fluid_plasma.plasma object, describing the different species
  in the quasineutral plasma.  
* phi0: function handle of x0,y0 for the electric potential
* ni0, ui0: function handles of x0,y0 for the density and velocity of
  ions  
* ne0, ue0: cell arrays with function handles of x0,y0 for the density
  and velocity of electrons, in the same order as given in plasma.
 
OUTPUT:
* Of: structure containing the following output fields:
    - HI,GI,UI,NI,PHI: arrays with the energy and flux integrals,
      velocity, density of ions and electric potential at the given
      points 
    - HE,GE,UE,NE: cell arrays, with arrays for each electron species in
      plasma, in the same order, for the energy and flux integrals, and
      for the velocity and density 
* I: structure containing all effective inputs, after adding any
  defaults to absent variables. 
%----------------------------------------------------------------------
Author: Mario Merino
Date: 20170326

-----Update------
Judit Nuez
Date: 20181020
%----------------------------------------------------------------------
%}
function [Of,I] = flow_solver(varargin)
 
%% Parse input
p = inputParser; 

% Plasma object
p.addParameter('plasma',fluid_plasma.plasma,@(x)isa(x,'fluid_plasma.plasma'));

% Find out how many species does plasma have
p.KeepUnmatched = true;
p.parse(varargin{:});
n_electrons = p.Results.plasma.n_electrons;
p.KeepUnmatched = false;

% Continue the validation     
p.addParameter('B',1,@isnumeric);

p.addParameter('X0',0,@isnumeric);
p.addParameter('Y0',0,@isnumeric);
p.addParameter('B0',1,@isnumeric);

p.addParameter('phi0',0,@isnumeric);
p.addParameter('ni0',@default_n,@(x)isa(x,'function_handle'));
p.addParameter('ui0',@default_ui,@(x)isa(x,'function_handle'));
    default_ncell(1:n_electrons) = {@default_n};
p.addParameter('ne0',default_ncell,@iscell);
    default_ucell(1:n_electrons) = {@default_ui};
p.addParameter('ue0',default_ucell,@iscell);
p.addParameter('Nlib',[10000 10000],@isnumeric);

p.parse(varargin{:}); % check all, and assign defaults to p.Results as needed.
 
% Place all input to the program in structure 'input' (to be given as last output)
I = p.Results;

% Clear unneeded variables
clear p default_ncell default_ucell

% Interpolation Library
[LIB] = fumagno.interpolation.library (I);

% Interpolation
[OI]   = fumagno.interpolation.interpolation (LIB,I); 

% Postprocessing
        Of.PHI = OI.PHI-OI.PHI(1,:,1);
        n = default_n(I.X0, I.Y0);
        Of.NI  = OI.NI.*n;
        Of.UI  = OI.UI.*n;
        Of.NE  = Of.NI; %QN
        Of.UE  = OI.UE.*n;
        Of.HE  = OI.HE; Of.GE  = OI.GE; 
        Of.HI  = OI.HI; Of.GI  = OI.GI;
       
        
        display ('I-FUMAGNO Completed')
end



%----------------------------------------------------------------------
%----------------------------------------------------------------------
%----------------------------------------------------------------------
%  
function n = default_n(x0,y0)
%Default plasma density in case no input is provided. Use a Gaussian
% profile with a peak density of 1 in the circle with R = 1
    n = exp(-(x0.^2+y0.^2)*3*log(10));
    n(x0.^2+y0.^2>1)=0;
end

function u = default_u(x0,y0)
%Default velocity in case no input is provided. Use u = 1 in the circle
%with R=1

    u = x0.*0 + 1;
    u(x0.^2+y0.^2>1) = 0;
end
