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
* O: structure containing the following output fields:
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
%----------------------------------------------------------------------
%}
function [O,I] = flow_solver(varargin)
 
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

p.addParameter('phi0',@(x0,y0)0,@(x)isa(x,'function_handle'));
p.addParameter('ni0',@default_n,@iscell);
p.addParameter('ui0',@default_u,@iscell);
    default_ncell(1:n_electrons) = {@default_n};
p.addParameter('ne0',default_ncell,@iscell);
    default_ucell(1:n_electrons) = {@default_u};
p.addParameter('ue0',default_ucell,@iscell);

p.parse(varargin{:}); % check all, and assign defaults to p.Results as needed.
 
% Place all input to the program in structure 'input' (to be given as last output)
I = p.Results;

% Clear unneeded variables
clear p default_ncell default_ucell
 
%% Compute H, G conserved quantities at each point
O.HI = fumagno.equations.Heq(I.plasma.ions,I.ui0(I.X0,I.Y0),I.ni0(I.X0,I.Y0),I.phi0(I.X0,I.Y0));
O.GI = fumagno.equations.Geq(I.ui0(I.X0,I.Y0),I.ni0(I.X0,I.Y0),I.B0);
for i_electrons = 1:n_electrons
    O.HE{i_electrons} = fumagno.equations.Heq(I.plasma.electrons{i_electrons},I.ue0{i_electrons}(I.X0,I.Y0),I.ne0{i_electrons}(I.X0,I.Y0),I.phi0(I.X0,I.Y0));
    O.GE{i_electrons} = fumagno.equations.Geq(I.ue0{i_electrons}(I.X0,I.Y0),I.ne0{i_electrons}(I.X0,I.Y0),I.B0);
end
 
%% Prepare interpolation libraries for guess (adequate guess for near-sonic cold ions + isothermal electrons)
q = I.plasma.ions.q; 
T0 = 0; % lumped temperature of all electron species
for i_electrons = 1:n_electrons 
    T0 = T0 + I.plasma.electrons{i_electrons}.T(I.ne0{i_electrons}(0,0)); 
end 
T0 = max(1,T0);
M0 = max(1,I.ui0(0,0))/I.plasma.cs(max(1,I.ni0(0,0))); % estimate of Mach number
uratio = linspace(1,20,10000); % velocity ratio interpolation library. !!! May need to change the limits
Bratio = uratio.*exp(-0.5*M0^2*(uratio.^2-1)); % magnetic field ratio interpolation library 
phidiff = T0*log(Bratio./uratio)/q; % potential difference interpolation library

%% Solve for U, N, PHI
O.PHI = I.B.*0; % Allocate
O.NI = I.B.*0;    
O.UI = I.B.*0;
for i_electrons = 1:n_electrons
    O.NE{i_electrons} = I.B.*0;    
    O.UE{i_electrons} = I.B.*0;
end
for j = 1:numel(I.B) % for each point
    if O.HI(j) == 0 % current streamline is empty, we are outside of MN
        O.NI(j) = NaN;
        O.UI(j) = NaN;
        O.PHI(j) = NaN;
        for i_electrons = 1:n_electrons
            O.NE{i_electrons}(j) = NaN;
            O.UE{i_electrons}(j) = NaN;
        end
        continue;
    end
    if I.B(j) == I.B0(j) % We are at the initial plane so return initial values
        O.NI(j) = I.ui0(I.X0(j),I.Y0(j));
        O.UI(j) = I.ui0(I.X0(j),I.Y0(j));
        O.PHI(j) = I.ui0(I.X0(j),I.Y0(j));
        for i_electrons = 1:n_electrons
            O.NE{i_electrons}(j) = I.ne0{i_electrons}(I.X0(j),I.Y0(j));
            O.UE{i_electrons}(j) = I.ue0{i_electrons}(I.X0(j),I.Y0(j));
        end
        continue;
    end        
    % Prepare a sensible initial guess for the supersonic branch
    ui = I.ui0(I.X0(j),I.Y0(j))*interp1(Bratio,uratio,I.B(j)/I.B0(j));
    phi = I.phi0(I.X0(j),I.Y0(j)) + interp1(Bratio,phidiff,I.B(j)/I.B0(j)); 
    ue = cell(1,n_electrons);
    for i_electrons = 1:n_electrons
        ue{i_electrons} = I.ue0{i_electrons}(I.X0(j),I.Y0(j))*interp1(Bratio,uratio,I.B(j)/I.B0(j));        
    end         
    % Put H and G in a form that can be used to call equationsystem
    for i_electrons = 1:n_electrons
        Hetemp{i_electrons} = O.HE{i_electrons}(j);
        Getemp{i_electrons} = O.GE{i_electrons}(j);
    end
    % Solve quasineutrality equation to determine phi, then the rest
    O.PHI(j) = fzero(@(phi)quasineutrality(phi,ui,ue,I.plasma,I.B(j),O.HI(j),O.GI(j),Hetemp,Getemp),phi); % fzero is shamefully slow and unflexible. ??? consider doing my own solver
    O.UI(j) = fzero(@(x)O.HI(j)-fumagno.equations.Heq(I.plasma.ions,x,O.GI(j)*I.B(j)/x,O.PHI(j)),ui);
    O.NI(j) = O.GI(j)*I.B(j)/O.UI(j);
    for i_electrons = 1:n_electrons
        O.UE{i_electrons}(j) = fzero(@(x)Hetemp{i_electrons} - ...
            fumagno.equations.Heq(I.plasma.electrons{i_electrons},x,Getemp{i_electrons}*I.B(j)/x,O.PHI(j)),ue{i_electrons});
        O.NE{i_electrons}(j) = Getemp{i_electrons}*I.B(j)/O.UE{i_electrons}(j);
    end 
end 

end % end main function

%----------------------------------------------------------------------
%----------------------------------------------------------------------
%----------------------------------------------------------------------
 
function n = default_n(x0,y0)
% Default plasma density in case no input is provided. Use a Gaussian
% profile with a peak density of 1 in the circle with R = 1
    n = exp(-(x0.^2+y0.^2)*3*log(10));
    n(x0.^2+y0.^2>1)=0;
end
function u = default_u(x0,y0)
% Default plasma density in case no input is provided. Use 
% u = 1 in the circle with R = 1
    u = x0.*0;
    u(x0.^2+y0.^2<=1) = 1;
end 
function err = quasineutrality(phi,ui,ue,plasma,B,Hi,Gi,He,Ge)
% Error equation used to solve for phi
    % ions
    ui = fzero(@(x)Hi-fumagno.equations.Heq(plasma.ions,x,Gi*B/x,phi),ui);
    ni = Gi*B/ui;    
    err = plasma.ions.q*ni; % quasineutrality error
    % electrons
    for i_electrons = 1:plasma.n_electrons
        % Get ne
        ue{i_electrons} = fzero(@(x)He{i_electrons}-fumagno.equations.Heq(plasma.electrons{i_electrons},x,Ge{i_electrons}*B/x,phi),ue{i_electrons});
        ne = Ge{i_electrons}*B/ue{i_electrons};
        % Add to quasineutrality error    
        err = err + plasma.electrons{i_electrons}.q*ne;
    end
end
