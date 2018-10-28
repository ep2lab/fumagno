%{
This is the function that finds the solution for each magnetic lilne by 
interpolation. It is important to notice, that it will be needed to do some
scaltions to adapt the problem to boundary conditions

INPUT: 

*I: Strucuture I that contains:
       - B,B0: arrays with magnetic field at each input point (B), and the
          corresponding initial condition point (B0)
       - X0,Y0: arrays with corresponding initial conditions. This and the
          previous input can be obtained one from the other with x0y0_direct or
          x0y0_inverse 
       - plasma: a fluid_plasma.plasma object, describing the different species
          in the quasineutral plasma.  
       - phi0: function handle of x0,y0 for the electric potential
       - ni0, ui0: function handles of x0,y0 for the density and velocity of
          ions  
       - ne0, ue0: cell arrays with function handles of x0,y0 for the density
          and velocity of electrons, in the same order as given in plasma.

*LIB: is a structure divided into 2 sub-structures (one for the divergent
 part and the other for the convergent part of the MN)
 Each substructure contains a library with the following variables:
  
    - ui,ue,ni,ne,phi,B: arrays with the velocity, plasma density,
    potencial drop and magnetic field.

OUTPUT:

*O: - HI,GI,UI,NI,PHI: arrays with the energy and flux integrals,
      velocity, density of ions and electric potential at the given
      points 

    - HE,GE,UE,NE: cell arrays, with arrays for each electron species in
      plasma, in the same order, for the energy and flux integrals, and
      for the velocity and density 

%----------------------------------------------------------------------
Author: Judit Nuez
Date: 20181020

%----------------------------------------------------------------------
%}

function [O] =interpolation (LIB,I)

[Bth,nth] = max(I.B(1,end,:));
[nr nl np] = size (I.B);
np = length(I.B(1,1,:));

        PHI = zeros(nr,nl,np);
        NI  = zeros(nr,nl,np);
        UI  = zeros(nr,nl,np);
        NE  = zeros(nr,nl,np);
        UE  = zeros(nr,nl,np);
        HE  = zeros(nr,nl,np);
        HI  = zeros(nr,nl,np);
        GE  = zeros(nr,nl,np);
        GI  = zeros(nr,nl,np);

        
display ('---------------------------------------------')
display ('I-FUMAGNO | Running the interpolation...')
display ('---------------------------------------------')

% Divergent MN 
if nth == 1

        for ir =  1:nr
            for il = 1:nl
            
            Bratio = I.B(ir,il,nth)/Bth;

            PHI(ir,il,:) =  interp1 (LIB.divergent.B(1:end-1)*Bratio,LIB.divergent.phi(1:end-1),I.B(ir,il,:));
            NI(ir,il,:)  =  interp1 (LIB.divergent.B(1:end-1)*Bratio,LIB.divergent.ni(1:end-1),I.B(ir,il,:));
            UI(ir,il,:)  =  interp1 (LIB.divergent.B(1:end-1)*Bratio,LIB.divergent.ui(1:end-1),I.B(ir,il,:));
            NE(ir,il,:)  =  interp1 (LIB.divergent.B(1:end-1)*Bratio,LIB.divergent.ne(1:end-1),I.B(ir,il,:));
            UE(ir,il,:)  =  interp1 (LIB.divergent.B(1:end-1)*Bratio,LIB.divergent.ue(1:end-1),I.B(ir,il,:));

            end
        end 
        
% Convergent MN  
elseif nth == length(I.B(1,1,:))
   
         for ir =  1:nr
            for il = 1:nl
                               
            Bratio = I.B(ir,il,nth)/Bth;
            
            PHI(ir,il,:) =  interp1 (LIB.convergent.B*Bratio,LIB.convergent.phi,I.B(ir,il,:));
            NI(ir,il,:)  =  interp1 (LIB.convergent.B*Bratio,LIB.convergent.ni,I.B(ir,il,:));
            UI(ir,il,:)  =  interp1 (LIB.convergent.B*Bratio,LIB.convergent.ui,I.B(ir,il,:));
            NE(ir,il,:)  =  interp1 (LIB.convergent.B*Bratio,LIB.convergent.ne,I.B(ir,il,:));
            UE(ir,il,:)  =  interp1 (LIB.convergent.B*Bratio,LIB.convergent.ue,I.B(ir,il,:));
            
            end 
         end 
% Convergent - divergent MN        
else 
       for il =  1:nl
          
            for ir = 1:nr
                
            Bratio = I.B(ir,il,nth)/Bth;
            
            PHI(ir,il,1:nth) =  interp1 (LIB.convergent.B*Bratio,LIB.convergent.phi,I.B(ir,il,1:nth));
            NI(ir,il,1:nth)  =  interp1 (LIB.convergent.B*Bratio,LIB.convergent.ni,I.B(ir,il,1:nth));
            UI(ir,il,1:nth)  =  interp1 (LIB.convergent.B*Bratio,LIB.convergent.ui,I.B(ir,il,1:nth));
            NE(ir,il,1:nth)  =  interp1 (LIB.convergent.B*Bratio,LIB.convergent.ne,I.B(ir,il,1:nth));
            UE(ir,il,1:nth)  =  interp1 (LIB.convergent.B*Bratio,LIB.convergent.ue,I.B(ir,il,1:nth));
            
            PHI(ir,il,nth:end) =  interp1 (LIB.divergent.B(1:end-1)*Bratio,LIB.divergent.phi(1:end-1),I.B(ir,il,nth:end));
            NI(ir,il,nth:end)  =  interp1 (LIB.divergent.B(1:end-1)*Bratio,LIB.divergent.ni(1:end-1),I.B(ir,il,nth:end));
            UI(ir,il,nth:end)  =  interp1 (LIB.divergent.B(1:end-1)*Bratio,LIB.divergent.ui(1:end-1),I.B(ir,il,nth:end));
            NE(ir,il,nth:end)  =  interp1 (LIB.divergent.B(1:end-1)*Bratio,LIB.divergent.ne(1:end-1),I.B(ir,il,nth:end));
            UE(ir,il,nth:end)  =  interp1 (LIB.divergent.B(1:end-1)*Bratio,LIB.divergent.ue(1:end-1),I.B(ir,il,nth:end));
            
            end
       end
end

% Save in output structure
 O.PHI = PHI; O.UI = UI; O.UE = UE; O.NI = NI; O.NE = NE;    
 
 O.GE = NE.*UE./I.B*I.B(ir,il,nth);
 O.GI = NI.*UI./I.B*I.B(ir,il,nth);
 O.HE = 0.5*UE.^2+PHI;
 O.HI = 0.5*UI.^2+PHI;
 
 
end 