%{
This is the function that interpolates the solution for each magnetic line.
The interpolation library has been previously computed using akiles2d
 
INPUT: (can be name:value pairs or a structure) 
* solution: structure that contains the interpolation library
* I: structure containing all effective inputs, after adding any
  defaults to absent variables. 

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
* Fum_Kin_solution: structure with the solution of the kinetic poblem in
the FMIL
%----------------------------------------------------------------------
Author: Judit Nuez
Date: 20181028

%----------------------------------------------------------------------
%}

function  [Fum_Kin_solution]  = kinetic_interp (solution, I)

Tad = (solution.electrons.Tz(1)+solution.electrons.Tr(1)+solution.electrons.Ttheta(1))./3; 
h_F = sqrt(I.B0./I.B);
      
   % define the indexes
[r, c1, c2] = size(h_F);
ll = {'','1','2','4'}; % Label for the whole population and each subpopulation
ee = {'electrons','ions'};
   
   % Allocate
n       = zeros(r,c1,c2);
u       = zeros(r,c1,c2);
Tz      = zeros(r,c1,c2);
Tr      = zeros(r,c1,c2);
Ttheta  = zeros(r,c1,c2);
qzz     = zeros(r,c1,c2);
qzr     = zeros(r,c1,c2);
qztheta = zeros(r,c1,c2);  
phi     = zeros(r,c1,c2);
  
 
 

 for ie = 1:length(ee) 
     
     for il = 1:length(ll)
       
            l = ll{il};
            e = ee{ie};
         
                                
                h_K = solution.h'; % random h used to solve akiles
                h_K = h_K(1:end-1);
          
                phi_K = solution.phi';
                phi_K = phi_K(1:end-1);
                phi(:,:,:) = interp1(h_K,phi_K,h_F(:,:,:))/Tad;
             
                % Plasma density
                n_K =solution.(e).(['n',l])';
                n_K = n_K(1:end-1);
                n(:,:,:) = interp1(h_K,n_K,h_F(:,:,:)); 
                
                u_K =solution.(e).(['u',l])';
                u_K = u_K(1:end-1);
                u(:,:,:) = interp1(h_K,u_K,h_F(:,:,:)); 
                
                
for ir = 1:r
  for ic1 = 1:c1
    for ic2 = 1:c2    
        
                % Axial temperature
                Tz_K =solution.(e).(['Tz',l])';
                Tz_K = Tz_K(1:end-1);
                Tz(ir,ic1,ic2) = interp1(h_K,Tz_K,h_F(ir,ic1,ic2))/Tad; 

                % Radial temperature
                Tr_K =solution.(e).(['Tr',l])';
                Tr_K = Tr_K(1:end-1);
                Tr(ir,ic1,ic2) = interp1(h_K,Tr_K,h_F(ir,ic1,ic2))/Tad;    

                % Azimutal temperature
                Ttheta_K =solution.(e).(['Ttheta',l])';
                Ttheta_K = Ttheta_K(1:end-1);
                Ttheta(ir,ic1,ic2) = interp1(h_K,Ttheta_K,h_F(ir,ic1,ic2))/Tad; 
    end          
  end
end


                % Axial flux of axial heat
                qzz_K =solution.(e).(['qzz',l])';
                qzz_K = qzz_K(1:end-1);
                qzz(:,:,:) = interp1(h_K,qzz_K,h_F(:,:,:))*(Tad^(-3/2)); 

                % Axial flux of radial heat
                qzr_K =solution.(e).(['qzr',l])';
                qzr_K = qzr_K(1:end-1);
                qzr(:,:,:) = interp1(h_K,qzr_K,h_F(:,:,:))*(Tad^(-3/2)); 

                % Axial flux of azimutal heat
                qztheta_K =solution.(e).(['qztheta',l])';
                qztheta_K = qztheta_K(1:end-1);
                qztheta(:,:,:) = interp1(h_K,qztheta_K,h_F(:,:,:))*(Tad^(-3/2)); 
     
    

    %Save the solution in a structure         
    x = I.phi0-phi;
    Fum_Kin_solution.h        = h_F; 
    Fum_Kin_solution.ne00p    = solution.ne00p;
    Fum_Kin_solution.npoints  = solution.npoints;
    Fum_Kin_solution.phi      = phi-x;
    Fum_Kin_solution.r        = solution.r;
    Fum_Kin_solution.errorfcn = solution.errorfcn;
    
    Fum_Kin_solution.(e).(['n',l])       = n.* I.ni0(I.X0,I.Y0);
    Fum_Kin_solution.(e).(['u',l])       = u.* I.ui0(I.X0,I.Y0);
    Fum_Kin_solution.(e).(['Tz',l])      = Tz;
    Fum_Kin_solution.(e).(['Tr',l])      = Tr;
    Fum_Kin_solution.(e).(['Ttheta',l])  = Ttheta;
    Fum_Kin_solution.(e).(['qzz',l])     = qzz;
    Fum_Kin_solution.(e).(['qzr',l])     = qzr;
    Fum_Kin_solution.(e).(['qztheta',l]) = qztheta;
    
     end 
 end  
  
 end 