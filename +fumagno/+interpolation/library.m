%{
This is the function that creates an interpolation library by solving the 
fluid equations for a random vector of plasma densities. There is a
different interpolation library for the convergent and the divergent part.

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
 
OUTPUT:
*LIB: is a structure divided into 2 sub-structures (one for the divergent
 part and the other for the convergent part of the MN)
 Each substructure contains a library with the following variables:
  
    - ui,ue,ni,ne,phi,B: arrays with the velocity, plasma density,
    potencial drop and magnetic field.

%----------------------------------------------------------------------
Author: Judit Nuez
Date: 20181020

%----------------------------------------------------------------------
%}
function [LIB] =library (I)

display ('---------------------------------------------')
display ('I-FUMAGNO | Creating interpolation library...')
display ('---------------------------------------------')

[Bth,nth] = max(I.B(1,end,:));

HI0 = fumagno.equations.Heq(I.plasma.ions,1,1,0);
HE0 = fumagno.equations.Heq(I.plasma.electrons{1},1,1,0);
GI0 = 1/Bth;
GE0 = 1/Bth;
c   = 1;
if nth == 1 %divergent MN
        N = I.Nlib(1);
        LIB.divergent.B  = zeros(N,1);
        LIB.divergent.ui = zeros(N,1);
        LIB.divergent.ue = zeros(N,1);

        LIB.divergent.ni = linspace (1,0,N);
        LIB.divergent.ne = LIB.divergent.ni;
        LIB.divergent.phi = (HE0 - I.plasma.electrons{1}.h(LIB.divergent.ni))./I.plasma.electrons{1}.q;
                
        for i=1:N-1
            a = linspace (0,1,N);
            b = downsample(a,N/10);
            
            if abs(a(i) - b(c))<1e-10
              perct = i / N*100;
              disp([num2str(round(perct,0)) '% completed'])  
                    
                      if c>0 && c<10
                        c = c+1;
                      end

            elseif i == N-1
              perct =100;
              disp([num2str(round(perct,0)) '% completed'])  
            
            end
            
             if i==1
                    ui_guess = 1;    
             else 
                    ui_guess = LIB.divergent.ui(i-1);
             end
             LIB.divergent.ui(i)   = fzero(@(x)HI0-fumagno.equations.Heq(I.plasma.ions,x,...
             LIB.divergent.ni(i),LIB.divergent.phi(i)),ui_guess);
             LIB.divergent.B(i) = LIB.divergent.ni(i)*LIB.divergent.ui(i)/GI0;   
             LIB.divergent.ue(i)   = GE0*LIB.divergent.B(i)/LIB.divergent.ne(i);
        end


        LIB.divergent.ui(end) = Inf;
        LIB.divergent.ue(end) = Inf;
        LIB.divergent.B(end)= 0;
        
        LIB.convergent.ui  = [];
        LIB.convergent.ue  = [];
        LIB.convergent.ni  = [];
        LIB.convergent.ne  = [];
        LIB.convergent.phi = [];
        LIB.convergent.B   = [];
        
      
elseif nth == length(I.B(1,1,:)) % convergent MN
        N = I.Nlib(2);
        LIB.convergent.B  = zeros(N,1);
        LIB.convergent.ui = zeros(N,1);
        LIB.convergent.ue = zeros(N,1);

        LIB.convergent.ni = linspace (1.5,1,N);
        LIB.convergent.ne = LIB.convergent.ni;
        LIB.convergent.phi = (HE0 - I.plasma.electrons{1}.h(LIB.convergent.ni))./I.plasma.electrons{1}.q;
                 
        for i=N:-1:2
            a = linspace (1,0,N);
            b = downsample(a,N/10);
            
            if abs(a(N-i+1) - b(c))<1e-10
              perct = 100-i / N*100;
              disp([num2str(round(perct,0)) '% completed'])  
                      if c>0 && c<10
                         c = c+1;
                      end 
            elseif i == 2
              perct =100;
              disp([num2str(round(perct,0)) '% completed'])  
            end    
            
            if i==N
               LIB.convergent.ui(i) = 1;   
               LIB.convergent.B(i) =  LIB.convergent.ni(i)* LIB.convergent.ui(i)/GI0;   
               LIB.convergent.ue(i)   = GE0* LIB.convergent.B(i)/ LIB.convergent.ne(i);
            elseif i==N-1
               ui_guess = 0.99;
               LIB.convergent.ui(i)   = fzero(@(x)HI0-fumagno.equations.Heq(I.plasma.ions,x, ...
               LIB.convergent.ni(i), LIB.convergent.phi(i)),ui_guess);
               LIB.convergent.B(i) =  LIB.convergent.ni(i)* LIB.convergent.ui(i)/GI0;   
               LIB.convergent.ue(i)   = GE0* LIB.convergent.B(i)/LIB.convergent.ne(i);
            else
               ui_guess = LIB.convergent.ui(i+1);    
               LIB.convergent.ui(i)   = fzero(@(x)HI0-fumagno.equations.Heq(I.plasma.ions,x,...
               LIB.convergent.ni(i), LIB.convergent.phi(i)),ui_guess);
               LIB.convergent.B(i) =  LIB.convergent.ni(i)* LIB.convergent.ui(i)/GI0;   
               LIB.convergent.ue(i)   = GE0* LIB.convergent.B(i)/ LIB.convergent.ne(i);
            end 

        end

        LIB.convergent.ue(end) = LIB.convergent.ui(end);
        LIB.convergent.B(end)= LIB.convergent.B(N-1)*0.8;
        
        LIB.divergent.ui  = [];
        LIB.divergent.ue  = [];
        LIB.divergent.ni  = [];
        LIB.divergent.ne  = [];
        LIB.divergent.phi = [];
        LIB.divergent.B   = [];
        

else % convergent-divergent MN
        
        % divergent
        N = I.Nlib(1);
        LIB.divergent.B  = zeros(N,1);
        LIB.divergent.ui = zeros(N,1);
        LIB.divergent.ue = zeros(N,1);

        LIB.divergent.ni = linspace (1,0,N);
        LIB.divergent.ne = LIB.divergent.ni;
        LIB.divergent.phi = (HE0 - I.plasma.electrons{1}.h(LIB.divergent.ni))./I.plasma.electrons{1}.q;
                
        for i=1:N-1
            a = linspace (0,1,N);
            b = downsample(a,N/10);
            
            if abs(a(i) - b(c))<1e-10
              perct = i / N*100 /2;
              disp([num2str(round(perct,0)) '% completed'])  
                  if c>0 && c<10
                     c = c+1;
                  end 
                  
             elseif i == N-1
              perct =50;
              disp([num2str(round(perct,0)) '% completed'])  
            end
             
             if i==1
                    ui_guess = 1;    
             else 
                    ui_guess = LIB.divergent.ui(i-1);
             end
             LIB.divergent.ui(i)   = fzero(@(x)HI0-fumagno.equations.Heq(I.plasma.ions,x,...
             LIB.divergent.ni(i),LIB.divergent.phi(i)),ui_guess);
             LIB.divergent.B(i) = LIB.divergent.ni(i)*LIB.divergent.ui(i)/GI0;   
             LIB.divergent.ue(i)   = GE0*LIB.divergent.B(i)/LIB.divergent.ne(i);
        end
        
        LIB.divergent.ui(end) = Inf;
        LIB.divergent.ue(end) = Inf;
        LIB.divergent.B(end)= 0;
        
        % convergent
        N = I.Nlib(2);
        LIB.convergent.B  = zeros(N,1);
        LIB.convergent.ui = zeros(N,1);
        LIB.convergent.ue = zeros(N,1);

        LIB.convergent.ni = linspace (1.5,1,N);
        LIB.convergent.ne = LIB.convergent.ni;
        LIB.convergent.phi = (HE0 - I.plasma.electrons{1}.h(LIB.convergent.ni))./I.plasma.electrons{1}.q;
    
        c = 1;            
    for i=N:-1:2
        
            a = linspace (1,0,N);
            b = downsample(a,N/10);
            
            if abs(a(N-i+1) - b(c))<1e-10
              perct = 50+(100-i/N*100)/2;
              disp([num2str(round(perct,0)) '% completed'])  
                      if c>0 && c<10
                         c = c+1;
                      end 
            elseif i == 2
              perct =100;
              disp([num2str(round(perct,0)) '% completed'])  
            end    
            
            if i==N
               LIB.convergent.ui(i) = 1;   
               LIB.convergent.B(i) =  LIB.convergent.ni(i)* LIB.convergent.ui(i)/GI0;   
               LIB.convergent.ue(i)   = GE0* LIB.convergent.B(i)/ LIB.convergent.ne(i);
            elseif i==N-1
               ui_guess = 0.9;
               LIB.convergent.ui(i)   = fzero(@(x)HI0-fumagno.equations.Heq(I.plasma.ions,x, ...
               LIB.convergent.ni(i), LIB.convergent.phi(i)),ui_guess);
               LIB.convergent.B(i) =  LIB.convergent.ni(i)* LIB.convergent.ui(i)/GI0;   
               LIB.convergent.ue(i)   = GE0* LIB.convergent.B(i)/LIB.convergent.ne(i);
               
            else
               ui_guess = LIB.convergent.ui(i+1);    
               LIB.convergent.ui(i)   = fzero(@(x)HI0-fumagno.equations.Heq(I.plasma.ions,x,...
               LIB.convergent.ni(i), LIB.convergent.phi(i)),ui_guess);
               LIB.convergent.B(i) =  LIB.convergent.ni(i)* LIB.convergent.ui(i)/GI0;   
               LIB.convergent.ue(i)   = GE0* LIB.convergent.B(i)/ LIB.convergent.ne(i);
            end 

    end

        LIB.convergent.ui(end) = 1;
        LIB.convergent.ue(end) = 1;
        LIB.convergent.B(end)=LIB.convergent.B(N-1)*0.8;
    
end 

end 