%{
Computes H of a given species on given points

INPUT: 
* species: a fluid_plasma.species object with the mass, charge,
  and h function of the current plasma species
* u,n,phi: arrays with the velocity, density and electric potential at
  the points of calculation

OUTPUT:
* H: the energy integral of the species at the points of calculation

%----------------------------------------------------------------------
Author: Mario Merino
Date: 20170313
%----------------------------------------------------------------------
%}
function H = Heq(species,u,n,phi)

%% Compute H
H = species.m*u.^2/2 + species.h(n) + species.q*phi;  
