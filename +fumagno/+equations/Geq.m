%{
Computes G of a given species on given points

INPUT: 
* u,n,B: arrays with the velocity, density and magnetic field at the
  points of calculation

OUTPUT:
* G: the flux integral of the species at the points of calculation

%----------------------------------------------------------------------
Author: Mario Merino
Date: 20170313
%----------------------------------------------------------------------
%}
function G = Geq(u,n,B)

%% Compute H
G = u.*n./B;
