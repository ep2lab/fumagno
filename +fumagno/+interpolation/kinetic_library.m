%{
This function provides an interpolation library
 
INPUT: (can be name:value pairs or a structure) 
* userdata: inputs needed by akiles2d solver
 
OUTPUT:
* solution: structure that contains the interpolation library, which is the
solution of akiles2d for a random vector of magtenic tube radius
%----------------------------------------------------------------------
Author: Judit Nuez
Date: 20181028
%----------------------------------------------------------------------
%}

function [solution] = kinetic_library (userdata, path)

[data,solution] = akiles2d.akiles2d(fullfile(path),userdata); 

end 
