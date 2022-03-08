clear 
close all
clc

mex -v -O -largeArrayDims -output computeGREEN_f90_mexed integration_fortran.f90 
%mex -v -O -largeArrayDims -output computeGREEN_f90_mexed integration_fortran.F90 % IN CASE THE ABOVE does not work, change the extension from .f90 to .F90
