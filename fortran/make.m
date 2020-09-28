clear 
close all
clc

mex -v -O -largeArrayDims -output computeGREEN_f90_mexed integration_fortran.f90 
