function [index] = intrpindx1(x,xprime)

%{
Function takes a grid of points representing an independent variable x, and
a value of the independent variable (not necessarily on the grid) at which
bilinear interpolation is to take place. Function returns the index i for
which:
 
    x_i <= x' <= x_i+1

NOTE: 
- function assumes that the grid data increases with increasing index
- xprime is a scalar 
%}

% could do a bunch of super-efficient code by assuming the independent
% variable data is sorted, etc.; but a simple find() will also do the trick

A = x-xprime;
index = find(A>0,1)-1;
