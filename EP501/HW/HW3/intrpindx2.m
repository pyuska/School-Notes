function [indexi,indexj] = intrpindx2(x,y,xprime,yprime)

%{
Function finds indices i,j for which

    x_i <= x' <= x_i+1    and    y_j <= y' <= y_j+1

given grids of independent variables x,y.

NOTE: 
- function assumes that the grid data increases with increasing index
- xprime and yprime are scalars
%}

A = x-xprime;
indexi = find(A>0,1)-1;

B = y-yprime;
indexj = find(A>0,1)-1;