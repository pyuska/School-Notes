function [a] = linlsqr(n,xdat,Ydat)

%{
Function generates coefficients for a least-squares fit of
linear-coefficient polynomial of order n. 

x and Y are measured data.  x is assumed to be exact, Y is assumed to 
include some noise. x and Y are assumed to be the same size, or at least 
conformable.

DEPENDENCIES: 
Gauss_elim.m from course repo: ~/EP501/linear_algebra/
backsub.m from same directory
%}

% N = length(Ydat);

% generate coefficient matrix
coeff = zeros(n+1); % initialize

% calculate sum i=1 > N of (x_i)^m for 0<=m<=2n outside of main for loop
% for efficiency. only have to calculate once this way, then can populate
% coeff matrix

xsum = zeros(2*n+1,1); % initialize
% xsum(1) = x^0; xsum(n) = x^(n-1)

for k = 0:2*n
    xsum(k+1) = sum(xdat.^(k));
end

b = zeros(n+1,1); % initialize solution vector

% main for loop, populate coefficient matrix and solution vector
for i = 1:n+1
    for j = 1:n+1
        coeff(i,j) = xsum(i+j-1);
    end
    b(i) = sum(Ydat.*xdat.^(i-1));
end

[Amod,~] = Gauss_elim(coeff,b);
a = backsub(Amod);