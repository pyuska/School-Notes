function [chisqr] = rcss(n,Ydat,sigy,poly)

%{
Function calculates the reduced chi-squared statistic (a goodness-of-fit
statistic) for a polynomial of order n. The polynomial should already be 
evaluated at every point where Ydat is defined.

NOTE: sigy, Ydat, and poly must all be either 
%}

% degrees of freedom of fit
nu = length(Ydat)-(n+1);

% produces garbage answer if any of 3 array inputs do not have identical
% dimensions to the others, so warn user and exit
if sum(size(sigy)==size(poly)) ~= 2
    error('Error: Array inputs to rcss.m are not equal in size or have mismatched dimensions. Transpose inputs as necessary to ensure identical dimensions.')
end

% otherwise compute desired value
chisqr = 1/nu*sum(1./(sigy.^2).*(Ydat-poly).^2);