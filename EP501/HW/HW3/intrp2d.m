function [finterp] = intrp2d(x,y,xprime,yprime,f)

%{
Find list of indices i,j in independent variable grids x,y such that

    x_i,n <= x_n' <= x_i+1,n    and    y_j,n <= y_n' <= y_j+1,n

for lists of points to interpolate at xprime,yprime.

Then, interpolate values of f at (xprime,yprime) using these indices.
Function returns the interpolated values of f that correspond to the xprime
and yprime grids.

NOTE:
- function assumes that the grid data increases with increasing index
- xprime and yprime are vectors
- dimensions of f are assumed to be lengths of x/y prime vectors

DEPENDENCIES: 
Gauss_elim.m from course repo: ~/EP501/linear_algebra/
backsub.m from same directory
%}

% find indices corresponding to points to interpolate at
indexlisti = zeros(length(xprime),1);
indexlistj = zeros(length(yprime),1);

% populate index vectors

% for every value in the xprime/yprime vectors, subtract the 
for a = 1:length(xprime)
    A = x-xprime(a);
    if isempty(find(A==0,1))
        indexlisti(a) = find(A>0,1)-1;
    else
        indexlisti(a) = find(A==0);
    end % if - isempty
end

for b = 1:length(yprime)
    B = y-yprime(b);
    if isempty(find(B==0,1))
        indexlistj(b) = find(B>0,1)-1;
    else
        indexlistj(b) = find(B==0);
    end % if - isempty
end

%{
Linear interpolation of a bivariate function f(x,y) is given by

    z = f(x,y) = a + bx + cy + dxy

substituting values for f,x, and y for the 4 points surrounding the point
to interpolate at allows for a,b,c, and d to be solved for using
Gauss elimination.
%}

% if xprime(n) == 1 || xprime(n) == length(xprime) (and same for yprime),
% then the point in question is on the edge of the domain, in which case 1D
% bilinear interpolation can be used: intrpindx1.m + some extra computation

% initialize grid of interplated points
finterp = zeros(length(xprime),length(yprime));

% modify f to add row/col of zeros at edge of domain to hotfix out of
% bounds error
fwork = cat(1,f,zeros(1,size(f,2)));
fwork = cat(2,fwork,zeros(size(fwork,1),1));

for m = 1:length(xprime)
    % define x_i, x_i+1
    x1 = x(indexlisti(m));
    
    % handle interpolation at domain boundaries
    if indexlisti(m)+1 > length(x)
        x2 = 0;
    else
        x2 = x(indexlisti(m)+1);
    end
    
    for n = 1:length(yprime)
        % define y_j, y_j+1
        y1 = y(indexlistj(n));
        
        % handle interpolation at domain boundaries
        if indexlistj(n)+1 > length(y)
            y2 = 0;
        else
            y2 = y(indexlistj(n)+1);
        end
        
        % i feel like there should be a cleaner way to implement this with
        % loops or SOMETHING but it's just not coming to me right now, so i'm
        % gonna hardcode it.
        
        % coefficient matrix fro Gauss_elim fn
        gaussMat = ...
            [1 x1 y1 x1*y1;...
            1 x1 y2 x1*y2;...
            1 x2 y1 x2*y1;...
            1 x2 y2 x2*y2];
        
        % b vector for Gauss_elim fn
        b = [fwork(indexlisti(m),indexlistj(n));...
             fwork(indexlisti(m),indexlistj(n)+1);...
             fwork(indexlisti(m)+1,indexlistj(n));...
             fwork(indexlisti(m)+1,indexlistj(n)+1)];
        [Amod,ord] = Gauss_elim(gaussMat,b);
        
        % solution vector
        abcd = backsub(Amod(ord,:));
        finterp(m,n) = abcd(1) + abcd(2)*xprime(m) + abcd(3)*yprime(n) + ...
            abcd(4)*xprime(m)*yprime(n);
    end % for - n
end % for - m
        
        
end % function

