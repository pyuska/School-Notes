function [root,it,success]=newton_exact_polyvec(P,Pprime,x0,maxit,tol,verbose)

% root=newton_exact(f,fprime)
% 
% Modified from newton_exact.m code from course repo. Finds a set of roots 
% corresponding to the polynomial P (input as a vector of coefficients in 
% descending order). This input format is much easier to modify for
% repeated root-finding iterations after polynomial deflation.

% Dependencies: polyeval.m

%% Error checking of input
narginchk(3,6);   %check for correct number of inputs to function
if (nargin<4)
    maxit=100;       %maximum number of iterations allowed
end %if
if (nargin<5)
    tol=1e-6;        %how close to zero we need to get to cease iterations
end %if
if (nargin<6)
    verbose=false;
end %if


%% Make sure we don't start at an inflection point with zero derivative
if (abs(polyeval(Pprime,x0))<tol)
    warning(' Attempting to start Newton iterations near an inflection point, you may wish to restart with a different guess...');
    x0=x0+1;   %bump the guess a ways off of initial value to see if we can get anything sensible
end %if


%% Newton iterations
it=1;
root=x0;
fval=polyeval(P,root);
converged=false;
while(~converged && it<=maxit)
    derivative=polyeval(Pprime,root);
    if (abs(derivative)<100*tol)    %this (inflection point) will end up kicking the root really far away...
        converged=false;
        warning(' Derivative close to zero, terminating iterations with failed convergence... ');
        break;
    else
        root=root-fval/derivative;    % update root estimate
        fval=polyeval(P,root);                  % see how far off we are from zero...
        if (verbose)
            fprintf(' iteration: %d; root:  %f + %f i; function value: %f, derivative:  %f \n',it,real(root),imag(root),fval,derivative);
        end %if
        it=it+1;
        converged=abs(fval)<tol;
    end %if
end %while
it=it-1;

if (~converged)
    warning('Used max number of iterations, or derivative near zero...')
    success=false;
else
    success=true;
end %if

end %function