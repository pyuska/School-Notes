function [x,nit]=sor(x0,A,b,tol,omega,verbose)

%{
Successive over-relaxation based on Jacobi script from EP501 repo
%}

%% Check the inputs
narginchk(3,6);
if nargin<4
    tol=1e-6;
end %if
if nargin<6
    verbose=false;
end %if


%% Setup iterations
maxit=200;    %max number of iterations
n=size(A,1);  %system size
residual=10*ones(n,1);
difftot=1e3+tol;   %max sure we enter iterations
x=x0;


%% Perform iterations
it=1;
while(difftot>tol && it<=maxit)
    difftotprev=difftot;
    resprev=residual;
    xprev=x;
%     fprintf('xprev = %.4f\n',xprev)
    for i=1:n
        residual(i)=b(i);
%         fprintf('residual(i) = %.4f\n',b(i))
        for j=1:n
            if j>(i-1)
%                 disp('j>(i-1)')
%                 disp('residual(i)=residual(i)-A(i,j)*xprev(j)')
%                 fprintf('i = %d; j = %d; A(i,j) = %.4f; xprev(j) = %.4f\n',i,j,A(i,j),xprev(j))
%                 disp(A(i,j)*xprev(j))
%                 fprintf('residual(i) = %.4f - (%.4f) = %.4f\n\n',residual(i),A(i,j)*xprev(j),residual(i)-A(i,j)*xprev(j))
                residual(i)=residual(i)-A(i,j)*xprev(j);
            else
%                 disp('j<=(i-1)')
%                 disp('residual(i)=residual(i)-A(i,j)*xprev(j)-A(i,j)*x(j)')
%                 fprintf('i = %d; j = %d; A(i,j) = %.4f; xprev(j) = %.4f; ',i,j,A(i,j),xprev(j))
%                 disp(A(i,j)*xprev(j))
%                 fprintf('x(j) = %.4f\n',x(j))
%                 disp(A(i,j)*x(j))
%                 fprintf('residual(i) = %.4f - (%.4f) - (%.4f) = %.4f\n\n',residual(i),A(i,j)*xprev(j),A(i,j)*x(j),residual(i)-A(i,j)*xprev(j)-A(i,j)*x(j))
                residual(i)=residual(i)-A(i,j)*x(j);
            end % if - j
        end %for
        x(i)=xprev(i)+omega*residual(i)/A(i,i);
%         disp('x(i)=xprev(i)+omega*residual(i)/A(i,i)')
%         fprintf('i = %d; x(i) = %.4f; xprev(i) = %.4f; omega = %.4f; residual(i) = %.4f; A(i,i) = %.4f\n\n',i,x(i),xprev(i),omega,residual(i),A(i,i))
    end %for
    difftot=sum(abs(residual-resprev));
    
    if (verbose)
        fprintf('x= ');
        for i=1:n
            fprintf('%f   ',x(i));
        end %for
        fprintf('\n');
        fprintf('it=%d; difftot = %e\n',it,difftot);
    end %if
    
    if (difftot>difftotprev)
        error('Solution appears to be diverging, check diagonal dominance...')
    end %if
%     fprintf('##################### END ITERATION %d #####################\n\n',it)
    it=it+1;
end %while

nit=it-1;
if (n==maxit)
    warning('Maximum iterations reached. Solution may not have converged fully...')
end %if

end %function
