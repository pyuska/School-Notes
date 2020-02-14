%% Vanilla forward elimination, single RHS
% nref=6;                %system size for larger reference problem
% Aref=randn(nref,nref);    %augmented matrix containing RHS of system of equations, in practice you'd want to check conditioning...
% bref=randn(nref,1);    %RHS
function [Awork] = forwardelim(A,b)
nref=length(A);
Aref = A;
bref = b;

%note that the elimination procedure coded below modifies the matrix B
% Awork=cat(2,Aref,bref);          %This is our working version of the matrix used to perform elimination (it will be modified)
Awork=cat(2,Aref,bref);
nans = size(Awork,2);

for pvtrw = 1:(nref-1)
    for elmrw = (pvtrw+1):nref % elimination row, below pivot row
        fact = Awork(elmrw,pvtrw)/Awork(pvtrw,pvtrw); % elimination multiplier
        for elem = pvtrw:nans
            Awork(elmrw,elem) = Awork(elmrw,elem)-fact*Awork(pvtrw,elem);
        end % for - elem
    end % for - elmrw
end % for - pvtrw
end % function
% disp('elim([Aref,bref]) = ');
% disp(Awork);