%% Doolittle LU Factorization

function [L,U] = doolittle(A)

% simple elimination, multiple RHS

% Q = [2 1 -1;-3 -1 2;-2 1 2]; % simple solved test matrix
% q = [8;-11;-3]; % simple solved test matrix
nref=length(A);
U = A;
L = zeros(nref);

%note that the elimination procedure coded below modifies the matrix B
% U = cat(2,Aref,bref);          %This is our working version of the matrix used to perform elimination (pvtrw.e. it will be modified)
% U = cat(2,Aref,bref,bref2,bref3);
nans = size(U,2); % number of b vectors concatenated to A

for pvtrw = 1:(nref-1)
    for elmrw = (pvtrw+1):nref % row below pivot row
        fact = U(elmrw,pvtrw)/U(pvtrw,pvtrw); % first nonzero element of pivot row
        for elem = pvtrw:nans
            U(elmrw,elem) = U(elmrw,elem)-fact*U(pvtrw,elem);
            if elem > pvtrw
                continue
            elseif elem == pvtrw
                L(elmrw,elem) = fact;
                L(elem,elem) = 1;
            end % if
        end % for - elem
    end % for - elmrw
end % for - pvtrw
L(nref,nref) = 1; % manually set last element to 1 to finish L matrix

end % function
% disp('elim([Aref,bref]) = ');
% disp(U);
% disp(L);