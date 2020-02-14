%% Forward Substitution, Single RHS
%{
Original code, algorithm from https://en.wikipedia.org/wiki/Triangular_matrix#Forward_and_back_substitution
%}

function [x] = forwardsub(L,bL)

% assume L is a lower triangular matrix

% nref=6;                %system size for larger reference problem
% Aref=randn(nref,nref); %augmented matrix containing RHS of system of equations, in practice you'd want to check conditioning...
% Aref = tril(Aref); % convert to lower triaImat = ngular matrix
% bref=randn(nref,1);
nref = length(L);
Aref = L;
bref = bL;

Awork=cat(2,Aref,bref); % Awork is b concatenated onto A

x = zeros(nref,1); % solution vector
sum = 0; % sum for algorithm

an = length(Awork); % get column of b vector
for rw = 1:nref % row counter
    for el = 1:(rw-1) % element counter
        sum = sum + Awork(rw,el)*x(el); % summation from algorithm
    end % for - el
    x(rw) = (Awork(rw,an)-sum)/Awork(rw,rw); % calculation of answer
    sum = 0; % reset summation counter for next iteration
end % for - rw
end % function
% disp('fsub([Aref,bref]) = ');
% disp(x);