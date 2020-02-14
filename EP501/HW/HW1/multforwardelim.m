%% Simple Elimination with Multiple RHS Vectors

function [Awork] = multforwardelim(A,varargin)

% limited to 3 RHS matrices/vectors

% begin with empty matrices, overwrite if specified in input
bref = [];
bref2 = [];
bref3 = [];
outmsg = 'elim([Aref]) = ';

switch nargin
    case 1
        Aref = A;
    case 2 % only A and b specified, but no constraints on dimensions
        bref = b;
        outmsg = 'elim([Aref,bref]) = ';
    case 3
        bref = b;
        bref2 = b2;
        outmsg = 'elim([Aref,bref,bref2]) = ';
    case 4
        bref = b;
        bref2 = b2;
        bref3 = b3;
        outmsg = 'elim([Aref,bref,bref2,bref3]) = ';
    otherwise
        ermsg = 'Need 1-4 inputs.';
        error(ermsg)
end % switch - nargin

nref=size(A,1);
Aref = A;

%note that the elimination procedure coded below modifies the matrix B
% Awork=cat(2,Aref,bref);          %This is our working version of the matrix used to perform elimination (pvtrw.e. it will be modified)
Awork=cat(2,Aref,bref,bref2,bref3);
nans = size(Awork,2); % number of b vectors concatenated to A

for pvtrw = 1:(nref-1)
    for elmrw = (pvtrw+1):nref % row below pivot row
        fact = Awork(elmrw,pvtrw)/Awork(pvtrw,pvtrw); % first nonzero element of pivot row
        for elem = pvtrw:nans
            Awork(elmrw,elem) = Awork(elmrw,elem)-fact*Awork(pvtrw,elem);
        end % for - elem
    end % for - elmrw
end % for - pvtrw

disp(outmsg);
disp(Awork);

% backsub(<A1, A2, or A3>) gives the answer it should, so this section of
% the code SHOULD be working correctly?
% A1 = Awork(:,1:9);
% A2 = cat(2,Awork(:,1:8),Awork(:,10));
% A3 = cat(2,Awork(:,1:8),Awork(:,11));

end % function