%% Gauss-Jordan Elimination

function [Awork] = gjelim(A,varargin)

%{
Simple forward elimination for multiple RHS, with added back-elimination to
convert A -> I
%}

% limited to 3 RHS matrices/vectors, but they have arbitrary dimensions
ninputs = length(varargin);

% begin with empty matrices, overwrite if specified in input
bref = [];
bref2 = [];
bref3 = [];
outmsg = 'Gauss-Jordan([Aref]) = ';

switch nargin
    case 2 % only A and b specified, but no constraints on dimensions
        bref = varargin{ninputs};
        outmsg = 'Gauss-Jordan([Aref,bref]) = ';
    case 3
        bref = varargin{ninputs-1};
        bref2 = varargin{ninputs};
        outmsg = 'Gauss-Jordan([Aref,bref,bref2]) = ';
    case 4
        bref = varargin{ninputs-2};
        bref2 = varargin{ninputs-1};
        bref3 = varargin{ninputs};
        outmsg = 'Gauss-Jordan([Aref,bref,bref2,bref3]) = ';
    otherwise
        ermsg = 'Need 1-4 inputs.';
        err(ermsg)
end % switch - nargin

% Q = [2 1 -1;-3 -1 2;-2 1 2]; % simple solved test matrix
% q = [8;-11;-3]; % simple solved test matrix
Aref = A;
nref=length(Aref);

%note that the elimination procedure coded below modifies the matrix B
% Awork=cat(2,Aref,bref);          %This is our working version of the matrix used to perform elimination (pvtrw.e. it will be modified)
Awork=cat(2,Aref,bref,bref2,bref3);
nans = size(Awork,2); % number of columns in Awork

% forward elimination
for pvtrw = 1:(nref-1)
%     disp(Awork) % debug statement
    for elmrw = (pvtrw+1):nref % row below pivot row
        fact = Awork(elmrw,pvtrw)/Awork(pvtrw,pvtrw); % elimination multiplier
%         fprintf('R%d - (%.2f)R%d -> R%d\n',elmrw,fact,pvtrw,elmrw) % debug statement
        for elem = pvtrw:nans
            Awork(elmrw,elem) = Awork(elmrw,elem)-fact*Awork(pvtrw,elem);
        end % for - elem
    end % for - elmrw
end % for - pvtrw
% fprintf('Forward elim done.\n') % debug statement
% halfway = Awork; % debug statement to check result of forward elimination

% backwards elimination to turn into diag matrix
for pvtrw = nref:-1:1
%     disp(Awork) % debug statement
    for elmrw = (pvtrw-1):-1:1
%         fprintf('Element being eliminated: (%d,%d)\n',elmrw,pvtrw) % debug statement
        fact = Awork(elmrw,pvtrw)/Awork(pvtrw,pvtrw); % pvtrw is value on diagonal, elmrw is value in same column
%         fprintf('R%d - (%.2f)R%d -> R%d\n',elmrw,fact,pvtrw,elmrw) % debug statement
        for elem = pvtrw:nans % ignore the rows that have already been eliminated
            Awork(elmrw,elem) = Awork(elmrw,elem)-fact*Awork(pvtrw,elem); %
        end % for - elem
    end % for - elmrw
end % for - pvtrw

% retrace diagonal and divide corresponding row by diagonal value to
% transform A into identity matrix
for i = 1:nref % iterate rows
    for j = nref+1:nans % iterate RHS columns
        Awork(i,j) = Awork(i,j)/Awork(i,i); % divide answer
    end % for - j
    Awork(i,i) = 1; % set diagonal to 1
end % for - i

% disp(outmsg);
% disp(Awork);
Awork = Awork(:,9:16); % use case: calculation of inverse
% dif = Awork(:,9:16)-inv(A);
% dif = Awork(:,nref+1)-(Aref\bref2);
% disp(dif)

end % function