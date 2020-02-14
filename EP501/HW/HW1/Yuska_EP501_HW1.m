%{
EP501 HW1
Paul Yuska

RUN THIS SCRIPT.

Main script.
Assumes test .mat files and instructor-provided functions are on the search
path.
%}

clc
fprintf('PAUL YUSKA - EP501 - HW1 - 10/3/19\n\n')

%% Load variables into workspace if they are not already there
if exist('A','var')==0 || exist('b','var')==0 || exist('bref2','var')==0 || exist('bref3','var')==0
    load testproblem.mat
end % if - exist

if exist('L','var')==0 || exist('bL','var')==0
    load lowertriang_testproblem.mat
end % if - exist

if exist('Ait','var')==0 || exist('bit','var')==0
    load iterative_testproblem.mat
end % if - exist

%% Problem 1

% Solve the test problem from the repository using your function and the 
% back-substitution function in the course repository. Show that the answer
% is the same as the provided Gauss_elim function and MATLAB's backslash
% '\' operator.

ans_s = backsub(forwardelim(A,b)); % call student-written function
[Amod,ord] = Gauss_elim(A,b); % call provided Gauss elim. function
ans_g = backsub(Amod(ord,:)); 
ans_m = A\b; % MATLAB built-in function

M = cat(2,ans_s,ans_g,ans_m);
col_titles = {'Student_Code','Provided_Gauss','MATLAB_Backslash'};

fprintf('\n')
disp('Problem 1: Parts A & B')
fprintf('Results of forward elimination -> back substitution\n\n')
disp(array2table(M,'VariableNames',col_titles))

disp('#######################################################')
fprintf('\n')

ans_s = forwardsub(L,bL); % call student-written function
ans_m = L\bL; % MATLAB built-in function
M = cat(2,ans_s,ans_m);
col_titles = {'Student_Code','MATLAB_Backslash'};

disp('Problem 1: Parts C & D')
fprintf('Results of forward substitution\n\n')
disp(array2table(M,'VariableNames',col_titles))

disp('#######################################################')
fprintf('\n')

%% Problem 2

% Calculate inverse of A using Gauss-Jordan elimination

Imat = eye(length(A)); % create identity matrix for calc. of inverse
ans_s = gjelim(A,Imat);
ans_m = inv(A);

disp('Problem 2')
fprintf('Calculation of Inverse of A\n\n')
disp('Student Code')
disp(ans_s)
disp('MATLAB Solution')
disp(ans_m)


disp('#######################################################')
fprintf('\n')

%% Problem 3

% Solve test problems using LU decomposition
[L1,U1] = doolittle(A); % perform LU decomposition

bp = forwardsub(L1,b);
ans_sb = backsub(cat(2,U1,bp)); % answer for RHS b
bp2 = forwardsub(L1,b2);
ans_sb2 = backsub(cat(2,U1,bp2)); % answer for RHS bref2
bp3 = forwardsub(L1,b3);
ans_sb3 = backsub(cat(2,U1,bp3)); % answer for RHS bref3

disp('Problem 3: Part C')
fprintf('Solving Multiple RHS''s via LU Decomposition\n\n')

M = cat(2,ans_sb,ans_sb2,ans_sb3);
col_titles = {'b_Solution','b2_Solution','b3_Solution'};
disp(array2table(M,'VariableNames',col_titles))
M = cat(2,A\b,A\b2,A\b3);
col_titles = {'MATLAB_b_Sol','MATLAB_b2_Sol','MATLAB_b3_Sol'};
disp(array2table(M,'VariableNames',col_titles))

disp('#######################################################')
fprintf('\n')

% Calculate inverse of A using Doolittle LU factorization

Imat = eye(length(A)); % set up identity matrix
invA = zeros(length(A));

for k = 1:length(A)
    bprime = forwardsub(L1,Imat(:,k));
    invA(:,k) = backsub(cat(2,U1,bprime)); % k'th column of inverse matrix
end % for - s

disp('Problem 3: Part D')
fprintf('Calculation of Inverse of A via LU Decomposition\n\n')

disp('Student Solution')
disp(invA)
disp('Matlab inv() Solution')
disp(inv(A))

disp('#######################################################')
fprintf('\n')

%% Problem 4
x0 = zeros(length(Ait),1);
ans_s = sor(x0,Ait,bit,1e-6,1);
ans_m = Ait\bit;

disp('Problem 4')
fprintf('Solving Iteratively via Successive Over-Relaxation\n\n')

M = cat(2,ans_s,ans_m);
col_titles = {'Student_Solution','MATLAB_Solution'};
disp(array2table(M,'VariableNames',col_titles))

disp('#######################################################')
fprintf('\n')

%% Problem 5
disp('Problem 5')
fprintf('Calculating Determinants via Gaussian Elimination\n\n')

ans_s = gaussdet(A);
ans_m = det(A);

disp('Student Solution')
disp(ans_s)
disp('Matlab det() Solution')
disp(ans_m)

disp('#######################################################')
fprintf('\n')