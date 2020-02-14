function [poly,rem] = polydiv(pcoeffs,N)
%{
for a given polynomial

P_n = a_n x^n + a_n-1 x^n-1 + ... + a_2 x^2 + a_1 x + a_0 

with a coefficients vector pcoeffs = [a_n a_n-1 ... a_1 a_0] 
and a given divisor N such that P(N) = 0, determine the deflated 
polynomial Q_n-1 (x) and the remainder R
%}
qcoeffs = zeros(1,length(pcoeffs));

% determine the coefficients of the deflated polynomial
for i = 1:length(pcoeffs)
    if i == 1
        qcoeffs(i) = pcoeffs(i);
    else
        qcoeffs(i) = pcoeffs(i) + qcoeffs(i-1)*N;
    end % if - i
end % for - i

poly = qcoeffs(1:(end-1));
rem = qcoeffs(end);
end