function res = polyeval(coef,x)
%{
assumes that the function is a polynomial of degree n in the form 

P_n(x) = a_n x^n + a_n-1 x^n-1 + ... + a_2 x^2 + a_1 x + a_0

that can be represented by a vector of size n of polynomial coefficients

coef = [a_n a_n-1 ... a_2 a_1 a_0]  .

Returns the value of the polynomial at the supplied x-value.
%}

res = 0;
len = length(coef);

for i = 1:len
    res = res + coef(i)*x^(len-i);
end % for - i

end % function