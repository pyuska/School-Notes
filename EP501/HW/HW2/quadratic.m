function [root1,root2] = quadratic(coeffs)

% function finds roots of quadratic specified by coefficient vector
% 'coeffs' = [a b c] --> P_2(x) = ax^2 + bx + c

a = coeffs(1);
b = coeffs(2);
c = coeffs(3);

root1 = (-b + sqrt(b^2 - 4*a*c))/(2*a);
root2 = (-b - sqrt(b^2 - 4*a*c))/(2*a);

end % function