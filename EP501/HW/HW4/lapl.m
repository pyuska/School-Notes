% HW4 Laplacian MATLAB Implelementation

close all

Q = 1;
a = 1;
eps = 8.854e-12;

nPts = 30;

x = linspace(-3*a,3*a,nPts);
y = linspace(-3*a,3*a,nPts);
z = linspace(-3*a,3*a,nPts);
laplac = zeros(nPts,nPts,nPts);

[X,Y,Z] = meshgrid(x,y,z);

% magnitude < a
for i = 1:nPts
    for j = 1:nPts
        for k = 1:nPts
            mag = sqrt(x(i)^2 + y(j)^2 + z(k)^2);
            
            if mag < a
                laplac(i,j,k) = (3*Q)/(4*a^3*eps*pi);
            else
                laplac(i,j,k) = 0;
            end
        end
    end
end

s = slice(X,Y,Z,laplac,0,0,0);
s(1).EdgeColor = 'none';
s(2).EdgeColor = 'none';
s(3).EdgeColor = 'none';
colorbar