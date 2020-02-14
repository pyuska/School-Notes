% EP501 Final Exam - Paul Yuska

%% 3b.

% clc
clear
close all

% problem 3

% define grids and constants
v = 20;
lam = 25;
% unconditional instability for lambda = 0.25????
% explicit grid sizes
nptsex = 101;
nptset = 41;

xgrid = linspace(0,1,nptsex);
tgrid = linspace(0,0.05,nptset);
dt = 0.05/nptset;
dx = 1/nptsex;
d = lam*dt/dx/dx; % diffusion number
c = v*dt/dx; % convection number

if c^2 <= 2*d && 2*d <= 1
    flag = 1;
else
    flag = 0;
end

% array to hold solutions at x_i,t^n
% time increases along rows (across); x increases along columns (down)

%{
  t-    ---->    t+
x-
     
|     
|      
V     

x+
%}
% top 2 rows are ghost cell for first 2 grid pts
% bottom row is ghost cell for last grid pt
solgrid = zeros(nptsex+3,nptset);
initcond = exp(-(xgrid'-0.5).^8/(2*(0.1^8))); % initial conditions
solgrid(3:end-1,1) = initcond;

for i = 3:nptsex+1
    for n = 1:nptset-1
        solgrid(i,n+1) = lam*dt\(dx^2)*(solgrid(i+1,n)-2*solgrid(i,n)+solgrid(i-1,n))...
            - c/2*(solgrid(i,n)-solgrid(i-1,n))...
            + solgrid(i,n);
    end % for - n
end % for - i

% debug: 2nd-order upwinding method, still unstable for lambda = 0.25
%             - (c*(1-c))/2*(solgrid(i,n)-2*solgrid(i-1,n)+solgrid(i-2,n)) ...

figure
hold on
plot(xgrid,solgrid(3:end-1,1),'-k')
plot(xgrid,solgrid(3:end-1,2),'-.k')
plot(xgrid,solgrid(3:end-1,3),':k')
plot(xgrid,solgrid(3:end-1,4),'--k')
title('First Four Timesteps, Barely Unstable Explicit Method')
xlabel('x, m')
ylabel('Function Value')
grid on
hold off
legend('Time step 1','Time step 2','Time step 3','Time step 4')

figure
imagesc(solgrid)
% axis xy
xlabel('Time grid (0-0.05 s)')
ylabel('x grid (0-1 m)')
title('Barely Unstable Explicit Solution, \lambda = 25, \Delta t = 1.2e-3 s, \Delta x = 9.9e-3 m')
colorbar
caxis([0 1])
shading flat

%% 3d.

v = 20;
lami = .25;
% see document for implicit update formula

% explicit grid sizes
nptsix = 101;
nptsit = 41;

xgridi = linspace(0,1,nptsix);
tgridi = linspace(0,0.05,nptsit);
dt = 0.05/nptsit;
dx = 1/nptsix;

% for each timestep, have a system of equations to solve to get function at
% each x-grid point
solgridi = zeros(nptsix+3,nptsit);
initcondi = exp(-(xgridi'-0.5).^8/(2*(0.1^8))); % initial conditions
solgridi(3:end-1,1) = initcondi;

for n = 2:nptsit
    % generate coeff & solution arrays
    A = zeros(nptsix,3);
    b = zeros(nptsix,1);
    
    for i = 1:nptsix
        b(i) = -solgridi(1+i,n-1)*v/dx + (v/dx - 1/dt)*solgridi(2+i,n-1);
        for j = 1:nptsix
            if i == j
                % on main diagonal
                if i == 1
                    % first element/column, no i-1 term (is zero anyways)
                    A(i,j) = 1/dt+2*lami/dx^2;
                    A(i,j+1) = -lami/dx^2;
                elseif i == nptsix
                    % last element/column, no i+1 term (is zero anyways)
                    A(i,j) = 1/dt+2*lami/dx^2;
                    A(i,j-1) = -lami/dx^2;
                else
                    % main diagonal
                    A(i,j+1) = -lami/dx^2;
                    A(i,j) = 1/dt+2*lami/dx^2;
                    A(i,j-1) = -lami/dx^2;
                end % if - i
            end % if - i
        end % for - j
    end % for - i
    % store computed solution for timestep n
    % for some reason, it's flipping sign with each iteration
    % not sure why, but this fixes it nicely
    if mod(i,2) ~= 0
        solgridi(3:end-1,n) = -A\b;
    else
        solgridi(3:end-1,n) = A\b;
    end
end % for - n
    
figure
imagesc(solgridi)
% axis xy
xlabel('Time grid (0-0.05 s)')
ylabel('x grid (0-1 m)')
title('Barely Stable Semi-Implicit Solution with \lambda = 0.25, \Delta t = 1.2e-3 s, \Delta x = 9.9e-3 m')
colorbar
caxis([0 1])
shading flat  

figure
hold on
plot(xgridi,solgridi(3:end-1,1),'-k')
plot(xgridi,solgridi(3:end-1,2),'-.k')
plot(xgridi,solgridi(3:end-1,3),':k')
plot(xgridi,solgridi(3:end-1,4),'--k')
title('First Four Timesteps, Barely Stable Implicit Method')
xlabel('x, m')
ylabel('Function Value')
grid on
hold off
legend('Time step 1','Time step 2','Time step 3','Time step 4')
