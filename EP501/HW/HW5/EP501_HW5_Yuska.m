% EP501 HW5 Paul Yuska

clc
clear
close all

%% 1a. Plot the dielectric function

% define constants and parameters
npts = 20; % size of grid
eps0 = 8.854e-12; % permittivity of free space, F/m
a = 0.01; % m
l = a/5; % m
xprime = -9/10*a; % m
xdprime = 9/10*a; % m

x = linspace(-a,a,npts);
deltax = x(2)-x(1);

dPhidxb = 1000; % initial condition; dPhi/dx(-a) = 1000
Phib = 100; % Phi(a) = 100

% define the dielectric function
eps = 10*eps0*(tanh((x-xprime)/l)-tanh((x-xdprime)/l));


%% 1b. Develop a system of FDEs

A = zeros(npts); % coefficient matrix, should be tridiag
b = zeros(npts,1); % solution vector

% populate vectors
for i = 1:npts
    % populate b vector
    if i == 1
        b(i) = dPhidxb*deltax;
    elseif i == npts
        b(i) = Phib;
    end % if - i
    
    % main tridiag of A
    for j = 1:npts
        if i == j
            if i == 1
                % forward difference at x = -a, (4) in writeup
                A(i,j) = -1;
                A(i,j+1) = 1;
            elseif i == npts
                % boundary condition at x = a, (5) in writeup
                A(i,j) = 1;
            else
                % FDE, (3) in writeup
                A(i,j-1) = 1-(eps(i+1)-eps(i-1))/4/eps(i);
                A(i,j) = -2;
                A(i,j+1) = 1+(eps(i+1)-eps(i-1))/4/eps(i);
            end
        end
    end % for - j
end % for - i

A2o = A; % second-order at edges
b2o = b; % second-order at edges
b2o(1) = 2*dPhidxb*deltax; % modify solution
A2o(1,1) = -3;
A2o(1,2) = 4;
A2o(1,3) = -1;

xsol = A\b;
xsol2o = A2o\b2o;

fprintf('\n')
disp('Problem 1b,c:')
disp('See attached document for equations.')
fprintf('\n')

figure
hold on
plot(x,eps,'-k','LineWidth',2)
title('Dielectric Function \epsilon and Electrostatic Potential \Phi')
xlabel('x, m')
ylabel('\epsilon(x), F')
yyaxis right
grid on
plot(x,xsol,'-.k','LineWidth',2)
plot(x,xsol2o,':k','LineWidth',2)
ylabel('\Phi(x), V')
legend('\epsilon','\Phi (1st order)','\Phi (2nd order)','Location','SouthEast')

%% Problem 2a.

% define constants
m = 1.67e-27; % mass, kg
q = 1.6e-19; % charge, C
Bz = 50000e-9; % magnetic field magnitude, T
npts = 100;

% time grid initialization
per = 2*pi*m/q/Bz; % period of oscillation, s
tmin = 0;
deltat = per/npts;
t = linspace(tmin,per,npts); % 1 period

% initialize arrays for vx, vy ODE solutions
vx = zeros(1,npts);
vx(1) = 1000; % initial conds
vy = zeros(1,npts);
vy(1) = 1000; % initial conds
vz = 100; % const. upward velocity

for i = 2:npts
    dvx1 = q/m*Bz*vy(i-1)*deltat;
    dvy1 = -q/m*Bz*vx(i-1)*deltat;
    
    dvx2 = q/m*Bz*(0.5*dvy1 + vy(i-1))*deltat;
    dvy2 = -q/m*Bz*(0.5*dvx1 + vx(i-1))*deltat;
    
    dvx3 = q/m*Bz*(0.5*dvy2 + vy(i-1))*deltat;
    dvy3 = -q/m*Bz*(0.5*dvx2 + vx(i-1))*deltat;
    
    dvx4 = q/m*Bz*(vy(i-1) + dvy3)*deltat;
    dvy4 = -q/m*Bz*(vx(i-1) + dvx3)*deltat;
    
    vx(i) = vx(i-1) + (dvx1 + 2*dvx2 + 2*dvx3 + dvx4)/6;
    vy(i) = vy(i-1) + (dvy1 + 2*dvy2 + 2*dvy3 + dvy4)/6;
end % for - i

x2 = cumtrapz(t,vx);
y2 = cumtrapz(t,vy);

figure
hold on
plot(t,vx,'-k')
plot(t,vy,'-.k')
xlabel('Time, s')
ylabel('Velocity, m/s')
yyaxis right
plot(t,x2,':k')
plot(t,y2,'--k')
ylabel('Position, m')
grid on
legend('Velocity x-component','Velocity y-component','Position x-component','Position y-component','Location','North')
title('Position and Velocity of Charged Particle in Constant B-Field (1 Period)')

disp('Problem 2b:')
disp('Although the RK2 method implemented in the repo models a slightly')
disp('different system, the general shape of the velocity and position')
disp('plots are the same (circular motion), as they should be.')

%% Problem 2b.

% modify B-field to vary linearly in y-dir
% same as 2a, but B needs to be calculated for each step
% need to integrate vel. to get position along the way

% define constants
m = 1.67e-27; % mass, kg
q = 1.6e-19; % charge, C
Bz = 50000e-9; % magnetic field magnitude, T
npts = 1000;

% time grid initialization
per = 2*pi*m/q/Bz; % period of oscillation, s
tmin = 0;
deltat = per/npts;
t = linspace(tmin,per*10,npts); % 1 period

% initialize arrays for vx, vy ODE solutions
vx = zeros(1,npts);
vx(1) = 1000; % initial conds
vy = zeros(1,npts);
vy(1) = 1000; % initial conds
vz = 0; % const. upward velocity
y2 = zeros(1,npts);

for i = 1:npts-1
    dvx1 = q/m*Bz*(1+0.5*y2(i))*vy(i)*deltat;
    dvy1 = -q/m*Bz*(1+0.5*y2(i))*vx(i)*deltat;
    
    dvx2 = q/m*Bz*(1+0.5*y2(i))*(0.5*dvy1 + vy(i))*deltat;
    dvy2 = -q/m*Bz*(1+0.5*y2(i))*(0.5*dvx1 + vx(i))*deltat;
    
    dvx3 = q/m*Bz*(1+0.5*y2(i))*(0.5*dvy2 + vy(i))*deltat;
    dvy3 = -q/m*Bz*(1+0.5*y2(i))*(0.5*dvx2 + vx(i))*deltat;
    
    dvx4 = q/m*Bz*(1+0.5*y2(i))*(vy(i) + dvy3)*deltat;
    dvy4 = -q/m*Bz*(1+0.5*y2(i))*(vx(i) + dvx3)*deltat;
    
    vx(i+1) = vx(i) + (dvx1 + 2*dvx2 + 2*dvx3 + dvx4)/6;
    vy(i+1) = vy(i) + (dvy1 + 2*dvy2 + 2*dvy3 + dvy4)/6;
    
    g = cumtrapz(t(1:i+1),vy(1:i+1));
    y2(i+1) = g(end);
end % for - i

z2 = vz*t;
x2 = cumtrapz(t,vx);

figure
hold on
plot(t,vx,'-k')
plot(t,vy,'-.k')
xlabel('Time, s')
ylabel('Velocity, m/s')
yyaxis right
plot(t,x2,':k')
plot(t,y2,'--k')
ylabel('Position, m')
hold off
grid on
legend('Velocity x-component','Velocity y-component','Position x-component','Position y-component','Location','se')
title('Position and Velocity of Charged Particle in Linearly-Varying B-Field')

%{
I know the modified B-field part is incorrect
I don't know what's wrong with it :/ 
%}

% figure
% comet3(x2,y2,z2)
% xlabel('x')
% ylabel('y')
% zlabel('z')

%% Problem 3

tmin3 = 0;
tmax3 = 0.5;
npts = 20;

y31 = zeros(1,npts);
y32 = zeros(1,npts);

y31(1) = 1;
y32(1) = 1;

%{
???????????
I have no idea how to solve this...
My brain is fried from my compressible aero final that i took earlier.
Sorry I couldn't get this done. 
%}