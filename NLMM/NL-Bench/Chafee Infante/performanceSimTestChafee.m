clear all
close all
clc 

%% Load benchmark
% Chafee Infante
% xdot = A*x + B*u(t) - fcubic(x)

k = 1000; % big model
% L = 5;
L = 1;

tic
[E,A,B,C] = ChafeeInfante_POD_Matrices(k,L); %E=I
toc

%% High Performance Test and Comparison of the different implementations
x0 = [1:1:k]';

tic
fcub1 = fcubic(x0);
toc

tic
fcub2 = fcubic2(x0);
toc

tic
f = fcubic3(k);
fcub3 = f(x0);
toc

%% High Performance Simulation Test and Comparison 
% Definition of simulation time
t0 = 0;
tEnd = 5;
dt = 0.01;
tSim = t0:dt:tEnd;

% Definition of the simulation input
u = @(t) 0.5*(cos(pi*t)+1);
% u = @(t) 1;

%% Implementation 2
x0 = sparse(k,1);
xdot = @(t,x,A,B) A*x + B*u(t) - fcubic2(x);

Jac_xdot = @(t,x,A,B) A - ChafeeInfante_Jac2(x);
optionsFOM = odeset('RelTol',1e-8,'AbsTol',1e-10,'Jac',Jac_xdot);

tic
[~, x] = ode15s(xdot, tSim, x0, optionsFOM, A, B);
y = C*x.';
time.simulation.original2 = toc;

%% Implementation 3
f = fcubic3(k);
fJac = ChafeeInfante_Jac3(k);

x0 = sparse(k,1);
xdot = @(t,x,A,B) A*x + B*u(t) - f(x);

Jac_xdot = @(t,x,A,B) A - fJac(x);
optionsFOM = odeset('RelTol',1e-8,'AbsTol',1e-10,'Jac',Jac_xdot);

tic
[~, x] = ode15s(xdot, tSim, x0, optionsFOM, A, B);
y = C*x.';
time.simulation.original3 = toc;

%% Combination: 2+3
f = fcubic3(k);

x0 = sparse(k,1);
xdot = @(t,x,A,B) A*x + B*u(t) - f(x);

Jac_xdot = @(t,x,A,B) A - ChafeeInfante_Jac2(x);
optionsFOM = odeset('RelTol',1e-8,'AbsTol',1e-10,'Jac',Jac_xdot);

tic
[~, x] = ode15s(xdot, tSim, x0, optionsFOM, A, B);
y = C*x.';
time.simulation.original4 = toc;
