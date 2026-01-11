clear all
close all
clc 

%% Load benchmark
% Chafee Infante
% xdot = A*x + B*u(t) - fcubic(x)

k = 10000;
% L = 5;
L = 1;

[E,A,B,C] = ChafeeInfante_POD_Matrices(k,L); %E=I

%% Simulations
% Definition of simulation time
t0 = 0;
tEnd = 5;
dt = 0.01;
tSim = t0:dt:tEnd;

% Definition of the simulation input
u = @(t) 0.5*(cos(pi*t)+1);
% u = @(t) 1;

f = fcubic3(k);
fJac = ChafeeInfante_Jac3(k);

% x0 = zeros(length(A),1);
x0 = sparse(k,1);
xdot = @(t,x,A,B) A*x + B*u(t) - f(x);

Jac_xdot = @(t,x,A,B) A - fJac(x);
optionsFOM = odeset('RelTol',1e-8,'AbsTol',1e-10,'Jac',Jac_xdot);

tic
[~, x] = ode15s(xdot, tSim, x0, optionsFOM, A, B);
y = C*x.';
time.simulation.original = toc;

% Output over space
discPoints = 0:L/k:(L-L/k);
figure(1)
surf(discPoints(1:5:end),tSim,x(:,1:5:end));
xlabel('x (m)')
ylabel('Time (s)')
zlabel('v(x,y)')

% Output over time
figure(2)
plot(tSim,y,'b');
grid on
legend('Full-order')
xlabel('Time(s)')
ylabel('v_1')
