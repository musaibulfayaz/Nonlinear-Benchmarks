clear all; close all; clc;

%% simulation parameter
k = 400;
nu = 0.01;
L = 5;

t0 = 0;
tEnd = 2;
dt = 0.01;
tSim = t0:dt:tEnd;

% excitation signals
u = @(t) 0.5*(cos(2*pi*t/10)+1);
% u = @(t) exp(-t);
% u = @(t) 1;
% u = @(t) 2*sin(pi*t);

%% Simulation of pure nonlinear Burger model (with inherent quadratic nonlinearity)
% TODO

% x0 = zeros(length(A),1);
% xdot = @(t,x,A,N,B) A*x + x.*x + N*x*u(t) + B*u(t);
% Jac_xdot = @(t,x,A,N,B) A + 2*x + N*u(t);
% optionsFOM = odeset('RelTol',1e-8,'AbsTol',1e-10,'Jac',Jac_xdot);
% tic
% [~, x] = ode15s(xdot, tSim, x0, optionsFOM, A, N, B);
% y = C*x.';
% time.simulation.original = toc;
% 
% figure;
% plot(tSim, y);
% xlabel('Time (s)');
% ylabel('Output y');

%% Simulation of the QB-version of Burger model
[E,A,H,N,B,C,v0] = Burgers_Matrices_nonzero_ic(k,nu); % E=I

nQB = size(A,1)

% Function and Jacobian of rhs of SISO QB-FOM (E=I)
xdotQB = @(t,x,A,H,N,B) A*x + H*(kron(x,x)) + (N*x + B)*u(t);
Jac_xdot = @(t,x,A,H,N,B) A + 2*H*(kron(x,speye(size(A,1)))) + N*u(t); 
% optionsFOM = odeset('RelTol',1e-3,'AbsTol',1e-6,'Jac',Jac_xdot);
optionsFOM = odeset('Jac',Jac_xdot);

x0QB = zeros(length(A),1);

tic
[~, xQB] = ode15s(xdotQB, tSim, x0QB, optionsFOM, A, H, N, B);
simTimeFOMQB = toc
xtQB = xQB + ones(size(xQB,1),1)*v0.';

discPoints = 0:L/k:(L-L/k);
intx = 15;
intt = 5;

figure;
surf(discPoints(1:intx:end),tSim(1:intt:end),xtQB(1:intt:end,1:intx:end));
xlabel('x (m)')
ylabel('Time (s)')
zlabel('v(x,y)')

%% Simulation of the bilinear version of Burger model
% TODO
% [A,N,B,C] = Burgers_Equation(k,nu);  % E=I
% 
% nBil = size(A,1)
% 
% N = reshape(N,[nBil,nBil]);
% 
% % Function and Jacobian of rhs of SISO Bil-FOM (E=I)
% xdotBil = @(t,x,A,N,B) A*x + (N*x + B)*u(t);
% Jac_xdotBil = @(t,x,A,N,B) A + N*u(t); 
% optionsFOMBil = odeset('RelTol',1e-6,'AbsTol',1e-9,'Jac',Jac_xdotBil);
% 
% x0 = zeros(length(A),1);
% xdot = @(t,x,A,H,N,B) A*x + H*kron(x,x) + N*x*u(t) + B*u(t);
% Jac_xdot = @(t,x,A,H,N,B) A + 2*H*(kron(x,speye(size(A,1)))) + N*u(t);
% optionsFOM = odeset('RelTol',1e-8,'AbsTol',1e-10,'Jac',Jac_xdot);
% tic
% [~, xt] = ode15s(xdot, tSim, x0, optionsFOM, A, H, N, B);
% original = toc
% x = xt + ones(size(xt,1),1)*v0.';
% 
% discPoints = 0:L/k:(L-L/k);
% figure(1)
% surf(discPoints(1:5:end),tSim,x(:,1:5:end));
% xlabel('x (m)')
% ylabel('Time (s)')
% zlabel('v(x,y)')
