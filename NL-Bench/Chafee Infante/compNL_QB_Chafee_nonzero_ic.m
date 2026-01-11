clear all; close all; clc;

%==========================================================================
%==========================================================================
% FIXME - TODO - Not completely done yet!
%==========================================================================
%==========================================================================

%% simulation parameter
k = 200;
L = 5;

t0 = 0;
tEnd = 2;
dt = 0.01;
tSim = t0:dt:tEnd;

% excitation signals
u = @(t) 0.5*(cos(pi*t)+1);
% u = @(t) 1;

%% Simulation of the QB-version of Chafee-Infante model
[E,A,H,N,B,C,v0QB] = ChafeeInfante_Matrices_nonzero_ic(k,L); %E=I

nQB = size(A,1)

% Function and Jacobian of rhs of SISO QB-FOM (E=I)
xdotQB = @(t,x,A,H,N,B) A*x + H*(kron(x,x)) + (N*x + B)*u(t);
Jac_xdot = @(t,x,A,H,N,B) A + 2*H*(kron(x,speye(size(A,1)))) + N*u(t); 
optionsFOM = odeset('RelTol',1e-8,'AbsTol',1e-10,'Jac',Jac_xdot);

x0QB = zeros(2*k,1);

tic
[~, xQB] = ode15s(xdotQB, tSim, x0QB, optionsFOM, A, H, N, B);
simTimeFOMQB = toc
xtQB = xQB + ones(size(xQB,1),1)*v0QB.';
yQB = C*xtQB.';

discPoints = 0:L/(2*k):(L-L/(2*k));
intx = 15;
intt = 5;

figure;
surf(discPoints(1:intx:end),tSim(1:intt:end),xtQB(1:intt:end,1:intx:end));
title('QB-FOM')
xlabel('x (m)')
ylabel('Time (s)')
zlabel('v(x,y)')

%% Simulation of pure nonlinear Chafee-Infante model (combination: 2+3)
[E,A,B,C] = ChafeeInfante_POD_Matrices(k,L); %E=I

nNL = size(A,1)

v0NL = v0QB(1:nNL);
% v0NL = v0QB(k+1:end);

f = fcubic3(k);

x0 = sparse(k,1);
xdot = @(t,x,A,B) A*x + B*u(t) - f(x);

Jac_xdot = @(t,x,A,B) A - ChafeeInfante_Jac2(x);
optionsFOM = odeset('RelTol',1e-8,'AbsTol',1e-10,'Jac',Jac_xdot);

tic
[~, xNL] = ode15s(xdot, tSim, x0, optionsFOM, A, B);
simTimeFOMNL = toc
xtNL = xNL + ones(size(xNL,1),1)*v0NL.';
yNL = C*xtNL.';

discPoints = 0:L/k:(L-L/k);
intx = 15;
intt = 5;

figure;
surf(discPoints(1:intx:end),tSim(1:intt:end),xtNL(1:intt:end,1:intx:end));
title('NL-FOM')
xlabel('x (m)')
ylabel('Time (s)')
zlabel('v(x,y)')

figure;
plot(tSim,yNL);
hold on;
plot(tSim,yQB,'--');
xlabel('Time (s)');
ylabel('Output y');
legend('NL-FOM','QB-FOM');
