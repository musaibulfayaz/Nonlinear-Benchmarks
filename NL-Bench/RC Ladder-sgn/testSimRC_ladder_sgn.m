%% initialize and prepare
clear
clc
close all

% RC-ladder initialisation
N = 500; 
[f,B,C] = RC_ladder_sgn(N);

n = size(B,1);

% x = randn(n,1);
% tic
% fJacLoop = RC_ladder_Jac(x,N,1);
% toc
% 
% tic
% fJacNoLoop = RC_ladder_Jac(x,N,0);
% toc

%% simulation parameter
dt = 1e-2;
tEnd = 10;
tSim = 0:dt:tEnd;

% excitation signals
% u = @(t) 0*(t>=0);
% u = @(t) 0*(t<4) + 3*(t>=4);
% u = @(t) 0*(t<0.3) + 1*(t>=0.3);
% u = @(t) 2*exp(-t);
u = @(t) exp(-t);
% u = @(t) 30*exp(0.1*t);
% u = @(t) t.^2;
% u = @(t) 10*sin(3*exp(0.1*t)*t);
% u = @(t) 10*sin(t.^3);
% u = @(t) 10*3*exp(0.1*t);

% u = @(t) (cos(2*pi*t*10)-6)*2;
 
% f1 = 10; f2 = 1000;
% u = @(t) sin(2*pi*f1*t)+sin(2*pi*f2*t);
% u = @(t) 5*sin(2*pi*f1*t);
% u = @(t) sin(2*pi*f1*t).*sin(2*pi*f2*t);

%% Simulation of nonlinear model
% x0 = zeros(n,1);
x0 = 0.1*ones(n,1);

xdotNonlin = @(t,x) f(x) + B*u(t);

% optionsJac = odeset('RelTol',1e-6,'AbsTol',1e-9,'Jacobian',@(t,x) RC_ladder_Jac(x,N));
% optionsNoJac = odeset('RelTol',1e-6,'AbsTol',1e-9);

optionsJac = odeset('Jacobian',@(t,x) RC_ladder_sgn_Jac(x,N));
optionsNoJac = [];


%-- ODE15s (ODE23s) for comparison
tic
[~,xJac] = ode15s(xdotNonlin,tSim,x0,optionsJac);
yJac = C*xJac';
simTime_yJac = toc

tic
[~,xNoJac] = ode15s(xdotNonlin,tSim,x0,optionsNoJac);
yNoJac = C*xNoJac';
simTime_yNoJac = toc

figure;
plot(tSim, yJac);
hold on;
plot(tSim, yNoJac, 'r--');
xlabel('Time (s)');
ylabel('Output y');

%-- Implicit Euler
tic
[~,xiEJac] = implicitEuler(xdotNonlin,tSim,x0,optionsJac);
yiEJac = C*xiEJac';
simTime_xiEJac = toc

tic
[~,xiENoJac] = implicitEuler(xdotNonlin,tSim,x0,optionsNoJac);
yiENoJac = C*xiENoJac';
simTime_xiENoJac = toc

plot(tSim, yiEJac, 'b--');
plot(tSim, yiENoJac, 'm--');
legend('ode15s Jac','ode15s NoJac','iEuler Jac','iEuler NoJac');
