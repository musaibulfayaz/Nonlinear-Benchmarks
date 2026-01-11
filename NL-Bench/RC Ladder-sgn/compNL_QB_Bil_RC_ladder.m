clear all; close all; clc;

%% simulation parameter
k = 200;

dt = 1e-2;
tEnd = 2;
tSim = 0:dt:tEnd;

% excitation signals
% u = @(t) 0*(t>=0);
% u = @(t) 0*(t<4) + 3*(t>=4);
% u = @(t) 0*(t<0.3) + 1*(t>=0.3);
% u = @(t) 2*exp(-t);
% u = @(t) exp(-t);
% u = @(t) 30*exp(0.1*t);
u = @(t) t.^2;
% u = @(t) 10*sin(3*exp(0.1*t)*t);
% u = @(t) 10*sin(t.^3);
% u = @(t) 10*3*exp(0.1*t);

% u = @(t) (cos(2*pi*t*10)-6)*2;
 
% f1 = 10; f2 = 1000;
% u = @(t) sin(2*pi*f1*t)+sin(2*pi*f2*t);
% u = @(t) 5*sin(2*pi*f1*t);
% u = @(t) sin(2*pi*f1*t).*sin(2*pi*f2*t);

%% Simulation of pure nonlinear RC-ladder model
[f,B,C] = RC_ladder(k); % E=I
% xdot = RC_ladder(k); % E=I

B = zeros(k,1);
B(1) = 1;

C = zeros(1,k);
C(1) = 1;

nNL = k

xdotNL = @(t,x) f(x) + B*u(t);
% xdotNL = @(t,x) xdot(x,u(t));
optionsFOMNL = odeset('RelTol',1e-6,'AbsTol',1e-9,'Jacobian',@(t,x) RC_ladder_Jac(x,k));
% optionsFOMNL = odeset('RelTol',1e-6,'AbsTol',1e-9);

% x0NL = 0.1*ones(nNL,1);
x0NL = zeros(nNL,1);
tic
[~,xNL] = ode15s(xdotNL,tSim,x0NL,optionsFOMNL);
simTimeFOMNL = toc
yNL = C*xNL';

figure;
plot(tSim,yNL);
xlabel('Time (s)');
ylabel('Output y');

%% Simulation of the QB-version of RC-ladder
[E,A,H,N,~,B,C] = RC_Matrices(k); % E=I

nQB = size(A,1)

% Function and Jacobian of rhs of SISO QB-FOM (E=I)
xdotQB = @(t,x,A,H,N,B) A*x + H*(kron(x,x)) + (N*x + B)*u(t);
Jac_xdotQB = @(t,x,A,H,N,B) A + 2*H*(kron(x,speye(size(A,1)))) + N*u(t); 
% optionsFOMQB = odeset('RelTol',1e-6,'AbsTol',1e-9,'Jac',Jac_xdotQB);
optionsFOMQB = odeset('RelTol',1e-6,'AbsTol',1e-9);

x0QB = zeros(nQB,1);
tic
[~, xQB] = ode15s(xdotQB, tSim, x0QB, optionsFOMQB, A, H, N, B);
simTimeFOMQB = toc
yQB = C*xQB.';

hold on;
plot(tSim,yQB,'--');

%% Simulation of the bilinear version of RC-ladder
[A,N,B,C] = RC_Circuit(k); % E=I

nBil = size(A,1)

N = reshape(N,[nBil,nBil]);

% Function and Jacobian of rhs of SISO Bil-FOM (E=I)
xdotBil = @(t,x,A,N,B) A*x + (N*x + B)*u(t);
Jac_xdotBil = @(t,x,A,N,B) A + N*u(t); 
optionsFOMBil = odeset('RelTol',1e-6,'AbsTol',1e-9,'Jac',Jac_xdotBil);
% optionsFOMBil = odeset('RelTol',1e-6,'AbsTol',1e-9);

x0Bil = zeros(nBil,1);
tic
[~, xBil] = ode15s(xdotBil, tSim, x0Bil, optionsFOMBil, A, N, B);
simTimeFOMBil = toc
yBil = C*xBil.';

hold on;
plot(tSim,yBil,'--');
legend('NL-FOM','QB-FOM','Bil-FOM');

