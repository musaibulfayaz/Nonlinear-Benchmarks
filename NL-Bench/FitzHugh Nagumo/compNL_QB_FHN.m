clear all; close all; clc;

%% simulation parameter
k = 50;
l=1;

t0 = 0;
tEnd = 5;
dt = 0.01;
tSim = t0:dt:tEnd;

% excitation signals
i0 = @(t) 5*10^4*t.^3.*exp(-15*t);
% i0 = @(t) 50*(sin(2*pi*t)-2);
u2 = @(t) 1;
u = @(t) [i0(t); u2(t)];

%% Simulation of pure nonlinear FHN-model
[E,A,B,C] = Fitz_Matrices_POD(k,l); % E~=I !!

% correct procedure!
A = E\A; %% Since E is diagonal, so these operations are very cheap and efficient. 
B = E\B; 
E = E\E; % E=I

Ed = E(1:2*k,1:2*k);
Ad = A(1:2*k,1:2*k);
Bd = B(1:2*k,:);
Cd = C(:,1:2*k);

% wrong procedure!
% [E,A,B,C] = Fitz_Matrices_POD(k,l); % E~=I !!
% Ed2 = E(1:2*k,1:2*k);
% Ad2 = Ed2\A(1:2*k,1:2*k);
% Bd2 = Ed2\B(1:2*k,:);
% Cd2 = C(:,1:2*k);

nNL = size(Ad,1)

xdotNL = @(t,x,A,B) A*x + Fitz_NL(x) +  B*u(t);
options = odeset('RelTol',1e-10,'AbsTol',1e-12);
% options = odeset('RelTol',1e-10,'AbsTol',1e-12,'Mass',Ed);

x0NL = zeros(2*k,1);

tic
[~,xNL] = ode15s(xdotNL,tSim,x0NL,options,Ad,Bd);
simTimeFOMNL = toc
yNL = Cd*xNL';

figure;
plot(tSim,yNL);

%% Simulation of the QB-version of FHN-model
[E,A,H,N,B,C] = FitzHughNagumo_Matrices(k);
N = [N,zeros(size(N))];

nQB = size(A,1)

% Function and Jacobian of rhs of MIMO QB-FOM (E~=I)
xdotQB = @(t,x,A,H,N,B) A*x + H*(kron(x,x)) + N*kron(u(t),x) + B*u(t);
Jac_xdot = @(t,x,A,H,N,B) A + 2*H*(kron(x,speye(size(A,1)))) + N*kron(u(t),speye(size(A,1))); 
optionsFOM = odeset('RelTol',1e-10,'AbsTol',1e-12,'Jac',Jac_xdot,'Mass',E);

x0QB = zeros(3*k,1);

tic
[~, xQB] = ode15s(xdotQB, tSim, x0QB, optionsFOM, A, H, N, B);
simTimeFOMQB = toc
yQB = C*xQB.';

hold on;
plot(tSim,yQB,'--');
xlabel('Time (s)');
ylabel('Outputs y');
legend('NL-FOM y1', 'NL-FOM y2', 'QB-FOM y1', 'QB-FOM y2');

