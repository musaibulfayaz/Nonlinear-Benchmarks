clear all; close all; clc;

%% simulation parameter
k = 50;
nu = 0.05;
% nu = 0.01;
% nu = 0.05;

t0 = 0;
tEnd = 10;
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

% for i=1:k
%     T(:,i) = H(:,i+(i-1)*k);
% end

% A2 = spdiags([-101*ones(k,1);0],0,k,k);
% A2 = spdiags([-101*ones(k,1);0],0,k,k) + spdiags(50.50*ones(k,1),1,k,k);
% A2 = spdiags([-101*ones(k,1);0],0,k,k) + spdiags(50.50*ones(k,1),-1,k,k);
A2 = spdiags([-101*ones(k,1);0],0,k,k) + spdiags(50.50*ones(k,1),1,k,k) + spdiags(50.50*ones(k,1),-1,k,k);


xDotNL = @(t,x) A*x + T*x.^2 + (N*x + B)*u_train(t);
xDotNL = @(t,x) A*x + A2*x.^2 + (N*x + B)*u_train(t);
JacxDotNL = @(t,x) A + 2*T*x + N*u_train(t);

%% Simulation of the QB-version of Burger model
[E,A,H,N,B,C] = Burgers_Matrices_zero_ic(k,nu); % E=I

nQB = size(A,1)

% Function and Jacobian of rhs of SISO QB-FOM (E=I)
xdotQB = @(t,x,A,H,N,B) A*x + H*(kron(x,x)) + (N*x + B)*u(t);
Jac_xdot = @(t,x,A,H,N,B) A + 2*H*(kron(x,speye(size(A,1)))) + N*u(t); 
optionsFOM = odeset('RelTol',1e-6,'AbsTol',1e-9,'Jac',Jac_xdot);

x0QB = zeros(nQB,1);

tic
[~, xQB] = ode15s(xdotQB, tSim, x0QB, optionsFOM, A, H, N, B);
simTimeFOMQB = toc;
yQB = C*xQB.';

figure;
plot(tSim, yQB);
xlabel('Time (s)');
ylabel('Output y');

%% Simulation of the bilinear version of Burger model
[A,N,B,C] = Burgers_Equation(k,nu);  % E=I

nBil = size(A,1)

N = reshape(N,[nBil,nBil]);

% Function and Jacobian of rhs of SISO Bil-FOM (E=I)
xdotBil = @(t,x,A,N,B) A*x + (N*x + B)*u(t);
Jac_xdotBil = @(t,x,A,N,B) A + N*u(t); 
optionsFOMBil = odeset('RelTol',1e-6,'AbsTol',1e-9,'Jac',Jac_xdotBil);

x0Bil = zeros(nBil,1);
tic
[~, xBil] = ode15s(xdotBil, tSim, x0Bil, optionsFOMBil, A, N, B);
simTimeFOMBil = toc
yBil = C*xBil.';

hold on;
plot(tSim,yBil,'--');
legend('QB-FOM','Bil-FOM');
