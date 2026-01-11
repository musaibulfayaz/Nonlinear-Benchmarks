clear
clc
close all

%% Prepare Model
%Load Model
sys_name = 'E100';

%% Load nonlinear model
[M, E, K, F, f_2nd, B, C, G, x0, v0] = ElectrostatBeam100;
n_2nd = size(M,1);

%% Convert model from second order to first order
E_1st = [speye(n_2nd,n_2nd) sparse(n_2nd,n_2nd); sparse(n_2nd,n_2nd) M];
A_1st = [sparse(n_2nd,n_2nd) speye(n_2nd,n_2nd); -K -E];
f_1st = @(x) [sparse(n_2nd,1); f_2nd(x)];
B_1st = [sparse(n_2nd,1); B];

n_1st = size(A_1st,1);

%% Simulation 
% Simulation parameters
dt = 0.1; %Step width
tEnd = 10;
tSim = 0:dt:tEnd;

%Excitation signal
u = @(t) 1*(t>=0);

% Initial state vector
x0_2nd = x0;
x0_1st = [x0; v0];

fTest2 = f_2nd(x0_2nd);
fTest1 = f_1st(x0_1st);

% Simulation of the nonlinear system
xdot = @(t,x) A_1st*x + f_1st(x) + B_1st*u(t);
           
options = odeset('Mass',E_1st,'RelTol',1e-6);
             
[t_nonlin,x_nonlin] = ode15s(xdot,tSim,x0_1st,options);
y_nonlin = C*x_nonlin.';


