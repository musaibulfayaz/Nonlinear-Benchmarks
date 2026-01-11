%% initialize and prepare
clear
clc
close all

% nonlinear model
n = 10; % original order of RC-ladder model 
%random vector
xVal=randn(n,1);
%symbolic array
name='x';
x = symbolic_array(name,n);
x = transpose(x);

symb=1;
[f,B,C,Ja] = RC_ladder(n,x,symb);
symb=0;
[~,~,~,Ja_Val]=RC_ladder(n,xVal,symb);
E = speye(n);

%% Validation of Jacobian
% symbolic Jacobian by Matlab's jacobian()
J  = SymJacobian(f,x);

% Spy
figure; spy(J); title('Jacobian calculated by Matlabs jacobian');
figure; spy(Ja_Val); title('analytical jacobian');
% Difference
JVal=subs(J,x,xVal);
Ja_subs=subs(Ja,x,xVal);
Diff=double(JVal-Ja_subs);  %Ja_subs: analytical jacobian calculated symbolic
Diff_Val=double(JVal-Ja_Val); %Ja_Val: analytical jacobian calculated with values
%Diff=double(JVal-Ja);

% simulation parameter
dt = 1e-2;
tEnd = 20;
tSim = 0:dt:tEnd;

% excitation signals
% u = @(t) 0*(t>=0);
% u = @(t) 0*(t<4) + 3*(t>=4);
% u = @(t) 0*(t<0.3) + 1*(t>=0.3);
u = @(t) 2*exp(-t);
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
x0 = zeros(n,1);
% x0 = 0.1*ones(n,1);

tic
xdotNonlin = @(t,x) f(x) + B*u(t);
[~,x_nonlin] = ode15s(xdotNonlin,tSim,x0,struct('RelTol',1e-6));
y_nonlin = C*x_nonlin.';
simTimey_nonlin = toc

figure;
plot(tSim, y_nonlin);
xlabel('Time (s)');
ylabel('Output y');
