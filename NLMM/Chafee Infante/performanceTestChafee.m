clear all
close all
clc 

%% Load benchmark
% Chafee Infante
% xdot = A*x + B*u(t) - fcubic(x)

k = 50000000; % huge model
% k = 5000; % big model
% L = 5;
L = 1;

% [~,A,B,C] = ChafeeInfante_POD_Matrices(k,L); %last too long with the
% current implementation

%% High Performance Test and Comparison of the different implementations
x0 = [1:1:k]';

tic
for i=1:100
    fcub1 = fcubic(x0);
end
toc

tic
for i=1:100
    fcub2 = fcubic2(x0);
end
toc

tic
f = fcubic3(k);
for i=1:100
    fcub3 = f(x0);
end
toc