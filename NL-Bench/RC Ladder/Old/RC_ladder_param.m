function [f_param,B,C] = RC_ladder_param(N)
%
% System: x_dot@(x,u)
%         y@(x)
%
% System: x_dot = f(x) + B*u;
%         y = C x
%
% Input: * N: order of original model
%          diodeparam: Parameter of the diode (Compare Shockley Equation: i(u) = i_s(exp(u*diodeparam)-1)) 

g = @(param,x) exp(param*x)+x-1;

A1 = spdiags(ones(N-1,1),-1,N,N)-speye(N);
A2 = spdiags([ones(N-1,1);0],0,N,N)-spdiags(ones(N,1),1,N,N);

% xdot = @(x,u) g(A1*x)-g(A2*x) + [u;sparse(N-1,1)];
% y = @(x) x(1);

%----------------------------- NEW ----------------------------------------
f_param = @(param,x) g(param,A1*x)-g(param,A2*x);

C = zeros(1,N);
C(1) = 1;
B = zeros(N,1);
B(1) = 1;

end

