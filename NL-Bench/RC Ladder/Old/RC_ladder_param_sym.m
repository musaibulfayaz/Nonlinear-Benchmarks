function [f,B,C] = RC_ladder_param_sym(N)
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

A1 = full(spdiags(ones(N-1,1),-1,N,N)-speye(N));
A1(1,1)=0;
A2 = full(spdiags([ones(N-1,1);0],0,N,N)-spdiags(ones(N,1),1,N,N));
A0 = sparse(N,N);
A0(1,1)=1; 
% xdot = @(x,u) g(A1*x)-g(A2*x) + [u;sparse(N-1,1)];
% y = @(x) x(1);

%----------------------------- NEW ----------------------------------------
f_temp = @(param,x) -g(param,A0*x)+ g(param,A1*x)-g(param,A2*x);
syms f(param,x) param x

x = symbolic_array('x',N);
f = sym(f_temp(param,x));

C = zeros(1,N);
C(1) = 1;
B = zeros(N,1);
B(1) = 1;

end

