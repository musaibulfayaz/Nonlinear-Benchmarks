function [f,B,C,fJac] = RC_ladder(N,x,symb)
%
% System: x_dot@(x,u)
%         y@(x)
%
% System: x_dot = f(x) + B*u;
%         y = C x
%
%         xdot = @(x,u) g(A1*x)-g(A2*x) + [u;sparse(N-1,1)];
%         y = @(x) x(1);
%
% Input: * N: order of original model

g = @(x) exp(40.0*x)+x-1;

A1 = full(spdiags(ones(N-1,1),-1,N,N)-speye(N));
A1(1,1)=0;
A2 = full(spdiags([ones(N-1,1);0],0,N,N)-spdiags(ones(N,1),1,N,N));
A0 = sparse(N,N);
A0(1,1)=1; 

%----------------------------- NEW ---------------------------------------%
f = @(x) -g(A0*x)+g(A1*x)-g(A2*x);

C = zeros(1,N);
C(1) = 1;
B = zeros(N,1);
B(1) = 1;

if nargout > 3
    fJac = RC_ladder_Jac(x,A1,A2,N,symb);
end
end

%%
%%------------------------ AUXILIARY FUNCTIONS --------------------------%%

function fJac = RC_ladder_Jac(x,A1,A2,N,symb)
%analytical calculation of jacobi matrix
e  = @(x) exp(40.0*x);
e1 = e(A1*x);
e2 = e(A2*x);

if symb==1
    fJac=sym(zeros(N,N));
    Jcurr=sym(zeros(N,N));
else
    fJac=sparse(N,N);
    Jcurr=sparse(N,N);
end
for i=1:1:N
    for j=1:1:N
        Jcurr(j,i)=A1(j,i)*e1(j,1)-A2(j,i)*e2(j,1);
    end
end
fJac = A1-A2 + 40*Jcurr;

end
