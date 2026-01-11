function [ J ] = Jacobi( A1,A2,e,N,x,symb)
%analytical calculation of jacobi matrix

if symb==1
    J=sym(zeros(N,N));
    Jcurr=sym(zeros(N,N));
    e1=diag(e(A1*x));
    e2=diag(e(A2*x));
else
    J=sparse(N,N);
    Jcurr=sparse(N,N);
    e1=spdiags(e(A1*x),0,N,N);
    e2=spdiags(e(A2*x),0,N,N);
end
    

J_p=transpose(A1)*e1-transpose(A2)*e2;

J=A1-A2+40*J_p;

end

