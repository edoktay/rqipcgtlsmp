function [omega,its] = pcgtls(f,R,sigma,maxit,precw)
% Preconditioned Conjugate Gradient to solve Total Least Squares problem
% (eye(n)-c^2inv(R')inv(R))omega = f in a user-defined precision. 
%
% Inputs
% f     : rhs
% R     : nxn upper triangular matrix
% sigma : squared constant c^2
% maxit : maximum number of iterations for the loop
% precw : working precision (1 - half, 2 - single, 3 - double, 4 - quad). Default value is 3.
%
% Outputs
% omega : solution vector
% its   : number of iterations performed

n = size(f,1);

if nargin > 5
    switch precw 
        case 3
            fprintf('**** Working precision is double.\n')
            f = double(f);   
            omega = double(zeros(n,1));
            p = double(zeros(n,1));
            s = double(zeros(n,1));
            q = double(zeros(n,1));
            R = double(R);
        case 2
            fprintf('**** Working precision is single.\n')
            f = single(f);
            omega = single(zeros(n,1));
            p = single(zeros(n,1));
            s = single(zeros(n,1));
            q = single(zeros(n,1));
            R = single(R);
        case 1
            fprintf('**** Working precision is half.\n')
            fp.format = 'h';
            chop([],fp);
            f = chop(f);
            omega = chop(zeros(n,1));
            p = chop(zeros(n,1));
            s = chop(zeros(n,1));
            q = chop(zeros(n,1));
            R = chop(R);
        case 4
            fprintf('**** Working precision is quad.\n')
            mp.Digits(34);
            f = mp(double(f),34);
            omega = mp(double(zeros(n,1)),34);
            p = mp(double(zeros(n,1)),34);
            s = mp(double(zeros(n,1)),34);
            q = mp(double(zeros(n,1)),34);
            R = mp(double(R),34);
    end
else
    fprintf('**** Working precision is double.\n')
    f = double(f);   
    omega = double(zeros(n,1));
    p = double(zeros(n,1));
    s = double(zeros(n,1));
    q = double(zeros(n,1));
    R = double(R);
end

p = R'\f;
s = p;
eta = norm(s)^2;

its = 0;
delta = 1;

while its<=maxit && delta ~= 0 
    its = its+1;
    q = R\p;
    delta = norm(p)^2-sigma*norm(q)^2;

    if delta == 0 
        break
    end

    alpha = eta/delta;
    omega = omega + alpha*q;
    q = R'\q;
    s = s-alpha*(p-single(sigma)*q);
    eta_old = eta;
    eta = norm(s)^2;
    beta = eta/eta_old;
    p = s + beta*p;
  
end

end