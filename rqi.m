function [x, its, errsigma, errx] = rqi(A,b,typec,precw,precc,precrqi)
% Multiprecision Rayleigh Quotient Iteration to solve Total Least Squares  
% problem with PCGSTL with R-factor of A from the Householder QR factorization as preconditioner. The algorithm uses one  
% step of inverse iteration in the beginning. The algorithm uses up to
% three different precision. Uniform precision with fp64 corresponds to the
% RQI-PCGTLS algorithm.
% For error calculation, quadruple precision is used. 
%
% Inputs
% A        : input matrix
% b        : rhs
% typec    : type of cholesky factorization (1 - cholesky with two-sided diagonal scaling, 2 - householder qr)
% precw    : working precision for PCGTLS (1 - half, 2 - single, 3 - double). Default is 3.
% precc    : cholesky precision in PCGTLS (1 - half, 2 - single, 3 - double). Default is 3.
% precrqi  : working precision for RQI (1 - half, 2 - single, 3 - double). Default is 3.
%
% Outputs
% x        : solution vector
% its      : number of iterations performed
% errsigma : error in sigma
% errx     : error in x

if precrqi == 3
    fprintf('**** RQI precision is double.\n')
    urqi = 'd';
    A = double(A);
    b = double(b);   
elseif precrqi == 2
    fprintf('**** RQI precision is single.\n')
    urqi = 's';
    A = single(A);
    b = single(b);
elseif precrqi == 1
    fprintf('**** RQI precision is half.\n')
    urqi = 'h';fp.format = uws;
    chop([],fp);
    A = chop(A);
    b = chop(b);
else
     error('Error: precrqi should be 1, 2, or 3.\n')  
end

Ab = [A,b];
[~,sigma1,v] = svds(double(Ab),1,'smallest');
zeta = v(end);
xtrue = v(1:end-1)./(-zeta);

% one step of inverse iteration for convergence
xls = A\b;
x = xls;
r = b-A*x;
sigma = r'*r/(1+x'*x);

% Cholesky/Householder QR factorization for the precondonditioner
[~,n] = size(A);

if precc == 3
    fprintf('**** Cholesky precision is double.\n')
    ucc = 'd';
elseif precc == 2
    fprintf('**** Cholesky precision is single.\n')
    ucc = 's';
elseif precc == 1
    fprintf('**** Cholesky precision is half.\n')
    ucc = 'h';fp.format = ucc;
    chop([],fp);
else
     error('Error: precc should be 1, 2, or 3.\n')    
end

if typec == 1
    scale.flag = 1; % 2 sided diag scaling (if 1)
    scale.cluf = 0; % custom lu fact (if 1)
    scale.pert = 1; %diagonal perturbation for low prec Chol (if 1)
    scale.theta = 0.1; % a number between (0,1]
    [uh,~,~,xmax] = float_params(ucc);
    [Ah,R] = spd_diag_scale(chop(chop(A)'*chop(A)));
    C = R;
    c = scale.pert;
    Ah = Ah+(c*uh*eye(n));
    mu = (scale.theta)*xmax;
    Ah = mu*Ah;
    U = chol_lp(Ah,ucc); 
    R = (1/sqrt(mu))*double(U)*diag(1./diag(C));

elseif typec == 2
    [~,R] = house_qr_lp(A,0,ucc);

else
     error('Error: typec should be 1 or 2.\n')      

end

v = R'\x;
u = R\v;
x = x+sigma*u;

% for termination criteria inside the loop 
rho = norm(r)^2/(norm(x)^2+1);
fk = -A'*r - rho*x;
gk = -b'*r+rho;
gamma = sqrt((norm(fk)^2 + gk^2)/(norm(x)^2+1));

its=0;
flag = 1;

errsigma = [mp(double(mp(double(abs(sqrt(mp(double(sigma),34))-mp(double(sigma1),34))),34)/mp(double(abs(sigma1)),34)),34)];
errx = [mp(double(mp(double(norm(mp(double(x),34)-mp(double(xtrue),34))),34)/mp(double(norm(mp(double(xtrue),34))),34)),34)];

while flag == 1 
    its=its+1;
    r = b-A*x;
    sigma = (r'*r)/(1+x'*x);
    f = -A'*r- sigma*x;
    g = -b'*r+sigma;

    [omega,~] = pcgtls(-f,R,sigma,its+1,precw);

    switch urqi
        case 'd'
            omega = double(omega);
        case 's'
            omega = single(omega);
        case 'h'
            omega = chop(omega);
    end
    
    z = x+omega;
    beta = (z'*f-g)/(z'*x+1);
    
    [u,~] = pcgtls(x,R,sigma,its+1,precw);

    switch urqi
        case 'd'
            u = double(u);
        case 's'
            u = single(u);
        case 'h'
            u = chop(u);
    end

    x = z+beta*u;

    % error computation
    errsigma = [errsigma;mp(double(mp(double(abs(sqrt(mp(double(sigma),34))-mp(double(sigma1),34))),34)/mp(double(abs(sigma1)),34)),34)];
    errx = [errx;mp(double(mp(double(norm(mp(double(x),34)-mp(double(xtrue),34))),34)/mp(double(norm(mp(double(xtrue),34))),34)),34)];

    % termination criteria
    gamma_old = gamma;
    gamma = sqrt((norm(f)^2 + g^2)/(norm(x)^2+1));
    if gamma>gamma_old
        flag = 0;
    end
end

% plot errors
figure
semilogy(1:numel(errx),errx,'--')
hold on
semilogy(1:numel(errsigma),errsigma)
hold off

xlim([1 numel(errx)])
atm = get(gca,'xticklabels');
xlab = [];
m = str2double(atm);
num = 1;
for i = 1:numel(m)
    if ceil(m(i)) == m(i)
        xlab(num) = m(i);
        num = num + 1;
    end
end
set(gca,'xticklabels',xlab);
set(gca,'xtick',xlab);
xlabel({'Rayleigh-Quotient iteration'},'Interpreter','latex');

set(gca,'FontSize',14)
a = get(gca,'Children');
set(a,'LineWidth',1);
set(a,'MarkerSize',10);
h = legend('Error of x','Error of sigma');
set(h,'Interpreter','latex');
if precw==precc && precc==precrqi
    title('RQI-PCGTLS','Interpreter','latex')
else
    title('RQI-PCGTLS-MP','Interpreter','latex')
end

end