%% Bjorck matrix test
close all
clear all

rng(1);
m = 30;
n = 15;
tol = 0.05;

Y = RandOrthMat(m);
Z = RandOrthMat(n);

for k = 1:n
	d(k) = 2^(-k+1);
end
D = diag(d);
D0 = [D;zeros(m-n,n)];

A_tilde = Y*D0*Z';

for k = 1:n
	x(k) = 1/k;
end

x = x';

b_tilde = A_tilde*x;

% random perturbations
r = tol*rand(m,1);
A = A_tilde;
b = b_tilde+r;
A = A+ 0.05*rand(m,n); 

precw = 2; precc = 1; typec = 2;  precrqi = 3;
[x1, its1, errsigma1, errx1] = rqi(A,b,typec,3,3,3);
figure
[x, its, errsigma, errx] = rqi(A,b,typec,precw,precc,precrqi);

% plot errors
figure
semilogy(1:numel(errx1),errx1,'b--')
hold on
semilogy(1:numel(errsigma1),errsigma1,'b')
hold on
semilogy(1:numel(errx),errx,'r--')
hold on
semilogy(1:numel(errsigma),errsigma,'r')
hold off

% Ensure only integers labeled on x axis
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
set(a,'LineWidth',1.5);
set(a,'MarkerSize',10);
h = legend('rerrx','rerrs','rerrx\_mp','rerrs\_mp');
set(h,'Interpreter','latex');
title('($u,u_p,u_c$) = (double,single,half)','Interpreter','latex')

%% Random matrix test
close all
clear all

rng(0);
m = 100;
n = 60;
tol = 1e-6;

A_tilde = rand(m,n);
b_tilde = ones(m,1);

% random perturbations
E = tol*rand(m,n);
r = tol*rand(m,1);
A = A_tilde + E;
b = b_tilde+r;

precw = 2; precc = 1; typec = 2;  precrqi = 3;
[x1, its1, errsigma1, errx1] = rqi(A,b,typec,3,3,3);
figure
[x, its, errsigma, errx] = rqi(A,b,typec,precw,precc,precrqi);

% plot errors
figure
semilogy(1:numel(errx1),errx1,'b--')
hold on
semilogy(1:numel(errsigma1),errsigma1,'b')
hold on
semilogy(1:numel(errx),errx,'r--')
hold on
semilogy(1:numel(errsigma),errsigma,'r')
hold off

% Ensure only integers labeled on x axis
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
set(a,'LineWidth',1.5);
set(a,'MarkerSize',10);
h = legend('rerrx','rerrs','rerrx\_mp','rerrs\_mp');
set(h,'Interpreter','latex');
title('($u,u_p,u_c$) = (double,single,half)','Interpreter','latex')

%% Delta matrix test
close all
clear all

rng(1);
m = 9;
n = 4;
tol = 1e-1;
sigma = 1e-2;

A_tilde = zeros(m,n);
A_tilde(1,1)=sigma;
A_tilde(3,2) = sigma;
A_tilde(7,3)= 1;
A_tilde(9,4)= 1;
b_tilde = ones(m,1);

% random perturbations
E1 = -1+2*rand(m,n);
r1 = -1+2*rand(m,1);
E = tol*E1.*A_tilde;
r = tol*r1.*b_tilde;
A = A_tilde + E;
b = b_tilde+r;

precw = 2; precc = 1; typec = 2;  precrqi = 3;
[x1, its1, errsigma1, errx1] = rqi(A,b,typec,3,3,3);
figure
[x, its, errsigma, errx] = rqi(A,b,typec,precw,precc,precrqi);

% plot errors
figure
semilogy(1:numel(errx1),errx1,'b--')
hold on
semilogy(1:numel(errsigma1),errsigma1,'b')
hold on
semilogy(1:numel(errx),errx,'r--')
hold on
semilogy(1:numel(errsigma),errsigma,'r')
hold off

% Ensure only integers labeled on x axis
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
set(a,'LineWidth',1.5);
set(a,'MarkerSize',10);
h = legend('rerrx','rerrs','rerrx\_mp','rerrs\_mp');
set(h,'Interpreter','latex');
title('($u,u_p,u_c$) = (double,single,half)','Interpreter','latex')

%% Vanhuffel matrix test
clear all
close all

rng(1);
m = 100;
tol = 1e-6;

A_tilde = -1.*ones(m,m-2);
for i = 1:m-2
    A_tilde(i,i) = m-1;
end
b_tilde = -1.*ones(m,1);
b_tilde(m-1,1) = m-1;

% random perturbations
E = tol*rand(m,m-2);
r = tol*rand(m,1);
A = A_tilde + E;
b = b_tilde+r;

precw = 2; precc = 1; typec = 2;  precrqi = 3;

[x1, its1, errsigma1, errx1] = rqi(A,b,typec,3,3,3);
figure
[x, its, errsigma, errx] = rqi(A,b,typec,precw,precc,precrqi);

% plot errors
figure
semilogy(1:numel(errx1),errx1,'b--')
hold on
semilogy(1:numel(errsigma1),errsigma1,'b')
hold on
semilogy(1:numel(errx),errx,'r--')
hold on
semilogy(1:numel(errsigma),errsigma,'r')
hold off

% Ensure only integers labeled on x axis
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
set(a,'LineWidth',1.5);
set(a,'MarkerSize',10);
h = legend('rerrx','rerrs','rerrx\_mp','rerrs\_mp');
set(h,'Interpreter','latex');
title('($u,u_p,u_c$) = (double,single,half)','Interpreter','latex')

%% Toeplitz matrix test
clear all
close all

rng(1);
m = 100;
omega = 2;
gamma = 0.001;
beta = 1.25;
tol = 1e-6;

t1 = zeros(m,1); 
ee = exp(omega^2/(2*beta^2));
t1(1)= 1/(sqrt(2*pi*beta^2)*ee);
A_tilde = toeplitz(t1);
A_tilde = A_tilde(:,1:m-2*omega);
b_tilde = ones(m,1);

% random perturbations
E = randn(m,m-2*omega);
E = E*gamma*norm(A_tilde)/norm(E);
r = randn(m,1);
r = r*gamma*norm(b_tilde)/norm(r);
A = A_tilde + E;
b = b_tilde+r;

precw = 2; precc = 1; typec = 2;  precrqi = 3;
[x1, its1, errsigma1, errx1] = rqi(A,b,typec,3,3,3);
figure
[x, its, errsigma, errx] = rqi(A,b,typec,precw,precc,precrqi);

% plot errors
figure
semilogy(1:numel(errx1),errx1,'b--')
hold on
semilogy(1:numel(errsigma1),errsigma1,'b')
hold on
semilogy(1:numel(errx),errx,'r--')
hold on
semilogy(1:numel(errsigma),errsigma,'r')
hold off

% Ensure only integers labeled on x axis
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
set(a,'LineWidth',1.5);
set(a,'MarkerSize',10);
h = legend('rerrx','rerrs','rerrx\_mp','rerrs\_mp');
set(h,'Interpreter','latex');
title('($u,u_p,u_c$) = (double,single,half)','Interpreter','latex')