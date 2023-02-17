%% Unconstrainted problem
clear all; close all; clc;
% problem 1: f(x) = (1/2)x'Px - q'x 
n = 10; rng('default');
% Generate a positive definite P
P = randn(n); P = (P+P')/2; P = P+ (abs(min(eig(P)))+0.01)*eye(n); q = randn(n,1);
% define a function
fun = @(x) 0.5*x'*P*x - q'*x ;
x0 = randn(n,1);
options = optimoptions(@fminunc,'Display','iter','PlotFcn','optimplotfval');
[x,fval] = fminunc(fun,x0,options)
z = P\q; % analytical solution
disp("Numerical and analytical solutions")
[x z]

%%
% problem 2: f(x) = sum(xlog(x))
clear all; close all;
n = 5;
fun = @(x) sum((x.*log(x)))
x0 = abs(randn(n,1));
options = optimoptions(@fminunc,'Display','iter','PlotFcn','optimplotfval');
[x,fval] = fminunc(fun,x0,options)

%%
% problem 3: f(x) = x1^2 + x1x2 + 1.5x2^2 -2log(x1) - log(x2)
clear all; close all;
x0 = [1;2];
fun = @(x) (x(1)^2 + x(1).*x(2) + 1.5*x(2)^2 - 2*log(x(1)) - log(x(2)))
options = optimoptions(@fminunc,'Display','iter','PlotFcn','optimplotfval');
[x,fval] = fminunc(fun,x0,options)

%%
% problem 4: f(x) = x1^2 - x1x2 + 2x2^2 - 2x1 + exp(x1+x2) 
clear all; close all;
fun = @(x) x(1)^2 - x(1)*x(2) + 2*x(2)^2 - 2*x(1) + exp(x(1) + x(2));
x0 = [5;10];
options = optimoptions(@fminunc,'Display','iter','PlotFcn','optimplotfval');
[x,fval] = fminunc(fun,x0,options)

%%
% problem 5: soft max loss in logistic regression
clear all; close all;
% generate y(-1,1) and x randomly
rng(1)
N = 200; n = 2;
y = (randi([-1 1],N,1)==1)*2-1;
x = [ones(N,1) randn(N,n)];
x0 = randn(n+1,1);
% define loss function
fun = @(m) 1 / N * sum(log(1 + exp(-y.*(x*m))))
options = optimoptions(@fminunc,'Display','iter','PlotFcn','optimplotfval');
[m,fval] = fminunc(fun,x0,options)
% plot data point and seperating line
ind = find(y == 1); indc = find(y == -1);
plot(x(ind,2),x(ind,3),'ob'); hold on; plot(x(indc,2),x(indc,3),'sr');
m1 = -5:0.001:5;
m2 = -(m(1)/m(3))*m1 - m(2)/m(3);
plot(m1,m2)
%%
% problem 5.2: create dataset that can be seperated by a line
clear all; close all;
% create dataset
N = 200; n = 2;
y = [ones(N/2,1); -ones(N/2,1)];
a = 0;
b = 1;
x1 = (b-a).*rand(100,1) + a;
c = 1;
d = 2;
x2 = (d-c).*rand(100,1) + c;
x = [x1;x2]
x = [ones(200,1) x randn(N,1)];
x0 = randn(n+1,1);
% define loss and fit
fun = @(m) 1 / N * sum(log(1 + exp(-y.*(x*m))))
options = optimoptions(@fminunc,'Display','iter','PlotFcn','optimplotfval');
[m,fval] = fminunc(fun,x0,options)
% plot
ind = find(y == 1); indc = find(y == -1);
plot(x(ind,2),x(ind,3),'ob'); hold on; plot(x(indc,2),x(indc,3),'sr');
m1 = -0:0.01:2;
m2 = -(m(1)/m(3)) - m(2)/m(3)*m1;
plot(m1,m2); hold on
%%
% problem 6: f(x) is rosenbrock function
clear all; close all;
fun = @(x) 100*(x(2) - x(1)^2)^2 + (1-x(1))^2;
x0 = [5;10];
options = optimoptions(@fminunc,'Display','iter','PlotFcn','optimplotfval');
[x,fval] = fminunc(fun,x0,options)

%% Non linear
% problem 1: 
clc; clear all; close all;
% generate dataset
rng('default') ;
a = 2; b = 1; c = 2; d = -4 ;
N = 200; x = linspace(-4,8,N) ; x = x';
fx = a*exp( -((x-b)/c).^2 ) + d ; 
y = fx+ 0.3*randn(N,1); 
plot(x,y,'o','linewidth',2); grid on; hold on;
% fit loss
fun = @(m1) norm(y - (m1(1) * exp(-(x - m1(2)).^2 / m1(3)^2) + m1(4)))^2
x0 = [1;1;1;1];
options = optimoptions('lsqnonlin','Display','iter');
[m1,resnorm,residual,exitflag,output] = lsqnonlin(fun,x0,[],[],options);
fx_hat = m1(1)*exp( -((x-m1(2))/m1(3)).^2 ) + m1(4) ; 
plot(x,fx_hat,'o','linewidth',2); grid on; hold on;

%%
% problem 2: sigmoid function
clc; clear all; close all;
rng('default') ;
k = 2;
b = 1;
N = 200; x = linspace(-4,8,N) ; x = x';
fx = k ./ (1 + exp(-b'*x)) ; 
y = fx+ 0.03*randn(N,1); 
plot(x,y,'o','linewidth',2); grid on; hold on;
fun = @(m1) (y - m1(1) ./ (1 + exp(-m1(2)'*x)))
x0 = [2;1];
options = optimoptions('lsqnonlin','Display','iter');
[m1,resnorm,residual,exitflag,output] = lsqnonlin(fun,x0,[],[],options);
fx_hat = m1(1) ./ (1 + exp(-m1(2)'*x)); 
plot(x,fx_hat,'o','linewidth',2); grid on; hold on;
%% LP
% problem 1: min c'x 
%            s.t. sum(x) <= 1, x >= 0
clear all; close all;
n = 5;
f = randn(n,1)
A = [ones(1,n); -eye(n)];
b = [1; zeros(n,1)];
x = linprog(f,A,b)
[f x]
%% problem 2: min c'x
%             s.t. l <= x <= u
clear all; close all;
n = 5;
rng(1);
f = randn(n,1);
A = [eye(n); -eye(n)];
b = [2*ones(n,1); -ones(n,1)];
x = linprog(f,A,b)
[f x f'*x*ones(n,1)]
%% ex3 problem 3: min c'x
%                 s.t. norm(x)_infinity <= 1
clear all; close all;
n = 5;
rng(1);
f = randn(n,1);
A = [eye(n)];
b = [ones(n,1)];
lb = -ones(n,1);
ub = ones(n,1);
x = linprog(f,[],[],[],[],lb,ub)
[f x f'*x*ones(n,1)]
%% ex4 problem 4: min c'x
%                 s.t. 0 <= x_1 <= x_2 <= ... <= x_n <= 1
clear all; close all;
rng(38);
n = 5;
f = randn(n,1);
f = [ 1.78862847,  0.43650985,  0.09649747, -1.8634927 , -0.2773882 ]; f = f'
b = [zeros(n,1);1];
A = -eye(n) + [zeros(1,n); eye(n-1) zeros(n-1,1)]; 
A = [A;zeros(1,n-1),1];
x = linprog(f,A,b);
[f x] 
%% problem 8: infinity norm
% min norm(Ax-y)_infinity
m = 4; n = 3;
rng(1)
f = [zeros(n,1);1];
a = randn(m,n); y = randn(m,1); y = [-0.3224172 , -0.38405435,  1.13376944, -1.09989127]';
a = [ 1.62434536, -0.61175641, -0.52817175;
        -1.07296862,  0.86540763, -2.3015387 ;
         1.74481176, -0.7612069 ,  0.3190391 ;
        -0.24937038,  1.46210794, -2.06014071]
A = [a, -ones(m,1); -a, -ones(m,1)];
b = [y; -y];
x = linprog(f, A, b);
a
[f x]
%% Quadratic Program
% problem 1: min (1/2)x'Px-q'x
% P is positive definite matrix
clear all; 
n = 2;
p = randn(n,1);
P = p'*p;
P = P + min(abs(eig(P))+0.01)*eye(n);
q = randn(2,1);
x = quadprog(P, -q)
analyatic = inv(P)*q
%% problem 2: min norm(Ax - y)_2
clear all;
rng(1)
m = 10; n = 3;
A = randn(m,n); y = randn(m,1);
P = 2*A'*A; q = -2*y'*A;
sol_un = quadprog(P,q)
analytic_sol = (inv(A'*A))*A'*y
%% problem 2.1: add constraint norm(x)_1 <= alpha
alpha = 0.1
h = alpha*[ones(n,1); ones(n,1)];
G = [eye(n); -eye(n)];
sol_con1 = quadprog(P,q,G,h)
%% problem 2.2: add constraint l <= x <= u
l = -0.1*ones(n,1);
u = 1*ones(n,1);
h = [u; -l];
G = [eye(n); -eye(n)]
sol_con2 = quadprog(P,q,G,h)
%% problem 3: Soft margin SVM
clear all ; load data-svm
N =100;
n = 2;
P = [eye(n), zeros(n,N+1); zeros(N+1,n), zeros(N+1,N+1)];
q = [zeros(n+1,1); 4*ones(N,1)];
h = [-ones(N,1); zeros(N,1)];
G = [-y.*x, -y, -eye(N); zeros(N,n+1), -eye(N)];
sol = quadprog(P, q, G,h);
ind = find(y == 1); indc = find(y == -1);
plot(x(ind,1),x(ind,2),'ob'); hold on; plot(x(indc,1),x(indc,2),'sr');
x1 = -1:0.1:3.5;
x2 = -x1*sol(1)/sol(2) - sol(3) / sol(2);
plot(x1,x2); hold on