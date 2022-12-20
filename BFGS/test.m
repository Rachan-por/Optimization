clear all;
x0 = [10;100];
eta = 0.01; mu = 10^(-4); epsilon = 10^(-4); itmax = 100;
[xmin, fmin, Xk, Fk, Gk, Lk, nF, nG, IFLAG] = BFGS(@Rosenbrock, x0, epsilon, mu, eta, itmax);
IterationK = (1:1:length(Fk)-1)';
Xk1 = Xk(1,1:end-1)';
Xk2 = Xk(2,1:end-1)';
Fk =  Fk(1:end-1)';
nF = nF';
nG = nG';
T = table(IterationK, Xk1, Xk2, Fk, nF, nG)
