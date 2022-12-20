clear all ; close all;
epsilon = 1e-4;
IterationMax = 100;
x0 = [8;5];
[xsolution, Xk, Fk, Jk, IFLAG, IterationUsed] = newton(@FunctionName, x0, epsilon, IterationMax);
IterationK = (0:1:IterationUsed)';
x1 = Xk(1,:)';
x2 = Xk(2,:)';
f1 = Fk(1,:)';
f2 = Fk(2,:)';
T = table(IterationK, x1, x2, f1, f2)
