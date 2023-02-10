clear all;
for i = 0.01:0.01:2;
    [x_min, f_min, IFLAG, IFunc, Ak, Bk, Xk1, Xk2] = golden(0,i,10^-4,100);
    y = 40*(x_min-1)^3-12*cos(3*x_min) ; % compute gradient of f at x = minimizer
    if abs(y) <= 10^(-5) ; % if gradient of f at x = minimizer converges to zero, then terminate this loop.
        break
    end
end
% Print the result
IterationK = (0:1:IFunc-2)';
Xk1 = Xk1';
Xk2 = Xk2';
Ak =  Ak';
Bk = Bk';
T = table(IterationK, Xk1, Xk2, Ak, Bk)
