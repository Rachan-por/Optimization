function [x_min, f_min, Xk, Fk, Gk, nF, nG, nH, IFLAG] = Newton(FunctionName, x0, epsilon, ep_rel, ep_abs, itmax);

%==========================================================================
% This Newton function solve n-dimensional optimization problem using
% modified Newton's method, where the line-search subprogram is solved
% using golden section method.
% FunctionName is the function we want to optimize.
% x0 is an initial point.
% epsilon is a parameter in stopping criterion.
% ep_real and ep_abs are the parameter used in stopping criterion of line search.
% itmax is the maximum iterations allowed.
% x_min is the estmate of the minimizer of f(x).
% f_min is f(x_min).
% Xk, Fk and Gk are the array containing xk, f(xk) and g(xk) respectively.
% nF, nG and nH are the numbers of evaluation of the objective function,
% gradient and Hessian matrix respectively.
% IFLAG = -999 indicates the algorithm fail to converge.
% IFLAG = 1 indicates the algorithm can converge.
%==========================================================================

    Xk = [x0];
    [f, g, h] = FunctionName(x0, 3);
    Fk = [f];
    Gk = [g];
    x_min = 0; f_min = 0; IFLAG = -999;
    nF = 1; nG = 1; nH = 1;
    for i = 1:itmax;
        s = -inv(h)*Gk(:,end);
        if s'*Gk(:,end) > 0;
            s = -Gk(:,end);
        end
        lambda = 0;
        b = 5;
        a = 0;
        t = (sqrt(5) - 1)/2; 
        x1 = b-(b-a)*t;
        x2 = a+(b-a)*t;
        f1 = FunctionName(Xk(:,end) + s*x1, 1); f2 = FunctionName(Xk(:,end) + s*x2,1);
        while lambda == 0;
            if f1 < f2;
                b = x2; a=a; x2 = x1; x1 = a + (1-t)*(b-a);
                f2 = f1;
                f1 = FunctionName(Xk(:,end) + s*x1, 1);
            else ;
                a = x1; b =b; x1 = x2; x2 = b - (1-t)*(b-a);
                f1 = f2;
                f2 = FunctionName(Xk(:,end) + s*x2,1);
            end
            x = Xk(:,end) + (x1+x2)/2 * s;
            [f, g] = FunctionName(x, 2);
            if abs(s'*g) <= abs(s'*Gk(:,end))*ep_rel + ep_abs;
                lambda = ( a + b ) / 2;
                break
            end
        end
        Xk =[Xk x];
        [f, g, h] = FunctionName(x, 3);
        nF = nF + 1; nG = nG +1; nH = nH +1;
        Fk = [Fk f];
        Gk = [Gk g];
        if abs(Fk(end-1)-Fk(end)) <= epsilon ;
            x_min = Xk(:,end);
            f_min = Fk(end);
            IFLAG = 1;
            break
        end
        end
end




