function [f, gradient] = Rosenbrock(x, options);
    x1 = x(1);
    x2 = x(2);
    if options == 1;
        f = 100*(x2-x1^2)^2 + (1-x1)^2;
        gradient = 0;
    elseif options == 2;
        f = 100*(x2-x1^2)^2 + (1-x1)^2;
        gradient = [400*x1^3 - 400*x1*x2-2+2*x1; 200*(x2-x1^2)];
    end
end
