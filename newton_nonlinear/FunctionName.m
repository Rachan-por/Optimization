function [f,j] = FunctionName(x)
    f = [(x(1)^2+x(2)^2 - 1); (5*x(1)^2-x(2) - 2)]; % f(x)
    j = [2*x(1), 2*x(2); 10*x(1), -1]; % Jacobian of f(x)
end
