function [f, gradient, Hessian] = FunctionName(x, options);
    if options == 1;
        f = 100*(x(2)-x(1)^2)^2 + (1-x(1))^2;
    elseif options == 2;
        f = 100*(x(2)-x(1)^2)^2 + (1-x(1))^2;
        gradient = [400*x(1)^3 - 400*x(1)*x(2)-2+2*x(1); 200*(x(2)-x(1)^2)];
    elseif options == 3;
        f = 100*(x(2)-x(1)^2)^2 + (1-x(1))^2;
        gradient = [400*x(1)^3 - 400*x(1)*x(2)-2+2*x(1); 200*(x(2)-x(1)^2)];
        Hessian = [1200*x(1)^2-400*x(2)+2, -400*x(1); -400*x(1), 200];
    end
end


       