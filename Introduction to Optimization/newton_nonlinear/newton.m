function [xsolution, Xk, Fk, Jk, IFLAG, IterationUsed] = newton(FunctionName, x0, epsilon, IterationMax)

%===============================================================
% Function newton is used to solve non-linear equation.
% FunctionName is the non-linear function.
% x0 is an initial point.
% epsilon is used in termination criterion.
% IterationMax is a number of permissible iterations.
% xsolution is the numberical solution obtained from newton
% Xk is an array for the iterations.
% Fk is an array of F(x)
% Jk is an array of Jacobian matrix.
% IFLAG = 1 indicates the algorithm fail to converge.
% IFLAG = 0 indicates the algorithm can converge.
% IterationUsed is a number of iteration used in compution the solution.
%===============================================================

    Xk = [x0];
    Fk = [];
    Jk = [];
    IFLAG = 0;
    for i = 1:IterationMax; 
        [f, j] = FunctionName(Xk(:,end)); % Find f(x) and jacobian of f(x)
        Fk = [Fk f]; % store Fk
        Jk = [Jk j]; % store Jk
        s = -(j^-1)*f; % compute s
        xk = Xk(:,end) + s; % update equation
        Xk = [Xk, xk]; % store xk
        if norm(Xk(:,i+1)- Xk(:,i)) < epsilon ; %terminate criterion
            break % terminate if this condition is true
        end
    end
    xsolution = Xk(:,end); % get solution of x
    [f, j] = FunctionName(Xk(:,end));
    Fk = [Fk f];
    Jk = [Jk j];
    IterationUsed = i;
    if IterationUsed == IterationMax ; % if this algorithm fail to converge return IFLAG = 1
        IFLAG = 1; 
    end
end