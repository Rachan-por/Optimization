function [xmin, fmin, Xk, Fk, Gk, Lk, nF, nG, IFLAG, nReset] = CG(x0, epsilon, mu, eta, itmax, option);

%==========================================================================
% CG is a function for solving the n-dimensional minimization problem
% using conjugate gradient methods with an inexact line search.
% x0 is an initial point.
% epsilon is a parameter in stopping criterion.
% mu and eta are the parameter used in the line search criteria.
% itmax is the maximum iterations allowed.
% if option = 1, then Fletcher-Reeves formula is used and if option = 2, then
% Polak-Ribi´ere formula is used.
% xmin is the estmate of the minimizer of f(x).
% fmin is f(x_min).
% Xk, Fk and Gk are the array containing xk, f(xk) and g(xk) respectively.
% Lk is an array containing lambda in each iteration.
% nF and nG are the numbers of evaluation of the function and gradient.
% IFLAG = 99 indicate this algorithm is fail to converge.
% IFLAG = 1 indicate this algorithm can converge.
% nReset is an integer vector recording the history of resetting sk in each iteration 
% where each element is set to be 0 when restart is not used, 1 when restart is used 
% because the angle between sk and −gk is too large, and 2 when sk
% does not have descent property
%==========================================================================


    % define variable that will be used
    Xk = [x0]; Fk = []; Gk = []; Lk = [];
    nF = [1]; nG = [1]; 
    IFLAG = 0; nReset = []; %IFLAG = 0 means this algorithm fails to converge
    [f, g] = Rosenbrock(x0, 2); Fk(end+1) = f; Gk(:,end+1) = g;
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    s = -g; %at k = 0; s = -gradient of f(x0)
    beta = 0.5; alpha = 1.5;xmin = [0;0];fmin =0; % define parameter beta, alpha
    for i = 1:itmax;
        % Line search
        nfunc = 0; ngrad = 0;
        lambda = 1
        f = Rosenbrock(Xk(:,end) + lambda*s,1);nfunc = nfunc+1;
        while f > Fk(end) + mu*lambda*s'Gk(:,end)  % Check Armijo's condition
            lambda = beta*lambda; %Update lambda
            f = Rosenbrock(Xk(:,end) + lambda*s, 1);nfunc = nfunc+1; 
        end
        [f,g] = Rosenbrock(Xk(:,end) + lambda*s, 2);
        nfunc = nfunc+1; ngrad = ngrad +1;
        while abs(s'*g) > -eta*s'*Gk(:,end); % Check strong's Wolf condition
            if s'*g>0; %check for using golden method
                [lambda, n_f, n_g] = golden(Xk(:,end),lambda,s,Gk(:,end),eta);
                nfunc = nfunc + n_f; ngrad = ngrad + n_g;
                break
            else
                lambda = alpha*lambda;% update lambda
                [f,g] = Rosenbrock(Xk(:,end) + lambda*s, 2);
                nfunc = nfunc+1; ngrad = ngrad +1;
            end
        end 
        % end of line search
        
        % update equation
        Xk(:,end+1) = Xk(:,end) + lambda*s;
        %%%%%%%%%%%%%%%%%%
        Lk(end+1) = lambda;
        [f,g] = Rosenbrock(Xk(:,end), 2) ; nfunc = nfunc +1; ngrad = ngrad+1;
        Fk(end+1) = f; % f(k+1)
        Gk(:, end+1) = g; % g(k+1)
        % termination condition
        if abs(Fk(:,end)) < epsilon ; % for Rosenbrock function is lowest at 0.
            xmin = Xk(:, end-1)
            fmin = Fk(:, end-1)
            IFLAG = 1 % IFLAG = 1 means this algorithm can converge
            nF(end+1) = nfunc; nG(end+1) = ngrad;
            break
        end
        if option == 1; % Fletcher-Reeves formula
            b = Gk(:,end)'*Gk(:,end) / (Gk(:,end-1)'*Gk(:,end-1));
        elseif option == 2;% Polak-Ribiere formula
            b = Gk(:,end)'*(Gk(:,end)-Gk(:,end-1)) / (Gk(:,end-1)'*Gk(:,end-1));
        end
        % watch-dog scheme
        if acos(s'*(-Gk(:,end-1)) / (norm(s)*norm(Gk(:,end-1)))) > 85*pi/180 ; % angle is too large
            nReset(end+1) = 1;
            s = -Gk(:,end); % update search direction
        elseif s'*Gk(:,end-1) >= 0 ; % sk doesn't have descent property.
            nReset(end+1) = 2;
            s = -Gk(:,end); % update search direction
        else 
            nReset(end+1) = 0;
            s = -Gk(:,end) + b*s % update search direction
        end
        nF(end+1) = nfunc; nG(end+1) = ngrad;
    end
end
    % Golden method 
    function [lambda, n_f, n_g] = golden(x, b, s,g,eta);
        t = (sqrt(5) - 1)/2; a= 0; 
        lambda = b;
        Xk1 = b-(b-a)*t; 
        Xk2 = a+(b-a)*t;
        n_f = 0;n_g=0;
        fx1 = Rosenbrock(x+Xk1*s,1);n_f=n_f+1;
        fx2 = Rosenbrock(x+Xk2*s,1);n_f=n_f+1;
        for i = 1:100;
            if fx2 >= fx1 ;
                b = Xk2; a = a; Xk2 = Xk1; Xk1  =a+(1-t)*(b-a);
                fx2 = fx1; 
                [fx1,g1] = Rosenbrock(x+Xk1*s,2);n_f=n_f+1;n_g=n_g+1;
            elseif fx1 > fx2;
                a = Xk1; b = b; Xk1 = Xk2; Xk2=b-(1-t)*(b-a);
                fx1 = fx2; 
                [fx2,g1] = Rosenbrock(x+Xk2*s,2);n_f=n_f+1;n_g=n_g+1;
            end
            
            if abs(s'*g1) <= -eta*s'*g;
                lambda = (Xk1+Xk2)/2;
                break
            end
        end
    end