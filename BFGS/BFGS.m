function [xmin, fmin, Xk, Fk, Gk, Lk, nF, nG, IFLAG] = BFGS(FcnName, x0, epsilon, mu, eta, itmax)

%==========================================================================
% BFGS function for solving the n-dimensional minimization
% problem using the BFGS method with the inexact line-search.
% FunctionName is the function we want to optimize.
% x0 is an initial point.
% epsilon is a parameter in stopping criterion.
% mu and eta are the parameter used in the line search criteria.
% itmax is the maximum iterations allowed.
% xmin is the estmate of the minimizer of f(x).
% fmin is f(x_min).
% Xk, Fk and Gk are the array containing xk, f(xk) and g(xk) respectively.
% Lk is an array containing lambda in each iteration.
% nF and nG are the numbers of evaluation of the function and gradient.
% IFLAG = 99 indicate this algorithm is fail to converge.
% IFLAG = 1 indicate this algorithm can converge.
%==========================================================================


    % Generate all variable that we will use
    nF = []; nG =[];nf=0;ng=0
    Lk = [];Fk = [];Gk = [];Xk = [x0];f=0;g=0;lambda=0;n1=0,n2=0
    B = [1 0; 0 1];
    IFLAG = 99; % IFLAG = 99 means this function cannot find the minimizer.
    [f,g] = FcnName(Xk(:,end), 2); nf = nf + 1; ng = ng + 1;
    Fk(end+1) = f; Gk(:,end+1) = g;
    % iteration = itmax times
    for i = 1:itmax ;
        s = inv(B) * (-g); % compute search dirction
        [lambda, n1, n2] = linesearch(Xk(:,end), s, f, g, mu, eta, FcnName) % using line search to find lambda_min
        Lk(end+1) = lambda; % store lambda in Lk
        nf = nf+n1; ng = ng+n2; % number of evaluation used
        x = Xk(:,end) + lambda*s; % compute xk+1
        Xk(:,end+1) = x; %store xk
        % compute function and gradient to store in Fk and Gk respectively
        [f,g] = FcnName(Xk(:,end), 2); nf = nf + 1; ng = ng + 1; 
        Fk(end+1) = f; Gk(:,end+1) = g;
        % compute Bk+1
        del_g = Gk(:,end)-Gk(:,end-1); del_x = Xk(:,end)-Xk(:,end-1);
        B = B + ((del_g)*(del_g)')/((del_g)'*del_x) - B*(del_x*(del_x)')*B / (del_x'*B*del_x)
        % store the number of evaluation used to nF and nG
        nF(end+1)=nf;nf=0;nG(end+1)=ng;ng=0
        % condition to terminate the iteration
        if norm(Gk(:,end)) < epsilon ;
            xmin = Xk(:,end);
            fmin = Fk(end);
            IFLAG = 1;
            break
        end
    end
end
% LineSearch
function [lambda, n_F, n_G] = linesearch(x, s, f, g, mu, eta, FcnName);
        % set lambda = 1, set alpha>1 and 0<beta<1
        lambda = 1; alpha = 1.5; beta = 0.5;f1=0;g1=0;n_F=0;n_G=0;
        [f1, g1] = FcnName(x+lambda*s,2); n_F = 1; n_G =1;
        % Armijo's condition
        % if lambda cannot satisfy this condition, set lambda = beta*lambda
        while f1 > f + mu*s'*g*lambda;
            lambda = lambda*beta
            f1 = FcnName(x+lambda*s,1); n_F = n_F+1
        end
        [f1, g1] = FcnName(x+lambda*s,2); n_F = n_F+1; n_G =n_G+1;

        % Strong Wolf's condition
        % If lambda cannot satisfy this condition
        % first check if s'*g1>0, then use golden method to find lambda_min
        % if s'*g1<0 set lambda = alpha*lambda 
        % and go back to check strong Wolf's condition
        if abs(s'*g1)>-eta*s'*g;
            while s'*g1 <0;
                lambda = lambda*alpha;
                [f1, g1] = FcnName(x+lambda*s,2); n_F = n_F+1; n_G =n_G+1;
            end
            [lambda,n_f] = golden(x, 0, lambda, s, FcnName,g,eta);
            n_F = n_F + n_f
        end
end 
    % Golden method 
    function [lambda, n_f] = golden(x, a, b, s, FcnName,g,eta);
        t = (sqrt(5) - 1)/2; 
        Xk1 = b-(b-a)*t; % Generate Xk1 to store x1_i and assign x1_0 = b-(b-a)*t
        Xk2 = a+(b-a)*t; % Generate Xk2 to store x2_i and assign x2_0 = a+(b-a)*t
        n_f = 0
        fx1 = FcnName(x+Xk1*s,1);n_f=n_f+1;
        fx2 = FcnName(x+Xk2*s,1);n_f=n_f+1;
        for i = 1:100;
            if fx2 > fx1 ;
                b = Xk2; a = a; Xk2 = Xk1; Xk1  =a+(1-t)*(b-a);
                fx2 = fx1; [fx1,g1] = FcnName(x+Xk1*s,2);n_f=n_f+1;
            elseif fx1 > fx2;
                a = Xk1; b = b; Xk1 = Xk2; Xk2=b-(1-t)*(b-a);
                fx1 = fx2; [fx2,g1] = FcnName(x+Xk2*s,2);n_f=n_f+1;
            end
            
            if abs(s'*g1) <= -eta*s'*g;
                lambda = (Xk1+Xk2)/2;
                break
            end
        end
    end




    
