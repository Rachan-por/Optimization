function [x_min, f_min, IFLAG, IFunc, Ak, Bk, Xk1, Xk2] = golden(a, b, epsilon, itmax)

%==========================================================================
% To solve one dimensional optimization problem using golden section method
% a anb b are the interval that brackets x_min.
% a is the lower bound, and b is the upper bound.
% epsilon is the parameter used in termination criterion.
% itmax is the maximum number of the iteration allowed.
% x_min is the estimate of the minimizer of f(x).
% f_min is f(x_min).
% IFLAG = 0 if the search is successful.
% IFLAG = -999 if the search is fail.
% IFunc is the total number of function evaluations used.
% Ak, Bk, Xk1, Xk2 are the vector storing the value of ai, bi, x1i, x2i.
%==========================================================================


    Ak = [a]; % Generate Ak to store a_i and assign a0 = a
    Bk = [b]; % Generate Bk to store b_i and assign b0 = b
    t = (sqrt(5) - 1)/2; 
    Xk1 = [b-(b-a)*t]; % Generate Xk1 to store x1_i and assign x1_0 = b-(b-a)*t
    Xk2 = [a+(b-a)*t]; % Generate Xk2 to store x2_i and assign x2_0 = a+(b-a)*t
    IFunc = 0; % create an IFunc value
    IFLAG = -999; % the search isn't sucessful.
    fx1 = compute_f(Xk1(end));  
    fx2 = compute_f(Xk2(end));
    IFunc = IFunc + 2;
    for i = 1:1:itmax ;
        % use nested function below the main function 
        % to compute a function that we want to find minimizer
        %fx1 = compute_f(Xk1(end));  
        %fx2 = compute_f(Xk2(end));
        % Implement Golden section algorithm
        if fx2 > fx1 ;
            Bk(end+1) = Xk2(end);
            Ak(end+1) = Ak(end);
            Xk2(end+1) = Xk1(end);
            Xk1(end+1) = Ak(end) + (1-t)*(Bk(end)-Ak(end-1));
            fx2 = fx1;
            fx1 = compute_f(Xk1(end));
            
        else
            Bk(end+1) = Bk(end);
            Ak(end+1) = Xk1(end);
            Xk1(end+1) = Xk2(end);
            Xk2(end+1) = Bk(end-1) - (1-t)*(Bk(end-1)-Ak(end));
            fx1 = fx2;
            fx2 = compute_f(Xk2(end));
        end
        IFunc = IFunc + 1; % The function was evaluated 2 times.
        if abs(Xk1(end)-Xk2(end)) <= epsilon; % if abs(fx1-fx2) < epsilon, then terminate the loop
            x_min = (Xk1(end)+Xk2(end))/2; % estimate x_min by using the mean of x1 and x2
            f_min = compute_f(x_min); % compute f(x_min)
            IFLAG = 0; % the search is sucessful.
            break % teminate this loop
        end
    end
    function y = compute_f(Xk) % create a function that we wanted to find a minimizer 
        y = 10*(Xk-1)^4 - 4*sin(3*Xk);% compute f(x)
    end 
end






