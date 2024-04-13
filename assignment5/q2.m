% Bisectioning algorithm - 
a = -10^6; % Start of interval
b = 10^6; % End of interval, adjust based on your knowledge of the function
tol = 10^-6; % Tolerance for the root's accuracy
maxIter = 1000; % Maximum number of iterations

root = bisectionMethod(a, b, tol, maxIter);
fprintf('Root found at: %f\n', root);

function root = bisectionMethod(a, b, tol, maxIter)
    % Define the function whose root we are trying to find
    f = @(x) ((exp(x / 200) * (sin((pi * x) / 500) + 2)) - 25);
    
    % Check if the initial interval is valid
    if f(a) * f(b) >= 0
        error('f(a) and f(b) must have different signs');
    end
    
    % Initialize the number of iterations
    iter = 0;
    
    % Main bisection algorithm
    while (b - a) / 2 > tol
        % Increment iteration count
        iter = iter + 1;
        
        % Prevent infinite loop
        if iter > maxIter
            error('Maximum iterations exceeded');
        end
        
        % Find midpoint
        c = (a + b) / 2;
        
        % Check if we have found the root or if the root is in the left or right half
        if f(c) == 0
            a = c;
            break; % Exact root found
        elseif f(a) * f(c) < 0
            b = c; % Root is in left half
        else
            a = c; % Root is in right half
        end
    end
    
    % Return the approximate root
    root = (a + b) / 2;
end