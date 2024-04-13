% Define the original function
y_func = @(x) sin(3.5 * pi * x) + 30 ./ (x + 10) + (x.^3) .* exp(-x ./ 4) + (0.3 * x).^5;

% Interval [0, 10]
a = 0;
b = 10;

% Start with a guess for the number of elements
n_elements = 2;

% Initialize the error
err = inf;

% Error threshold
err_thresh = 1;

% Maximum number of iterations to prevent infinite loop
max_iter = 1000;
iter = 0;

while err > err_thresh && iter < max_iter
    iter = iter + 1;
    dx = (b - a) / n_elements;
    x_nodes = a:dx:b; % Nodes
    % Compute the y values for the line segments
    psi = y_func(x_nodes);
    
    % Calculate the error using numerical integration
    % The interp1 function creates a function for psi using the values
    % computed at the nodal points.
    err_func = @(x) (interp1(x_nodes, psi, x, 'linear', 'extrap') - y_func(x)).^2;
    err = integral(err_func, a, b);
    
    if err > err_thresh
        n_elements = n_elements + 1;
            % Increase the number of elements step-wise.
    end
end

% Plotting
x_fine = linspace(a, b, 1000);
y_exact = y_func(x_fine);
y_approx = interp1(x_nodes, psi, x_fine, 'linear');

figure;
plot(x_fine, y_exact, 'b-', 'LineWidth', 2);
hold on;
plot(x_fine, y_approx, 'r--', 'LineWidth', 2);
legend('Exact Function', 'Approximated Function', 'Location', 'Best');
title('Function Approximation using Linear Finite Elements');
xlabel('x');
ylabel('y');
grid on;

n_basis = 0;
% Output the final number of elements and the error
fprintf('Final number of basis functions: %d\n', n_elements);
fprintf('Final error: %f\n', err);

numCols = 3;
numRows = ceil(length(psi') / numCols);
paddedArray = [psi'; nan(numRows*numCols - length(psi'), 1)]; % Pad with NaNs
psir = reshape(paddedArray, [numRows, numCols]);
tableArray = array2table(psir);
disp(tableArray);
