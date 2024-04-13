% Define the function y
y = @(x) sin(7*pi*x/2) + 30./(x+10) + exp(-x/4).*x.^3 + (3*x/10).^5;

% Define the integration limits
a = 0;
b = 10;

% Define the number of basis functions (line segments)
num_basis_functions = 2:100; 

% Initialize the minimum number of basis functions to a large value
min_basis_functions = 1000000;

% Compute the error for each number of basis functions
for i = 1:length(num_basis_functions)
    % Define the interval
    x_values = linspace(a, b, num_basis_functions(i)); % Equally spaced points
    
    % Compute the y values for the line segments
    y_values = y(x_values);
    
    % Create the approximate function psi with line segments
    psi = @(x) interp1(x_values, y_values, x, 'linear');
    
    % Compute the error
    integrand = @(x) (y(x) - psi(x)).^2;
    error = integral(integrand, a, b);
    
    % Check if the error is less than 1 and update the minimum basis functions
    if error < 1
        min_basis_functions = num_basis_functions(i);
        break; % Exit the loop if minimum basis functions found
    end
end

% Display the minimum number of basis functions and the corresponding error
disp(['Minimum Number of Basis Functions for Error Less Than 1: ' num2str(min_basis_functions)]);
disp(['Error for Minimum Number of Basis Functions: ' num2str(error)]);

% Define the interval with the minimum number of basis functions
x_values_min = linspace(a, b, min_basis_functions); % Equally spaced points

% Compute the y values for the line segments with the minimum basis functions
y_values_min = y(x_values_min);

% Create the approximate function psi with the minimum number of basis functions
psi_min = @(x) interp1(x_values_min, y_values_min, x, 'linear');

% Plot the original function y and the approximate function psi with the minimum number of basis functions
figure;
hold on;
fplot(y, [a, b], 'b', 'LineWidth', 1.5); % Original function y in blue
fplot(psi_min, [a, b], 'r--', 'LineWidth', 1.5); % Approximate function psi with minimum basis functions in red dashed line

% Initialize arrays to store intersection points
intersection_x = [];
intersection_y = [];

% Plot the connections for each basis function
for j = 1:length(x_values_min)-1
    % Calculate the coordinates of the current connection line
    x1 = [x_values_min(j+1), x_values_min(j)]; % x-values for the current line
    y1 = [0, y(x_values_min(j))]; % y-values for the current line
    
    % Plot the connection line
    plot(x1, y1, 'g'); 
    
    % Calculate the coordinates of the next connection line
    x2 = [x_values_min(j), x_values_min(j+1)]; % x-values for the next line
    y2 = [0, y(x_values_min(j+1))]; % y-values for the next line
    
    % Plot the next connection line
    plot(x2, y2, 'g'); 
    
    % Find the intersection point between the current and next connection lines
    [xi, yi] = polyxpoly(x1, y1, x2, y2);
    
    % Append the arrays
    intersection_x = [intersection_x, xi];
    intersection_y = [intersection_y, yi];
end

% Store the x and y coordinates of the intersection points in a new array called 'centres'
centres = [intersection_x; intersection_y]';

% Plot the intersection points
plot(intersection_x, intersection_y, 'ro', 'MarkerSize', 8); % Red circle marker for centres
ylabel('y');
xlabel('x');
title('y and psi with Minimum Basis Functions');
legend('y', 'psi', 'Line Segments', 'Location', 'best'); 
hold off;