% Define the function
f = @(n, r) n.^2 - 4.*n + r.^2 - 6.*r + 18;

% Create a grid of values for n and r
n_values = linspace(-10, 10, 100); % Adjust the range as needed
r_values = linspace(-10, 10, 100); % Adjust the range as needed
[n, r] = meshgrid(n_values, r_values);

% Evaluate the function at each point in the grid
z = f(n, r);

% Plot the surface
figure;
surf(n, r, z);
xlabel('n');
ylabel('r');
zlabel('f(n, r)');
title('Plot of f(n, r) = n^2 - 4n + r^2 - 6r + 18');
