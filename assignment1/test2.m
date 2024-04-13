% Define the system of equations as anonymous functions
equations = @(x) [2*x(1) - 4 + x(3)*(4/x(1)^2); 2*x(2) - 6 + x(3); x(2) - 4/x(1)];

% Initial guess
initial_guess = [1, 1, 1];

% Solve numerically
numerical_solution = fsolve(equations, initial_guess);

disp('Numerical Solution:');
disp(['n = ', num2str(numerical_solution(1))]);
disp(['r = ', num2str(numerical_solution(2))]);
disp(['lambda = ', num2str(numerical_solution(3))]);
