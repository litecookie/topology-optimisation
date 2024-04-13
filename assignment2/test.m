% Define the symbolic variable
syms x;

% Define the function
f = exp(3*x + 1);

% Define the value(s) you want to substitute into the function
% For example, if you want to substitute x = 2
x_value = 1.5;

% Substitute the value into the function
f_substituted = subs(f, x, x_value);

% Display the result
disp(double(f_substituted))