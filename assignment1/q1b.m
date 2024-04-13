% Define the objective function
f = @(n, r) n.^2 - 4*n + r.^2 - 6*r + 18;

% Defining curve limits.
max_n = 10;
min_n = -10;
max_r = 10;
min_r = -10;

% Define the range for n and r variables
nRange = linspace(min_n, max_n, 20);
rRange = linspace(min_r, max_r, 20);

% Create a meshgrid for n and r
[n, r] = meshgrid(nRange, rRange);

% Evaluate the objective function at each point in the meshgrid
z = f(n, r);


% ------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Defining the functions for performing differentiation operation
syms n2 r2 lambda;

fx = n2^2 - 4*n2 + r2^2 - 6*r2 + 18;

% Calculate the partial derivatives with respect to n and r
df_dn = diff(fx, n2);
df_dr = diff(fx, r2);

mini_n = solve(df_dn == 0, n2);
mini_r = solve(df_dr == 0, r2);

count = 0;

% Verify if the variables solve the constraint equations
if mini_r > 4
    count = count+1;
end

if mini_n > 4/mini_r
    count = count + 1;
end

% Calculations shown in notes.
calculated_n = 1.47;
calculated_r = 2.72;

if count > 0
    % Shown in calculation as well.
    disp("Constraint where r - 4/n <=0 are inactive")
    fminima = subs(fx, [n2, r2], [calculated_n, calculated_r]);
    disp("Minimum value of f");
    disp(fminima);
    disp("At n and r values:");
    disp(calculated_n);
    disp(calculated_r);
end
% ------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Plot contours
figure;
contour(n, r, z, 20);           % 20 contour lines are plotted
hold on;

% Plotting the constraints
fimplicit(@(n, r) r - 4, 'r-', 'LineWidth', 2);  % Curve r - 4 = 0
fimplicit(@(n, r) n - 4./r, 'b-', 'LineWidth', 2);  % Curve n = 10/r
plot(1.47, 2.72, 'ro', 'MarkerSize', 10);

xlabel('n');
ylabel('r');
title('Contour Plot of f(n, r) = n^2 - 4n + r^2 - 6r + 18');
colorbar;  % Display colorbar

hold off;