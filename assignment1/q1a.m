% Objective function
f = @(n, r) n.^2 - 4*n + r.^2 - 6*r + 18;

% Define curve limits.
max_n = 10;
min_n = -10;
max_r = 10;
min_r = -10;

% Meshgrid (n , r) and corresponding equations
nRange = linspace(min_n, max_n, 20);
rRange = linspace(min_r, max_r, 20);
[n, r] = meshgrid(nRange, rRange);
z = f(n, r);

% ------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Defining the functions for calculating minima.
syms n2 r2;

fx = n2^2 - 4*n2 + r2^2 - 6*r2 + 18;

df_dn = diff(fx, n2);
df_dr = diff(fx, r2);

% Calculate the minimum of n and r by equation the first derivative to zero
mini_n = solve(df_dn == 0, n2);
mini_r = solve(df_dr == 0, r2);

% Verify if the variables solve the constraint equations
count = 0;
if mini_r > 4
    count = count+1;
end
if mini_n > 10/mini_r
    count = count + 1;
end

if count == 0
    % Since count is zero, both the constraints are inactive.
    disp("Both constraints are inactive")
    fminima = subs(fx, [n2, r2], [mini_n, mini_r]);
    
    disp("Minimum value of f");
    disp(fminima);
    disp("At n and r values:");
    disp(mini_n);
    disp(mini_r);
end
% ------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Plot contours
figure;
contour(n, r, z, 20);  
hold on;

% Plot the curves r - 4 = 0 and n = 10/r
fimplicit(@(n, r) r - 4, 'r-', 'LineWidth', 2);  % Curve r - 4 = 0
fimplicit(@(n, r) n - 10./r, 'b-', 'LineWidth', 2);  % Curve n = 10/r
plot(mini_n, mini_r, 'ro', 'MarkerSize', 10);

xlabel('n');
ylabel('r');
title('Contour Plot of f(n, r) = n^2 - 4n + r^2 - 6r + 18');
colorbar;

hold off;