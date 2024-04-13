% -------------------------------------------------------------------------

% Matlab code for Problem Set 6 - Question 4
% Von Mises stress - Adjusted & Relaxed.

% -------------------------------------------------------------------------

% Assuming dx, dy, and rho are already in your workspace

length = 51;
X = zeros(length, length);
Y = zeros(length, length);

for i=1:length
    for j = 1:length
        X(i, j) = j;
        Y(i, j) = -i;
    end
end

nlength = 50;
nx = zeros(nlength, nlength);
ny = zeros(nlength, nlength);

for i=1:nlength
    for j = 1:nlength
        nx(i, j) = ( X(i, j) + X(i, j+1) )*0.5;
        ny(i, j) = ( Y(i, j) + Y(i+1, j) )*0.5;
    end
end

% Initialize von_mises arrays
von_mises = zeros(nlength,nlength);
adjusted_von_mises = zeros(nlength,nlength);
gamma = 0.2;
yield_strength = 325; % (MPa)
relaxed_von_mises = zeros(nlength,nlength);

E = 70*10^3; % Young's modulus (MPa)
nu = 0.3; % Poisson's ratio
penal = 3; % Penalisation factor

for i=1:nlength
    for j=1:nlength
        
        % element geometry
        x_e = [X(i+1, j) X(i+1, j+1) X(i, j+1) X(i, j)];
        y_e = [Y(i+1, j) Y(i+1, j+1) Y(i, j+1) Y(i, j)];
        d = [dx(i+1, j) dy(i+1, j) dx(i+1, j+1) dy(i+1, j+1) dx(i, j+1) dy(i, j+1) dx(i, j) dy(i, j)];
        area = 1; % area of element (mm^2)
        
        x = mean(x_e);
        y = mean(y_e);
        
        % matrix of elastic constants
        D = (E*((rho(i, j))^penal)/(1-nu^2)) * [1 nu 0; nu 1 0; 0 0 (1-nu)/2]; % (Pa)

        % ------------------------------------------------------------------------
        % Von Misses Stresses computation
        H = (1 / area) * [(y - y_e(4)), 0, -(y - y_e(4)), 0, (y - y_e(1)), 0 , -(y - y_e(1)), 0;
        0, (x - x_e(2)), 0, -(x - x_e(1)), 0, (x - x_e(1)), 0, -(x - x_e(2));
        (x - x_e(2)), (y - y_e(4)), -(x - x_e(1)), -(y - y_e(4)), (x - x_e(1)), (y - y_e(1)),-(x - x_e(2)), -(y - y_e(1))];
        
        % Stress & Strain from lecture notes
        strain = H*d';
        stress = D*strain;
        
        % Utilising the average of computed four nodal stresses to
        % calculate von mises stress at elemental center.
        
        von_mises(i, j) = sqrt(stress(1)^2 - (stress(1) * stress(2)) + stress(2)^2 + (3 * stress(3)^2));
    end
end

p = 1:500;
pnorm_values = zeros(size(p));

for i = 1:500
    pnorm_values(i) = (sum(sum((max(0, von_mises)).^p(i))))^(1/p(i));
end

figure;
plot(p, pnorm_values, 'b'); % Blue color for P-norm vs p plot
hold on;

hold off; % Release the hold on the plot
xlabel('p');
ylabel('p norm');
title('p-norm vs p');
