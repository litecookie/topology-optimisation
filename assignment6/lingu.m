% Load data from the MATLAB file
load('problem_set_06_data.mat');

% Define material properties (Young's modulus, yield strength, Poisson's ratio)
E = 70e9; % Young's modulus in Pa (70 GPa)
Sy = 325e6; % Yield strength in Pa (325 MPa)
nu = 0.3; % Poisson's ratio

% Define relaxation factor gamma
gamma = 0.2;

% Determine the coordinates of element centers (assuming uniform grid)
nelx = size(dx, 2) - 1; % number of elements in X direction
nely = size(dy, 1) - 1; % Number of elements in Y direction

X = linspace(0.5, nelx-0.5, nelx);
Y = linspace(0.5, nely-0.5, nely);

% Initialize arrays to store von Mises stresses at element centers
sigma_VM_adjusted = zeros(nely, nelx);
sigma_VM_relaxed = zeros(nely, nelx);

% Initialize von_mises arrays
von_mises = zeros(nelx, nely);

% Calculate stress components and von Mises stresses at element centers
for i = 1:nely % Skip boundary elements
    for j = 1:nelx % Skip boundary elements
        x_e = [i-1, i, i, i-1];
        y_e = [j-1, j-1, j, j];
        
        % Calculate strain components using displacements and strain-displacement relationship
        D = (E *(rho(i, j)^0.3) / (1 - nu^2)) * [1 nu 0; nu 1 0; 0 0 0.5 - 0.5 * nu];
        H = (1 / area) * [(Y(j) - y_e(4)), 0, -(Y(j) - y_e(4)), 0, (Y(j) - y_e(1)), 0, -(Y(j) - y_e(1)), 0;
        0, (X(i) - x_e(2)), 0, -(X(i) - x_e(1)), 0, (X(i) - x_e(1)), 0, -(X(i) - x_e(2));
        (X(i) - x_e(2)), (Y(j) - y_e(4)), -(X(i) - x_e(1)), -(Y(j) - y_e(4)), (X(i) - x_e(1)), (Y(j) - y_e(1)), -(X(i) - x_e(2)), -(Y(j) - y_e(1))];
        
        % Calculate stress components using strain components and material properties
        
    end
end

% Plot contours of adjusted von Mises stress
figure;
contourf(sigma_VM_adjusted);
colorbar;
title('Adjusted Von Mises Stress at Element Centers');
xlabel('Element X');
ylabel('Element Y');

% Plot contours of relaxed von Mises stress
figure;
contourf(sigma_VM_relaxed);
colorbar;
title('Relaxed Von Mises Stress at Element Centers');
xlabel('Element X');
ylabel('Element Y');