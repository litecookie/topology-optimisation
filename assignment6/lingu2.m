% Load data from the MATLAB file
load('problem_set_06_data.mat');

% Define material properties (Young's modulus, yield strength, Poisson's ratio)
E = 70; % Young's modulus in Pa (70 GPa)
Sy = 325e6; % Yield strength in Pa (325 MPa)
nu = 0.3; % Poisson's ratio

% Define relaxation factor gamma and penal value
gamma = 0.2;
penal = 3;

% Determine the coordinates of element centers (assuming uniform grid)
nelx = size(dx, 2) - 1; % number of elements in X direction
nely = size(dy, 1) - 1; % Number of elements in Y direction

X = linspace(0.5, nelx-0.5, nelx);
Y = linspace(0.5, nely-0.5, nely);

% Initialize arrays to store von Mises stresses and their adjusted/relaxed versions
von_mises = zeros(nely, nelx);
adjusted_von_mises = zeros(nely, nelx);
relaxed_von_mises = zeros(nely, nelx);

% Calculate stress components and von Mises stresses at element centers
for i = 1:nely % Skip boundary elements
    for j = 1:nelx % Skip boundary elements
        x_e = [i-1, i, i, i-1];
        y_e = 1 .* [j-1, j-1, j, j];
        
        % Define displacement vector from precomputed tables
        d = [dx(i, j), dy(i, j), dx(i, j+1), dy(i, j+1), dx(i+1, j+1), dy(i+1, j+1), dx(i+1, j), dy(i+1, j)];

        % area of the element
        area = 1;

        % Calculate strain components using displacements and strain-displacement relationship
        D = (E * (rho(i, j)^penal) / (1 - nu^2)) * [1 nu 0; nu 1 0; 0 0 0.5 - 0.5 * nu];
        H = (1 / area) * [(Y(j) - y_e(4)), 0, -(Y(j) - y_e(4)), 0, (Y(j) - y_e(1)), 0, -(Y(j) - y_e(1)), 0;
            0, (X(i) - x_e(2)), 0, -(X(i) - x_e(1)), 0, (X(i) - x_e(1)), 0, -(X(i) - x_e(2));
            (X(i) - x_e(2)), (Y(j) - y_e(4)), -(X(i) - x_e(1)), -(Y(j) - y_e(4)), (X(i) - x_e(1)), (Y(j) - y_e(1)), -(X(i) - x_e(2)), -(Y(j) - y_e(1))];

        % Calculate strain components using displacements and strain-displacement relationship
        strain = H * d';
        stress = D * strain;

        % Calculate von Mises stress
        von_mises(i, j) = sqrt(stress(1)^2 - (stress(1) * stress(2)) + stress(2)^2 + (3 * stress(3)^2));
        
        % Calculate adjusted von Mises stress
        adjusted_von_mises(i, j) = von_mises(i, j) / ((rho(i, j))^penal);
        
        % Calculate relaxed von Mises stress
        relaxed_von_mises(i, j) = von_mises(i, j) / (((rho(i, j))^penal) * (1 - gamma + (gamma / rho(i, j))));
    end
end

% Plot contours of von Mises stress
figure(1);
contourf(von_mises);
colorbar;
title('Von Mises Stress at Element Centers');
xlabel('Element X');
ylabel('Element Y');

% Plot contours of von Mises stress
figure(2);
contourf(adjusted_von_mises);
colorbar;
title('Adjusted Von Mises Stress');
xlabel('Element X');
ylabel('Element Y');

% Plot contours of von Mises stress
figure(3);
contourf(relaxed_von_mises);
colorbar;
title('Relaxed Von Mises Stress');
xlabel('Element X');
ylabel('Element Y');
