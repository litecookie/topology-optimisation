% -------------------------------------------------------------------------

% Matlab code for Problem Set 6 - Question 3
% Strain energy and gradient of compliance.

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
strain_energy = zeros(nlength,nlength);
adjusted_von_mises = zeros(nlength,nlength);
gamma = 0.2;
relaxed_von_mises = zeros(nlength,nlength);
gradient = zeros(nlength,nlength);

E = 70*10^9; % Young's modulus (Pa)
nu = 0.3; % Poisson's ratio
penal = 3; % Penalisation factor

total_strain_energy = 0;
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
        % Calculate total strain energy for every element with deformation.
        volume = 1; % Assuming that the deformation is minimal.
        strain_energy(i, j) = (0.5).*strain'*stress;

        total_strain_energy = total_strain_energy + strain_energy(i, j);

        % calculate the local stiffness matrix from current Gauss point
        D_star = D / E;
        K_star = E * H' * D_star * H;

        gradient(i, j) = -penal * (rho(i, j)^(penal-1)) * d * K_star * d';
    end
end

disp(total_strain_energy);

% Plotting the contour of rho
figure; % Opens a new figure window
contourLevels = linspace(min(gradient(:)), max(gradient(:)), 400); % Adjust 20 to increase or decrease the number of levels
contourf(nx, ny, gradient, contourLevels); % Filled contour plot
colorbar; % Shows the color scale
title('Gradient of Compliance Contours');
xlabel('X Coordinate');
ylabel('Y Coordinate');