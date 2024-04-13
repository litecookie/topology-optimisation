% -------------------------------------------------------------------------

% Matlab code for Problem Set 3 - Question 1

% -------------------------------------------------------------------------

% clear data space
clear;
% close figure windows
close all;
% element geometry
x_e = [0 2 2 0]; y_e = [-1 -1 2 2]; % nodal locations from bottom left node first (m)
area = 6; % area of element (m^2)

% initialize a few variables
K = zeros(8); % total stiffness matrix
K_star = zeros(8); % current term of stiffness matrix
f = zeros(8,1); % nodal force vector
d = zeros(8,1); % nodal displacement vector

t = 0.05; % thickness of element - not used for 2D case
E = 70e9; % Young's modulus (Pa)
nu = 0.3; % Poisson’s ratio
% matrix of elastic constants
D = (E/(1-nu^2)) * [1 nu 0; nu 1 0; 0 0 (1-nu)/2]; % (Pa)

% set up Gauss quadrature
xi = zeros(2);
eta = zeros(2);
xi(1) = - 1 / sqrt(3); xi(2) = 1 / sqrt(3); % Gauss points for x-direction
eta(1) = - 1 / sqrt(3); eta(2) = 1 / sqrt(3); % Gauss points for y-direction
w(1) = 1; w(2) = 1; % weights for Gauss quadrature
J_det = 1.5; % determinant of Jacobian of transformation to (xi, eta)
x = [1+xi(1) 1+xi(2)]; y = [(1+3*eta(1))/2 (1+3*eta(2))/2]; % physical locations of Gauss points (m)

for i = 1:2
    for j = 1:2
        % value of H at current Gauss point
        H = (1 / area) * [(y(j) - y_e(4)), 0, -(y(j) - y_e(4)), 0 (y(j) - y_e(1)), 0 , -(y(j) - y_e(1)), 0; 0, (x(i) - x_e(2)), 0, -(x(i) - x_e(1)), 0, (x(i) - x_e(1)), 0, -(x(i) - x_e(2));(x(i) - x_e(2)), (y(j) - y_e(4)), -(x(i) - x_e(1)), -(y(j) - y_e(4)), (x(i) - x_e(1)), (y(j) - y_e(1)), (x(i) - x_e(2)), -(y(j) - y_e(1))];
        % contribution to stiffness matrix from current Gauss point
        K_star = w(i) * w(j) * J_det * H' * D * H;
        % current total stiffness matrix
        K = K + K_star;
    end
end

K = t.*K;

disp(K);