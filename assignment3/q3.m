% -------------------------------------------------------------------------

% Matlab code for Problem Set 3 - Question 3

% -------------------------------------------------------------------------

% clear data space
clear;
% close figure windows
close all;

E = 70e9; % Young's modulus (Pa)
nu = 0.3; % Poisson’s ratio
% matrix of elastic constants
D = (E/(1-nu^2)) * [1 nu 0; nu 1 0; 0 0 (1-nu)/2]; % (Pa)

% Displacements at the nodal locations
d = 1e-3.*[5 5 -5 10 -10 15 10 10];

% Coordinates of vertices of rectangle element declared.
x1 = 0;
y1 = -1;
x2 = 2;
y2 = -1;
x3 = 2;
y3 = 2;
x4 = 0;
y4 = 2;

% Center displacement calculation
cx = 1;
cy = 0.5;
A = 3*2;

Nf = zeros(4);
Nf(1) = (1/A)*(cx - x2)*(cy - y4);
Nf(2) = -(1/A)*(cx - x1)*(cy - y4);
Nf(3) = (1/A)*(cx - x1)*(cy - y1);
Nf(4) = -(1/A)*(cx - x2)*(cy - y1);
        
N = [Nf(1) 0 Nf(2) 0 Nf(3) 0 Nf(4) 0; 0 Nf(1) 0 Nf(2) 0 Nf(3) 0 Nf(4)];
  
% Computed displacement values for the given center point.
v = N*d';

H = (1/A)*[(cy - y4) 0 -(cy - y4) 0 (cy - y1) 0 -(cy - y1) 0; 0 (cx - x4) 0 -(cx - x1) 0 (cx - x1) 0 -(cx - x1); (cx - x2) (cy - y4) -(cx - x1) -(cy - y4) (cx - x1) (cy - y1) -(cx - x2) -(cy - y1)];

% Strain function calculated.
strain = H*d';

stress = D*strain;

disp(stress);