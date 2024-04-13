% -------------------------------------------------------------------------

% Matlab code for Problem Set 3 - Question 2

% -------------------------------------------------------------------------

% clear data space
clear;
% close figure windows
close all;

% Coordinates of vertices of rectangle element declared.

x1 = 0;
y1 = -1;
x2 = 2;
y2 = -1;
x3 = 2;
y3 = 2;
x4 = 0;
y4 = 2;

% Variables declared to plot the rectangle element in the graph

lx = [x1 x2 x3 x4 x1];
ly = [y1 y2 y3 y4 y1];

% Displacements at the nodal locations

d = (1e-3).*[10 10 5 5 -5 10 -10 15];

% Declare the variables to plot the contour.

n = 20;

jx = linspace(0, 2, n);
jy = linspace(-1, 2, n);
[nx, ny] = meshgrid(jx, jy);

A = 3*2;

% U and V vectors of displacements initialised.
u_inter = zeros(n, n);
v_inter = zeros(n, n);

for j = 1:n
    for i = 1:n
        
        % Iteration to initialise Nf matrix - 
        % substituted for each co-ordinate in the graph
        % Compute the interpolated linear displacement
        % Nfc matrix initialised from lecture equations.

        Nf = zeros(4);
        Nf(1) = (1/A)*(jx(i) - x2)*(jy(j) - y4);
        Nf(2) = -(1/A)*(jx(i) - x1)*(jy(j) - y4);
        Nf(3) = (1/A)*(jx(i) - x1)*(jy(j) - y1);
        Nf(4) = -(1/A)*(jx(i) - x2)*(jy(j) - y1);
        
        N = [Nf(1) 0 Nf(2) 0 Nf(3) 0 Nf(4) 0; 0 Nf(1) 0 Nf(2) 0 Nf(3) 0 Nf(4)];
        
        u{i, j} = N*d';
                
        u_inter(i, j) = u{i, j}(1);
        v_inter(i, j) = u{i, j}(2);
    end
end

% Center displacement calculation

cx = 1;
cy = 0.5;

Nf = zeros(4);
Nf(1) = (1/A)*(cx - x2)*(cy - y4);
Nf(2) = -(1/A)*(cx - x1)*(cy - y4);
Nf(3) = (1/A)*(cx - x1)*(cy - y1);
Nf(4) = -(1/A)*(cx - x2)*(cy - y1);
        
N = [Nf(1) 0 Nf(2) 0 Nf(3) 0 Nf(4) 0; 0 Nf(1) 0 Nf(2) 0 Nf(3) 0 Nf(4)];
        
v = N*d';
  
c_displacement = sqrt(v(1)^2+v(2)^2);

disp(c_displacement);
% ------------------------------------------------------------------------

% ------------------------------------------------------------------------
% Contours for interpolated displacement


figure(1);
quiver(nx, ny, u_inter, v_inter, 'r', 'LineWidth', 1);
hold on;
plot(lx, ly, 'b', 'LineWidth',1);

xlabel('x');
ylabel('y');
title('Interpolated Displacement Vectors for given element');
hold off;