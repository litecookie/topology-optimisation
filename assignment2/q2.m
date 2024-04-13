% -------------------------------------------------------------------------

% Matlab code for Problem Set 2 - Question 2

% -------------------------------------------------------------------------

% Coordinates of vertices of rectangle element declared.

x1 = 1;
y1 = 1;
x2 = 3;
y2 = 1;
x3 = 3;
y3 = 2;
x4 = 1;
y4 = 2;

% Variables declared to plot the rectangle element in the graph

lx = [x1 x2 x3 x4 x1];
ly = [y1 y2 y3 y4 y1];

% Displacements of the rectangle element declared.

syms xa ya

ux = (3/8)*xa - (2/3)*(xa/(ya^2));
uy = (4/(xa^3)) - (ya/8);

n = 4;                              % Number of nodes in the element.
dux = zeros(n, 1);
duy = zeros(n, 1);

for i = 1:4
    dux(i) = subs(ux, [xa, ya], [lx(i), ly(i)]);
    duy(i) = subs(uy, [xa, ya], [lx(i), ly(i)]);
end

% Preallocate C for efficiency
d = zeros(length(dux) + length(duy), 1);

% Interleave elements
d(1:2:end) = dux;  % Assign elements of A to odd indices of C
d(2:2:end) = duy;  % Assign elements of B to even indices of C

% Declare the variables to plot the contour.

jx = linspace(1, 3, 30);
jy = linspace(1, 2, 30);
[nx, ny] = meshgrid(jx, jy);

l = length(jx);

A = 2*1;

% Interpolated displacements variable initialised.
dx = zeros(l, l);

% Exact analytical displacements variable initialised.
adx = zeros(l, l);

% Nf matrix initialised from lecture equations

syms x y;

Nf = cell(4);

Nf{1} = (1/A)*(x - x2)*(y - y4);
Nf{2} = -(1/A)*(x - x1)*(y - y4);
Nf{3} = (1/A)*(x - x1)*(y - y1);
Nf{4} = -(1/A)*(x - x2)*(y - y1);

for i = 1:l
    for j = 1:l
        
        % Iteration to initialise Nf matrix - 
        % substituted for each co-ordinate in the graph
        % Compute the interpolated linear displacement
        
        Nfc = zeros(4);
        for ck = 1:4
            Nfc(ck) = subs(Nf{ck}, [x, y], [jx(i), jy(j)]);
        end

        N = [Nfc(1) 0 Nfc(2) 0 Nfc(3) 0 Nfc(4) 0; 0 Nfc(1) 0 Nfc(2) 0 Nfc(3) 0 Nfc(4)];
        
        u{i, j} = N*d;
        
        dx(i,j) = sqrt( (u{i,j}(1))^2 + (u{i,j}(2))^2 );

        % Compute the exact analytical displacement

        adxi = subs(ux, [xa, ya], [jx(i), jy(j)]);
        adyi = subs(uy, [xa, ya], [jx(i), jy(j)]);

        adx(i,j) = sqrt( (adxi)^2 + (adyi)^2 );
        
    end
end

% ------------------------------------------------------------------------

% ------------------------------------------------------------------------
% Contours for interpolated displacement

figure(1);

contourf(ny, nx, dx);
clabel(contourf(ny, nx, dx));
colorbar;
hold on; 

xlabel('x');
ylabel('y');
title('Interpolated Displacement Contour for given element');

% ------------------------------------------------------------------------
% Contours for exact displacement

figure(2);

contourf(ny, nx, adx);
clabel(contourf(ny, nx, adx));
colorbar;
hold on; 

xlabel('x');
ylabel('y');
title('Exact Analytical Displacement Contour for given element');
