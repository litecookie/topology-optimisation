% My student number is 1010188967
% Corresponding A - 9 
% Corresponding B - 6
% Corresponding C - 7

% Node location for element starting with the top left node as node 1

x_e = 0.01.*[0 400 400 0 0];
y_e = 0.01.*[0 0 -100 -100 0];

% Node location for hole starting with the top left node as node 1

A = 9;
B = 6;
C = 7;

x_h2 = 150 + 50 + 5*A;
x_h = [150 x_h2 x_h2 150 150];

y_h1 = -100 + (20 + 2*B + 40 - B);
y_h2 = -100 + (40 - B);
y_h = [y_h1 y_h1 y_h2 y_h2 y_h1];

load_x = 400;
load_y = (-100+10*C);

E = 71.7*10^9;

scale = 1;
px1 = 150; px2 = 245;
py1 = -34; py2 = -66;
ulx1 = (px1+1)/scale; ulx2 = (px2+1)/scale;
uly1 = (py1+1)/scale; uly2 = (py2+1)/scale;

nelx = 160;
nely = 40;
elx = 400/2.5;
ely = 30/2.5;
n2 = (nely+1)*(elx-1) +ely;