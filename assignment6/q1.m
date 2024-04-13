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

% Plotting the contour of rho
figure; % Opens a new figure window
contourf(nx, ny, rho, 'LineColor', 'none'); % Filled contour plot
colorbar; % Shows the color scale
title('Rho Contours');
xlabel('X Coordinate');
ylabel('Y Coordinate');
