% Element Geometry
x_e = [0 2 2 0]; y_e = [-1 -1 2 2 ]; % Nodal Locations (m)
area = 6; % Area of Element (m^2)

% Displacements at the corner points (m)
for i = 1:4
   u_corners = [0.005 -0.005 -0.010 0.010];
   v_corners = [0.005 0.010 0.015 0.010];
end

% Generate a meshgrid for the region
[x, y] = meshgrid(linspace(0, 2, 25), linspace(-1, 2, 25));

% Shape functions
N1 = @(x, y) (1/6) .* (x - x_e(2)) .* (y - y_e(4));
N2 = @(x, y) (-1/6) .* (x - x_e(1)) .* (y - y_e(4));
N3 = @(x, y) (1/6) .* (x - x_e(1)) .* (y - y_e(1));
N4 = @(x, y) (-1/6) .* (x - x_e(2)) .* (y - y_e(1));

% Interpolate displacements using Shape functions
u_interpolated = N1(x, y) * u_corners(1) + N2(x, y) * u_corners(2) + N3(x, y) * u_corners(3) + N4(x, y) * u_corners(4);
v_interpolated = N1(x, y) * v_corners(1) + N2(x, y) * v_corners(2) + N3(x, y) * v_corners(3) + N4(x, y) * v_corners(4);

% Calculate the total displacement
total_displacement_interpolated = sqrt(u_interpolated.^2 + v_interpolated.^2);

% Create the figure
figure;

% Add corner points to the plot for reference
scatter(x_e, y_e, 'ro', 'filled');
hold on;

% Connect corners with lines
plot([x_e, x_e(1)], [y_e, y_e(1)], 'k', 'LineWidth', 1.5);

% Plot displacement vectors throughout the element
quiver(x(:), y(:), u_interpolated(:), v_interpolated(:), 'b', 'LineWidth', 1);

hold off;

title('Vectors');
xlabel('x');
ylabel('y');
% Element Center Point (1.0,0.5)
center_x = 1.0;
center_y = 0.5;

% Find the index of the closest point in the meshgrid
[~, idx_x] = min(abs(x(1, :) - center_x));
[~, idx_y] = min(abs(y(:, 1) - center_y));

% Interpolated displacement components at (1.0,0.5)
u_at_center = u_interpolated(idx_y, idx_x);
v_at_center = v_interpolated(idx_y, idx_x);
total_displacement_at_center = total_displacement_interpolated(idx_y, idx_x);

% Display the values
disp(['u_interpolated at ', (num2str(center_x)), ',', (num2str(center_y)), ' (center) = ', num2str(u_at_center)]);
disp(['v_interpolated at ', (num2str(center_x)), ',', (num2str(center_y)), ' (center) = ', num2str(v_at_center)]);
disp(['Total Displacement at ', num2str(center_x), ',', num2str(center_y), ' (center) = ', num2str(total_displacement_at_center)]);