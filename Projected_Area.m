% Parameters
p = 200; % Diameter of the paraboloid
a = 15; % Half-length of the square plate side
h = 120; % Height of the square plate above the bottom of the paraboloid

% Focus of the paraboloid
f = p / 2;

% Discretize the square plate edges
num_points = 100; % Number of points on each edge
t = linspace(-a, a, num_points);

% Generate the points on the edges of the square plate
x_plate = [t, a*ones(1,num_points), fliplr(t), -a*ones(1,num_points)];
y_plate = [-a*ones(1,num_points), t, a*ones(1,num_points), fliplr(t)];
z_plate = h * ones(size(x_plate));

% Focus point
focus = [0, 0, f];

% Initialize arrays for projected points
x_proj = zeros(size(x_plate));
y_proj = zeros(size(x_plate));
z_proj = zeros(size(x_plate));

% Project each point onto the paraboloid
for i = 1:length(x_plate)
    % Direction vector from the point on the plate to the focus
    direction = focus - [x_plate(i), y_plate(i), z_plate(i)];
    direction = direction / norm(direction); % Normalize the direction vector
    
    % Find the intersection of the line with the paraboloid
    syms t;
    x_t = x_plate(i) + t * direction(1);
    y_t = y_plate(i) + t * direction(2);
    z_t = z_plate(i) + t * direction(3);
    paraboloid_eq = z_t == (x_t^2 + y_t^2) / (4 * f);
    t_sol = solve(paraboloid_eq, t);
    t_sol = double(t_sol);
    
    % Select the correct intersection point (t > 0)
    t_proj = t_sol(t_sol > 0);
    
    % Calculate the projected point coordinates
    x_proj(i) = x_plate(i) + t_proj * direction(1);
    y_proj(i) = y_plate(i) + t_proj * direction(2);
    z_proj(i) = (x_proj(i)^2 + y_proj(i)^2) / (4 * f);
end

% Calculate area of the projected shape using polyarea
area = polyarea(x_proj, y_proj);

% Display the result
fprintf('The projected area on the paraboloid is: %.2f\n', area);

% Plot for visualization
figure;
hold on;
% Plot the paraboloid
[X, Y] = meshgrid(linspace(-p, p, 100), linspace(-p, p, 100));
Z = (X.^2 + Y.^2) / (4 * f);
surf(X, Y, Z, 'EdgeColor', 'none', 'FaceAlpha', 0.5);

% Plot the square plate
plot3(x_plate, y_plate, z_plate, 'r', 'LineWidth', 2);

% Plot the projection points
plot3(x_proj, y_proj, z_proj, 'bo', 'MarkerSize', 5, 'MarkerFaceColor', 'b');

% Plot settings
xlabel('X');
ylabel('Y');
zlabel('Z');
title('Projection of Square Plate onto Paraboloid');
axis equal;
grid on;
hold off;
