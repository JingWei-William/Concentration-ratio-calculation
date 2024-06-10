% Parameters
p = 200; % Diameter of the paraboloid
a = 15; % Half-length of the square plate side
f = p / 2; % Focus of the paraboloid

% Range of heights
h_values = linspace(100, 200, 101); % Adjust the range and number of values as needed

% Array to store the calculated areas
areas = zeros(size(h_values));

% Discretize the square plate edges
num_points = 100; % Number of points on each edge
t = linspace(-a, a, num_points);

% Generate the points on the edges of the square plate
x_plate = [t, a*ones(1,num_points), fliplr(t), -a*ones(1,num_points)];
y_plate = [-a*ones(1,num_points), t, a*ones(1,num_points), fliplr(t)];

% Loop over different heights
for j = 1:length(h_values)
    h = h_values(j);
    z_plate = h * ones(size(x_plate));

    % Initialize arrays for projected points
    x_proj = zeros(size(x_plate));
    y_proj = zeros(size(x_plate));
    z_proj = zeros(size(x_plate));

    % Project each point onto the paraboloid
    for i = 1:length(x_plate)
        % Direction vector from the point on the plate to the focus
        direction = [0, 0, f] - [x_plate(i), y_plate(i), z_plate(i)];
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
    areas(j) = polyarea(x_proj, y_proj);
end

% Display the results
for j = 1:length(h_values)
    fprintf('Height h = %.2f: Projected area = %.2f\n', h_values(j), areas(j));
end

% Plot the results
figure;
plot(h_values, areas, 'LineWidth', 2);
xlabel('Height h');
ylabel('Projected Area');
title('Projected Area on Paraboloid vs. Height of Square Plate');
grid on;
