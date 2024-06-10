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
parfor j = 1:length(h_values)
    h = h_values(j);
    z_plate = h * ones(size(x_plate));

    % Compute the direction vectors for each height
    directions = [zeros(size(x_plate)); zeros(size(y_plate)); f * ones(size(x_plate))] - [x_plate; y_plate; z_plate];
    directions = directions ./ vecnorm(directions);

    % Project each point onto the paraboloid
    t_proj = zeros(size(x_plate));
    
    for i = 1:length(x_plate)
        % Initial guess for t using quadratic approximation
        t_guess = (4*f*(h - z_plate(i))) / (1 + (x_plate(i)^2 + y_plate(i)^2) / (4*f*h));
        
        % Refine the guess using Newton-Raphson method
        t_proj(i) = t_guess;
        for k = 1:5
            x_t = x_plate(i) + t_proj(i) * directions(1,i);
            y_t = y_plate(i) + t_proj(i) * directions(2,i);
            z_t = z_plate(i) + t_proj(i) * directions(3,i);
            f_val = z_t - (x_t^2 + y_t^2) / (4 * f);
            df_dt = directions(3,i) - (2*x_t*directions(1,i) + 2*y_t*directions(2,i)) / (4*f);
            t_proj(i) = t_proj(i) - f_val / df_dt;
        end
    end
    
    % Calculate the projected point coordinates
    x_proj = x_plate + t_proj .* directions(1,:);
    y_proj = y_plate + t_proj .* directions(2,:);
    
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
