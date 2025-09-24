% Define time steps
timestep = 0:0.001:1; 

% Compute linear motion (assumes linearmotion is a function)
position = linearmotion(timestep); 
disp(position)

% Compute signal based on position and time
sig = signal(position, timestep);  % store in 'sig' instead of 'signal'

% Plot linear motion vs time
figure;
plot(timestep, position, 'b', 'DisplayName', 'Linear Motion'); 
hold on;  % Retain current plot when adding new plots

% Plot signal vs time
plot(timestep, sig, 'r', 'DisplayName', 'Signal');

% Add labels, title, grid, and legend
xlabel('Time (s)');
ylabel('Value');
title('Linear Motion and Signal over Time');
grid on;
legend show;
