%% OCT SIMULATION DATA: AVERAGES FROM PAPERS WE FOUND. 
% System spec -------------------------------------------
velocity = 0.3;
wavelength = 850e-9;
waveNumber = (2 * pi) / wavelength;
fs = (3e7) / wavelength; % Sampling frequency (Hz)

%% INPUT VARIABLES
% Dependent Variables ------------------------------------
T = 1e-5; % Total sampling time duration (seconds)
dt = 1e-9;
%r = 0:0.00001:0.1;
r = 0; % fixed
t = 0:dt:T;
%t = 0; % fixed
% Noise ---------------------------------------------
noise_amplitude = 1e-14; % Adjust for desired noise strength
continuous_amp_noise = noise_amplitude * randn(size(t));
noiseZ_mag = 1e-14; %1e-10
continuous_Z_noise = noiseZ_mag * randn(size(t));
% Amplitude --------------------------------------------
I0 = 2;
absorbtion_u = 5;
amp0 = @(r) I0 * exp(-absorbtion_u * r);
ampT = @(r) amp0(r) + continuous_amp_noise;
%% MOTION SELECTION

motionType = 'sinusoidal';
Z0 = 0
switch motionType
    case 'linear'
        Z = @(t) Z0 + velocity*t;

    case 'quadratic'
        Z = @(t) Z0 + 50e-4 * (velocity/(2*T^(1))) * t.^2;
    case 'cubic'
        Z = @(t) Z0 + 100e-3*(velocity/(3*T^(2))) * t.^3;
    case 'quartic'
        Z = @(t) Z0 + 10*(velocity/(4*T^(3))) * t.^4;
    case 'polynomial(n=20)'
        Z = @(t) Z0 + 2e-2 * (velocity / (20 * T^(19))) * t.^20;
    case 'sinusoidal'
        f = 1e5; 
        A = velocity/(2*pi*f); 
        Z = @(t) Z0 + A*100e-4*sin(2*pi*f*t);
    case 'sawtooth'
        f = 1e6; 
        A = velocity/(2*pi*f);
        Z = @(t) Z0 + A*sawtooth(2*pi*f*t); 
    case 'pulse'
        f = 1e5;
        A = velocity/(2*pi*f);
        Z = @(t) Z0 + A*square(2*pi*f*t);
    case 'ramp'
        f = 1e6;
        A = velocity/(2*pi*f);
        Z = @(t) Z0 + A*sawtooth(2*pi*f*t); 
end
%%  CONVOLUTION WITH RECTANGULAR SAMPLING FUNCTION
position = @(t) Z(t) + continuous_Z_noise;
veloc = gradient(position(t), dt);
signal = ampT(r) .* (exp(1i*waveNumber*position(t)));




windowWidth = 1000;    % 16, 32, 100, 1000: NUM SAMPLES

rectWindow = ones(1, windowWidth) / windowWidth;
convSignal = conv(signal, rectWindow, 'same');  
%convSignal = conv(signal, rectWindow);

phase = angle(convSignal);
phase = unwrap(phase);                          
dPdt = gradient(phase, dt);
velocity_calc = dPdt / waveNumber;
%%  CONTROL

orig_phase_info = imag(signal);
orig_intensity_info = real(signal);
orig_phase = angle(signal);
orig_dPdt = gradient(orig_phase, dt);
orig_velocity_calc = orig_dPdt / waveNumber;
%% Plotting, Data Visualization
% ------------------ Existing figure ------------------
figure;
sgtitle('Signal Comparison: Original vs Convolved');
grid on;
hold on;

% ------------------ 1. Position vs. Time -------------------
subplot(4,2,1);
plot(t, position(t), 'b');
title(['Position vs. Time (' motionType ')']);
xlabel('Time (s)');
ylabel('Position');
xlim([0 T]);

% ------------------ 2. Position (duplicate for alignment) --
subplot(4,2,2);
plot(t, convSignal);
title('Convolved Phase vs. Time');
xlabel('Time (s)');
ylabel('Position');
xlim([0 T]);

% ------------------ 3. Original Signal ---------------------
subplot(4,2,3);
plot(t, orig_intensity_info);
title('Original Signal vs. Time');
xlabel('Time (s)');
ylabel('Signal');

% ------------------ 4. Convolved Signal --------------------
subplot(4,2,4);
plot(t, position(t));
title('Convolved Signal vs. Time');
xlabel('Time (s)');
ylabel('Signal');

% ------------------ 5. Original Phase ----------------------
subplot(4,2,5);
plot(t, orig_phase);
title('Original Phase vs. Time');
xlabel('Time (s)');
ylabel('Phase');

% ------------------ 6. Convolved Phase ---------------------
subplot(4,2,6);
plot(t, phase);
title('Convolved Phase vs. Time');
xlabel('Time (s)');
ylabel('Phase');

% ------------------ 7. Original Velocity -------------------
subplot(4,2,7);
plot(t, orig_velocity_calc);
title('Original Velocity vs. Time');
xlabel('Time (s)');
ylabel('Velocity');
ylim([-velocity*6 velocity*6]);

% ------------------ 8. Convolved Velocity ------------------
subplot(4,2,8);
plot(t, velocity_calc);
title('Convolved Velocity vs. Time');
xlabel('Time (s)');
ylabel('Velocity');
ylim([-velocity*6 velocity*6]);

sgtitle('Motion and Signal Analysis');

% ------------------ 9. Overlap of Original and Convolved Velocity ---------
% Create a new axes spanning the figure width
% optionally, use same figure or create new
figure;
subplot(1,1,1); % single axes

plot(t, orig_velocity_calc, 'b', 'LineWidth', 1.5); hold on;
plot(t, velocity_calc, 'r', 'LineWidth', 1.5); hold on;
plot(t, veloc, 'g', 'LineWidth', 1.5); 
xlabel('Time (s)');
ylabel('Velocity');
title('Comparison: Original vs Convolved Velocity');
legend({'Original Velocity', 'Convolved Velocity', 'Calculated Velocity'}, 'Location', 'best');
grid on;