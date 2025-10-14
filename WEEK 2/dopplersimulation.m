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
r = 0:0.00001:0.1;
% r = 0; % fixed
t = 0:dt:T;
%t = 0; % fixed
% Noise ---------------------------------------------
noise_amplitude = 1e-3; % Adjust for desired noise strength
continuous_amp_noise = noise_amplitude * randn(size(t));
noiseZ_mag = 1e-10; %1e-10
continuous_Z_noise = noiseZ_mag * randn(size(t));
% Amplitude --------------------------------------------
I0 = 2;
absorbtion_u = 5;
amp0 = @(r) I0 * exp(-absorbtion_u * r);
ampT = @(r) amp0(r) + continuous_amp_noise;
%% MOTION SELECTION

motionType = 'cubic';
Z0 = 0
switch motionType
    case 'linear'
        Z = @(t) Z0 + velocity*t;

    case 'quadratic'
        % Degree 2 polynomial
        Z = @(t) Z0 + 50 * (velocity/(2*T^(1))) * t.^2;
        
    case 'cubic'
        % Degree 3 polynomial
        Z = @(t) Z0 + 75*(velocity/(3*T^(2))) * t.^3;
        
    case 'quartic'
        % Degree 4 polynomial
        Z = @(t) Z0 + 10*(velocity/(4*T^(3))) * t.^4;
        
    case 'quintic'
        % Degree 5 polynomial
        Z = @(t) Z0 + 10*(velocity/(5*T^(4))) * t.^5;
    case 'sinusoidal'
       
        f = 1e5; % Hz
        A = velocity/(2*pi*f); 
        Z = @(t) Z0 + A*sin(2*pi*f*t);

    case 'sawtooth'
        
        f = 1e6; % Hz
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

signal = ampT(r) .* (exp(1i*waveNumber*position(t)));



pulse_train = 1 * (mod(0:dt:T, 2) <= 1e-5);
windowWidth = 100;    % 16, 32, 100, 1000: NUM SAMPLES

rectWindow = ones(1, windowWidth) / windowWidth;
convSignal = conv(signal, rectWindow, 'same');  

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
plot(t, phase_info);
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
plot(t, phase_info);
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
plot(t, velocity_calc, 'r', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Velocity');
title('Comparison: Original vs Convolved Velocity');
legend({'Original Velocity', 'Convolved Velocity'}, 'Location', 'best');
grid on;