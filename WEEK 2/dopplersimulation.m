% System spec -------------------------------------------
velocity = 0.3;
wavelength = 850e-9;
waveNumber = (2 * pi) / wavelength;
fs = (3e7) / wavelength; % Sampling frequency (Hz)
   
% Dependent Variables ------------------------------------
T = 1e-5; % Total sampling time duration (seconds)
dt = 1e-9;
%r = 0:0.00001:0.1;
r = 0; % fixed
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

motionType = 'quadratic';
Z0 = 0
switch motionType
    case 'linear'
        Z = @(t) Z0 + velocity*t;

    case 'quadratic'
        % Degree 2 polynomial
        Z = @(t) Z0 + 10 * (velocity/(2*T^(1))) * t.^2;
        
    case 'cubic'
        % Degree 3 polynomial
        Z = @(t) Z0 + 25*(velocity/(3*T^(2))) * t.^3;
        
    case 'quartic'
        % Degree 4 polynomial
        Z = @(t) Z0 + 10*(velocity/(4*T^(3))) * t.^4;
        
    case 'quintic'
        % Degree 5 polynomial
        Z = @(t) Z0 + 10*(velocity/(5*T^(4))) * t.^5;
    case 'sinusoidal'
        % z(t) = A sin(ωt), choose ω and A so max vel ≈ velocity
        f = 1e5; % Hz
        A = velocity/(2*pi*f); 
        Z = @(t) Z0 + A*sin(2*pi*f*t);

    case 'sawtooth'
        % Proper sawtooth: linear rise, sharp drop
        f = 1e5; % Hz
        A = velocity/(2*pi*f);
        Z = @(t) Z0 + A*sawtooth(2*pi*f*t); % width=1 → classic sawtooth

    case 'pulse'
        f = 1e5;
        A = velocity/(2*pi*f);
        Z = @(t) Z0 + A*square(2*pi*f*t);

    case 'ramp'
        % Symmetric triangular wave
        f = 1e5;
        A = velocity/(2*pi*f);
        Z = @(t) Z0 + A*sawtooth(2*pi*f*t); % width=0.5 → triangle
end

position = @(t) Z(t) + continuous_Z_noise;

% Signal -----------------------------------------------
signal = ampT(r) .* (exp(1i*waveNumber*position(t)));
% pulse_train = (1/(round(T/dt))+1) * (mod(0:dt:T, 2) <= 1e-4);
pulse_train = 1 * (mod(0:dt:T, 2) <= 1e-5);
% Convolve with the rectangular sampling function
%convSignal = ifft(conv(fft(signal), fft(ones(1,round(T/dt))), 'same'));
convSignal = conv(signal, ones(1,round(T/dt)), 'same')
%convSignal = ifft(fft(signal) .* fft(ones(1,round(T/dt)+1)));
%windowSize = round(T/dt) + 1;  
%rectWindow = ones(1, windowSize);  
%convSignal = conv(signal, rectWindow, 'same');  
%convSignal = convSignal / sum(rectWindow); 
% window_width_sec = 1e-9; 
% 
% window_width_samples = round(window_width_sec / dt);
% 
% averaging_window = ones(1, window_width_samples) / window_width_samples;
% 
% averaged_position = conv(signal, averaging_window, 'same');


%convSignal = signal .* pulse_train

% window_width_sec = 1e-5;
% window_width_samples = round(window_width_sec / dt);
% averaging_window = ones(1, window_width_samples) / window_width_samples;
% 
% averaged_position = conv(signal, averaging_window, 'same');
% 
% time = 0:dt:T; 
% convSignal= cumtrapz(signal, time);

% % 
% phase_info = imag(signal);
% intensity_info = real(signal);
% phase = angle(signal);

% % 
phase_info = imag(convSignal);
intensity_info = real(convSignal);
phase = angle(convSignal);


%phase_info = imag(averaged_position);
%intensity_info = real(averaged_position);
%phase = angle(averaged_position);

%phase_info = imag(averaged_position);
%intensity_info = real(averaged_position);
%phase = angle(averaged_position);


%phase = unwrap(angle(averaged_position));

dPdt = gradient(phase, dt);
velocity_calc = dPdt / waveNumber;

% Plots -----------------------------------------------
figure;
grid on;
hold on;

% ------------------ Position Info -------------------
% title(['Position vs. time (' motionType ')']);
% xlabel('time (s)');
% ylabel('Position');
% plot(t, position(t), 'b');
% xlim([0 T]); % ensure one full window



% ---------------------- Signal Info -----------------


% title('Signal vs. time');
% xlabel('time(s)');
% ylabel('Signal');
% %plot(t, intensity_info);
% plot(t, phase_info);


% ---------------------- Phase change ----------------
% 
% title('Phase vs. time');
% xlabel('time(s)');
% ylabel('Phase');
% plot (t, angle(signal));

% ----------------------Velocity calculation ----------
% 
% title('Velocity vs. time');
% xlabel('time(s)');
% ylabel('Velocity');
% ylim([-velocity*6 velocity*6])
% plot(t, velocity_calc);


figure; % Open a new figure window

% ------------------ 1. Position vs. Time -------------------
subplot(2,2,1); % Top-left
plot(t, position(t), 'b');
title(['Position vs. Time (' motionType ')']);
xlabel('Time (s)');
ylabel('Position');
xlim([0 T]);

% ------------------ 2. Signal vs. Time ---------------------
subplot(2,2,2); % Top-right
plot(t, intensity_info); % or intensity_info if desired
title('Signal vs. Time');
xlabel('Time (s)');
ylabel('Signal');

% ------------------ 3. Phase vs. Time ----------------------
subplot(2,2,3); % Bottom-left
plot(t, angle(convSignal));
title('Phase vs. Time');
xlabel('Time (s)');
ylabel('Phase');

% ------------------ 4. Velocity vs. Time -------------------
subplot(2,2,4); % Bottom-right
plot(t, velocity_calc);
title('Velocity vs. Time');
xlabel('Time (s)');
ylabel('Velocity');
ylim([-velocity*6 velocity*6]);

% ------------------ Optional Formatting -------------------
sgtitle('Motion and Signal Analysis'); % Add a global title for all subplots