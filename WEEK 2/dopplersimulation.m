
timestep = 0:0.001:2*pi; 
sig = signal(0, timestep);

phase_unwrapped = unwrap(angle(sig));


phase_change = diff(phase_unwrapped);

dt = timestep(2) - timestep(1);
t_phase = (timestep(1:end-1) + timestep(2:end)) / 2;


f_d_inst = phase_change ./ (2*pi*dt);    

lambda = [];        
f0 = 77e9;           

c = 3e8;            

if isempty(lambda)
    lambda = c / f0;
end


velocity = (f_d_inst .* lambda) / 2;  




figure('Name','Doppler & Velocity','NumberTitle','off','Position',[100 100 900 500]);

subplot(2,1,1);
plot(t_phase, f_d_inst, 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Doppler Shift f_D (Hz)');
title('Instantaneous Doppler Frequency');
grid on;

subplot(2,1,2);
plot(t_phase, velocity, 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Radial Velocity v (m/s)');
title(['Estimated Radial Velocity ( \lambda = ' num2str(lambda) ' m )']);
grid on;


fprintf('Doppler frequency (Hz): mean = %.3f, std = %.3f\n', mean(f_d_inst), std(f_d_inst));
disp(mean(velocity))
