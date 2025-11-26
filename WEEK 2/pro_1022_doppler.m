%% --- Full Doppler Velocity Pipeline with N-frame Window Averaging ---
clear
load('matlab_motion_.5v_5_1.mat');  % 1024x1024x256

%% --- Parameters ---
N = 10         % number of frames per averaging window
T = 0.2;       % time between successive A-scans (s)
f0 = 455e-9;      % system frequency (Hz)
theta = 0;     % angle in radians for velocity calculation
n = 5;         % parameter for velocity scaling

[rows, cols, frames] = size(im_data);

%% --- Step 1: Hilbert transform ---
dop = hilbert(im_data);
I = real(dop);
Q = imag(dop);

%% --- Step 2: Compute average phase shift over N-frame windows ---
% avgPhaseShift = [];  % initialize
% 
% k = 1;
% while k + N - 1 <= frames
%     numerator = zeros(rows, cols);
%     denominator = zeros(rows, cols);
% 
%     % Sum over N-1 consecutive pairs within this window
%     for i = 0:(N-2)
%         numerator = numerator + ( Q(:,:,k+i) .* I(:,:,k+i+1) - I(:,:,k+i) .* Q(:,:,k+i+1) );
%         denominator = denominator + ( I(:,:,k+i) .* I(:,:,k+i+1) + Q(:,:,k+i) .* Q(:,:,k+i+1) );
%     end
% 
%     avgPhaseShift(:,:,end+1) = -atan(numerator ./ denominator);  % append to array
%     k = k + N;  % jump by N frames for the next window
% end


avgPhaseShift = [];  % initialize

for k = 1:(frames - N + 1)  % sliding window from frame 1 to frame (frames-N+1)
    numerator = zeros(rows, cols);
    denominator = zeros(rows, cols);
    
    % Sum over N-1 consecutive pairs within this window
    for i = 0:(N-2)
        numerator = numerator + ( Q(:,:,k+i) .* I(:,:,k+i+1) - I(:,:,k+i) .* Q(:,:,k+i+1) );
        denominator = denominator + ( I(:,:,k+i) .* I(:,:,k+i+1) + Q(:,:,k+i) .* Q(:,:,k+i+1) );
    end
    
    avgPhaseShift(:,:,end+1) = -atan(numerator ./ denominator);  % append to array
end


%% --- Step 3: Doppler frequency ---
fd = avgPhaseShift / (2*pi*T);

%% --- Step 4: Velocity ---
velocity = (f0 ./ (2*1.33*cos(theta))) .* fd;

%% --- Step 5: Average velocity per frame/window ---
D_doppler = squeeze(mean(mean(velocity,1),2)) * 10e9;  % 1D array, one value per N-frame window

%% --- Step 6: Generate PZT reference ---


%% --- Step 7: Time vector ---
time_t = (1:length(D_doppler)) * T * N;  % account for N-frame windowing

%% --- Step 8: Plot Doppler vs PZT ---
figure;
d_const = ones(1, length(D_doppler))
plot(time_t, D_doppler, 'b', 'LineWidth', 1.5);
plot(time_t(1:length(d_const)),circshift(0.5*1000*(d_const)*2.8/75,6));
xlabel('time (s)');
ylabel('velocity (nm/s)');
legend('doppler estimated', 'pzt');
title('Average Doppler Velocity vs PZT Reference');
grid on;

%%
% N_size=1024;
% x=zeros(N_size);x0=1:1:N_size;
% for i=1:1:N_size
%     x(:,i)=x0;
% end
% y=x';
% %%
% data_input=im_data;
% N_window=256;
% data1d=abs(fft(permute(data_input(512,512,1:N_window),[3 1 2])));
% 
% N0=4;[max_fft,peak_index]=max(data1d(N0:N_window/2));
% 
% data3d_fd=ifft(data_input(:,:,1+N_window*0:N_window*1),N_window,3);
% 
% 
% 
% 
% %%
% roi_top=128;roi_bottom=1024-128;
% roi_left=128;roi_right=1024-128;
% 
% 
% %%
% N_interval=100;% interval between data points used to calculate Doppler phase shift
% 
% data3d_fd_half=data3d_fd;
% 
% data3d_fd_half(:,:,1:2)=0;data3d_fd_half(:,:,N_window/2+1:N_window)=0;
% 
% data3d_fd_half=circshift(data3d_fd_half,[0 0 -(peak_index+N0-1-1)]);
% 
% data3d_fd_half=fft(data3d_fd_half,[],3);
% 
% data_doppler_3D=angle(data3d_fd_half(:,:,1+N_interval:N_window).*conj(data3d_fd_half(:,:,1:N_window-N_interval)));
% 
% 
% data_doppler_3D(data_doppler_3D>0.9*pi)=data_doppler_3D(data_doppler_3D>0.9*pi)-pi;
% 
% data_doppler_3D(data_doppler_3D<-0.9*pi)=data_doppler_3D(data_doppler_3D<-0.9*pi)+pi;
% 
% 
% 
% figure
% for i=1:1:size(N_window-N_interval,1)
%     imagesc(data_doppler_3D(:,:,i))
%     drawnow
% end
% 
% N=19;
% d_doppler=permute(data_doppler_3D(330,527,:),[3 1 2])*455/1.33/4/pi/N_interval/0.2;
% for i=1:1:floor(size(d_doppler, 1)/N/2)
%     for j=1:1:N
%         d_const((i-1)*N*2+j)=1;
%     end
% 
% 
%     for j=1:1:N
%         d_const((i-1)*N*2+j+N)=-1;
%     end
% end
% %%
% 
% time_t=(1:1:size(d_doppler, 1))*0.2;
% figure
% plot(time_t,(d_doppler))
% hold on
% plot(time_t(1:length(d_const)),circshift(0.5*1000*(d_const)*2.8/75,6))
% xlabel('time')
% ylabel('velocity(nm/s)')
% legend ('doppler estimated','pzt')
