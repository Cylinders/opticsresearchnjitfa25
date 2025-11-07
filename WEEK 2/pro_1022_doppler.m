clear
load('matlab_motion_.5v_5_1.mat');


%%
N_size=1024;
x=zeros(N_size);x0=1:1:N_size;
for i=1:1:N_size
    x(:,i)=x0;
end
y=x';
%%
data_input=im_data;
N_window=256;
data1d=abs(fft(permute(data_input(512,512,1:N_window),[3 1 2])));
figure
plot(data1d(3:64))
N0=4;[max_fft,peak_index]=max(data1d(N0:N_window/2));

data3d_fd=ifft(data_input(:,:,1+N_window*0:N_window*1),N_window,3);




%%
roi_top=128;roi_bottom=1024-128;
roi_left=128;roi_right=1024-128;


%%
N_interval=100;% interval between data points used to calculate Doppler phase shift
data3d_fd_half=data3d_fd;
data3d_fd_half(:,:,1:2)=0;data3d_fd_half(:,:,N_window/2+1:N_window)=0;
data3d_fd_half=circshift(data3d_fd_half,[0 0 -(peak_index+N0-1-1)]);
data3d_fd_half=fft(data3d_fd_half,[],3);
data_doppler_3D=angle(data3d_fd_half(:,:,1+N_interval:N_window).*conj(data3d_fd_half(:,:,1:N_window-N_interval)));


 data_doppler_3D(data_doppler_3D>0.9*pi)=data_doppler_3D(data_doppler_3D>0.9*pi)-pi;
 data_doppler_3D(data_doppler_3D<-0.9*pi)=data_doppler_3D(data_doppler_3D<-0.9*pi)+pi;



figure
for i=1:1:size(N_window-N_interval,1)
    imagesc(data_doppler_3D(:,:,i))
    drawnow
end

N=19;
d_doppler=permute(data_doppler_3D(330,527,:),[3 1 2])*455/1.33/4/pi/N_interval/0.2;
for i=1:1:floor(size(d_doppler, 1)/N/2)
    for j=1:1:N
        d_const((i-1)*N*2+j)=1;
    end


    for j=1:1:N
        d_const((i-1)*N*2+j+N)=-1;
    end
end
%%

time_t=(1:1:size(d_doppler, 1))*0.2;
figure
plot(time_t,(d_doppler))
hold on
plot(time_t(1:length(d_const)),circshift(0.5*1000*(d_const)*2.8/75,6))
xlabel('time')
ylabel('velocity(nm/s)')
legend ('doppler estimated','pzt')
