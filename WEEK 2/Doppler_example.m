%
% clear;
load('matlab_mirror_.15_40_.25.mat');
im_data_pro=im_data;% 3D measurement

%% determine x,y modulation and obtain M-OCPM image 
% The method was described in detail in the following paper
%X. Liu, R. Bhakta, E. Kryvorutsky et al., "Imaging cells and nanoparticles using modulated optically 
% computed phase microscopy," Scientific Reports 15, 3157 (2025).
peak_index=4;N0=4;
size0=size(im_data_pro);
N_window=size0(3);
data1d=abs(ifft(im_data_pro(304,592,:)));
N0=4;[max_fft,peak_index]=max(data1d(N0:N_window/2))
data3d_fd=ifft(im_data_pro(:,:,1:N_window),N_window,3);

data_3d_plane=data3d_fd(:,:,peak_index+N0-1);
data_3d_plane_fd=fft2(data_3d_plane,1024,1024); 

im_int_fd_y=fft(data_3d_plane,1024,1);
Y_1d=mean(abs(im_int_fd_y'));
[v_y,N_y]=max(Y_1d);

im_int_fd_x=fft(data_3d_plane,1024,2);
X_1d=mean(abs(im_int_fd_x));
[v_x,N_x]=max(X_1d);

data_3d_plane_fd(800:900,:)=0;
data_3d_plane_fd_shift=circshift(data_3d_plane_fd,[-N_y -N_x]);% frequency shifting or demodulation
im_output=ifft2(data_3d_plane_fd_shift,1024,1024);% convert to spatial domain
%% M-DPM: Modulated Doppler Microscopy

im_data_k=ifft(im_data_pro,[],3);% Fourier transform
im_data_k(:,:,1:2)=0;
im_data_k(:,:,65:128)=0;% eliminate half of the freqnecy components of the real measurement

im_data_complex=im_data_k;


im_data_kx=fft(im_data_complex,[],2);
im_data_kx_ky=fft(im_data_kx,[],1);
im_data_kx_ky(1,:,:)=0;
im_data_kx_ky(:,1,:)=0;

im_data_kx_ky(800:900,:,:)=0;
im_data_kx_ky_shift=circshift(im_data_kx_ky,[-(N_y) -(N_x) -6]);

im_data_shift=ifft(im_data_kx_ky_shift,[],1);
im_data_shift=ifft(im_data_shift,[],2);

im_f=fft(im_data_shift,[],3);
% convert the signal back to time domain and im_f is a complex array from which Doppler signal can be extracted

Doppler0=angle(im_f(:,:,3:128).*conj(im_f(:,:,1:126)));% calculate Doppler signal
% skip one sampling point for Doppler calculation

%%

time_t=(1:1:126)*0.2;% sampling interval is 0.2s
time_pzt=1:1:100;
v_pzt=zeros(size(time_pzt)); v_pzt(21:70)=0.10; v_pzt_c=cumsum(v_pzt);

dis_doppler=455*cumsum(permute(mean(mean(Doppler0(286-0:286+0,584-0:584+0,:,1),1),2),[3 1 2 4]))/4/pi/1.33/2;
% Convert Doppler phase to displacement
dis_pzt=1000*v_pzt_c*1.5*2.8/75;% piezo displacement

%%
doppler_raw=(mean((Doppler0(:,:,30:1:80,1)),3));% average to improve signal to noise ratio
doppler_raw=455*(doppler_raw)/4/pi/1.33/0.4;% convert Doppler phase shift to velocity
im_std=sqrt(mean((455*Doppler0(:,:,30:1:80,1)/4/pi/1.33/0.4).^2,3)-doppler_raw.^2);
im_std=im_std;mask0=zeros(1024);mask0(im_std>75)=0;mask0(im_std<=75)=1;
% use the temporal variation of Doppler signal to create a mast for display
V_max=25;
map_r=zeros(256,1);map_r(157:256)=(1:100)/100;map_g=zeros(256,1);map_g(1:100)=(100:-1:1)/100;map_b=zeros(256,1);map_c=[map_r,map_g,map_b];
% create a colormap to display Doppler signal
%%
figure('Position',[256 256 1200 800])
ax1=subplot(2,3,1), imagesc(Doppler0(:,:,50)), title('Single frame doppler');
ax2=subplot(2,3,2), imagesc(doppler_raw), title('Averaged doppler');
ax3=subplot(2,3,3), imagesc(doppler_raw.*mask0,[-V_max V_max]),colormap(ax3, map_c), title('Thresholded,color-coded doppler');
ax4=subplot(2,3,4), plot(time_t,dis_doppler,'k','LineWidth',2), hold on, plot(time_pzt*0.25,dis_pzt,'r--','LineWidth',2)
xlabel('time(s)','FontSize',10),ylabel('displacement(nm)','FontSize',10), title('velocity tracking');
ax5=subplot(2,3,5),imagesc(abs(im_output)),title('Amplitude Image');
ax6=subplot(2,3,6),imagesc(angle(im_output)),title('Phase Image');
