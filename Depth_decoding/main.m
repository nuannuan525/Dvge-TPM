%% get the time traces of all active neurons by CaImAn
% clear all;close all;clc;
% nam =['Im_L.tif'];
% run(fullfile(pwd, 'CaImAn-MATLAB-master', 'demo_script_main.m'));
% save('data_after_CNMF_L.mat','Yr','A_or','C_or','b2','f2','Cn','options','mask_cell','Df');
% clear all;close all;clc;
% nam =['Im_R.tif'];
% run(fullfile(pwd, 'CaImAn-MATLAB-master', 'demo_script_main.m'));
% save('data_after_CNMF_R.mat','Yr','A_or','C_or','b2','f2','Cn','options','mask_cell','Df');

%% load parameters

clear all;close all;clc;
load('data_after_CNMF_L.mat','Yr','A_or','C_or','b2','f2','Cn','options','mask_cell','Df');
[c_L,NeurTT_L] = preprocess_components_L(mask_cell,C_or,15,80,Df);
load('data_after_CNMF_R.mat','Yr','A_or','C_or','b2','f2','Cn','options','mask_cell','Df');
[c_R,NeurTT_R] = preprocess_components_R(mask_cell,C_or,5,40,Df);

load(['psf_parameters.mat'],'psf_11','psf_22'); % parameters of psf

noise_image=32393;  
Im_L = tifread('Im_L.tif')-noise_image; Im_L(Im_L<0)=0;
Im_R = tifread('Im_R.tif')-noise_image; Im_L(Im_L<0)=0;

%% proj_pairs

[select_pairs_all,select_pairs_unique,NeurTT_1,ROI,psf_radio_new,distance] = proj_pairs(NeurTT_L,NeurTT_R,c_L,c_R,Im_L,Im_R,psf_11,psf_22); % select projection pairs
ROI_final=select_pairs_unique(:,6)>0;
final_combinations=select_pairs_unique(ROI_final,:);
NeurTT=NeurTT_1(ROI_final,:);

%% show depth and Ga traces
[nx,ny,n_cell]=size(c_L);

depth_A=zeros(nx,ny);
Img_A=sum(c_L(:,:,final_combinations(:,1)),3);

for hh=1:size(final_combinations,1)
    depth_A(c_L(:,:,final_combinations(hh,1))>0)=max(34-final_combinations(hh,6),1);
end

figure(101);
imagesc(depth_A);axis image; colormap(jet);
xticklabels '';yticklabels '';
for hh=1:size(final_combinations,1)
    mask_L1=(c_L(:,:,final_combinations(hh,1)))>0;
    boundaries_L1 = bwboundaries(mask_L1);
    boundary = boundaries_L1{1};
    hold on;plot(boundary(:,2), boundary(:,1), 'k', 'LineWidth', 0.5); % draw boundary
end

figure(102);
frame=1.67/2; % Hz
t=[0:1:size(NeurTT,2)-1].*1/frame;
N_Profile=1:1:size(NeurTT,1);
h=subplot(121);
set(h,'position',[0.1 0.15,0.64 0.8]);
imagesc(t,N_Profile,NeurTT);ax = gca; 
ax.FontSize = 14;
ylabel('Profile number','Fontsize',15);
xlabel('Times (s)','Fontsize',15);
h = colorbar; 
h.Position = [0.77, 0.22, 0.03, 0.65];
caxis([0 8])

%% Ga_flu traces
n_pc=34;% depth decoding range of Dvge-TPM
frame=1.67/2; % Hz
t=[0:1:size(NeurTT,2)-1].*1/frame;
depths=max(n_pc-final_combinations(:,6),1);
[depths, sort_idx] = sort(depths);
signals = NeurTT(sort_idx,:);

num_neurons = size(NeurTT,1);
num_time_points =size(NeurTT,2); 

x=t;

figure;
set(gcf, 'Color', 'w');
hold on;
cmap = colormap('jet'); 
depth_norm = mat2gray(depths,[0 n_pc]); 

for ii=1:num_neurons
   base_color = interp1(linspace(0, 1, size(cmap, 1)), cmap, depth_norm(ii));
   y=mat2gray(signals(ii,:),[min(signals(ii,:)) 2.7])-(ii-1)*0.9;
    plot(x,y, 'LineWidth', 1.5,'color',base_color);
    hold on;
    plot(x,y,'k','LineWidth',0.5);
end
ylim([-66 2]);axis off;
