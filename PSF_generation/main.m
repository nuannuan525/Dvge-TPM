%% code to generate grad beam with mutil focus
%   Summary of this function goes here
%   basic imaging condition setting
%   wavelength of laser : lambda (um)
%   the focus lenghth : z (um)
%   GradF_phase_1 and GradF_phase_2 contains the focal length、 phase modulator and Binary matrix of mutil foci

clear all;close all;clc;
load('GradF_phase.mat');
nmcoeff(1,1) = 2; 
nmcoeff(2,1) = 0;
ImgStim = FourierOptics();
N=size(GradF_phase_1{3},2);
phase_sum_1=zeros(N,N);
for ii=1:size(GradF_phase_1{1},2)
    nmcoeff(3,1) =GradF_phase_1{1}(ii);  % set coefficient 
    ImgStim = setZernike_nmcoeff(ImgStim,(nmcoeff));
    phase=ImgStim.ObjPl.phase;
    phase_Binary=squeeze(GradF_phase_1{3}(ii,:,:));
    phase_sum_1=(phase-ones(N,N).*GradF_phase_1{2}(ii)).*phase_Binary+phase_sum_1;
end

phase_sum_2=zeros(N,N);
for ii=1:size(GradF_phase_2{1},2)
    nmcoeff(3,1) =GradF_phase_2{1}(ii);  % set coefficient 
    ImgStim = setZernike_nmcoeff(ImgStim,(nmcoeff));
    phase=ImgStim.ObjPl.phase;
    phase_Binary=squeeze(GradF_phase_2{3}(ii,:,:));
    phase_sum_2=(phase-ones(N,N).*GradF_phase_1{2}(ii)).*phase_Binary+phase_sum_2;
end

%% %% show the lightfield
zvect=ImgStim.Img_plane.zlist;
xvect=ImgStim.Img_plane.xx;

ImgStim.ObjPl.phase=phase_sum_1;
ImgStim.ObjPl.LightField = ImgStim.ObjPl.Pupil.*squeeze(sum(GradF_phase_1{3})).*exp(1i.*ImgStim.ObjPl.phase);
[ImgStim,lightfield_1] = ImgStim.getLightField() ;    % complex amplitude
I1=squeeze(lightfield_1).*conj(squeeze(lightfield_1));  % intensity

ImgStim.ObjPl.phase=phase_sum_2;
ImgStim.ObjPl.LightField = ImgStim.ObjPl.Pupil.*squeeze(sum(GradF_phase_2{3})).*exp(1i.*ImgStim.ObjPl.phase);
[ImgStim,lightfield_2] = ImgStim.getLightField() ;    % complex amplitude
I2=squeeze(lightfield_2).*conj(squeeze(lightfield_2));  % intensity

figure(1);
subplot 211;plot(zvect,I1./max(I1));title(['PSF_1']);title(['Intensity']);xlabel('zvect (\mum)');ylabel('xvect (\mum)');ylim([0 1.05])
subplot 212;plot(zvect,I2./max(I2));title(['PSF_2']);title(['Intensity']);xlabel('zvect (\mum)');ylabel('xvect (\mum)');ylim([0 1.05])




