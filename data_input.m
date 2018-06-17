%% clear window
clc; close all; 
addpath (genpath ('v2-1 data'));

%% load in V2 data
for i=1:length(psth)
    if isnan(psth(i,1))
        psth(i,1)=0;
    end
end

mov = loadimfile ('v2-1 data/V2Data05/V2Data5/DGrat/r0345/onegrat-apr-64.index60.pixel.imsm');
figure;
subplot(2,2,1);imagesc(mov(:,:,100));
subplot(2,2,2);imagesc(mov(:,:,600));
mov_compressed = imresize3(mov, [32 32 8000]); 
% psth=imresize(psth,[800,1],'bilinear');
psth_normalised = psth./ max(psth);
subplot(2,2,3:4);plot(psth_normalised);

figure; 
subplot(2,2,1);imagesc(mov_compressed(1:20,1:20,100));
title('input stimuli at time bin 150');
subplot(2,2,2);imagesc(mov_compressed(1:20,1:20,670));
title('input stimuli at time bin 670');
subplot(2,2,3:4);plot(psth_normalised);
xlabel('time bins'); ylabel('spike probability');
title('spike response');

mov_medium= mov_compressed(:,:,1:4000);
psth_medium  = psth_normalised(1:4000,1); 

mov_short= mov_compressed(:,:,1:1000);
psth_short  = psth_normalised(1:1000,1); 

%% dissect auditory data
% [Ndim, Tsample]=size(stimulus_compressed);
[Ndim, Tsample]=size(tstimulus);
pixel=zeros(16,20,Tsample);
Nlags=20;
for i=Nlags:Tsample
    pixel(:,:,i) = ...
        tstimulus(1:Ndim,i-Nlags+1:i);
%         stimulus_compressed(1:Ndim,i-Nlags+1:i)
        
end

figure; 
subplot(2,2,1); imagesc(pixel(:,:,1000));
% title('Time-Lagged Stimuli Input', 'FontSize', 20);
colormap(colormap(jet(256)));
xlabel('Time(sec)', 'FontSize', 20)
ylabel('Freq (Hz)', 'FontSize', 20)
subplot(2,2,2); imagesc(pixel(:,:,2000));
colormap(colormap(jet(256)));
xlabel('Time(sec)', 'FontSize', 20)
ylabel('Freq (Hz)', 'FontSize', 20)
subplot(2,2,3); imagesc(pixel(:,:,3000));
colormap(colormap(jet(256)));
xlabel('Time(sec)', 'FontSize', 20)
ylabel('Freq (Hz)', 'FontSize', 20)
subplot(2,2,4); imagesc(pixel(:,:,4000));
colormap(colormap(jet(256)));
xlabel('Time(sec)', 'FontSize', 20)
ylabel('Freq (Hz)', 'FontSize', 20)

psth= response_compressed';

%% gabor filter
%{
RGB = mov(:,:,1500);
I = rgb2gray(RGB);
wavelength = 4;
orientation = 90;
[mag,phase] = imgaborfilt(mov,wavelength,orientation);

figure;
subplot(1,3,1);
imshow(I);
title('Original Image');
subplot(1,3,2);
imshow(mag,[])
title('Gabor magnitude');
subplot(1,3,3);
imshow(phase,[]);
title('Gabor phase');
%}
