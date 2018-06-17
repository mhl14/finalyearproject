%% Add Functions and Data folders to search path
addpath(genpath('FinalYearProject')); 
clc; close all; 

%% Reshape parameters
a=A_mean;
ydim = 16;
xdim = 16;
Nlags=xdim; Ndim=ydim*xdim;
h=h_mean; 
J= reshape(J_mean,Ndim,Ndim); 
start = 0; 
step = 0.0194363;
nfft = 128;
window = 1224;
overlap = 0.3;
freq_compress = 4;
corre_coef = zeros(1,5);
self_coef = zeros(1,5);
% step = 0.0194363:(0.194363/2-0.0194363)/10:0.194363/2;
% window = 1224:(12240/2-1224)/10:12240/2;

%% Run looop 
% for t=1:1:11
    %% Load and prepare spectrogram of testing stimulus
toelist = ...
    'response/concat_chan_19_20_electrode_8_3/10 reps/ss001m_497_17_s1_toe.txt';
stop = 60; %30 for fast, 60 for normal, 90 for slow
[~, psth, spec]= plot_raster_SMI2(toelist, start, stop, step, nfft, window, overlap);
time_compress = 1; %5.74 for fast, 11.5 for normal, 17.2 for slow
[tstimulus, tresponse, tresponse_nor, ~, Nsample, tNdim]=compression(spec, psth', freq_compress, time_compress); 
close all 

% Fit MNE model with trained parameters to get testing results 
[~, Tsample]=size(tresponse);

tNsamples_compressed = Nsample - (Nlags-1); 
tNdimtotal = tNdim*Nlags;
stim_ = zeros(tNdimtotal, tNsamples_compressed);
for i=1:Nlags
    stim_(tNdim*(i-1)+1:tNdim*i,:) = ...
        tstimulus(:,i:tNsamples_compressed+i-1);
end
stim_=stim_';

for b = 1:1:Tsample
     pSpike = 1./(1+exp(a+stim_*h'+sum(stim_.*(stim_*J),2)));  % Nsamples x 1
end

pSpike_nor = pSpike ./max(pSpike);

figure; subplot(2,1,1); 
newT = start:((stop-start)/length(tstimulus)):stop;
freqs=0:22050/(nfft/2):22050;
imagesc(newT,freqs,tstimulus);axis xy;
colormap(colormap(jet(256)));axis([start stop 0 15000]);
xlabel('time(s)'); 
ylabel('frequency (Hz)');
title('real stimulus');

subplot(2,1,2); hold on; 
xbins = start:((stop-start)/length(tresponse)):stop;
plot(xbins(1,Nlags+1:end),tresponse_nor(:,Nlags:end),'r'); 
xbins2 = start:((stop-start)/length(pSpike_nor)):stop;
plot(xbins2(1,2:end), pSpike_nor,'b'); 
axis([start stop 0 1]);
legend('real spike', 'estimated spike');
xlabel('time(s)'); ylabel('spike probability');
title('real sampled and tested spike response');
%}
%% convolution of both data
% convcoef = corrcoef(tresponse(:,Nlags:end), pSpike')
convcoef_nor1 = corrcoef(tresponse_nor(:,Nlags:end), pSpike_nor');
self_conv1 = self_cc_10reps(toelist, stop, step);
corre_coef(1,1) = convcoef_nor1(1,2);
self_coef(1,1) = self_conv1(1,2);

%% Load and prepare spectrogram of testing stimulus
toelist = ...
    'response/concat_chan_19_20_electrode_8_3/30 reps/ss001m_497_17_s1_2fast_toe.txt';
stop = 30; %30 for fast, 60 for normal, 90 for slow
[~, psth, spec]= plot_raster_SMI2(toelist, start, stop, step, nfft, window, overlap);
time_compress = 0.5; %5.74 for fast, 11.5 for normal, 17.2 for slow
[tstimulus, tresponse, tresponse_nor, ~, Nsample, tNdim]=compression(spec, psth', freq_compress, time_compress); 
close all 

% Fit MNE model with trained parameters to get testing results 
[~, Tsample]=size(tresponse);

tNsamples_compressed = Nsample - (Nlags-1); 
tNdimtotal = tNdim*Nlags;
stim_ = zeros(tNdimtotal, tNsamples_compressed);
for i=1:Nlags
    stim_(tNdim*(i-1)+1:tNdim*i,:) = ...
        tstimulus(:,i:tNsamples_compressed+i-1);
end
stim_=stim_';

for b = 1:1:Tsample
     pSpike2 = 1./(1+exp(a+stim_*h'+sum(stim_.*(stim_*J),2)));  % Nsamples x 1
end

pSpike_nor2 = pSpike2 ./max(pSpike2);

figure; subplot(2,1,1); 
newT = start:((stop-start)/length(tstimulus)):stop;
freqs=0:22050/(nfft/2):22050;
imagesc(newT,freqs,tstimulus);axis xy;
colormap(colormap(jet(256)));axis([start stop 0 15000]);
xlabel('time(s)'); 
ylabel('frequency (Hz)');
title('real stimulus');

subplot(2,1,2); hold on; 
xbins = start:((stop-start)/length(tresponse)):stop;
plot(xbins(1,Nlags+1:end),tresponse_nor(:,Nlags:end),'r'); 
xbins2 = start:((stop-start)/length(pSpike_nor2)):stop;
plot(xbins2(1,2:end), pSpike_nor2,'b'); 
axis([start stop 0 1]);
legend('real spike', 'estimated spike');
xlabel('time(s)'); ylabel('spike probability');
title('real sampled and tested spike response');
%}
%% convolution of both data
% convcoef = corrcoef(tresponse(:,Nlags:end), pSpike')
convcoef_nor2 = corrcoef(tresponse_nor(:,Nlags:end), pSpike_nor2');
self_conv2 = self_cc(toelist, stop, step);
corre_coef(1,2) = convcoef_nor2(1,2);
self_coef(1,2) = self_conv2(1,2);
%{
%% Load and prepare spectrogram of testing stimulus
toelist = ...
    'response/concat_chan_19_20_electrode_8_3/30 reps/ss001m_497_17_s1_shift_1_toe.txt';
stop = 60; %30 for fast, 60 for normal, 90 for slow
[~, psth, spec]= plot_raster_SMI2(toelist, start, stop, step, nfft, window, overlap);
time_compress = 1; %5.74 for fast, 11.5 for normal, 17.2 for slow
[tstimulus, tresponse, tresponse_nor, ~, Nsample, tNdim]=compression(spec, psth', freq_compress, time_compress); 
close all 

%% Fit MNE model with trained parameters to get testing results 
[~, Tsample]=size(tresponse);

tNsamples_compressed = Nsample - (Nlags-1); 
tNdimtotal = tNdim*Nlags;
stim_ = zeros(tNdimtotal, tNsamples_compressed);
for i=1:Nlags
    stim_(tNdim*(i-1)+1:tNdim*i,:) = ...
        tstimulus(:,i:tNsamples_compressed+i-1);
end
stim_=stim_';

for b = 1:1:Tsample
     pSpike3 = 1./(1+exp(a+stim_*h'+sum(stim_.*(stim_*J),2)));  % Nsamples x 1
end

pSpike_nor3 = pSpike3 ./max(pSpike3);

figure; subplot(2,1,1); 
newT = start:((stop-start)/length(tstimulus)):stop;
freqs=0:22050/(nfft/2):22050;
imagesc(newT,freqs,tstimulus);axis xy;
colormap(colormap(jet(256)));axis([start stop 0 15000]);
xlabel('time(s)'); 
ylabel('frequency (Hz)');
title('real stimulus');

subplot(2,1,2); hold on; 
xbins = start:((stop-start)/length(tresponse)):stop;
plot(xbins(1,Nlags+1:end),tresponse_nor(:,Nlags:end),'r'); 
xbins2 = start:((stop-start)/length(pSpike_nor3)):stop;
plot(xbins2(1,2:end), pSpike_nor3,'b'); 
axis([start stop 0 1]);
legend('real spike', 'estimated spike');
xlabel('time(s)'); ylabel('spike probability');
title('real sampled and tested spike response');
%}

%% convolution of both data
% convcoef = corrcoef(tresponse(:,Nlags:end), pSpike')
convcoef_nor3 = corrcoef(tresponse_nor(:,Nlags:end), pSpike_nor3');
self_conv3 = self_cc(toelist, stop, step);
corre_coef(1,3) = convcoef_nor3(1,2);
self_coef(1,3) = self_conv3(1,2);

%% Load and prepare spectrogram of testing stimulus
toelist = ...
    'response/concat_chan_19_20_electrode_8_3/30 reps/ss001m_497_17_s1_shift_6_toe.txt';
stop = 60; %30 for fast, 60 for normal, 90 for slow
[~, psth, spec]= plot_raster_SMI2(toelist, start, stop, step, nfft, window, overlap);
[tstimulus, tresponse, tresponse_nor, ~, Nsample, tNdim]=compression(spec, psth', freq_compress, time_compress); 
close all 

%% Fit MNE model with trained parameters to get testing results 
[~, Tsample]=size(tresponse);

tNsamples_compressed = Nsample - (Nlags-1); 
tNdimtotal = tNdim*Nlags;
stim_ = zeros(tNdimtotal, tNsamples_compressed);
for i=1:Nlags
    stim_(tNdim*(i-1)+1:tNdim*i,:) = ...
        tstimulus(:,i:tNsamples_compressed+i-1);
end
stim_=stim_';

for b = 1:1:Tsample
     pSpike4 = 1./(1+exp(a+stim_*h'+sum(stim_.*(stim_*J),2)));  % Nsamples x 1
end

pSpike_nor4 = pSpike4 ./max(pSpike4);

figure; subplot(2,1,1); 
newT = start:((stop-start)/length(tstimulus)):stop;
freqs=0:22050/(nfft/2):22050;
imagesc(newT,freqs,tstimulus);axis xy;
colormap(colormap(jet(256)));axis([start stop 0 15000]);
xlabel('time(s)'); 
ylabel('frequency (Hz)');
title('real stimulus');

subplot(2,1,2); hold on; 
xbins = start:((stop-start)/length(tresponse)):stop;
plot(xbins(1,Nlags+1:end),tresponse_nor(:,Nlags:end),'r'); 
xbins2 = start:((stop-start)/length(pSpike_nor4)):stop;
plot(xbins2(1,2:end), pSpike_nor4,'b'); 
axis([start stop 0 1]);
legend('real spike', 'estimated spike');
xlabel('time(s)'); ylabel('spike probability');
title('real sampled and tested spike response');
%}

%% convolution of both data
% convcoef = corrcoef(tresponse(:,Nlags:end), pSpike')
convcoef_nor4 = corrcoef(tresponse_nor(:,Nlags:end), pSpike_nor4');
self_conv4 = self_cc(toelist, stop, step);
corre_coef(1,4) = convcoef_nor4(1,2);
self_coef(1,4) = self_conv4(1,2);
%}
%% Load and prepare spectrogram of testing stimulus
toelist = ...
    'response/concat_chan_19_20_electrode_8_3/30 reps/ss001m_497_17_s1_slow_toe.txt';
stop = 90; %30 for fast, 60 for normal, 90 for slow
[~, psth, spec]= plot_raster_SMI2(toelist, start, stop, step, nfft, window, overlap);
time_compress = 2; %5.74 for fast, 11.5 for normal, 17.2 for slow
[tstimulus, tresponse, tresponse_nor, ~, Nsample, tNdim]=compression(spec, psth', freq_compress, time_compress); 
close all 

% Fit MNE model with trained parameters to get testing results 
[Tdim, Tsample]=size(tresponse);

tNsamples_compressed = Nsample - (Nlags-1); 
tNdimtotal = tNdim*Nlags;
stim_ = zeros(tNdimtotal, tNsamples_compressed);
for i=1:Nlags
    stim_(tNdim*(i-1)+1:tNdim*i,:) = ...
        tstimulus(:,i:tNsamples_compressed+i-1);
end
stim_=stim_';

for b = 1:1:Tsample
     pSpike5 = 1./(1+exp(a+stim_*h'+sum(stim_.*(stim_*J),2)));  % Nsamples x 1
end

pSpike_nor5 = pSpike5 ./max(pSpike5);

figure; subplot(2,1,1); 
newT = start:((stop-start)/length(tstimulus)):stop;
freqs=0:22050/(nfft/2):22050;
imagesc(newT,freqs,tstimulus);axis xy;
colormap(colormap(jet(256)));axis([start stop 0 15000]);
xlabel('time(s)'); 
ylabel('frequency (Hz)');
title('real stimulus');

subplot(2,1,2); hold on; 
xbins = start:((stop-start)/length(tresponse)):stop;
plot(xbins(1,Nlags+1:end),tresponse_nor(:,Nlags:end),'r'); 
xbins2 = start:((stop-start)/length(pSpike_nor5)):stop;
plot(xbins2(1,2:end), pSpike_nor5,'b'); 
axis([start stop 0 1]);
legend('real spike', 'estimated spike');
xlabel('time(s)'); ylabel('spike probability');
title('real sampled and tested spike response');
%}

%% convolution of both data
% convcoef = corrcoef(tresponse(:,Nlags:end), pSpike')
convcoef_nor5 = corrcoef(tresponse_nor(:,Nlags:end), pSpike_nor5');
self_conv5 = self_cc(toelist, stop, step);
corre_coef(1,5) = convcoef_nor5(1,2);
self_coef(1,5) = self_conv5(1,2);

% end 

%% visualise % data in bar chart
song = categorical({'m\_497\_17\_s1','m\_497\_17\_s1\_2fast','m\_497\_17\_s1\_shift\_1', ...
    'm\_497\_17\_s1\_shift\_6', 'm\_497\_17\_s1\_slow'});
figure; 
percdata = zeros(4,5);
for i=1:5
%     for j= 1:11
   percdata(1,i) = corre_coef(1,i).*100./self_coef(1,i);
%     end
end 
bar(song,percdata');
title('percentage correlation coefficient with different stimuli');
ylabel('% correlation coefficient');
% legend('19.4ms','27.2ms','35ms','42.8ms','50.5ms','58.3ms','66.1ms','73.9ms','81.6ms', '89.4ms','97.2ms');
