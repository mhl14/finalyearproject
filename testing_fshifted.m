%% Add Functions and Data folders to search path
addpath(genpath('FinalYearProject')); 
clc; close all; 

%% Reshape parameters
a=A_mean;

ydim = MNE_params.Ndim;
xdim = MNE_params.Nlags;
Nlags=xdim; Ndim=ydim*xdim;

h=h_mean; 

J= reshape(J_mean,Ndim,Ndim); 

%% Load and prepare spectrogram of testing stimulus
toelist = ...
    'response/concat_chan_23_24_electrode_7_4/10 reps/ss001m_497_17_s1_toe.txt';
start = 0; 
stop = 60; %30 for fast, 60 for normal, 90 for slow
step = 0.0194363;
nfft = 128;
window = 1224;
overlap = 0.3;
[ax, psth, spec]= plot_raster_SMI2(toelist, start, stop, step, nfft, window, overlap);
freq_compress = 4;
time_compress = 1; %5.74 for fast, 11.5 for normal, 17.2 for slow
[tstimulus, tresponse, tresponse_nor, Nsample, respNsample, tNdim]=compression(spec, psth', freq_compress, time_compress); 

%% Fit MNE model with trained parameters to get testing results 
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
     pSpike = 1./(1+exp(a+stim_*h'+sum(stim_.*(stim_*J),2)));  % Nsamples x 1
end

pSpike_nor = pSpike ./max(pSpike);

figure; subplot(2,1,1); imagesc(tstimulus); axis xy;
xlabel('time bins'); ylabel('frequency bins');
title('real stimulus');

subplot(2,1,2); hold on
plot(tresponse_nor(:,Nlags:end),'r'); plot(pSpike_nor,'b'); 
axis([0 Tsample 0 1]); legend('real spike', 'estimated spike');
xlabel('time bins'); ylabel('number of spikes');
title('real sampled and tested spike response');

% subplot(3,1,3); hold on
% plot(tresponse(:,Nlags:end),'r'); plot(pSpike,'b'); 
% axis([0 Tsample 0 1]); legend('real spike', 'estimated spike');
% xlabel('time bins'); ylabel('number of spikes');
% title('real sampled and tested spike response');

%% convolution of both data
% convcoef = corrcoef(tresponse(:,Nlags:end), pSpike')
convcoef_nor1 = corrcoef(tresponse_nor(:,Nlags:end), pSpike_nor')
[self_conv1] = self_cc_10reps(toelist, stop, step)

%% Load and prepare spectrogram of testing stimulus
toelist = ...
    'response/concat_chan_23_24_electrode_7_4/30 reps/ss001m_497_17_s1_2fast_toe.txt';
stop = 30; %30 for fast, 60 for normal, 90 for slow
[ax, psth, spec]= plot_raster_SMI2(toelist, start, stop, step, nfft, window, overlap);
time_compress = 0.5; %5.74 for fast, 11.5 for normal, 17.2 for slow
[tstimulus, tresponse, tresponse_nor, Nsample, respNsample, tNdim]=compression(spec, psth', freq_compress, time_compress); 

%% Fit MNE model with trained parameters to get testing results 
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
     pSpike2 = 1./(1+exp(a+stim_*h'+sum(stim_.*(stim_*J),2)));  % Nsamples x 1
end

pSpike_nor2 = pSpike2 ./max(pSpike2);

figure; subplot(2,1,1); imagesc(tstimulus); axis xy;
xlabel('time bins'); ylabel('frequency bins');
title('real stimulus');

subplot(2,1,2); hold on
plot(tresponse_nor(:,Nlags:end),'r'); plot(pSpike_nor2,'b'); 
axis([0 Tsample 0 1]); legend('real spike', 'estimated spike');
xlabel('time bins'); ylabel('number of spikes');
title('real sampled and tested spike response');

% subplot(3,1,3); hold on
% plot(tresponse(:,Nlags:end),'r'); plot(pSpike,'b'); 
% axis([0 Tsample 0 1]); legend('real spike', 'estimated spike');
% xlabel('time bins'); ylabel('number of spikes');
% title('real sampled and tested spike response');

%% convolution of both data
% convcoef = corrcoef(tresponse(:,Nlags:end), pSpike')
convcoef_nor2 = corrcoef(tresponse_nor(:,Nlags:end), pSpike_nor2')
[self_conv2] = self_cc(toelist, stop, step)

%% Load and prepare spectrogram of testing stimulus
toelist = ...
    'response/concat_chan_23_24_electrode_7_4/30 reps/ss001m_497_17_s1_shift_1_toe.txt';
stop = 60; %30 for fast, 60 for normal, 90 for slow
[ax, psth, spec]= plot_raster_SMI2(toelist, start, stop, step, nfft, window, overlap);
time_compress = 1; %5.74 for fast, 11.5 for normal, 17.2 for slow
[tstimulus, tresponse, tresponse_nor, Nsample, respNsample, tNdim]=compression(spec, psth', freq_compress, time_compress); 

%% Fit MNE model with trained parameters to get testing results 
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
     pSpike3 = 1./(1+exp(a+stim_*h'+sum(stim_.*(stim_*J),2)));  % Nsamples x 1
end

pSpike_nor3 = pSpike3 ./max(pSpike3);

figure; subplot(2,1,1); imagesc(tstimulus); axis xy;
xlabel('time bins'); ylabel('frequency bins');
title('real stimulus');

subplot(2,1,2); hold on
plot(tresponse_nor(:,Nlags:end),'r'); plot(pSpike_nor3,'b'); 
axis([0 Tsample 0 1]); legend('real spike', 'estimated spike');
xlabel('time bins'); ylabel('number of spikes');
title('real sampled and tested spike response');

% subplot(3,1,3); hold on
% plot(tresponse(:,Nlags:end),'r'); plot(pSpike,'b'); 
% axis([0 Tsample 0 1]); legend('real spike', 'estimated spike');
% xlabel('time bins'); ylabel('number of spikes');
% title('real sampled and tested spike response');

%% convolution of both data
% convcoef = corrcoef(tresponse(:,Nlags:end), pSpike')
convcoef_nor3 = corrcoef(tresponse_nor(:,Nlags:end), pSpike_nor3')
[self_conv3] = self_cc(toelist, stop, step)

%% Load and prepare spectrogram of testing stimulus
toelist = ...
    'response/concat_chan_23_24_electrode_7_4/30 reps/ss001m_497_17_s1_shift_6_toe.txt';
stop = 60; %30 for fast, 60 for normal, 90 for slow
[ax, psth, spec]= plot_raster_SMI2(toelist, start, stop, step, nfft, window, overlap);
[tstimulus, tresponse, tresponse_nor, Nsample, respNsample, tNdim]=compression(spec, psth', freq_compress, time_compress); 

%% Fit MNE model with trained parameters to get testing results 
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
     pSpike4 = 1./(1+exp(a+stim_*h'+sum(stim_.*(stim_*J),2)));  % Nsamples x 1
end

pSpike_nor4 = pSpike4 ./max(pSpike4);

figure; subplot(2,1,1); imagesc(tstimulus); axis xy;
xlabel('time bins'); ylabel('frequency bins');
title('real stimulus');

subplot(2,1,2); hold on
plot(tresponse_nor(:,Nlags:end),'r'); plot(pSpike_nor4,'b'); 
axis([0 Tsample 0 1]); legend('real spike', 'estimated spike');
xlabel('time bins'); ylabel('number of spikes');
title('real sampled and tested spike response');

% subplot(3,1,3); hold on
% plot(tresponse(:,Nlags:end),'r'); plot(pSpike,'b'); 
% axis([0 Tsample 0 1]); legend('real spike', 'estimated spike');
% xlabel('time bins'); ylabel('number of spikes');
% title('real sampled and tested spike response');

%% convolution of both data
% convcoef = corrcoef(tresponse(:,Nlags:end), pSpike')
convcoef_nor4 = corrcoef(tresponse_nor(:,Nlags:end), pSpike_nor3')
[self_conv4] = self_cc(toelist, stop, step)

%% Load and prepare spectrogram of testing stimulus
toelist = ...
    'response/concat_chan_23_24_electrode_7_4/30 reps/ss001m_497_17_s1_slow_toe.txt';
stop = 90; %30 for fast, 60 for normal, 90 for slow
[ax, psth, spec]= plot_raster_SMI2(toelist, start, stop, step, nfft, window, overlap);
time_compress = 2; %5.74 for fast, 11.5 for normal, 17.2 for slow
[tstimulus, tresponse, tresponse_nor, Nsample, respNsample, tNdim]=compression(spec, psth', freq_compress, time_compress); 

%% Fit MNE model with trained parameters to get testing results 
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

figure; subplot(2,1,1); imagesc(tstimulus); axis xy;
xlabel('time bins'); ylabel('frequency bins');
title('real stimulus');

subplot(2,1,2); hold on
plot(tresponse_nor(:,Nlags:end),'r'); plot(pSpike_nor5,'b'); 
axis([0 Tsample 0 1]); legend('real spike', 'estimated spike');
xlabel('time bins'); ylabel('number of spikes');
title('real sampled and tested spike response');

% subplot(3,1,3); hold on
% plot(tresponse(:,Nlags:end),'r'); plot(pSpike,'b'); 
% axis([0 Tsample 0 1]); legend('real spike', 'estimated spike');
% xlabel('time bins'); ylabel('number of spikes');
% title('real sampled and tested spike response');

%% convolution of both data
% convcoef = corrcoef(tresponse(:,Nlags:end), pSpike')
convcoef_nor5 = corrcoef(tresponse_nor(:,Nlags:end), pSpike_nor5')
[self_conv5] = self_cc(toelist, stop, step)
