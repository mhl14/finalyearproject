%% Add Functions and Data folders to search path
addpath(genpath('freq_shifted_song')); 
addpath(genpath('UROP')); 
clc; close all; 

%% Load and prepare spectrogram of stimulus and spike response
toelist = ...
    'response/concat_chan_37_38_electrode_23_29/ss001m_497_17_s1_slow_toe.txt';
start = 0; 
stop = 90; %30 for fast, 60 for normal, 90 for slow
figure; 
[ax, psth1, spec1]= plot_raster_SMI2(toelist, start, stop);

toelist = ...
    'response/concat_chan_41_42_electrode_28_25/ss001m_497_17_s1_m_462_s2prime_slow_toe.txt';
start = 0; 
stop = 90; %30 for fast, 60 for normal, 90 for slow
figure; 
[ax, psth2, spec2]= plot_raster_SMI2(toelist, start, stop);

toelist = ...
    'response/concat_chan_41_42_electrode_28_25/ss001m_462_s2prime_slow_toe.txt';
start = 0; 
stop = 90; %30 for fast, 60 for normal, 90 for slow
figure; 
[ax, psth3, spec3]= plot_raster_SMI2(toelist, start, stop);

toelist = ...
    'response/concat_chan_41_42_electrode_28_25/ss001165_s_12_slow_toe.txt';
start = 0; 
stop = 90; %30 for fast, 60 for normal, 90 for slow
figure; 
[ax, psth4, spec4]= plot_raster_SMI2(toelist, start, stop);

%{
%% Reshape and compress spectrogram of stimuli
stimulus = spec1(1:end, 1:61953); %sometimes for slow songs, spec(1:end, 1:61953);
[Ndim, Nsamples]=size(stimulus); 
stimulus = stimulus-repmat(mean(stimulus,2),[1, Nsamples]);
stimulus = stimulus./repmat(std(stimulus,0,2),[1, Nsamples]);
freq_compress = 4; %decrease compression for lower freq song
time_compress = 17.2; %5.74 for fast, 11.5 for normal, 17.2 for slow
%frequency bin: 64/4 = 16, ~= @1.4kHz
%time bin: 30/3600 ~= @16.7ms
Ndim_compressed = ceil(Ndim/freq_compress);
Nsamples_compressed = ceil(Nsamples/time_compress);
stimulus_compressed1 = imresize(stimulus, [Ndim_compressed Nsamples_compressed],'bilinear'); 
%{
stimulus = spec2(1:end, 1:61953); %sometimes for slow songs, spec(1:end, 1:61953);
[Ndim, Nsamples]=size(stimulus); 
stimulus = stimulus-repmat(mean(stimulus,2),[1, Nsamples]);
stimulus = stimulus./repmat(std(stimulus,0,2),[1, Nsamples]);
freq_compress = 4; %decrease compression for lower freq song
time_compress = 17.2; %5.74 for fast, 11.5 for normal, 17.2 for slow
%frequency bin: 64/4 = 16, ~= @1.4kHz
%time bin: 30/3600 ~= @16.7ms
Ndim_compressed = ceil(Ndim/freq_compress);
Nsamples_compressed = ceil(Nsamples/time_compress);
stimulus_compressed2 = imresize(stimulus, [Ndim_compressed Nsamples_compressed],'bilinear'); 

stimulus = spec3(1:end, 1:61953); %sometimes for slow songs, spec(1:end, 1:61953);
[Ndim, Nsamples]=size(stimulus); 
stimulus = stimulus-repmat(mean(stimulus,2),[1, Nsamples]);
stimulus = stimulus./repmat(std(stimulus,0,2),[1, Nsamples]);
freq_compress = 4; %decrease compression for lower freq song
time_compress = 17.2; %5.74 for fast, 11.5 for normal, 17.2 for slow
%frequency bin: 64/4 = 16, ~= @1.4kHz
%time bin: 30/3600 ~= @16.7ms
Ndim_compressed = ceil(Ndim/freq_compress);
Nsamples_compressed = ceil(Nsamples/time_compress);
stimulus_compressed3 = imresize(stimulus, [Ndim_compressed Nsamples_compressed],'bilinear'); 

stimulus = spec4(1:end, 1:61953); %sometimes for slow songs, spec(1:end, 1:61953);
[Ndim, Nsamples]=size(stimulus); 
stimulus = stimulus-repmat(mean(stimulus,2),[1, Nsamples]);
stimulus = stimulus./repmat(std(stimulus,0,2),[1, Nsamples]);
freq_compress = 4; %decrease compression for lower freq song
time_compress = 17.2; %5.74 for fast, 11.5 for normal, 17.2 for slow
%frequency bin: 64/4 = 16, ~= @1.4kHz
%time bin: 30/3600 ~= @16.7ms
Ndim_compressed = ceil(Ndim/freq_compress);
Nsamples_compressed = ceil(Nsamples/time_compress);
stimulus_compressed4 = imresize(stimulus, [Ndim_compressed Nsamples_compressed],'bilinear'); 
%}

%% Reshape and compress spike responses
response = psth1(1:61953,1)'; %psth(1:61953,1);
[respNdim, respNsamples]=size(response); 
respNsamples_compressed = ceil(respNsamples/time_compress);
response_compressed = imresize(response, [1 respNsamples_compressed],'bilinear');
response_compressed1 = response_compressed ./max(response_compressed); 
%{
response = psth2(1:61953,1)'; %psth(1:61953,1);
[respNdim, respNsamples]=size(response); 
respNsamples_compressed = ceil(respNsamples/time_compress);
response_compressed = imresize(response, [1 respNsamples_compressed],'bilinear');
response_compressed2 = response_compressed ./max(response_compressed); 

response = psth3(1:61953,1)'; %psth(1:61953,1);
[respNdim, respNsamples]=size(response); 
respNsamples_compressed = ceil(respNsamples/time_compress);
response_compressed = imresize(response, [1 respNsamples_compressed],'bilinear');
response_compressed3 = response_compressed ./max(response_compressed); 

response = psth4(1:61953,1)'; %psth(1:61953,1);
[respNdim, respNsamples]=size(response); 
respNsamples_compressed = ceil(respNsamples/time_compress);
response_compressed = imresize(response, [1 respNsamples_compressed],'bilinear');
response_compressed4 = response_compressed ./max(response_compressed); 
%}
%}

%% plot compressed stimuli and spike repsonse
%{
figure;
subplot(2,1,1); imagesc(stimulus_compressed1(:, 1:end)); axis xy;
xlabel('time bins'); ylabel('frequency bins');
colormap(colormap(jet(256)));
title('stimuli');
subplot(2,1,2); plot(response_compressed1); 
xlabel('time bins'); ylabel('number of spikes');
title('spike reponses');
%}
