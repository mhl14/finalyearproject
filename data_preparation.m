%% Add Functions and Data folders to search path
addpath(genpath('Functions')); 
addpath(genpath('Data')); 
clc; close all; 
clear all; 

%% Load and prepare spectrogram of stimulus
stimulus_path = 'Data/Stimuli/songs_workspace.mat'; 
stim_file = load(stimulus_path);

window = double(stim_file.window);
nfft = double(stim_file.nfft);
noverlap = double(stim_file.noverlap);

stimulus1 = double(stim_file.m_211_s_08_song_60sec_24kHz);
stimulus2 = double(stim_file.m_211_s_16_song_60sec_24kHz);
stimulus3 = double(stim_file.m_211_s_17_song_60sec_24kHz);
stimulus4 = double(stim_file.m_211_s_25_song_60sec_24kHz);

%if the song comes in wavefront file, we have to convert into spectrogram.
[stimulus1] = spectrogram(stimulus1,window,noverlap,nfft,24e3);
% Normalise and discard phase information (i.e.: remove imaginary component)
stimulus1 = 20*log10(stimulus1); %get power density spectrogram
% stimulus1 = stimulus1./max(stimulus1); %normalize 
% Normalise and discard phase information (i.e.: remove imaginary component)
stimulus1= real(stimulus1);
stimulus1 = stimulus1(2:end,:); %Remoce DC component, remove first dimension
% imagesc(stimulus1(2:64,1:end))

[stimulus2]= spectrogram(stimulus2,window,noverlap,nfft,24e3);
% Normalise and discard phase information (i.e.: remove imaginary component)
stimulus2 = 20*log10(stimulus2); %get power density spectrogram
% stimulus2 = stimulus2./max(stimulus2); %normalize 
% Normalise and discard phase information (i.e.: remove imaginary component)
stimulus2= real(stimulus2);
stimulus2 = stimulus2(2:end,:); %Remoce DC component, remove first dimension
% imagesc(stimulus2(2:128,3000:5000))

[stimulus3]= spectrogram(stimulus3,window,noverlap,nfft,24e3);
% Normalise and discard phase information (i.e.: remove imaginary component)
stimulus3 = 20*log10(stimulus3); %get power density spectrogram
% stimulus3 = stimulus3./max(stimulus3); %normalize 
% Normalise and discard phase information (i.e.: remove imaginary component)
stimulus3= real(stimulus3);
stimulus3 = stimulus3(2:end,:); %Remoce DC component, remove first dimension
% imagesc(stimulus3(2:128,3000:5000))

[stimulus4]= spectrogram(stimulus4,window,noverlap,nfft,24e3);
% Normalise and discard phase information (i.e.: remove imaginary component)
stimulus4 = 20*log10(stimulus4); %get power density spectrogram
% stimulus4 = stimulus4./max(stimulus4); %normalize 
% Normalise and discard phase information (i.e.: remove imaginary component)
stimulus4= real(stimulus4);
stimulus4 = stimulus4(2:end,:); %Remoce DC component, remove first dimension
% imagesc(stimulus4(2:128,3000:5000))

[stimulus5,Fs]= audioread('Data/wetransfer-8ade3c/m_462_s2prime.wav');
stimulus5 = resample(stimulus5,24000,Fs);
[stimulus5]= spectrogram(stimulus5,window,noverlap,nfft,24e3);
% Normalise and discard phase information (i.e.: remove imaginary component)
stimulus5 = 20*log10(stimulus5); %get power density spectrogram
% stimulus5 = stimulus5./max(stimulus5); %normalize 
% Normalise and discard phase information (i.e.: remove imaginary component)
stimulus5= real(stimulus5);
stimulus5 = stimulus5(2:end,:); %Remoce DC component, remove first dimension
% imagesc(stimulus5(2:128,3000:5000))

[stimulus6,Fs]= audioread('Data/wetransfer-8ade3c/m_497_17_s1.wav');
stimulus6 = resample(stimulus6,24000,Fs);
[stimulus6]= spectrogram(stimulus6,window,noverlap,nfft,24e3);
% Normalise and discard phase information (i.e.: remove imaginary component)
stimulus6 = 20*log10(stimulus6); %get power density spectrogram
% stimulus6 = stimulus6./max(stimulus6); %normalize 
% Normalize and discard phase information (i.e.: remove imaginary component)
stimulus6= real(stimulus6);
stimulus6 = stimulus6(2:end,:); %Remoce DC component, remove first dimension
% imagesc(stimulus6(2:128,3000:5000))

% Stimulus is a matrix with size NdimxNsamples
stimulus= [stimulus1 stimulus2 stimulus3 stimulus4 stimulus5 stimulus6];

%% Load and prepare response
% Point to spike trace
% Format is .mat file with binary response vector
response_path = 'Data/wetransfer-8ade3c/spike_vector.mat'; 
response_file = load(response_path);

response1 = double(response_file.spike_vector_m_211_s_08_24kHz');
response2 = double(response_file.spike_vector_m_211_s_16_24kHz');
response3 = double(response_file.spike_vector_m_211_s_17_24kHz');
response4 = double(response_file.spike_vector_m_211_s_25_24kHz');
response5 = double(response_file.spike_vector_m_462_s2prime_24kHz');
response6 = double(response_file.spike_vector_m_497_17_s1_24kHz');

% Respone is of size 1xNsamples
response = [response1 response2 response3 response4 response5 response6];
% imagesc(response(1,1:end))

%% plot stimuli and spike repsonse
%{
figure; 
subplot(2,1,1); imagesc(stimulus(:, 1:end)); axis xy;
xlabel('time bins'); ylabel('frequency bins');
title('stimuli');
subplot(2,1,2); plot(response); 
xlabel('time bins'); ylabel('number of spikes');
title('spike reponses');
%}

%% Reshape and compress spectrogram of stimuli
% Determine dimensions of spectrogram
% ydim = number of frequency bands
% xdim = number of time bins
[Ndim, Nsamples]=size(stimulus); 

% The stimulus units is changed to be z-score
% (zero mean, normalized by standard deviation),
% but could be put into your own unit system (or 
% maybe your stimulus file is already saved
% in the proper units).
stimulus = stimulus-repmat(mean(stimulus,2),[1, Nsamples]);
stimulus = stimulus./repmat(std(stimulus,0,2),[1, Nsamples]);

% Choosing compression for spectrogram
% Without compression, ultrasonic spectrogram is usually too large for analysis
freq_compress = 4; %decrease compression for lower freq song
time_compress = 5; %long song compress more

Ndim_compressed = ceil(Ndim/freq_compress);
Nsamples_compressed = ceil(Nsamples/time_compress);
%ceil rouns up element to nearest integers towards infinity 

% Using 'imresize' since spectrogram is just an image
%compress spectrogram 
stimulus_compressed = imresize(stimulus, [Ndim_compressed Nsamples_compressed],'bilinear'); 
%imagesc(stimulus_compressed);
% stimulus_size=size(stimulus_compressed);

%% Reshape and compress spike responses
% Compress number of time bins to be equal to copmressed stimulus
response_compressed = imresize(response, [1 Nsamples_compressed],'bilinear');

%Normalise response values to 0-1
response_compressed = response_compressed ./max(response_compressed); 
% response_compressed2 = (response_compressed - min(response_compressed))/...
%     ( max(response_compressed) - min(response_compressed));

% response_compressed = smooth(response_compressed);
% response_compressed = response_compressed.'; 
%}
%% plot compressed stimuli and spike repsonse
figure;
subplot(2,1,1); imagesc(stimulus_compressed(:, 1:end)); axis xy;
% x=linspace(0,18000); 
% xticks([0 20 40 60 80 100 120 140 160 180 200 220 240 260 280 300]);
xlabel('time bins'); ylabel('frequency bins');
title('stimuli');
subplot(2,1,2); plot(response_compressed); axis([0 Nsamples_compressed 0 1]);
% x=linspace(0,18000); 
% xticks([0 20 40 60 80 100 120 140 160 180 200 220 240 260 280 300]);
xlabel('time bins'); ylabel('number of spikes');
title('spike reponses');
