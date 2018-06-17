%% Add Functions and Data folders to search path
addpath(genpath('Functions')); 
addpath(genpath('Data')); 
addpath(genpath('results'));
clc; close all; 
% clear all;

%% Insert data
% load('results\5songs_training_MNE_results.mat');

%% Load and prepare spectrogram of testing stimulus
tstimulus= stimulus2;

%% Reshape and compress spectrogram of stimuli
% Determine dimensions of spectrogram
% ydim = number of frequency bands
% xdim = number of time bins
[tNdim, tNsamples]=size(tstimulus); 

% The stimulus units is changed to be z-score
% (zero mean, normalized by standard deviation),
% but could be put into your own unit system (or 
% maybe your stimulus file is already saved
% in the proper units).
tstimulus = tstimulus-repmat(mean(tstimulus,2),[1, tNsamples]);
tstimulus = tstimulus./repmat(std(tstimulus,0,2),[1, tNsamples]);

% Choosing compression for spectrogram
% Without compression, ultrasonic spectrogram is usually too large for analysis
freq_compress = 4; %decrease compression for lower freq song
time_compress = 6.25; %long song compress more

tNdim_compressed = ceil(tNdim/freq_compress);
tNsamples_compressed = ceil(tNsamples/time_compress);
%ceil rouns up element to nearest integers towards infinity 

% Using 'imresize' since spectrogram is just an image
%compress spectrogram 
tstimulus_compressed = imresize(tstimulus, [tNdim_compressed tNsamples_compressed],'bilinear'); 
tstimulus=tstimulus_compressed;

%% Load and prepare testing spike responses
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
tresponse = [response1 response2 response3 response4 response5 response6];
% imagesc(response(1,1:end))

%% Reshape and compress spike responses
[respNdim, respNsamples]=size(tresponse); 
respNsamples_compressed = ceil(respNsamples/time_compress);
% Compress number of time bins to be equal to copmressed stimulus
tresponse_compressed = imresize(tresponse, [1 respNsamples_compressed],'bilinear');

% %Normalise response values to 0-1
tresponse_compressed = tresponse_compressed ./max(tresponse_compressed);
% response_size=size(response_compressed); 
% response_compressed = smooth(response_compressed);
tresponse= tresponse_compressed; 

% response_compressed = smooth(response_compressed);
song = 2;
int=respNsamples_compressed/6;
tresponse= tresponse(1,int*(song-1)+1:int*song);

%% Reshape parameters
a=A_mean;

ydim = MNE_params.Ndim;
xdim = MNE_params.Nlags;
Nlags=xdim; Ndim=ydim*xdim;

h=h_mean; 

J= reshape(J_mean,Ndim,Ndim); 

%% Fit MNE model with trained parameters to get testing results 
[Tdim, Tsample]=size(tresponse);

tNsamples_compressed = tNsamples_compressed - (Nlags-1); 
Ndimtotal = tNdim_compressed*Nlags;
stim_ = zeros(Ndimtotal, tNsamples_compressed);
for i=1:Nlags
    stim_(tNdim_compressed*(i-1)+1:tNdim_compressed*i,:) = ...
        tstimulus(:,i:tNsamples_compressed+i-1);
end
stim_=stim_';

for b = 1:1:Tsample
     pSpike = 1./(1+exp(a+stim_*h'+sum(stim_.*(stim_*J),2)));  % Nsamples x 1
end

figure; subplot(2,1,1); imagesc(tstimulus); axis xy;
xlabel('time bins'); ylabel('frequency bins');
title('real stimulus');

subplot(2,1,2); hold on
plot(tresponse,'r'); plot(pSpike,'b'); 
axis([0 Tsample 0 1]); legend('real spike', 'estimated spike');
xlabel('time bins'); ylabel('number of spikes');
title('real sampled and tested spike response');

%% convolution of both data
convcoef = corrcoef(tresponse(:,Nlags:end), pSpike')
