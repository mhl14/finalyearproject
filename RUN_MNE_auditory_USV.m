  %% Add Functions folder to search path
addpath(genpath('Functions'))
addpath(genpath('Data'))
clc; close all; clear all; 

%% Load and prepare spectrogram of stimulus
stim_file = load('Data/Stimuli/USV_stim.mat'); 
stimulus = stim_file.stimulus;

%if the song comes in wavefront file, we have to convert into spectrogram.

% Normalise and discard phase information (i.e.: remove imaginary component)
% stimulus = 20*log10(stimulus);
% stimulus = stimulus./max(stimulus);
stimulus = real(stimulus);

% Determine dimensions of spectrogram
% ydim = number of frequency bands
% xdim = number of time bins
[Ndim, Nsamples]=size(stimulus); %64,503538

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
time_compress = 256; %long song compress more

Ndim_compressed = ceil(Ndim/freq_compress);
Nsamples_compressed = ceil(Nsamples/time_compress);

% Using 'imresize' since spectrogram is just an image
%compress spectrogram 
stimulus_compressed = imresize(stimulus, [Ndim_compressed Nsamples_compressed],'bilinear'); 
% Transpose because MNE function takes stimulus as NsamplesxNfreq
% stimulud_compressed = stimulus_compressed';
% imshow(stimulus_compressed);
stimulus_size=size(stimulus_compressed); %16x984

% stim_train = stimulus(1:round(length(stimulus)*.75),:);
% stim_test = stimulus(round(length(stimulus)*.75)+1:end,:);
% stimulus = stim_train;
% [Ndim, Nsamples]=size(stimulus');

%% MNE setup
MNE_params.Ndim = Ndim_compressed; %number of frequency bands in STRF 
MNE_params.Nlags = 16; %times the length of a bin gives the time length of STRF
% # time lags to use in spatiotemporal models (=1 is just a spatial model)
MNE_params.order = 2;   % order of MNE model to fit (1 or 2)
MNE_params.fittype = 0;   % 0 for regular fitting, 1 for random fitting
MNE_params.Njack = 10; % Number of jackknives to use 

%% Load and prepare spike response

% Point to spike trace
% Format is .mat file with binary response vector
response_path = 'Data/Responses/_binary_spike_trace_trode_1_cluster_no_15.mat';
% response_path = 'Algorithm_Functions/ chan_23_24.mat'; 
response_file = load(response_path);
response = double(response_file.binary_spike_trace);

% Average response over 10 repititions
% Experiment has 10 repititions of the same stimulus played continuously,
% so we assume that the response of the neuron can be averaged of those
% repetitions
num_rep = 10;
response = response(1:end-mod(length(response),num_rep));
response = reshape(response,[num_rep,floor(length(response)/num_rep)]); 
response = mean(response,1);

% Compress number of time bins to be equal to copmressed stimulus
response_compressed = imresize(response, [1 Nsamples_compressed],'bilinear');

% %Normalise response values to 0-1
response_compressed = (response_compressed - min(response_compressed))/...
    ( max(response_compressed) - min(response_compressed));
% response_compressed = smooth(response_compressed);
response_size=size(response_compressed); %1x984

% % % % response = response>0;
% threshold = 0.3;
% R = response;
% [i,j]=find(R<threshold);
% response(i,j) = 0; % an attempt at thresholding
% [i,j]=find(R>threshold);
% response(i,j) = 1; % an attempt at thresholding

% Split into train and test set
% resp_train = response(1:round(length(response)*.75),:);
% resp_test = response(round(length(response)*.75)+1:end,:);
% response = resp_train;

%% Fit MNE model
[J_mean, h_mean]  = MNE(stimulus_compressed, response_compressed, MNE_params);
[V,D] = eig(reshape(J_mean,MNE_params.Nlags*MNE_params.Ndim,MNE_params.Nlags*MNE_params.Ndim));

%% Plot results and misc
fig = figure;
subplot(3,3,1:3)
eigenvalues = diag(D);
[eigenvalues_sorted,index] = sort(eigenvalues);
plot(eigenvalues_sorted,'o');
xlabel('dimensions')
ylabel('eigenvalues')

ydim = MNE_params.Ndim;
xdim = MNE_params.Nlags;
colorbar

subplot(3,3,4)
eig_sorted_1 = V(:,index(1));
imagesc(reshape(eig_sorted_1,ydim,xdim))
axis xy
xlabel('time band')
ylabel('frequency band')
title('largest eigenvalue')

subplot(3,3,5)
eig_sorted_2=V(:,index(2));
imagesc(reshape(eig_sorted_2,ydim,xdim))
axis xy
xlabel('time band')
ylabel('frequency band')
title('second largest eigenvalue')

subplot(3,3,6)
eig_sorted_3=V(:,index(3));
imagesc(reshape(eig_sorted_3,ydim,xdim))
axis xy
xlabel('time band')
ylabel('frequency band')
title('third largest eigenvalue')

subplot(3,3,7)
eig_sorted_end=V(:,index(end));
imagesc(reshape(eig_sorted_end,ydim,xdim))
axis xy
xlabel('time band')
ylabel('frequency band')
title('smallest eigenvalue')

subplot(3,3,8)
eig_sorted_end_1=V(:,index(end-1));
imagesc(reshape(eig_sorted_end_1,ydim,xdim))
axis xy
xlabel('time band')
ylabel('frequency band')
title('second smallest eigenvalue')

subplot(3,3,9)
eig_sorted_end_2=V(:,index(end-2));
imagesc(reshape(eig_sorted_end_2,ydim,xdim))
axis xy
xlabel('time band')
ylabel('frequency band')
title('third smallest eigenvalue')

% figure
% title('Input Stimulus')
% for i = 1:16
%     subplot(4,4,i)
%     show_stim(stimulus,i,ydim,xdim)
%     colormap parula
% end

% figure;
% plot(response_compressed)
% title ('response');

% figure;
% plot(response_file.binned_spike_trace(1:section*32))

% eval('J_random')
% end
