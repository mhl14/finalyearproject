%% Add Functions and Data folders to search path
addpath(genpath('FinalYearProject')); 
clc; close all; 
clear all;

%% Load and prepare spectrogram of stimulus and spike response
start = 0; 
stop = 60; %30 for fast, 60 for normal, 90 for slow
step = 0.0194363;
nfft = 128;
window = 1224;
overlap = 0.3;

toelist = ...
    'response/concat_chan_22_23_electrode_10_7/ss001m_211_s_17_toe.txt';
figure; 
[~, psth1, spec1]= plot_raster_SMI2(toelist, start, stop, step, nfft, window, overlap);

toelist = ...
    'response/concat_chan_22_23_electrode_10_7/ss001m_462_S2prime_m_497_17_s1_toe.txt';
figure; 
[~, psth2, spec2]= plot_raster_SMI2(toelist, start, stop, step, nfft, window, overlap);

toelist = ...
    'response/concat_chan_22_23_electrode_10_7/ss001m_462_s2prime_toe.txt';
figure; 
[~, psth3, spec3]= plot_raster_SMI2(toelist, start, stop, step, nfft, window, overlap);

toelist = ...
    'response/concat_chan_22_23_electrode_10_7/ss001165_s_12_toe.txt';
figure; 
[~, psth4, spec4]= plot_raster_SMI2(toelist, start, stop, step, nfft, window, overlap);
%{
toelist = ...
    'response/concat_chan_22_23_electrode_10_7/ss001m_462_S2prime_min100ms_m_497_17_s1_toe.txt';
figure; 
[~, psth5, spec5]= plot_raster_SMI2(toelist, start, stop, step, nfft, window, overlap);

toelist = ...
    'response/concat_chan_22_23_electrode_10_7/ss001165_s_12_m_462_S2prime_m_497_17_s1_toe.txt';
figure; 
[~, psth6, spec6]= plot_raster_SMI2(toelist, start, stop, step, nfft, window, overlap);

toelist = ...
    'response/concat_chan_22_23_electrode_10_7/ss001165_s_12_toe.txt';
figure; 
[~, psth7, spec7]= plot_raster_SMI2(toelist, start, stop, step, nfft, window, overlap);

% start = 0; 
% stop = 30; %30 for fast, 60 for normal, 90 for slow
% step = 0.0194363;
% nfft = 128;
% window = 1224;
% overlap = 0.3;

toelist = ...
    'response/concat_chan_22_23_electrode_10_7/ss001m_165_s46_toe.txt';
figure; 
[~, psth8, spec8]= plot_raster_SMI2(toelist, start, stop, step, nfft, window, overlap);

toelist = ...
    'response/concat_chan_22_23_electrode_10_7/ss001m_165_s48_toe.txt';
figure; 
[~, psth9, spec9]= plot_raster_SMI2(toelist, start, stop, step, nfft, window, overlap);

toelist = ...
    'response/concat_chan_22_23_electrode_10_7/ss001m_165_s56_toe.txt';
figure; 
[~, psth10, spec10]= plot_raster_SMI2(toelist, start, stop, step, nfft, window, overlap);

toelist = ...
    'response/concat_chan_22_23_electrode_10_7/ss001m_165_s58_toe.txt';
figure; 
[~, psth11, spec11]= plot_raster_SMI2(toelist, start, stop, step, nfft, window, overlap);

toelist = ...
    'response/concat_chan_22_23_electrode_10_7/ss001m_165_s59_toe.txt';
figure; 
[~, psth12, spec12]= plot_raster_SMI2(toelist, start, stop, step, nfft, window, overlap);

% start = 0; 
% stop = 90; %30 for fast, 60 for normal, 90 for slow
% step = 0.0194363;
% nfft = 128;
% window = 1224;
% overlap = 0.3;

toelist = ...
    'response/concat_chan_22_23_electrode_10_7/ss001m_165_s60_toe.txt';
figure; 
[~, psth13, spec13]= plot_raster_SMI2(toelist, start, stop, step, nfft, window, overlap);

toelist = ...
    'response/concat_chan_22_23_electrode_10_7/ss001m_165_s63_toe.txt';
figure; 
[~, psth14, spec14]= plot_raster_SMI2(toelist, start, stop, step, nfft, window, overlap);

toelist = ...
    'response/concat_chan_22_23_electrode_10_7/ss001m_462_s2prime_toe.txt';
figure; 
[~, psth15, spec15]= plot_raster_SMI2(toelist, start, stop, step, nfft, window, overlap);

toelist = ...
    'response/concat_chan_22_23_electrode_10_7/ss001m_497_20_4_9_s1_toe.txt';
figure; 
[~, psth16, spec16]= plot_raster_SMI2(toelist, start, stop, step, nfft, window, overlap);
%}
%% Load spectrogram and spike responses 
% load('results\6songs_data_input.mat');

response = [psth1' psth2' psth3' psth4' ];% ...
%     psth5' psth6' psth7' psth8' ...
%      psth9' psth10' psth11' psth12' ...
%     psth13' psth14' psth15' psth16'

stimulus = [spec1 spec2 spec3 spec4 ]; % ...
%     spec5 spec6 spec7 spec8 ...
%      spec9 spec10 spec11 spec12 ...
%     spec13 spec14 spec15 spec16
    
%% MNE
freq_compress = 4;
time_compress = 4; %5.74 for fast, 11.5 for normal, 17.2 for slow
[stimulus_compressed, response_compressed]=compression(stimulus, response, freq_compress, time_compress);
[A_mean, J_mean, h_mean, V, D, MNE_params, eigenvalues_sorted] = RUN_MNE_auditory(stimulus_compressed, response_compressed);
