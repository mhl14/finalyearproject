% this file takes all songs and converts them in the spectrographic
% reprepsentation outputing the S,F,T,P matrices.

clear all; close all;

window = 128; nfft = 128; noverlap = window/2;


% cd('/Users/andreikozlov/Work/bird/Analysis/sortedData/st750/sortedData/comb_stims_run_25_51_and_28_34/songs');
cd('/Users/andreikozlov/work/bird/Analysis/sortedData/st750/sortedData/comb_run_0640_2032_2314_1651/control_birds_insample');


filenames = []; %initialize the variable

filenames = [filenames; dir('*.wav')]; %  indicate which files to select

for i = 1:length(filenames); %loop through all files
    
    infile = filenames(i).name;
    song = wavread(infile);
    song_r = resample(song,24000,44100);
    
            infile = strrep(filenames(i).name,'.wav',''); %remove '.wav' from the filename
 
 [S,F,T,P] = spectrogram(song_r,window,noverlap,nfft,24000);
 
varname_S = genvarname(['S_',infile]); varname_F = genvarname(['F_',infile]); varname_T = genvarname(['T_',infile]); varname_P = genvarname(['P_',infile]);
eval([varname_S '= S;']); eval([varname_F '= F;']); eval([varname_T '= T;']); eval([varname_P '= P;']);


end
 

save('song_Ps','P_m*');
clear all;

song_Ps = load('song_Ps.mat');
c = struct2cell(song_Ps);
stimulus = horzcat(c{:});
 
 