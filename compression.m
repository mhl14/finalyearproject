function [stimulus_compressed, response_compressed, response_compressed_nor, Nsamples_compressed, respNsamples_compressed, Ndim_compressed]=compression(stimulus, response, freq_compress, time_compress)
%% clear data
clc; 

%% 
[Ndim, Nsamples]=size(stimulus);
for i=1:Ndim
    for j=1:Nsamples
        if stimulus(i,j)== -Inf 
            stimulus(i,j)=0; 
        end
    end
end

%% Reshape and compress spectrogram of stimuli
% stimulus = spec1(1:end, 1:end); %sometimes for slow songs, spec(1:end, 1:61953); 
stimulus = stimulus-repmat(mean(stimulus,2),[1, Nsamples]);
stimulus = stimulus./repmat(std(stimulus,0,2),[1, Nsamples]);
% freq_compress = 4; %decrease compression for lower freq song
% time_compress = 11.5; %5.74 for fast, 11.5 for normal, 17.2 for slow
%frequency bin: 64/4 = 16, ~= @1.4kHz
%time bin: 30/3600 ~= @16.7ms
Ndim_compressed = ceil(Ndim/freq_compress);
Nsamples_compressed = ceil(Nsamples/time_compress);
stimulus_compressed = imresize(stimulus, [Ndim_compressed Nsamples_compressed],'bilinear'); 
%{
stimulus = spec2(1:end, 1:end); %sometimes for slow songs, spec(1:end, 1:61953);
[Ndim, Nsamples]=size(stimulus); 
stimulus = stimulus-repmat(mean(stimulus,2),[1, Nsamples]);
stimulus = stimulus./repmat(std(stimulus,0,2),[1, Nsamples]);
freq_compress = 4; %decrease compression for lower freq song
time_compress = 11.5; %5.74 for fast, 11.5 for normal, 17.2 for slow
%frequency bin: 64/4 = 16, ~= @1.4kHz
%time bin: 30/3600 ~= @16.7ms
Ndim_compressed = ceil(Ndim/freq_compress);
Nsamples_compressed = ceil(Nsamples/time_compress);
stimulus_compressed2 = imresize(stimulus, [Ndim_compressed Nsamples_compressed],'bilinear'); 

stimulus = spec3(1:end, 1:end); %sometimes for slow songs, spec(1:end, 1:61953);
[Ndim, Nsamples]=size(stimulus); 
stimulus = stimulus-repmat(mean(stimulus,2),[1, Nsamples]);
stimulus = stimulus./repmat(std(stimulus,0,2),[1, Nsamples]);
freq_compress = 4; %decrease compression for lower freq song
time_compress = 11.5; %5.74 for fast, 11.5 for normal, 17.2 for slow
%frequency bin: 64/4 = 16, ~= @1.4kHz
%time bin: 30/3600 ~= @16.7ms
Ndim_compressed = ceil(Ndim/freq_compress);
Nsamples_compressed = ceil(Nsamples/time_compress);
stimulus_compressed3 = imresize(stimulus, [Ndim_compressed Nsamples_compressed],'bilinear'); 

stimulus = spec4(1:end, 1:end); %sometimes for slow songs, spec(1:end, 1:61953);
[Ndim, Nsamples]=size(stimulus); 
stimulus = stimulus-repmat(mean(stimulus,2),[1, Nsamples]);
stimulus = stimulus./repmat(std(stimulus,0,2),[1, Nsamples]);
freq_compress = 4; %decrease compression for lower freq song
time_compress = 11.5; %5.74 for fast, 11.5 for normal, 17.2 for slow
%frequency bin: 64/4 = 16, ~= @1.4kHz
%time bin: 30/3600 ~= @16.7ms
Ndim_compressed = ceil(Ndim/freq_compress);
Nsamples_compressed = ceil(Nsamples/time_compress);
stimulus_compressed4 = imresize(stimulus, [Ndim_compressed Nsamples_compressed],'bilinear'); 

stimulus = spec5(1:end, 1:end); %sometimes for slow songs, spec(1:end, 1:61953);
[Ndim, Nsamples]=size(stimulus); 
stimulus = stimulus-repmat(mean(stimulus,2),[1, Nsamples]);
stimulus = stimulus./repmat(std(stimulus,0,2),[1, Nsamples]);
freq_compress = 4; %decrease compression for lower freq song
time_compress = 11.5; %5.74 for fast, 11.5 for normal, 17.2 for slow
%frequency bin: 64/4 = 16, ~= @1.4kHz
%time bin: 30/3600 ~= @16.7ms
Ndim_compressed = ceil(Ndim/freq_compress);
Nsamples_compressed = ceil(Nsamples/time_compress);
stimulus_compressed5 = imresize(stimulus, [Ndim_compressed Nsamples_compressed],'bilinear'); 

stimulus = spec6(1:end, 1:end); %sometimes for slow songs, spec(1:end, 1:61953);
[Ndim, Nsamples]=size(stimulus); 
stimulus = stimulus-repmat(mean(stimulus,2),[1, Nsamples]);
stimulus = stimulus./repmat(std(stimulus,0,2),[1, Nsamples]);
freq_compress = 4; %decrease compression for lower freq song
time_compress = 11.5; %5.74 for fast, 11.5 for normal, 17.2 for slow
%frequency bin: 64/4 = 16, ~= @1.4kHz
%time bin: 30/3600 ~= @16.7ms
Ndim_compressed = ceil(Ndim/freq_compress);
Nsamples_compressed = ceil(Nsamples/time_compress);
stimulus_compressed6 = imresize(stimulus, [Ndim_compressed Nsamples_compressed],'bilinear'); 

stimulus = spec7(1:end, 1:end); %sometimes for slow songs, spec(1:end, 1:61953);
[Ndim, Nsamples]=size(stimulus); 
stimulus = stimulus-repmat(mean(stimulus,2),[1, Nsamples]);
stimulus = stimulus./repmat(std(stimulus,0,2),[1, Nsamples]);
freq_compress = 4; %decrease compression for lower freq song
time_compress = 11.5; %5.74 for fast, 11.5 for normal, 17.2 for slow
%frequency bin: 64/4 = 16, ~= @1.4kHz
%time bin: 30/3600 ~= @16.7ms
Ndim_compressed = ceil(Ndim/freq_compress);
Nsamples_compressed = ceil(Nsamples/time_compress);
stimulus_compressed7 = imresize(stimulus, [Ndim_compressed Nsamples_compressed],'bilinear'); 

stimulus = spec8(1:end, 1:end); %sometimes for slow songs, spec(1:end, 1:61953);
[Ndim, Nsamples]=size(stimulus); 
stimulus = stimulus-repmat(mean(stimulus,2),[1, Nsamples]);
stimulus = stimulus./repmat(std(stimulus,0,2),[1, Nsamples]);
freq_compress = 4; %decrease compression for lower freq song
time_compress = 11.5; %5.74 for fast, 11.5 for normal, 17.2 for slow
%frequency bin: 64/4 = 16, ~= @1.4kHz
%time bin: 30/3600 ~= @16.7ms
Ndim_compressed = ceil(Ndim/freq_compress);
Nsamples_compressed = ceil(Nsamples/time_compress);
stimulus_compressed8 = imresize(stimulus, [Ndim_compressed Nsamples_compressed],'bilinear'); 

stimulus = spec9(1:end, 1:end); %sometimes for slow songs, spec(1:end, 1:61953);
[Ndim, Nsamples]=size(stimulus); 
stimulus = stimulus-repmat(mean(stimulus,2),[1, Nsamples]);
stimulus = stimulus./repmat(std(stimulus,0,2),[1, Nsamples]);
freq_compress = 4; %decrease compression for lower freq song
time_compress = 11.5; %5.74 for fast, 11.5 for normal, 17.2 for slow
%frequency bin: 64/4 = 16, ~= @1.4kHz
%time bin: 30/3600 ~= @16.7ms
Ndim_compressed = ceil(Ndim/freq_compress);
Nsamples_compressed = ceil(Nsamples/time_compress);
stimulus_compressed9 = imresize(stimulus, [Ndim_compressed Nsamples_compressed],'bilinear'); 

stimulus = spec10(1:end, 1:end); %sometimes for slow songs, spec(1:end, 1:61953);
[Ndim, Nsamples]=size(stimulus); 
stimulus = stimulus-repmat(mean(stimulus,2),[1, Nsamples]);
stimulus = stimulus./repmat(std(stimulus,0,2),[1, Nsamples]);
freq_compress = 4; %decrease compression for lower freq song
time_compress = 11.5; %5.74 for fast, 11.5 for normal, 17.2 for slow
%frequency bin: 64/4 = 16, ~= @1.4kHz
%time bin: 30/3600 ~= @16.7ms
Ndim_compressed = ceil(Ndim/freq_compress);
Nsamples_compressed = ceil(Nsamples/time_compress);
stimulus_compressed10 = imresize(stimulus, [Ndim_compressed Nsamples_compressed],'bilinear'); 

stimulus = spec11(1:end, 1:end); %sometimes for slow songs, spec(1:end, 1:61953);
[Ndim, Nsamples]=size(stimulus); 
stimulus = stimulus-repmat(mean(stimulus,2),[1, Nsamples]);
stimulus = stimulus./repmat(std(stimulus,0,2),[1, Nsamples]);
freq_compress = 4; %decrease compression for lower freq song
time_compress = 11.5; %5.74 for fast, 11.5 for normal, 17.2 for slow
%frequency bin: 64/4 = 16, ~= @1.4kHz
%time bin: 30/3600 ~= @16.7ms
Ndim_compressed = ceil(Ndim/freq_compress);
Nsamples_compressed = ceil(Nsamples/time_compress);
stimulus_compressed11 = imresize(stimulus, [Ndim_compressed Nsamples_compressed],'bilinear'); 

stimulus = spec12(1:end, 1:end); %sometimes for slow songs, spec(1:end, 1:61953);
[Ndim, Nsamples]=size(stimulus); 
stimulus = stimulus-repmat(mean(stimulus,2),[1, Nsamples]);
stimulus = stimulus./repmat(std(stimulus,0,2),[1, Nsamples]);
freq_compress = 4; %decrease compression for lower freq song
time_compress = 11.5; %5.74 for fast, 11.5 for normal, 17.2 for slow
%frequency bin: 64/4 = 16, ~= @1.4kHz
%time bin: 30/3600 ~= @16.7ms
Ndim_compressed = ceil(Ndim/freq_compress);
Nsamples_compressed = ceil(Nsamples/time_compress);
stimulus_compressed12 = imresize(stimulus, [Ndim_compressed Nsamples_compressed],'bilinear'); 

stimulus = spec13(1:end, 1:end); %sometimes for slow songs, spec(1:end, 1:61953);
[Ndim, Nsamples]=size(stimulus); 
stimulus = stimulus-repmat(mean(stimulus,2),[1, Nsamples]);
stimulus = stimulus./repmat(std(stimulus,0,2),[1, Nsamples]);
freq_compress = 4; %decrease compression for lower freq song
time_compress = 11.5; %5.74 for fast, 11.5 for normal, 17.2 for slow
%frequency bin: 64/4 = 16, ~= @1.4kHz
%time bin: 30/3600 ~= @16.7ms
Ndim_compressed = ceil(Ndim/freq_compress);
Nsamples_compressed = ceil(Nsamples/time_compress);
stimulus_compressed13 = imresize(stimulus, [Ndim_compressed Nsamples_compressed],'bilinear'); 

stimulus = spec14(1:end, 1:end); %sometimes for slow songs, spec(1:end, 1:61953);
[Ndim, Nsamples]=size(stimulus); 
stimulus = stimulus-repmat(mean(stimulus,2),[1, Nsamples]);
stimulus = stimulus./repmat(std(stimulus,0,2),[1, Nsamples]);
freq_compress = 4; %decrease compression for lower freq song
time_compress = 11.5; %5.74 for fast, 11.5 for normal, 17.2 for slow
%frequency bin: 64/4 = 16, ~= @1.4kHz
%time bin: 30/3600 ~= @16.7ms
Ndim_compressed = ceil(Ndim/freq_compress);
Nsamples_compressed = ceil(Nsamples/time_compress);
stimulus_compressed14 = imresize(stimulus, [Ndim_compressed Nsamples_compressed],'bilinear'); 

stimulus = spec15(1:end, 1:end); %sometimes for slow songs, spec(1:end, 1:61953);
[Ndim, Nsamples]=size(stimulus); 
stimulus = stimulus-repmat(mean(stimulus,2),[1, Nsamples]);
stimulus = stimulus./repmat(std(stimulus,0,2),[1, Nsamples]);
freq_compress = 4; %decrease compression for lower freq song
time_compress = 11.5; %5.74 for fast, 11.5 for normal, 17.2 for slow
%frequency bin: 64/4 = 16, ~= @1.4kHz
%time bin: 30/3600 ~= @16.7ms
Ndim_compressed = ceil(Ndim/freq_compress);
Nsamples_compressed = ceil(Nsamples/time_compress);
stimulus_compressed15 = imresize(stimulus, [Ndim_compressed Nsamples_compressed],'bilinear'); 

stimulus = spec16(1:end, 1:end); %sometimes for slow songs, spec(1:end, 1:61953);
[Ndim, Nsamples]=size(stimulus); 
stimulus = stimulus-repmat(mean(stimulus,2),[1, Nsamples]);
stimulus = stimulus./repmat(std(stimulus,0,2),[1, Nsamples]);
freq_compress = 4; %decrease compression for lower freq song
time_compress = 11.5; %5.74 for fast, 11.5 for normal, 17.2 for slow
%frequency bin: 64/4 = 16, ~= @1.4kHz
%time bin: 30/3600 ~= @16.7ms
Ndim_compressed = ceil(Ndim/freq_compress);
Nsamples_compressed = ceil(Nsamples/time_compress);
stimulus_compressed16 = imresize(stimulus, [Ndim_compressed Nsamples_compressed],'bilinear'); 
%}

%% Reshape and compress spike responses
% response = psth1'; %psth(1:61953,1);
% response = response';
[respNdim, respNsamples]=size(response); 
respNsamples_compressed = ceil(respNsamples/time_compress);
response_compressed = imresize(response, [1 respNsamples_compressed],'bilinear');
response_compressed_nor = response_compressed ./max(response_compressed); 
%{
response = psth2'; %psth(1:61953,1);
[respNdim, respNsamples]=size(response); 
respNsamples_compressed = ceil(respNsamples/time_compress);
response_compressed = imresize(response, [1 respNsamples_compressed],'bilinear');
response_compressed2 = response_compressed ./max(response_compressed); 

response = psth3'; %psth(1:61953,1);
[respNdim, respNsamples]=size(response); 
respNsamples_compressed = ceil(respNsamples/time_compress);
response_compressed = imresize(response, [1 respNsamples_compressed],'bilinear');
response_compressed3 = response_compressed ./max(response_compressed); 

response = psth4'; %psth(1:61953,1);
[respNdim, respNsamples]=size(response); 
respNsamples_compressed = ceil(respNsamples/time_compress);
response_compressed = imresize(response, [1 respNsamples_compressed],'bilinear');
response_compressed4 = response_compressed ./max(response_compressed); 

response = psth5'; %psth(1:61953,1);
[respNdim, respNsamples]=size(response); 
respNsamples_compressed = ceil(respNsamples/time_compress);
response_compressed = imresize(response, [1 respNsamples_compressed],'bilinear');
response_compressed5 = response_compressed ./max(response_compressed); 

response = psth6'; %psth(1:61953,1);
[respNdim, respNsamples]=size(response); 
respNsamples_compressed = ceil(respNsamples/time_compress);
response_compressed = imresize(response, [1 respNsamples_compressed],'bilinear');
response_compressed6 = response_compressed ./max(response_compressed); 

response = psth7'; %psth(1:61953,1);
[respNdim, respNsamples]=size(response); 
respNsamples_compressed = ceil(respNsamples/time_compress);
response_compressed = imresize(response, [1 respNsamples_compressed],'bilinear');
response_compressed7 = response_compressed ./max(response_compressed); 

response = psth8'; %psth(1:61953,1);
[respNdim, respNsamples]=size(response); 
respNsamples_compressed = ceil(respNsamples/time_compress);
response_compressed = imresize(response, [1 respNsamples_compressed],'bilinear');
response_compressed8 = response_compressed ./max(response_compressed); 

response = psth9'; %psth(1:61953,1);
[respNdim, respNsamples]=size(response); 
respNsamples_compressed = ceil(respNsamples/time_compress);
response_compressed = imresize(response, [1 respNsamples_compressed],'bilinear');
response_compressed9 = response_compressed ./max(response_compressed); 

response = psth10'; %psth(1:61953,1);
[respNdim, respNsamples]=size(response); 
respNsamples_compressed = ceil(respNsamples/time_compress);
response_compressed = imresize(response, [1 respNsamples_compressed],'bilinear');
response_compressed10 = response_compressed ./max(response_compressed); 

response = psth11'; %psth(1:61953,1);
[respNdim, respNsamples]=size(response); 
respNsamples_compressed = ceil(respNsamples/time_compress);
response_compressed = imresize(response, [1 respNsamples_compressed],'bilinear');
response_compressed11 = response_compressed ./max(response_compressed); 

response = psth12'; %psth(1:61953,1);
[respNdim, respNsamples]=size(response); 
respNsamples_compressed = ceil(respNsamples/time_compress);
response_compressed = imresize(response, [1 respNsamples_compressed],'bilinear');
response_compressed12 = response_compressed ./max(response_compressed); 

response = psth13'; %psth(1:61953,1);
[respNdim, respNsamples]=size(response); 
respNsamples_compressed = ceil(respNsamples/time_compress);
response_compressed = imresize(response, [1 respNsamples_compressed],'bilinear');
response_compressed13 = response_compressed ./max(response_compressed); 

response = psth14'; %psth(1:61953,1);
[respNdim, respNsamples]=size(response); 
respNsamples_compressed = ceil(respNsamples/time_compress);
response_compressed = imresize(response, [1 respNsamples_compressed],'bilinear');
response_compressed14 = response_compressed ./max(response_compressed); 

response = psth15'; %psth(1:61953,1);
[respNdim, respNsamples]=size(response); 
respNsamples_compressed = ceil(respNsamples/time_compress);
response_compressed = imresize(response, [1 respNsamples_compressed],'bilinear');
response_compressed15 = response_compressed ./max(response_compressed); 

response = psth16'; %psth(1:61953,1);
[respNdim, respNsamples]=size(response); 
respNsamples_compressed = ceil(respNsamples/time_compress);
response_compressed = imresize(response, [1 respNsamples_compressed],'bilinear');
response_compressed16 = response_compressed ./max(response_compressed); 
%}

%% plot compressed stimuli and spike repsonse

figure;
subplot(2,1,1); imagesc(stimulus_compressed(:, 1:end)); axis xy;
xlabel('time bins'); ylabel('frequency bins');
colormap(colormap(jet(256)));
title('stimuli');
subplot(2,1,2); plot(response_compressed); 
xlabel('time bins'); ylabel('number of spikes');
title('spike reponses');

end