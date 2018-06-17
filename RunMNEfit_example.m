% function RunMNEfit_example(response_path, stimulus)

% clear all; close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Open stimulus and arrange according to # of time lags %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% cd('/Users/andreikozlov/Work/bird/Analysis/sortedData/st750/songs');
% 
% load songs_workspace.mat
% 
% stimulus = [P_m_211_s_16_song_60sec_24kHz P_m_211_s_17_song_60sec_24kHz P_m_211_s_25_song_60sec_24kHz P_m_462_24kHz P_m_462_s2_4_9_24kHz];


% stim_file = load('D:\Sihao\stimulus.mat') %NFFT=128, overlap = NFFT/2
stimulus = stim_file.stimulus;
[num_dim, time_bins]=size(stimulus);

stimulus = stimulus(2:end,:);
stimulus = stimulus(:,1:floor(time_bins/10)); %Stimulus has 10 repetitions
[num_dim, time_bins]=size(stimulus);

stimulus = imresize(stimulus, [(num_dim)/4 time_bins/32],'bilinear');
size(stimulus);
section = length(stimulus);
stimulus = stimulus(:, 1:section);
stimulus = 20*log10(stimulus);
stimulus = real(stimulus);

[Ndim, Nsamples]=size(stimulus);
stimulus=stimulus';
% The stimulus units is changed to be z-score
% (zero mean, normalized by standard deviation),
% but could be put into your own unit system (or 
% maybe your stimulus file is already saved
% in the proper units).
 stimulus = stimulus-repmat(mean(stimulus),[Nsamples,1]);
stimulus = stimulus./repmat(std(stimulus),[Nsamples,1]);



cellnum = 750;
Nf = 16; %number of frequency bands in STRF
Nlags = 16; %times the length of a bin gives the time length of STRF
order   = 2;   % order of MNE model to fit
fittype = 0;   % 0 for regular fitting, 1 for random fitting
njack   = 4;   % # jackknives to run (also determines the size of each jackknives)
%Nlags   = 1;   % # time lags to use in spatiotemporal models (=1 is just a
%spatial model)

% redefine stimulus to include time lags
if Nlags>1
    Nsamples = Nsamples - (Nlags-1); %total length of stimulus minus 19 time bins
    Ndimtotal = Ndim*Nlags; %16x20
    stim_ = zeros(Nsamples,Ndimtotal);
    for i=1:Nlags
        stim_(:,Ndim*(i-1)+1:Ndim*i) = stimulus(i:Nsamples+i-1,:);
    end
else
    stim_ = stimulus;
end
%clear stimulus;
Nd = sqrt(Ndim); % # pixels per side of image


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Open responses and cut according to # of time lags %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% cd('/Users/andreikozlov/Work/bird/Analysis/sortedData/st750/sortedData/chan_19_20');
response_path = 'D:\Sihao\2016-04-05_16-58-04\_binary_spike_trace_trode_5_cluster_no_13.mat';
response_file = load(response_path);
response = double(response_file.binned_spike_trace);
%Average response over 10 repititions
response = response(1:end-mod(length(response),10));

response = reshape(response,[10,floor(length(response)/10)]); 
response = mean(response,1);

% %%%%%%%%%%Use if edges removed from stimuli%%%%%%%%%%
% fs = length(response)/1650 
% resp_cut=[];
% for i = 1:330
%     
%    start = 1 + round(fs*.5)+round(fs*5)*(i-i);
%    stop = start + round(fs*4);
%    
%    resp_temp = response(:,start:stop); 
%    resp_cut = [resp_cut resp_temp];
%   
% end
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% response = [spike_vector_m_211_s_16_song_60sec; spike_vector_m_211_s_17_song_60sec; spike_vector_m_211_s_25_song_60sec; spike_vector_m_462_s2prime; spike_vector_m_462_s2_4_9];

% response = response'; % this depends on how the response file was
%prepared:check the rows-columns

response = imresize(response, [1 time_bins/32],'bilinear');
size(response)
% 
% % response_odd_columns = response(:,1:2:end-1); %do this if the original number of columns is odd
% response_odd_columns = response(:,1:2:end); %if it is even
% response_even_columns = response(:,2:2:end);
% response_averaged = (response_odd_columns+response_even_columns)./2;
% response = response_averaged;
% 
% %and again:
% response_odd_columns = response(:,1:2:end-1);
% response_even_columns = response(:,2:2:end);
% response_averaged = response_odd_columns+response_even_columns;
% response = response_averaged;
% 
% %and again
% response_odd_columns = response(:,1:2:end-1); 
% %response_odd_columns = response(:,1:2:end); 
% response_even_columns = response(:,2:2:end);
% response_averaged = response_odd_columns+response_even_columns;
% response = response_averaged;
% 
% %and again
% response_odd_columns = response(:,1:2:end-1); 
% %response_odd_columns = response(:,1:2:end); 
% response_even_columns = response(:,2:2:end);
% response_averaged = response_odd_columns+response_even_columns;
% response = response_averaged;
% 
% %and again
% response_odd_columns = response(:,1:2:end-1); 
% %response_odd_columns = response(:,1:2:end); 
% response_even_columns = response(:,2:2:end);
% response_averaged = response_odd_columns+response_even_columns;
% response = response_averaged;
% 

response = response(:, 1:section);
response = response';
 
% % %Now we reshuffle response for control;
% [length_resp_1,length_resp_2]=size(response);
% indeks = randperm(length_resp_1);
% response = response(indeks);


% Now we normalize spike counts in each bin by the maximal value for the minimal models (for MID they are not
% normalized

response = response./max(response);
% response = response>0;
threshold = 0.5;
R = response;
[i,j]=find(R<threshold);
response(i,j) = 0; % an attempt at thresholding

% % response = response.^2;

resp = response;
resp_ = resp(Nlags:length(resp));



h=[];
J=[];
pfinal_all_jacks = [];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% break into training and test sets
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for jack = 1:njack; %loop over all njacks
    stim = stim_;
    resp = resp_;
Ntest = floor(length(stim)/njack);  % could be changed to do different # of jackknives
teststim = stim(1+(jack-1)*Ntest:jack*Ntest,:);
testresp = resp(1+(jack-1)*Ntest:jack*Ntest);
stim(1+(jack-1)*Ntest:jack*Ntest,:) = [];
resp(1+(jack-1)*Ntest:jack*Ntest) = [];


[Nsamples,Ndim] = size(stim);
psp = mean(mean(resp));   % spike probability
avg = (stim'*resp)/Nsamples;  % Ndim x Nrep
avg = mean(avg,2);  % Ndim x 1
avgs = [psp;avg];
if order>1
    avgsqrd = stim'*(repmat(resp,[1,Ndim]).*stim)/Nsamples;  % Ndim x Ndim
    avgsqrd = reshape(avgsqrd,[Ndim^2,1]);
    avgs = [avgs;avgsqrd];
end
% Initialize parameters
pstart = log(1/avgs(1)-1);
pstart(2:Ndim+1) = .001*(2*rand([1,Ndim])-1);
if order>1
    temp = .0005*(2*rand([Ndim,Ndim])-1);
    pstart(Ndim+2:length(pstart)+Ndim^2) = reshape((temp+temp'),[1,Ndim^2]);

    clear temp;
end


disp('Starting optimization');
tic
celltype = '';  % ignore this
%MNEfit(stim, resp, teststim, testresp, celltype, cellnum, jack, order, Nd, fittype);

%Instead of calling the MNEfit function, I just copy the script here below.



% Run conjugate gradient algorithm
pfinal = frprmn(pstart, @logloss, @dlogloss, stim, resp, teststim, testresp, order, avgs, Nd, fittype);

% 
% % % Save results
% if strcmp(celltype, 'RGC')==1 || strcmp(celltype, 'LGN')==1
%     if fittype==0
%         save([celltype '_' num2str(cellnum) '_MNE' num2str(order) '_jack_' num2str(jack) '.mat'],'pfinal');
%     else
%         save([celltype '_' num2str(cellnum) '_MNE' num2str(order) '_jack_' num2str(jack) '_random.mat'],'pfinal');
%     end
% else
%     if fittype==0
%         
% 
%           save(['st_' num2str(cellnum) '_m_211_s16_s17_s25_m_462_s2prime_m_462_s2_4_9_24kHz' '_Nlags' num2str(Nlags) '_nfft128_Nf16_bin21p6ms' '_jack_' num2str(jack) '_of_' num2str(njack) '.mat'],'pfinal');
%           
%           
%     else
%         %save(['ModelCell_' num2str(cellnum) '_MNE' num2str(order) '_jack_' num2str(jack) '_random.mat'],'pfinal');
%     end
% end

disp(['Optimization took ' num2str(toc/60) ' minutes']);

%And now plot the results and save the figures
%close all

h=[h;pfinal(2:Nlags*Nf+1)];
pfinal_all_jacks = [pfinal_all_jacks;pfinal];


 J=[J;pfinal(Nlags*Nf+2:end)];
[V,D] = eig(reshape(pfinal(Nlags*Nf+2:end),Nlags*Nf,Nlags*Nf));

fig=figure;
subplot(3,3,1:3)
eigenvalues = diag(D);
[eigenvalues_sorted,index] = sort(eigenvalues);
plot(eigenvalues_sorted,'o');

subplot(3,3,4)
eig_sorted_1 = V(:,index(1));
imagesc(reshape(eig_sorted_1,Nf,Nlags))
axis xy

subplot(3,3,5)
eig_sorted_2=V(:,index(2));
imagesc(reshape(eig_sorted_2,Nf,Nlags))
axis xy

subplot(3,3,6)
eig_sorted_3=V(:,index(3));
imagesc(reshape(eig_sorted_3,Nf,Nlags))
axis xy

subplot(3,3,7)
eig_sorted_end=V(:,index(end));
imagesc(reshape(eig_sorted_end,Nf,Nlags))
axis xy

subplot(3,3,8)
eig_sorted_end_1=V(:,index(end-1));
imagesc(reshape(eig_sorted_end_1,Nf,Nlags))
axis xy

subplot(3,3,9)
eig_sorted_end_2=V(:,index(end-2));
imagesc(reshape(eig_sorted_end_2,Nf,Nlags))

saveas(fig,[response_path '_thresh_' num2str(threshold) '_' num2str(section) 'samples' 'jack' num2str(jack) 'of' num2str(njack) 'jacks_' num2str(Nf) 'x' num2str(Nlags) '.png'])
close all;
% subplot(2,2,3)
% plot(diag(D),'o');


% figure_name = (['st_' num2str(cellnum) '_MNE' num2str(order) '_jack_' num2str(jack) 'Nlags' num2str(Nlags) '_freq_band_step_' num2str(freq_band) '_hJ']); 
% hgsave(figure_name); 
% 
% axis xy


%figure_name = (['st_' num2str(cellnum) '_m_211_s16_s17_s25_m_462_s2prime_m_462_s2_4_9_24kHz' '_Nlags' num2str(Nlags) '_nfft128_Nf16' '_jack_' num2str(jack) '_of_' num2str(njack)]);


% hgsave(figure_name); 
h_mean = mean(h);


end
J_mean = mean(J);
[V,D] = eig(reshape(J_mean,Nlags*Nf,Nlags*Nf));


fig = figure;
subplot(3,3,1:3)
eigenvalues = diag(D);
[eigenvalues_sorted,index] = sort(eigenvalues);
plot(eigenvalues_sorted,'o');

subplot(3,3,4)
eig_sorted_1 = V(:,index(1));
imagesc(reshape(eig_sorted_1,Nf,Nlags))
axis xy

subplot(3,3,5)
eig_sorted_2=V(:,index(2));
imagesc(reshape(eig_sorted_2,Nf,Nlags))
axis xy

subplot(3,3,6)
eig_sorted_3=V(:,index(3));
imagesc(reshape(eig_sorted_3,Nf,Nlags))
axis xy

subplot(3,3,7)
eig_sorted_end=V(:,index(end));
imagesc(reshape(eig_sorted_end,Nf,Nlags))
axis xy

subplot(3,3,8)
eig_sorted_end_1=V(:,index(end-1));
imagesc(reshape(eig_sorted_end_1,Nf,Nlags))
axis xy

subplot(3,3,9)
eig_sorted_end_2=V(:,index(end-2));
imagesc(reshape(eig_sorted_end_2,Nf,Nlags))
saveas(fig,[response_path '_thresh_' num2str(threshold) '_' num2str(section) 'samples' num2str(njack) 'jacks_' num2str(Nf) 'x' num2str(Nlags) '.png'])

save([response_path '_thresh_' num2str(threshold) '_' num2str(section) 'samples' num2str(njack) 'jacks_' num2str(Nf) 'x' num2str(Nlags) '.mat'],'J','h','pfinal_all_jacks')

figure;
imagesc(stimulus');
title('Input')

figure;
plot(response_file.spike_times_in_seconds,'*')
figure;
plot(response)
% figure;
% plot(response_file.binned_spike_trace(1:section*32))



%eval('J_random')
% end

