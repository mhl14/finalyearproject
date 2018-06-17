clc; close all; 

%% Load and prepare spectrogram of stimulus
stim_file = load('C:\Users\Edward\Documents\Imperial\UROP\stim.mat'); %NFFT=128, overlap = NFFT/2 %nonequispaced fast Fourier transform
stimulus = stim_file.stim_vectors;
[num_dim, time_bins]=size(stimulus); %900,450
num_repeat = 10;
stimulus = stimulus(2:end,:); %Remoce DC component, remove first dimension
stimulus = stimulus(:,1:floor(time_bins/num_repeat)); %Stimulus has 10 repetitions, so make divisble by 10
%floor rounds the element of X to the nearest integers
[num_dim, time_bins]=size(stimulus); %899,45
freq_compress = 4;
time_compress = 512;
stimulus = imresize(stimulus, [(num_dim)/freq_compress time_bins/time_compress],'bilinear'); %compress spectrogram 
%imresize compress "stimulus" into "(num_dim)/freq_compress" rows and
%"time_bins/time_compress" columns

size(stimulus) %225x16
section = length(stimulus); %225
stimulus = stimulus'; %transpose matrix
stimulus = stimulus(:, 1:section);
stimulus = double(stimulus); 
stimulus = 20*log10(stimulus); %get voltage values
stimulus = real(stimulus); %real values

[Ndim, Nsamples]=size(stimulus); %16,225
stimulus=stimulus'; %transpose matrix
% The stimulus units is changed to be z-score
% (zero mean, normalized by standard deviation),
% but could be put into your own unit system (or 
% maybe your stimulus file is already saved
% in the proper units).
stimulus = stimulus-repmat(mean(stimulus),[Nsamples,1]); 
%repmat creats a large matrix B consisting of Nsamples-by-1 tiling of
%copies of mean(stimulus) 
stimulus = stimulus./repmat(std(stimulus),[Nsamples,1]); %normalize?

stim_train = stimulus(1:round(length(stimulus)*.75),:); %first result
stim_test = stimulus(round(length(stimulus)*.75)+1:end,:); %the rest
stimulus = stim_train;
[Nsamples,Ndim]=size(stimulus); %169,16

%% MNE setup
% cellnum = 750;
Nf = Ndim; %number of frequency bands in Spectual-Temporal Receptive Field 
Nlags = 16; %times the length of a bin gives the time length of STRF, ie. total time/ length of a bin
order   = 2;   % order of MNE model to fit, second-order
fittype = 0;   % 0 for regular fitting, 1 for random fitting
njack   = 10;   % no. of jackknives to run (also determines the size of each jackknives)
%Nlags   = 1;   % no. of time lags to use in spatiotemporal models (=1 is just a
%spatial model)

% redefine stimulus to include time lags
if Nlags>1
    Nsamples = Nsamples - (Nlags-1); %total length of stimulus minus 19 time bins
    Ndimtotal = Ndim*Nlags; %16x20
    stim_ = zeros(Nsamples,Ndimtotal); %zero matrix
    for i=1:Nlags
        stim_(:,Ndim*(i-1)+1:Ndim*i) = stimulus(i:Nsamples+i-1,:);
    end
else
    stim_ = stimulus;
end

size(stim_); %154x16
%clear stimulus;
Nd = sqrt(Ndim); % no. of pixels per side of image, i.e 1 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Open responses and cut according to no. of time lags %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Point to spike trace
% Format is .mat file with binary response vector
response_path = 'C:\Users\Edward\Documents\Imperial\UROP\_binary_spike_trace_trode_1_cluster_no_15.mat';
response_file = load(response_path);
response = double(response_file.binned_spike_trace);
size(response); %1x5035386

%Average response over 10 repititions
response = response(1:end-mod(length(response),num_repeat));
%return 0 if both input are the same, otherwise return length(response)
response = reshape(response,[num_repeat,floor(length(response)/num_repeat)]); 
% calculates the length of the dimension represented by [ ], 
% such that the product of the dimensions equals NUMEL(response).
response = mean(response,1); %get the mean of response along dimension 1
size(response); %1x5035386

response = imresize(response, [1 time_compress],'bilinear'); %compress spectrogram 
%imresize compress "stimulus" into 1 rows and
%"section" columns
size(response) %1x512
response = response(:, 1:section);
response = response'; %transpose matrix 
size(response); %225x1

response = response./max(response);
% response = response>0;
% threshold = 0.3;
% R = response;
% [i,j]=find(R<threshold);
% response(i,j) = 0; % an attempt at thresholding

resp_train = response(1:round(length(response)*.75),:); %first result
resp_test = response(round(length(response)*.75)+1:end,:); %the rest
response = resp_train;
size(response); %169x1

resp = response;
resp_ = resp(Nlags:length(resp));
size(resp); %169x1

h=[];
J=[];
pfinal_all_jacks = [];
%% MNE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% break into training and test sets
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for jack = 1:njack; %loop over all njacks
    stim = stim_;
    size(stim); %154x16
    resp = resp_;
    size(resp); %154x1
    
    Ntest = floor(length(stim)/njack);  % could be changed to do different no. of jackknives, Ntest= 25
    teststim = stim(1+(jack-1)*Ntest:jack*Ntest,:);
    testresp = resp(1+(jack-1)*Ntest:jack*Ntest);
    stim(1+(jack-1)*Ntest:jack*Ntest,:) = [];
    resp(1+(jack-1)*Ntest:jack*Ntest) = [];
    
    [Nsamples,Ndim] = size(stim); %139,16
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
    tic %tic, by itself, saves the current time that TOC uses later to measure the time elapsed between the two.

    % Run conjugate gradient algorithm
    pfinal = frprmn(pstart, @logloss, @dlogloss, stim, resp, teststim, testresp, order, avgs, Nd, fittype);

    disp(['Optimization took ' num2str(toc/60) ' minutes']);

    %And now plot the results and save the figures

    h=[h;pfinal(2:Nlags*Nf+1)];
    size(h); %1x16
    pfinal_all_jacks = [pfinal_all_jacks;pfinal];
    size(pfinal_all_jacks); %1x273

    J=[J;pfinal(Nlags*Nf+2:end)];
    size(J); %1x256
    [V,D] = eig(reshape(pfinal(Nlags*Nf+2:end),Nlags*Nf,Nlags*Nf)); %get eigenvalues and eigenvectors
    size(V); %16x16
    size(D); %16x16
    close all
    
    fig=figure;
    subplot(3,3,1:3) %first 3 out of the 9
    %display eigenvalues
    eigenvalues = diag(D);
    [eigenvalues_sorted,index] = sort(eigenvalues); %sort in ascending or descending order
    plot(eigenvalues_sorted,'o');

    subplot(3,3,4) %4th out of the 9
    %eigenvector of first eigenvalue
    eig_sorted_1 = V(:,index(1));
    imagesc(reshape(eig_sorted_1,Nf,Nlags))
    axis xy

    subplot(3,3,5) %5th out of the 9
    %eigenvector of second eigenvalue
    eig_sorted_2=V(:,index(2));
    imagesc(reshape(eig_sorted_2,Nf,Nlags))
    axis xy

    subplot(3,3,6) %5th out of the 9
    %eigenvector of third eigenvalue
    eig_sorted_3=V(:,index(3));
    imagesc(reshape(eig_sorted_3,Nf,Nlags))
    axis xy

    subplot(3,3,7)  %7th out of the 9
    %eigenvector of last eigenvalue
    eig_sorted_end=V(:,index(end));
    imagesc(reshape(eig_sorted_end,Nf,Nlags))
    axis xy

    subplot(3,3,8) %8th out of the 9
    %eigenvector of second last eigenvalue
    eig_sorted_end_1=V(:,index(end-1));
    imagesc(reshape(eig_sorted_end_1,Nf,Nlags))
    axis xy

    subplot(3,3,9) %9th out of the 9
    %eigenvector of third last eigenvalue
    eig_sorted_end_2=V(:,index(end-2));
    imagesc(reshape(eig_sorted_end_2,Nf,Nlags))

    %saveas(fig,[response_path '_thresh_' num2str(threshold) '_' num2str(section) 'samples' 'jack' num2str(jack) 'of' num2str(njack) 'jacks_' num2str(Nf) 'x' num2str(Nlags) '.png'])
    %close all;
end

J_mean = mean(J); %get mean 
size(J_mean); %1x256
h_mean = mean(h); %get mean
size(h_mean); %1x16
[V,D] = eig(reshape(J_mean,Nlags*Nf,Nlags*Nf)); %get eigenvalues of means
size(V); %16x16
size(D); %16x16

close all
%% Plot results and misc
fig = figure;
subplot(3,3,1:3)
eigenvalues = diag(D);
[eigenvalues_sorted,index] = sort(eigenvalues);
plot(eigenvalues_sorted,'o');
title('eigenvalues')

subplot(3,3,4)
eig_sorted_1 = V(:,index(1));
size(eig_sorted_1); %16x1
size(Nf); %1x1
size(Nlags); %1x1
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

% saveas(fig,[response_path '_thresh_' num2str(threshold) '_' num2str(section) 'samples' num2str(njack) 'jacks_' num2str(Nf) 'x' num2str(Nlags) '.png'])

% save([response_path 'MNEeval_thresh_' num2str(threshold) '_' num2str(section) 'samples' num2str(njack) 'jacks_' num2str(Nf) 'x' num2str(Nlags) '.mat'],'J','h','pfinal_all_jacks','stim_test','resp_test','pfinal_all_jacks')

figure;
imagesc(stimulus'); %display image with scaled colors
title('Input stimulus')

figure;
plot(response_file.spike_times_in_seconds,'*')
title('spike time in seconds')

figure;
plot(response)
title('response')

% figure;
% plot(response_file.binned_spike_trace(1:section*32))


%eval('J_random')
% end


