%% Load and prepare spectrogram of stimulus
stim_file = load('C:\Users\Edward\Documents\Imperial\UROP\stim.mat'); %NFFT=128, overlap = NFFT/2
stimulus = stim_file.stimulus;
[num_dim, time_bins]=size(stimulus);
num_rep = 10;
stimulus = stimulus(2:end,:); %Remoce DC component
stimulus = stimulus(:,1:floor(time_bins/num_rep)); %Stimulus has 10 repetitions, so make divisble by 10
[num_dim, time_bins]=size(stimulus);
freq_compress = 4;
time_compress = 512;
stimulus = imresize(stimulus, [(num_dim)/freq_compress time_bins/time_compress],'bilinear'); %compress spectrogram 
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

stim_train = stimulus(1:round(length(stimulus)*.75),:);
stim_test = stimulus(round(length(stimulus)*.75)+1:end,:);
stimulus = stim_train;
[Nsamples,Ndim]=size(stimulus);

%% MNE setup
% cellnum = 750;
Nf = Ndim; %number of frequency bands in STRF 
Nlags = 16; %times the length of a bin gives the time length of STRF
order   = 2;   % order of MNE model to fit
fittype = 0;   % 0 for regular fitting, 1 for random fitting
njack   = 10;   % # jackknives to run (also determines the size of each jackknives)
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

% Point to spike trace
% Format is .mat file with binary response vector
response_path = 'C:\Users\Edward\Documents\Imperial\UROP\_binary_spike_trace_trode_1_cluster_no_15.mat';
response_file = load(response_path);
response = double(response_file.binned_spike_trace);
%Average response over 10 repititions
response = response(1:end-mod(length(response),num_rep));

response = reshape(response,[num_rep,floor(length(response)/num_rep)]); 
response = mean(response,1);

response = imresize(response, [1 time_bins/time_compress],'bilinear');
size(response)
response = response(:, 1:section);
response = response';
 



response = response./max(response);
% response = response>0;
%threshold = 0.3;
%R = response;
%[i,j]=find(R<threshold);
%response(i,j) = 0; % an attempt at thresholding

resp_train = response(1:round(length(response)*.75),:);
resp_test = response(round(length(response)*.75)+1:end,:);
response = resp_train;


resp = response;
resp_ = resp(Nlags:length(resp));



h=[];
J=[];
pfinal_all_jacks = [];
%% MNE
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


% Run conjugate gradient algorithm
pfinal = frprmn(pstart, @logloss, @dlogloss, stim, resp, teststim, testresp, order, avgs, Nd, fittype);

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



end
J_mean = mean(J);
h_mean = mean(h);
[V,D] = eig(reshape(J_mean,Nlags*Nf,Nlags*Nf));

%% Plot results and misc
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

save([response_path 'MNEeval_thresh_' num2str(threshold) '_' num2str(section) 'samples' num2str(njack) 'jacks_' num2str(Nf) 'x' num2str(Nlags) '.mat'],'J','h','pfinal_all_jacks','stim_test','resp_test','pfinal_all_jacks')

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


