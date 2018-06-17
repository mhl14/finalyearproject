%% clear window
clc; close all; 
addpath (genpath ('v2-1 data'));

%% New receptive field plot 
% MNE parameters setup
MNE_params.Ndim = 16; %number of frequency bands in STRF 
MNE_params.Nlags = 20; %times the length of a bin gives the time length of STRF
% # time lags to use in spatiotemporal models (=1 is just a spatial model)
MNE_params.order = 2;   % order of MNE model to fit (1 or 2)
MNE_params.fittype = 0;   % 0 for regular fitting, 1 for random fitting
MNE_params.Njack = 4; % Number of jackknives to use, 4 per songs 


% Fit MNE model
[V,D] = eig(reshape(J_mean,MNE_params.Nlags*MNE_params.Ndim,MNE_params.Nlags*MNE_params.Ndim));
% [V,D] = eig(J_mean);

% Get confidence intervals for eigenvalues

Nf = MNE_params.Ndim; %number of frequency bands in STRF
Nlags = MNE_params.Nlags ; %times the length of a bin gives the time length of STRF
J_reshaped = reshape(J_mean,Nlags*Nf,Nlags*Nf);

%now we reshuffle J and save it in J_r
J_r = zeros(size(J_reshaped)); % initialize square matrix of size J (currently 320x320)
J_r_alltogether = zeros(size(J_reshaped));

ind_r = randperm(length(J_mean)); % generate 320 random integers
J_r_alltogether = J_mean(ind_r); % permute J elements using those indices
J_r_alltogether = reshape(J_r_alltogether,Nlags*Nf,Nlags*Nf); %make it square
 
%here we symmetrize this matrix J_r_alltogether; here we scramble diagonal and non-diagonal elements together
J_r_alltogether = diag(diag(J_r_alltogether)) + triu(J_r_alltogether,1) + triu(J_r_alltogether,1)'; 

% now we permute diagonal and non-diagonal elements separately, for comparizon
ind_diag = randperm(size(J_reshaped,1)); %randomly permute 320 integers from 1 to 320
J_diagonal = diag(J_reshaped); %get the diagonal
J_r_diagonal = J_diagonal(ind_diag); %permute the diagonal of J_reshaped using those indices

% put zeros on the diagonal
J_reshaped_zerodiag = tril(J_reshaped,-1)+triu(J_reshaped,1) + diag(zeros(length(J_reshaped),1)); 
J_reshaped_zerodiag_1d = reshape(J_reshaped_zerodiag,[],1); % reshape as 320x1 vector
ind = randperm(size(J_reshaped_zerodiag_1d,1)); %get random integers for permutation indices
J_reshaped_zerodiag_1d_shuffled = J_reshaped_zerodiag_1d(ind); %permute the vector
%make it square
J_reshaped_zerodiag_1d_shuffled_square = reshape(J_reshaped_zerodiag_1d_shuffled,Nlags*Nf,Nlags*Nf); 
[row,col,v] = find(J_reshaped_zerodiag_1d_shuffled_square); %find nonzero elements, i.e. all nondiagonal elements
zeros_lenJ = zeros(length(J_reshaped),1); % we need zeros to fill in the matrix to keep its size 320x320
matrix_temp = vertcat(v,zeros_lenJ); % reshuffled offdiagonal elements and 320 zeroes, to keep the corect size
%zeros are in the last column, they will get away since we don't use the upper tiangle
matrix_temp = reshape(matrix_temp,length(J_reshaped),length(J_reshaped)); 
% reconstruct the final reshuffled matrix by putting the reshuffled lower 
% triangles (so it's square) and the reshuffled diagonal
J_r = tril(matrix_temp,-1)+tril(matrix_temp,-1)'+diag(J_r_diagonal); 

%take its eigenvalues of J_r_alltogether, and sort them
[V_r_alltogether,D_r_alltogether] = eig(J_r_alltogether);
eigenvalues_r_alltogether = diag(D_r_alltogether);
[eigenvalues_sorted_r_alltogether,indx_alltogether] = sort(eigenvalues_r_alltogether);
% figure; hold on
% plot(eigenvalues_sorted_r_alltogether,'b');

%take its eigenvalues of this J_r, and sort them
[V_r,D_r] = eig(J_r);
eigenvalues_r = diag(D_r);
[eigenvalues_sorted_r,indx] = sort(eigenvalues_r);
% plot(eigenvalues_sorted_r,'r');

% now we use a second method to control for statistical significance: we
% create 500 random gaussian matrices (symmetrical) with size and moments
% equal to those of J, and calculate a distribution of their eigenvalues.
% We use this distribution to detect eigenvalues of J that are significant
% i.e. those that lie outside.

%initialize to empty matrices
eig_max_all = [];
eig_min_all = [];
eig_all = [];

for i = 1:500 %we will use 500 random matrices   
    %here we generate a random matrix of size, mean and variance equal to those of J 
    J_random = mean(J_mean)+std(J_mean).*randn(size(J_mean)); 

     J_random = reshape(J_random,Nlags*Nf,Nlags*Nf); %here we make it square

     J_rand = diag(diag(J_random)) + triu(J_random,1) + triu(J_random,1)'; %here we symmetrize this matrix

     % get its eigenvalues
    [V_rand,D_rand] = eig(J_rand);
    eigenvalues_rand = diag(D_rand);
    %[eigenvalues_sorted_random,index] = sort(eigenvalues_rand);

    eig_max = max(eigenvalues_rand);
    eig_min = min(eigenvalues_rand);

    eig_max_all = [eig_max_all eig_max];
    eig_min_all = [eig_min_all eig_min];

    eig_all = [eig_all eigenvalues_rand];
end

eig_all = reshape(eig_all,[],1);
[eig_all_sorted,ix] = sort(eig_all);

%plot the eigenvalues of the reshuffled matrix J (J_r), and the reference
%lines for the min and max of the distribution of random generated eigenvalues
% plot(eigenvalues_sorted_r,'r');
% refline(0,min(eig_all))
% refline(0,max(eig_all))

max_eig_all = max(eig_all);
min_eig_all = min(eig_all);
%}

% Subplot results 
clc; close all; 
ydim = MNE_params.Ndim; 
xdim = MNE_params.Nlags;
nfft = 128;
step = 0.0194363;

figure;
% subplot(3,3,1)
% % imagesc(h_mean);
% newT = -20*step:step:0;
% freqs=0:22050/(nfft/2):22050;
% % imagesc(newT,freqs,reshape(h_mean,ydim,xdim));axis xy;
% imagesc(reshape(h_mean,ydim,xdim)); axis xy; 
% colormap(colormap(jet(256)));
% xlabel('Time(sec)')
% ylabel('Freq (kHz)')
% title('Linear Filter')

subplot(3,4,1:4)
% hold on
eigenvalues = diag(D);
[eigenvalues_sorted,index] = sort(eigenvalues);
plot(eigenvalues_sorted,'o');
xlabel('Feature', 'FontSize', 14)
ylabel('Magnitude', 'FontSize', 14)
title('Eigenspectrum of the Quadratic Filter J', 'FontSize', 20)
refline(0,min_eig_all)
refline(0,max_eig_all)
% hold off 

subplot(3,4,5)
eig_sorted_1 = V(:,index(1));
newT = -20*step:step:0;
freqs=0:(nfft/2000):22;
imagesc(newT,freqs,reshape(eig_sorted_1,ydim,xdim));axis xy;
colorbar('ticks',[-0.1,0,0.1]);
% imagesc(reshape(eig_sorted_1,ydim,xdim))
% axis xy
% xlabel('Time(sec)')
ylabel('Freq (kHz)', 'FontSize', 14)
xlabel('Time(sec)', 'FontSize', 14)
colormap(colormap(jet(256)));

subplot(3,4,6)
eig_sorted_2=V(:,index(2));
newT = -20*step:step:0;
freqs=0:(nfft/2000):22;
imagesc(newT,freqs,reshape(eig_sorted_2,ydim,xdim));axis xy;
colorbar('ticks',[-0.1,0,0.1]);
% imagesc(reshape(eig_sorted_2,ydim,xdim))
% axis xy
colormap(colormap(jet(256)));
ylabel('Freq (kHz)', 'FontSize', 14)
xlabel('Time(sec)', 'FontSize', 14)
% title('second smallest eigenvalue')
title('Excitatory Receptive Field Features', 'FontSize', 20)

subplot(3,4,7)
eig_sorted_3=V(:,index(3));
newT = -20*step:step:0;
freqs=0:(nfft/2000):22;
imagesc(newT,freqs,reshape(eig_sorted_3,ydim,xdim));axis xy;
colorbar('ticks',[-0.1,0,0.1]);
% imagesc(reshape(eig_sorted_3,ydim,xdim))
% axis xy
ylabel('Freq (kHz)', 'FontSize', 14)
xlabel('Time(sec)', 'FontSize', 14)
colormap(colormap(jet(256)));
% title('third smallest eigenvalue')


subplot(3,4,8)
eig_sorted_4=V(:,index(4));
newT = -20*step:step:0;
freqs=0:(nfft/2000):22;
imagesc(newT,freqs,reshape(eig_sorted_4,ydim,xdim));axis xy;
colorbar('ticks',[-0.1,0,0.1]);
% imagesc(reshape(eig_sorted_3,ydim,xdim))
% axis xy
ylabel('Freq (kHz)', 'FontSize', 14)
xlabel('Time(sec)', 'FontSize', 14)
colormap(colormap(jet(256)));
% title('third smallest eigenvalue')

subplot(3,4,9)
eig_sorted_end=V(:,index(end));
newT = -20*step:step:0;
freqs=0:(nfft/2000):22;
imagesc(newT,freqs,reshape(eig_sorted_end,ydim,xdim));axis xy;
colorbar('ticks',[-0.1,0,0.1]);
% imagesc(reshape(eig_sorted_end,ydim,xdim))
% axis xy
xlabel('Time(sec)', 'FontSize', 14)
ylabel('Freq (kHz)', 'FontSize', 14)
colormap(colormap(jet(256)));

subplot(3,4,10)
eig_sorted_end_1=V(:,index(end-1));
newT = -20*step:step:0;
freqs=0:(nfft/2000):22;
imagesc(newT,freqs,reshape(eig_sorted_end_1,ydim,xdim));axis xy;
colorbar('ticks',[-0.1,0,0.1]);
% imagesc(reshape(eig_sorted_end_1,ydim,xdim))
% axis xy
xlabel('Time(sec)', 'FontSize', 14)
ylabel('Freq (kHz)', 'FontSize', 14)
colormap(colormap(jet(256)));
% title('second largest eigenvalue')
title('Inhibitory Receptive Field Features', 'FontSize', 20)

subplot(3,4,11)
eig_sorted_end_2=V(:,index(end-2));
newT = -20*step:step:0;
freqs=0:(nfft/2000):22;
imagesc(newT,freqs,reshape(eig_sorted_end_2,ydim,xdim));axis xy;
colorbar('ticks',[-0.1,0,0.1]);
% imagesc(reshape(eig_sorted_end_2,ydim,xdim))
% axis xy
xlabel('Time(sec)', 'FontSize', 14)
ylabel('Freq (kHz)', 'FontSize', 14)
colormap(colormap(jet(256)));
% title('third largest eigenvalue')


subplot(3,4,12)
eig_sorted_end_3=V(:,index(end-3));
newT = -20*step:step:0;
freqs=0:(nfft/2000):22;
imagesc(newT,freqs,reshape(eig_sorted_end_3,ydim,xdim));axis xy;
colorbar('ticks',[-0.1,0,0.1]);
% imagesc(reshape(eig_sorted_end_2,ydim,xdim))
% axis xy
xlabel('Time(sec)', 'FontSize', 14)
ylabel('Freq (kHz)', 'FontSize', 14)
colormap(colormap(jet(256)));
% title('third largest eigenvalue')

%% New visual receptive field plot 
% MNE parameters setup
MNE_params.Ndim = 16; %number of frequency bands in STRF 
MNE_params.Nlags = 20; %times the length of a bin gives the time length of STRF
% # time lags to use in spatiotemporal models (=1 is just a spatial model)
MNE_params.order = 2;   % order of MNE model to fit (1 or 2)
MNE_params.fittype = 0;   % 0 for regular fitting, 1 for random fitting
MNE_params.Njack = 4; % Number of jackknives to use, 4 per songs 

% Fit MNE model
[V,D] = eig(reshape(J_mean,MNE_params.Nlags*MNE_params.Ndim,MNE_params.Nlags*MNE_params.Ndim));
% [V,D] = eig(J_mean);

% Get confidence intervals for eigenvalues

Nf = MNE_params.Ndim; %number of frequency bands in STRF
Nlags = MNE_params.Nlags ; %times the length of a bin gives the time length of STRF
J_reshaped = reshape(J_mean,Nlags*Nf,Nlags*Nf);

%now we reshuffle J and save it in J_r
J_r = zeros(size(J_reshaped)); % initialize square matrix of size J (currently 320x320)
J_r_alltogether = zeros(size(J_reshaped));

ind_r = randperm(length(J_mean)); % generate 320 random integers
J_r_alltogether = J_mean(ind_r); % permute J elements using those indices
J_r_alltogether = reshape(J_r_alltogether,Nlags*Nf,Nlags*Nf); %make it square
 
%here we symmetrize this matrix J_r_alltogether; here we scramble diagonal and non-diagonal elements together
J_r_alltogether = diag(diag(J_r_alltogether)) + triu(J_r_alltogether,1) + triu(J_r_alltogether,1)'; 

% now we permute diagonal and non-diagonal elements separately, for comparizon
ind_diag = randperm(size(J_reshaped,1)); %randomly permute 320 integers from 1 to 320
J_diagonal = diag(J_reshaped); %get the diagonal
J_r_diagonal = J_diagonal(ind_diag); %permute the diagonal of J_reshaped using those indices

% put zeros on the diagonal
J_reshaped_zerodiag = tril(J_reshaped,-1)+triu(J_reshaped,1) + diag(zeros(length(J_reshaped),1)); 
J_reshaped_zerodiag_1d = reshape(J_reshaped_zerodiag,[],1); % reshape as 320x1 vector
ind = randperm(size(J_reshaped_zerodiag_1d,1)); %get random integers for permutation indices
J_reshaped_zerodiag_1d_shuffled = J_reshaped_zerodiag_1d(ind); %permute the vector
%make it square
J_reshaped_zerodiag_1d_shuffled_square = reshape(J_reshaped_zerodiag_1d_shuffled,Nlags*Nf,Nlags*Nf); 
[row,col,v] = find(J_reshaped_zerodiag_1d_shuffled_square); %find nonzero elements, i.e. all nondiagonal elements
zeros_lenJ = zeros(length(J_reshaped),1); % we need zeros to fill in the matrix to keep its size 320x320
matrix_temp = vertcat(v,zeros_lenJ); % reshuffled offdiagonal elements and 320 zeroes, to keep the corect size
%zeros are in the last column, they will get away since we don't use the upper tiangle
matrix_temp = reshape(matrix_temp,length(J_reshaped),length(J_reshaped)); 
% reconstruct the final reshuffled matrix by putting the reshuffled lower 
% triangles (so it's square) and the reshuffled diagonal
J_r = tril(matrix_temp,-1)+tril(matrix_temp,-1)'+diag(J_r_diagonal); 

%take its eigenvalues of J_r_alltogether, and sort them
[V_r_alltogether,D_r_alltogether] = eig(J_r_alltogether);
eigenvalues_r_alltogether = diag(D_r_alltogether);
[eigenvalues_sorted_r_alltogether,indx_alltogether] = sort(eigenvalues_r_alltogether);
% figure; hold on
% plot(eigenvalues_sorted_r_alltogether,'b');

%take its eigenvalues of this J_r, and sort them
[V_r,D_r] = eig(J_r);
eigenvalues_r = diag(D_r);
[eigenvalues_sorted_r,indx] = sort(eigenvalues_r);
% plot(eigenvalues_sorted_r,'r');

% now we use a second method to control for statistical significance: we
% create 500 random gaussian matrices (symmetrical) with size and moments
% equal to those of J, and calculate a distribution of their eigenvalues.
% We use this distribution to detect eigenvalues of J that are significant
% i.e. those that lie outside.

%initialize to empty matrices
eig_max_all = [];
eig_min_all = [];
eig_all = [];

for i = 1:500 %we will use 500 random matrices   
    %here we generate a random matrix of size, mean and variance equal to those of J 
    J_random = mean(J_mean)+std(J_mean).*randn(size(J_mean)); 

     J_random = reshape(J_random,Nlags*Nf,Nlags*Nf); %here we make it square

     J_rand = diag(diag(J_random)) + triu(J_random,1) + triu(J_random,1)'; %here we symmetrize this matrix

     % get its eigenvalues
    [V_rand,D_rand] = eig(J_rand);
    eigenvalues_rand = diag(D_rand);
    %[eigenvalues_sorted_random,index] = sort(eigenvalues_rand);

    eig_max = max(eigenvalues_rand);
    eig_min = min(eigenvalues_rand);

    eig_max_all = [eig_max_all eig_max];
    eig_min_all = [eig_min_all eig_min];

    eig_all = [eig_all eigenvalues_rand];
end

eig_all = reshape(eig_all,[],1);
[eig_all_sorted,ix] = sort(eig_all);

%plot the eigenvalues of the reshuffled matrix J (J_r), and the reference
%lines for the min and max of the distribution of random generated eigenvalues
% plot(eigenvalues_sorted_r,'r');
% refline(0,min(eig_all))
% refline(0,max(eig_all))

max_eig_all = max(eig_all);
min_eig_all = min(eig_all);
%}

% Subplot results 
clc; close all; 
ydim = MNE_params.Ndim; 
xdim = MNE_params.Nlags;
nfft = 128;
step = 0.0194363;

figure;
% subplot(3,3,1)
% % imagesc(h_mean);
% newT = -20*step:step:0;
% freqs=0:22050/(nfft/2):22050;
% % imagesc(newT,freqs,reshape(h_mean,ydim,xdim));axis xy;
% imagesc(reshape(h_mean,ydim,xdim)); axis xy; 
% colormap(colormap(jet(256)));
% xlabel('Time(sec)')
% ylabel('Freq (kHz)')
% title('Linear Filter')

subplot(3,4,1:4)
% hold on
eigenvalues = diag(D);
[eigenvalues_sorted,index] = sort(eigenvalues);
plot(eigenvalues_sorted,'o');
xlabel('Feature', 'FontSize', 14)
ylabel('Magnitude', 'FontSize', 14)
title('Eigenspectrum of the Quadratic Filter J', 'FontSize', 20)
refline(0,min_eig_all)
refline(0,max_eig_all)
% hold off 

subplot(3,4,5)
eig_sorted_1 = V(:,index(1));
imagesc(reshape(eig_sorted_1,ydim,xdim));axis xy;
colorbar('ticks',[-0.1,0,0.1]);
% imagesc(reshape(eig_sorted_1,ydim,xdim))
% axis xy
% xlabel('Time(sec)')
% ylabel('Freq (kHz)', 'FontSize', 14)
% xlabel('Time(sec)', 'FontSize', 14)
colormap(colormap(jet(256)));

subplot(3,4,6)
eig_sorted_2=V(:,index(2));
imagesc(reshape(eig_sorted_2,ydim,xdim));axis xy;
colorbar('ticks',[-0.1,0,0.1]);
% imagesc(reshape(eig_sorted_2,ydim,xdim))
% axis xy
colormap(colormap(jet(256)));
% xlabel('Time(sec)', 'FontSize', 14)
% title('second smallest eigenvalue')
title('Excitatory Receptive Field Features', 'FontSize', 20)

subplot(3,4,7)
eig_sorted_3=V(:,index(3));
imagesc(reshape(eig_sorted_3,ydim,xdim));axis xy;
colorbar('ticks',[-0.1,0,0.1]);
% imagesc(reshape(eig_sorted_3,ydim,xdim))
% axis xy
% xlabel('Time(sec)', 'FontSize', 14)
colormap(colormap(jet(256)));
% title('third smallest eigenvalue')


subplot(3,4,8)
eig_sorted_4=V(:,index(4));
imagesc(reshape(eig_sorted_4,ydim,xdim));axis xy;
colorbar('ticks',[-0.1,0,0.1]);
% imagesc(reshape(eig_sorted_3,ydim,xdim))
% axis xy
% xlabel('Time(sec)', 'FontSize', 14)
colormap(colormap(jet(256)));
% title('third smallest eigenvalue')

subplot(3,4,9)
eig_sorted_end=V(:,index(end));
imagesc(reshape(eig_sorted_end,ydim,xdim));axis xy;
colorbar('ticks',[-0.1,0,0.1]);
% imagesc(reshape(eig_sorted_end,ydim,xdim))
% axis xy
% xlabel('Time(sec)', 'FontSize', 14)
% ylabel('Freq (kHz)', 'FontSize', 14)
colormap(colormap(jet(256)));

subplot(3,4,10)
eig_sorted_end_1=V(:,index(end-1));
imagesc(reshape(eig_sorted_end_1,ydim,xdim));axis xy;
colorbar('ticks',[-0.1,0,0.1]);
% imagesc(reshape(eig_sorted_end_1,ydim,xdim))
% axis xy
% xlabel('Time(sec)', 'FontSize', 14)
% ylabel('Freq (kHz)')
colormap(colormap(jet(256)));
% title('second largest eigenvalue')
title('Inhibitory Receptive Field Features', 'FontSize', 20)

subplot(3,4,11)
eig_sorted_end_2=V(:,index(end-2));
imagesc(reshape(eig_sorted_end_2,ydim,xdim));axis xy;
colorbar('ticks',[-0.1,0,0.1]);
% imagesc(reshape(eig_sorted_end_2,ydim,xdim))
% axis xy
% xlabel('Time(sec)', 'FontSize', 14)
% ylabel('Freq (kHz)')
colormap(colormap(jet(256)));
% title('third largest eigenvalue')


subplot(3,4,12)
eig_sorted_end_3=V(:,index(end-3));
% newT = -20*step:step:0;
% freqs=0:22050/(nfft/2):22050;
imagesc(reshape(eig_sorted_end_3,ydim,xdim));axis xy;
colorbar('ticks',[-0.1,0,0.1]);
% imagesc(reshape(eig_sorted_end_2,ydim,xdim))
% axis xy
% xlabel('Time(sec)', 'FontSize', 14)
% ylabel('Freq (kHz)')
colormap(colormap(jet(256)));
% title('third largest eigenvalue')

%% Invariance comparison 4 plot 
figure
subplot(2,2,1)
corre_coef = [corre_coef18_20; corre_coef5_12; corre_coef7_4; corre_coef8_3];
corre_coef= sort(corre_coef);
hold on
for i=1:4
    scatter(corre_coef(i,2),corre_coef(i,1));
end
xlabel('fast-warped','FontSize',20);
ylabel('normal speed','FontSize',20);
axis([-0.2 0.2 -0.2 0.2]);

subplot(2,2,2)
hold on
corre_coef = [corre_coef18_20; corre_coef5_12; corre_coef7_4; corre_coef8_3];
for i=1:4
    scatter(corre_coef(i,3),corre_coef(i,1));
end
xlabel('20% shifted up in frequency','FontSize',20);
ylabel('normal speed','FontSize',20);
axis([-0.2 0.2 -0.2 0.2]);

subplot(2,2,3)
hold on
corre_coef = [corre_coef18_20; corre_coef5_12; corre_coef7_4; corre_coef8_3];
for i=1:4
    scatter(corre_coef(i,4),corre_coef(i,1));
end
xlabel('20% shifted down in frequency','FontSize',20);
ylabel('normal speed','FontSize',20)
axis([-0.2 0.2 -0.2 0.2]);

subplot(2,2,4)
hold on
corre_coef = [corre_coef18_20; corre_coef5_12; corre_coef7_4; corre_coef8_3];
for i=1:4
    scatter(corre_coef(i,5),corre_coef(i,1));
end
xlabel('slow-warped','FontSize',20);
ylabel('normal speed','FontSize',20)
axis([-0.2 0.2 -0.2 0.2]);

%% plotting window size analysis
figure; 
percdata = zeros(11,5);
for i=1:5
    for j= 1:11
        percdata(j,i) = corre_coef(j,i).*100./self_coef(j,i);
    end
end 
hold on; 
plot(step, percdata(1:11,1), 'lineWidth',3 );
plot(step, percdata(1:11,2), 'lineWidth',3 );
% plot(step, percdata(1:11,3));
% plot(step, percdata(1:11,4));
plot(step, percdata(1:11,5), 'lineWidth',3 );
% y=0; plot(step, y, 'linewidth', 10);
grid on
xlabel('window size(sec)', 'FontSize', 28);
title('Window Size Analysis', 'FontSize', 28);
ylabel('% correlation coefficient', 'FontSize', 28);
legend({'normal speed','fast-warped','slow-warped'}, 'FontSize', 12);


%% testing prediction
% Reshape parameters
a=A_mean;
ydim = 16;
xdim = 20;
Nlags=xdim; Ndim=ydim*xdim;
h=h_mean; 
J= reshape(J_mean,Ndim,Ndim); 
start = 0; 
nfft = 128;
window = 1224;
overlap = 0.3;
freq_compress = 4;
corre_coef = zeros(11,5);
self_coef = zeros(11,5);
% step = 0.0194363:(0.194363/2-0.0194363)/10:0.194363/2;

toelist = ...
    'response/concat_chan_17_18_electrode_5_12/30 reps/ss001m_497_17_s1_slow_toe.txt';
stop = 90; %30 for fast, 60 for normal, 90 for slow
[~, psth, spec]= plot_raster_SMI2(toelist, start, stop, step(1,t), nfft, window, overlap);
time_compress = 2; %5.74 for fast, 11.5 for normal, 17.2 for slow
[tstimulus, tresponse, tresponse_nor, ~, Nsample, tNdim]=compression(spec, psth', freq_compress, time_compress); 
close all 
% tresponse = tresponse'; 
[Tdim, Tsample]=size(tresponse);
% [tNdim, Nsample]=size(tstimulus);
tNsamples_compressed = Nsample - (Nlags-1); 
tNdimtotal = tNdim*Nlags;
stim_ = zeros(tNdimtotal, tNsamples_compressed);
for i=1:Nlags
    stim_(tNdim*(i-1)+1:tNdim*i,:) = ...
        tstimulus(:,i:tNsamples_compressed+i-1);
end
stim_=stim_';

for b = 1:1:Tsample
     pSpike5 = 1./(1+exp(a+stim_*h'+sum(stim_.*(stim_*J),2)));  % Nsamples x 1
end

pSpike_nor5 = pSpike5 ./max(pSpike5);
% tresponse_nor = tresponse./ max(tresponse); 
hold on; 
xbins = start:((stop-start)/length(tresponse)):stop;
plot(xbins(1,Nlags+1:end),tresponse_nor(:,Nlags:end),'r', 'lineWidth',3 ); 
xbins2 = start:((stop-start)/length(pSpike_nor5)):stop;
plot(xbins2(1,2:end), pSpike_nor5,'b', 'lineWidth',3 ); 
axis([start stop 0 1]);
legend({'real spike', 'predicted spike'}, 'FontSize', 16);
xlabel('Time(sec)', 'FontSize', 28); ylabel('Spike probability', 'FontSize', 28);
title('Real Sampled and Predicted Spike Response Binned at 97.2ms', 'FontSize', 28);

%% discussion 3 spectrogram plot
% normal speed 
subplot(2,4,2:3)
nfft = 128;
window = 1224;
overlap = 0.3;
xmin=0;
xmax=60;
stimpath = 'freq_shifted_song/stimulus';
stimulus1 = 'response/concat_chan_17_18_electrode_5_12/10 reps/ss001m_497_17_s1_toe.txt';
nlap = round(window*overlap);
[stimfile, subjectID, ~, site, sort1, ~, nreps, ...
    nspikes, toes, alltoes,~] = readtoe_2(stimulus1);
stimfile = strrep(stimfile, '.wav', '');
figname = strrep(sprintf('subj:%s  site:%s  sort1:%s  stim:%s', ...
    subjectID,site,sort1,stimfile), '_', '\_');
fullstim = [stimpath '/' stimfile '.wav' ];
[Y,FS]=audioread(fullstim);
[~,~,T,P] = spectrogram(Y,window,nlap,nfft,FS, 'yaxis');

newT = xmin:((xmax-xmin)/length(T)):xmax;
freqs=0:22050/(nfft):22050;
if(length(newT)~=length(newT))
    newT=newT(1:length(T));
end
clim = [-200  -65];
spec = 20*log10(P);
spec = spec(2:end, 1:end);
imagesc(newT,freqs,spec,clim);  
axis xy;
colormap(colormap(jet(256)));

axis([xmin xmax 0 22050]); %11025
xlabel('Time(sec)', 'FontSize', 20);
ylabel('Freq (Hz)','FontSize', 20);
title('Normal Speed','FontSize', 24);

% fast-warped
subplot(2,4,5:6)
xmax=60;
stimpath = 'freq_shifted_song/stimulus';
% stimulus2 = 'response/concat_chan_17_18_electrode_5_12/30 reps/ss001m_497_17_s1_2fast_toe.txt';
stimulus2 = 'response/concat_chan_17_18_electrode_5_12/30 reps/ss001m_497_17_s1_shift_1_toe.txt';
nlap = round(window*overlap);
[stimfile, subjectID, ~, site, sort1, ~, nreps, ...
    nspikes, toes, alltoes,~] = readtoe_2(stimulus2);
stimfile = strrep(stimfile, '.wav', '');
figname = strrep(sprintf('subj:%s  site:%s  sort1:%s  stim:%s', ...
    subjectID,site,sort1,stimfile), '_', '\_');
fullstim = [stimpath '/' stimfile '.wav' ];
[Y,FS]=audioread(fullstim);
[~,~,T,P] = spectrogram(Y,window,nlap,nfft,FS, 'yaxis');

newT = xmin:((xmax-xmin)/length(T)):xmax;
freqs=0:22050/(nfft):22050;
if(length(newT)~=length(newT))
    newT=newT(1:length(T));
end
clim = [-200  -65];
spec = 20*log10(P);
spec = spec(2:end, 1:end);
imagesc(newT,freqs,spec,clim);  
axis xy;
colormap(colormap(jet(256)));

axis([xmin xmax 0 22050]); %11025
xlabel('Time(sec)','FontSize', 20);
ylabel('Freq (Hz)','FontSize', 20);
title('20% Shifted Up in Frequency','FontSize', 24);
% title('Fast-Warped','FontSize', 24);

% %% slow-warped
subplot(2,4,7:8)
xmax=60;
stimpath = 'freq_shifted_song/stimulus';
% stimulus3 = 'response/concat_chan_17_18_electrode_5_12/30 reps/ss001m_497_17_s1_slow_toe.txt';
stimulus3 = 'response/concat_chan_17_18_electrode_5_12/30 reps/ss001m_497_17_s1_shift_6_toe.txt';
nlap = round(window*overlap);
[stimfile, subjectID, ~, site, sort1, ~, nreps, ...
    nspikes, toes, alltoes,~] = readtoe_2(stimulus3);
stimfile = strrep(stimfile, '.wav', '');
figname = strrep(sprintf('subj:%s  site:%s  sort1:%s  stim:%s', ...
    subjectID,site,sort1,stimfile), '_', '\_');
fullstim = [stimpath '/' stimfile '.wav' ];
[Y,FS]=audioread(fullstim);
[~,~,T,P] = spectrogram(Y,window,nlap,nfft,FS, 'yaxis');

newT = xmin:((xmax-xmin)/length(T)):xmax;
freqs=0:22050/(nfft):22050;
if(length(newT)~=length(newT))
    newT=newT(1:length(T));
end
clim = [-200  -65];
spec = 20*log10(P);
spec = spec(2:end, 1:end);
imagesc(newT,freqs,spec,clim);  
axis xy;
colormap(colormap(jet(256)));

axis([xmin xmax 0 22050]); %11025
xlabel('Time(sec)','FontSize', 20);
ylabel('Freq (Hz)','FontSize', 20);
title('20% Shifted Down in Frequency','FontSize', 24);
% title('Slow-Warped','FontSize', 24);

%% dissect spectrogram plot
figure; 
% ydim = MNE_params.Ndim; 
% xdim = MNE_params.Nlags;
nfft = 128;
step = 0.0194363;

subplot(2,2,1); 
% title('Time-Lagged Stimuli Input', 'FontSize', 20);
colormap(colormap(jet(256)));
newT = -20*step:step:0;
freqs=0:(nfft/2000):22;
imagesc(newT,freqs,pixel(:,:,1000));axis xy;
xlabel('Time(sec)', 'FontSize', 20)
ylabel('Freq (kHz)', 'FontSize', 20)
title('Time-Lagged Stimulus at 9.74s', 'FontSize', 20)
colorbar('ticks',[-1,0,1]);
subplot(2,2,2); 
imagesc(newT,freqs,pixel(:,:,1001));axis xy;
colormap(colormap(jet(256)));
xlabel('Time(sec)', 'FontSize', 20)
ylabel('Freq (kHz)', 'FontSize', 20)
title('Time-Lagged Stimulus at 9.75s', 'FontSize', 20)
colorbar('ticks',[-1,0,1]);
subplot(2,2,3); 
imagesc(newT,freqs,pixel(:,:,1002));axis xy;
colormap(colormap(jet(256)));
xlabel('Time(sec)', 'FontSize', 20)
ylabel('Freq (kHz)', 'FontSize', 20)
title('Time-Lagged Stimulus at 9.76s', 'FontSize', 20)
colorbar('ticks',[-1,0,1]);
subplot(2,2,4); 
imagesc(newT,freqs,pixel(:,:,1003));axis xy;
colormap(colormap(jet(256)));
xlabel('Time(sec)', 'FontSize', 20)
ylabel('Freq (kHz)', 'FontSize', 20)
title('Time-Lagged Stimulus at 9.77s', 'FontSize', 20)
colorbar('ticks',[-1,0,1]);

%% MNE vs QC prediction plot
a=A_meanMNE;
ydim = 16;
xdim = 20;
Nlags=xdim; Ndim=ydim*xdim;
h=h_meanMNE; 
J= reshape(J_meanMNE,Ndim,Ndim); 
start = 0; 
step = 0.0194363;
nfft = 128;
window = 1224;
overlap = 0.3;
freq_compress = 4;

toelist = ...
    'response/concat_chan_19_20_electrode_8_3/10 reps/ss001m_497_17_s1_toe.txt';
stop = 60; %30 for fast, 60 for normal, 90 for slow
[~, psth, spec]= plot_raster_SMI2(toelist, start, stop, step, nfft, window, overlap);
time_compress = 1; %5.74 for fast, 11.5 for normal, 17.2 for slow
[tstimulus, tresponse, tresponse_nor, ~, Nsample, tNdim]=compression(spec, psth', freq_compress, time_compress); 
close all 
[~, Tsample]=size(tresponse);

tNsamples_compressed = Nsample - (Nlags-1); 
tNdimtotal = tNdim*Nlags;
stim_ = zeros(tNdimtotal, tNsamples_compressed);
for i=1:Nlags
    stim_(tNdim*(i-1)+1:tNdim*i,:) = ...
        tstimulus(:,i:tNsamples_compressed+i-1);
end
stim_=stim_';

for b = 1:1:Tsample
     pSpikeMNE = 1./(1+exp(a+stim_*h'+sum(stim_.*(stim_*J),2)));  % Nsamples x 1
end

pSpike_norMNE = pSpikeMNE ./max(pSpikeMNE);

a=A_meanQC;
h=h_meanQC; 
J= reshape(J_meanQC,Ndim,Ndim); 

for b = 1:1:Tsample
     pSpikeQC = 1./(1+exp(a+stim_*h'+sum(stim_.*(stim_*J),2)));  % Nsamples x 1
end

pSpike_norQC = pSpikeQC ./max(pSpikeQC);

hold on; 
xbins = start:((stop-start)/length(tresponse)):stop;
plot(xbins(1,Nlags+1:end),tresponse_nor(:,Nlags:end),'r', 'lineWidth',3 ); 
xbins2 = start:((stop-start)/length(pSpike_norMNE)):stop;
plot(xbins2(1,2:end), pSpike_norMNE,'b', 'lineWidth',3 ); 
xbins3 = start:((stop-start)/length(pSpike_norMNE)):stop;
plot(xbins3(1,2:end), pSpike_norQC,'m', 'lineWidth',3 ); 
axis([start stop 0 1]);
legend({'real spike', 'predicted spike by MNE','predicted spike by QC'}, 'FontSize', 16);
xlabel('Time(sec)', 'FontSize', 28); ylabel('Spike probability', 'FontSize', 28);
% title('Real Sampled and Predicted Spike Response Binned at 97.2ms', 'FontSize', 28);


