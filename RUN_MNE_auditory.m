function [A_mean, J_mean, h_mean, V, D, MNE_params, eigenvalues_sorted] = RUN_MNE_auditory(stimulus_compressed, response_compressed)
%% Add Functions folder to search path
addpath(genpath('UROP/Functions')); 
addpath(genpath('freq_shifted_song')); 
clc; close all; 

%% MNE parameters setup
MNE_params.Ndim = 16; %number of frequency bands in STRF 
MNE_params.Nlags = 20; %times the length of a bin gives the time length of STRF
% # time lags to use in spatiotemporal models (=1 is just a spatial model)
MNE_params.order = 2;   % order of MNE model to fit (1 or 2)
MNE_params.fittype = 0;   % 0 for regular fitting, 1 for random fitting
MNE_params.Njack = 4; % Number of jackknives to use, 4 per songs 

%% Fit MNE model
[A_mean, J_mean, h_mean]= MNE(stimulus_compressed, response_compressed, MNE_params);
[V,D] = eig(reshape(J_mean,MNE_params.Nlags*MNE_params.Ndim,MNE_params.Nlags*MNE_params.Ndim));

%% Get confidence intervals for eigenvalues

Nf = 16; %number of frequency bands in STRF
Nlags = 20; %times the length of a bin gives the time length of STRF
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

%% Subplot results 
clc; close all; 
ydim = MNE_params.Ndim; 
xdim = MNE_params.Nlags;

figure;
subplot(3,3,1)
imagesc(reshape(h_mean,ydim,xdim));
axis xy; 
colormap(colormap(jet(256)));
xlabel('time band')
ylabel('frequency band')
title('linear parameter h')

subplot(3,3,2:3)
hold on
eigenvalues = diag(D);
[eigenvalues_sorted,index] = sort(eigenvalues);
plot(eigenvalues_sorted,'o');
xlabel('Feature')
ylabel('Magnitude')
title('Eigenvalues sorted with confidence interval generated with noise')

plot(eigenvalues_sorted_r_alltogether,'b');
plot(eigenvalues_sorted_r,'r');
refline(0,min_eig_all)
refline(0,max_eig_all)
hold off 

subplot(3,3,4)
eig_sorted_1 = V(:,index(1));
imagesc(reshape(eig_sorted_1,ydim,xdim))
axis xy
colormap(colormap(jet(256)));
title('smallest eigenvalue')

subplot(3,3,5)
eig_sorted_2=V(:,index(2));
imagesc(reshape(eig_sorted_2,ydim,xdim))
axis xy
colormap(colormap(jet(256)));
title('second smallest eigenvalue')

subplot(3,3,6)
eig_sorted_3=V(:,index(3));
imagesc(reshape(eig_sorted_3,ydim,xdim))
axis xy
colormap(colormap(jet(256)));
title('third smallest eigenvalue')

subplot(3,3,7)
eig_sorted_end=V(:,index(end));
imagesc(reshape(eig_sorted_end,ydim,xdim))
axis xy
colormap(colormap(jet(256)));
title('largest eigenvalue')

subplot(3,3,8)
eig_sorted_end_1=V(:,index(end-1));
imagesc(reshape(eig_sorted_end_1,ydim,xdim))
axis xy
colormap(colormap(jet(256)));
title('second largest eigenvalue')

subplot(3,3,9)
eig_sorted_end_2=V(:,index(end-2));
imagesc(reshape(eig_sorted_end_2,ydim,xdim))
axis xy
colormap(colormap(jet(256)));
title('third largest eigenvalue')

%% Plot receptive fields of 8 most significant eigenvalues
%{
figure;
title ('receptive fields of 8 most positive and negative significant eigenvalues')
subplot(2,8,1)
eig_sorted_1 = V(:,index(1));
imagesc(reshape(eig_sorted_1,ydim,xdim))
axis xy
title('smallest')

subplot(2,8,2)
eig_sorted_2=V(:,index(2));
imagesc(reshape(eig_sorted_2,ydim,xdim))
axis xy
title('second smallest')

subplot(2,8,3)
eig_sorted_3=V(:,index(3));
imagesc(reshape(eig_sorted_3,ydim,xdim))
axis xy
title('third smallest')

subplot(2,8,4)
eig_sorted_4=V(:,index(4));
imagesc(reshape(eig_sorted_4,ydim,xdim))
axis xy
title('fourth smallest')

subplot(2,8,5)
eig_sorted_5=V(:,index(5));
imagesc(reshape(eig_sorted_5,ydim,xdim))
axis xy
title('fifth smallest')

subplot(2,8,6)
eig_sorted_6=V(:,index(6));
imagesc(reshape(eig_sorted_6,ydim,xdim))
axis xy
title('sixth smallest')

subplot(2,8,7)
eig_sorted_7=V(:,index(7));
imagesc(reshape(eig_sorted_7,ydim,xdim))
axis xy
title('seventh smallest')

subplot(2,8,8)
eig_sorted_8=V(:,index(8));
imagesc(reshape(eig_sorted_8,ydim,xdim))
axis xy
title('eighth smallest')

subplot(2,8,9)
eig_sorted_end=V(:,index(end));
imagesc(reshape(eig_sorted_end,ydim,xdim))
axis xy
title('largest')

subplot(2,8,10)
eig_sorted_end_1=V(:,index(end-1));
imagesc(reshape(eig_sorted_end_1,ydim,xdim))
axis xy
title('second largest')

subplot(2,8,11)
eig_sorted_end_2=V(:,index(end-2));
imagesc(reshape(eig_sorted_end_2,ydim,xdim))
axis xy
title('third largest')

subplot(2,8,12)
eig_sorted_end_3=V(:,index(end-3));
imagesc(reshape(eig_sorted_end_3,ydim,xdim))
axis xy
title('fourth largest')

subplot(2,8,13)
eig_sorted_end_4=V(:,index(end-4));
imagesc(reshape(eig_sorted_end_4,ydim,xdim))
axis xy
title('fifth largest')

subplot(2,8,14)
eig_sorted_end_5=V(:,index(end-5));
imagesc(reshape(eig_sorted_end_5,ydim,xdim))
axis xy
title('sixth largest')

subplot(2,8,15)
eig_sorted_end_6=V(:,index(end-6));
imagesc(reshape(eig_sorted_end_6,ydim,xdim))
axis xy
title('seventh largest')

subplot(2,8,16)
eig_sorted_end_7=V(:,index(end-7));
imagesc(reshape(eig_sorted_end_7,ydim,xdim))
axis xy
title('eighth largest')
%}

end 