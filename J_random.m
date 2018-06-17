

function [max_eig_all, min_eig_all] = J_random() %generates 500 random matrices with mean and var of J, then extracts the eigenvalues

clear all; close all;


% cd('/Users/andreikozlov/Work/bird/Analysis/sortedData/st750/sortedData/chan_19_20');


% cellnum = 1101;
% cellnum = 750;
Nf = 16; %number of frequency bands in STRF
Nlags = 16; %times the length of a bin gives the time length of STRF
njack   = 4;   % # jackknives to run (also determines the size of each jackknives)


% Now we determine the statistical significance


file = load('D:\Sihao\2016-04-04_18-22-02\_binary_spike_trace_trode_1_cluster_no_15.mat_thresh_0.4_157356samples4jacks_16x16_.mat') % do it for all four jacks

J = mean(file.J);
h = mean(file.h);
pfinal = mean(file.pfinal_all_jacks);

a1 = pfinal(1); h1=pfinal(2:Nlags*Nf+1); J1=pfinal(Nlags*Nf+2:end); clear pfinal;


J_reshaped = reshape(J,Nlags*Nf,Nlags*Nf);
 
[V,D] = eig(J_reshaped);

eigenvalues = diag(D);
[eigenvalues_sorted,index] = sort(eigenvalues);

%plot them
figure
plot(eigenvalues_sorted,'o');
hold on

%now we reshuffle J and save it in J_r
J_r = zeros(size(J_reshaped)); % initialize square matrix of size J (currently 320x320)
J_r_alltogether = zeros(size(J_reshaped));

ind_r = randperm(length(J)); % generate 320 random integers
J_r_alltogether = J(ind_r); % permute J elements using those indices
J_r_alltogether = reshape(J_r_alltogether,Nlags*Nf,Nlags*Nf); %make it square
 
J_r_alltogether = diag(diag(J_r_alltogether)) + triu(J_r_alltogether,1) + triu(J_r_alltogether,1)'; %here we symmetrize this matrix J_r_alltogether; here we scramble diagonal and non-diagonal elements together


% now we permute diagonal and non-diagonal elements separately, for comparizon

ind_diag = randperm(size(J_reshaped,1)); %randomly permute 320 integers from 1 to 320
J_diagonal = diag(J_reshaped); %get the diagonal
J_r_diagonal = J_diagonal(ind_diag); %permute the diagonal of J_reshaped using those indices



J_reshaped_zerodiag = tril(J_reshaped,-1)+triu(J_reshaped,1) + diag(zeros(length(J_reshaped),1)); % put zeros on the diagonal
J_reshaped_zerodiag_1d = reshape(J_reshaped_zerodiag,[],1); % reshape as 320x1 vector
ind = randperm(size(J_reshaped_zerodiag_1d,1)); %get random integers for permutation indices
J_reshaped_zerodiag_1d_shuffled = J_reshaped_zerodiag_1d(ind); %permute the vector
J_reshaped_zerodiag_1d_shuffled_square = reshape(J_reshaped_zerodiag_1d_shuffled,Nlags*Nf,Nlags*Nf); %make it square
[row,col,v] = find(J_reshaped_zerodiag_1d_shuffled_square); %find nonzero elements, i.e. all nondiagonal elements
zeros_lenJ = zeros(length(J_reshaped),1); % we need zeros to fill in the matrix to keep its size 320x320
matrix_temp = vertcat(v,zeros_lenJ); % reshuffled offdiagonal elements and 320 zeroes, to keep the corect size
matrix_temp = reshape(matrix_temp,length(J_reshaped),length(J_reshaped)); %zeros are in the last column, they will get away since we don't use the upper tiangle
J_r = tril(matrix_temp,-1)+tril(matrix_temp,-1)'+diag(J_r_diagonal); % reconstruct the final reshuffled matrix by putting the reshuffled lower triangles (so it's square) and the reshuffled diagonal


%take its eigenvalues of J_r_alltogether, and sort them

[V_r_alltogether,D_r_alltogether] = eig(J_r_alltogether);
eigenvalues_r_alltogether = diag(D_r_alltogether);
[eigenvalues_sorted_r_alltogether,indx_alltogether] = sort(eigenvalues_r_alltogether);
plot(eigenvalues_sorted_r_alltogether,'b');

%take its eigenvalues of this J_r, and sort them
[V_r,D_r] = eig(J_r);
eigenvalues_r = diag(D_r);
 [eigenvalues_sorted_r,indx] = sort(eigenvalues_r);
plot(eigenvalues_sorted_r,'r');


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
    
 J_random = mean(J)+std(J).*randn(size(J)); %here we generate a random matrix of size, mean and variance equal to those of J
 
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
refline(0,min(eig_all))
refline(0,max(eig_all))

max_eig_all = max(eig_all);
min_eig_all = min(eig_all);






