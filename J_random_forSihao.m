

function J_random_forSihao %generates 500 random matrices with mean and var of J, then extracts the eigenvalues

clear all; close all;


cd('/Volumes/SIHAO/MNE_for_Sihao/MNE/Algorithm_Functions');



Nf = 16; %number of frequency bands in STRF
Nlags = 20; %times the length of a bin gives the time length of STRF
njack   = 4;   % # jackknives to run (also determines the size of each jackknives)


% Now we determine the statistical significance


load st_750_m_211_s16_s17_s25_m_462_s2prime_m_462_s2_4_9_24kHz_Nlags20_nfft128_Nf16_bin21p6ms_jack_1_of_4.mat % do it for all four jacks

a1 = pfinal(1); h1=pfinal(2:Nlags*Nf+1); J1=pfinal(Nlags*Nf+2:end); clear pfinal;

load st_750_m_211_s16_s17_s25_m_462_s2prime_m_462_s2_4_9_24kHz_Nlags20_nfft128_Nf16_bin21p6ms_jack_2_of_4.mat 
a2 = pfinal(1); h2=pfinal(2:Nlags*Nf+1); J2=pfinal(Nlags*Nf+2:end); clear pfinal;

load st_750_m_211_s16_s17_s25_m_462_s2prime_m_462_s2_4_9_24kHz_Nlags20_nfft128_Nf16_bin21p6ms_jack_3_of_4.mat 

a3 = pfinal(1); h3=pfinal(2:Nlags*Nf+1); J3=pfinal(Nlags*Nf+2:end); clear pfinal;

load st_750_m_211_s16_s17_s25_m_462_s2prime_m_462_s2_4_9_24kHz_Nlags20_nfft128_Nf16_bin21p6ms_jack_4_of_4.mat
a4 = pfinal(1); h4=pfinal(2:Nlags*Nf+1); J4=pfinal(Nlags*Nf+2:end); clear pfinal;

a = mean(a1:a4); h = (h1+h2+h3+h4)/4; J = (J1+J2+J3+J4)/4;

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



J_reshaped_zerodiag = tril(J_reshaped,-1)+triu(J_reshaped,1) + diag(zeros(320,1)); % put zeros on the diagonal
J_reshaped_zerodiag_1d = reshape(J_reshaped_zerodiag,[],1); % reshape as 320x1 vector
ind = randperm(size(J_reshaped_zerodiag_1d,1)); %get random integers for permutation indices
J_reshaped_zerodiag_1d_shuffled = J_reshaped_zerodiag_1d(ind); %permute the vector
J_reshaped_zerodiag_1d_shuffled_square = reshape(J_reshaped_zerodiag_1d_shuffled,Nlags*Nf,Nlags*Nf); %make it square
[row,col,v] = find(J_reshaped_zerodiag_1d_shuffled_square); %find nonzero elements, i.e. all nondiagonal elements
zeros_320 = zeros(320,1); % we need zeros to fill in the matrix to keep its size 320x320
matrix_temp = vertcat(v,zeros_320); % reshuffled offdiagonal elements and 320 zeroes, to keep the corect size
matrix_temp = reshape(matrix_temp,320,320); %zeros are in the last column, they will get away since we don't use the upper tiangle
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


figure_name = ('st_750_m_211_s_08_24kHz_Nlags20_nfft128_Nf16_bin21p6ms_statistics');


hgsave(figure_name)



number_of_negative_eigenvalues_vis_J_reshuffled = length(find(eigenvalues_sorted<min(eigenvalues_r)))
number_of_positive_eigenvalues_vis_J_reshuffled = length(find(eigenvalues_sorted>max(eigenvalues_r)))

number_of_negative_eigenvalues_vis_J_random = length(find(eigenvalues_sorted<(min(eig_all))))
number_of_positive_eigenvalues_vis_J_random = length(find(eigenvalues_sorted>(max(eig_all))))

number_signif_eigenvals = [number_of_negative_eigenvalues_vis_J_reshuffled number_of_positive_eigenvalues_vis_J_reshuffled number_of_negative_eigenvalues_vis_J_random number_of_positive_eigenvalues_vis_J_random];




save('st_750_m_211_s_08_24kHz_Nlags20_nfft128_Nf16_bin21p6ms_statistics.mat', 'a', 'h','J','J_r','eig_all','number_signif_eigenvals');



% and now plot the features
figure
subplot(3,3,1)
imagesc(reshape(h,Nf,Nlags))
axis xy
 
subplot(3,3,2)
imagesc(reshape(J,Nlags*Nf,Nlags*Nf))
axis xy
axis square

subplot(3,3,3)
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
axis xy



figure_name = ('st_750_m_211_s_08_24kHz_Nlags20_nfft128_Nf16_bin21p6ms');



hgsave(figure_name)


% now plot as a separate figure 10 most negative eigenvectors
% Note that until st750 the positive and negative eigenvalues were
% interchanged by mistake (positive were ploted as negative and vice versa)

figure

subplot(2,5,1)
eig_sorted_1=V(:,index(1));
imagesc(reshape(eig_sorted_1,Nf,Nlags))
axis xy

subplot(2,5,2)
eig_sorted_2=V(:,index(2));
imagesc(reshape(eig_sorted_2,Nf,Nlags))
axis xy

subplot(2,5,3)
eig_sorted_3=V(:,index(3));
imagesc(reshape(eig_sorted_3,Nf,Nlags))
axis xy

subplot(2,5,4)
eig_sorted_4=V(:,index(4));
imagesc(reshape(eig_sorted_4,Nf,Nlags))
axis xy

subplot(2,5,5)
eig_sorted_5=V(:,index(5));
imagesc(reshape(eig_sorted_5,Nf,Nlags))
axis xy

subplot(2,5,6)
eig_sorted_6=V(:,index(6));
imagesc(reshape(eig_sorted_6,Nf,Nlags))
axis xy

subplot(2,5,7)
eig_sorted_7=V(:,index(7));
imagesc(reshape(eig_sorted_7,Nf,Nlags))
axis xy

subplot(2,5,8)
eig_sorted_8=V(:,index(8));
imagesc(reshape(eig_sorted_8,Nf,Nlags))
axis xy

subplot(2,5,9)
eig_sorted_9=V(:,index(9));
imagesc(reshape(eig_sorted_9,Nf,Nlags))
axis xy

subplot(2,5,10)
eig_sorted_10=V(:,index(10));
imagesc(reshape(eig_sorted_10,Nf,Nlags))
axis xy

figure_name_negative = [figure_name '_10mostnegative_eigenvectors'];
hgsave(figure_name_negative)

% now plot as a separate figure 10 most positive eigenvectors
figure

subplot(2,5,1)
eig_sorted_end=V(:,index(end));
imagesc(reshape(eig_sorted_end,Nf,Nlags))
axis xy

subplot(2,5,2)
eig_sorted_end_1=V(:,index(end-1));
imagesc(reshape(eig_sorted_end_1,Nf,Nlags))
axis xy

subplot(2,5,3)
eig_sorted_end_2=V(:,index(end-2));
imagesc(reshape(eig_sorted_end_2,Nf,Nlags))
axis xy

subplot(2,5,4)
eig_sorted_end_3=V(:,index(end-3));
imagesc(reshape(eig_sorted_end_3,Nf,Nlags))
axis xy

subplot(2,5,5)
eig_sorted_end_4=V(:,index(end-4));
imagesc(reshape(eig_sorted_end_4,Nf,Nlags))
axis xy

subplot(2,5,6)
eig_sorted_end_5=V(:,index(end-5));
imagesc(reshape(eig_sorted_end_5,Nf,Nlags))
axis xy

subplot(2,5,7)
eig_sorted_end_6=V(:,index(end-6));
imagesc(reshape(eig_sorted_end_6,Nf,Nlags))
axis xy

subplot(2,5,8)
eig_sorted_end_7=V(:,index(end-7));
imagesc(reshape(eig_sorted_end_7,Nf,Nlags))
axis xy

subplot(2,5,9)
eig_sorted_end_8=V(:,index(end-8));
imagesc(reshape(eig_sorted_end_8,Nf,Nlags))
axis xy

subplot(2,5,10)
eig_sorted_end_9=V(:,index(end-9));
imagesc(reshape(eig_sorted_end_9,Nf,Nlags))

figure_name_positive = [figure_name '_10mostpositive_eigenvectors'];
hgsave(figure_name_positive)




