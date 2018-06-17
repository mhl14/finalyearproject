% function [J_choose_eigs] = choose_eigs(J_square, Nlags, Nf,number_signif_eigenvals,chosen_eigenvalue)
function [J_choose_eigs] = choose_eigs(J_square, Nlags, Nf,number_signif_eigenvals)
[V,D]=eig(J_square);
eigenvalues = diag(D);
[eigenvalues_sorted,index] = sort(eigenvalues);

% get the number of significant eigenvals, negative and positive, saved
% in the corresponding mat file

    number_neg_eigenvals = number_signif_eigenvals(3);
    number_pos_eigenvals = number_signif_eigenvals(4);


% this part is to reassign values to chosen eigenvalues
% new_D = D;
% for i=1:6
% new_D(index(i),index(i)) = D(index(7),index(7));
% end
% 
% for i=0:5
%     new_D(index(end-i),index(end-i)) = D(index(end-6),index(end-6));
% end

% this part is to instead choose one or few significant eigenvalues, replacing noise (insignificant ones) by
% zeros. This is to generate predictions and determine correlation coefficients.

new_D = zeros(Nlags*Nf,Nlags*Nf);
%%%% here we put the significant positive eigenvalues from D into new_D 
% for i=0:((number_pos_eigenvals-1)+0) % because we start from the end and moving backwards thus first i must be zero to give the index(end); thus numpos-1 and not numpos as the # of + eigenvals
%     new_D(index(end-i),index(end-i)) = D(index(end-i),index(end-i));
% end

% j = chosen_eigenvalue; % choose eigenvalues one by one, e.g. j= 1 chooses the first most negative eigenvalue; j = end chooses the most positive one, etc. Here we pass it from corr_coeff in a loop
j2 = 2; j3 = 3; j5 = 5; j6 = 6;
% new_D(index(j),index(j)) = D(index(j),index(j));
new_D(index(j2),index(j2)) = D(index(j2),index(j2));
new_D(index(j3),index(j3)) = D(index(j3),index(j3));
% new_D(index(j5),index(j5)) = D(index(j5),index(j5));
% new_D(index(j6),index(j6)) = D(index(j6),index(j6));


% %%% choose them in a loop
% for i=1:(320/2)
% new_D(index(i),index(i)) = D(index(i),index(i));
% end


% %%% here we put the significant negative eigenvalues from D into new_D
% for i=1:number_neg_eigenvals+0
% new_D(index(i),index(i)) = D(index(i),index(i));
% end

% new_D(320,320)=eigenvalues_sorted(320);

J_choose_eigs = V*new_D*V';