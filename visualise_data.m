%% clear 
clc; %close all; 

%% visualise data
% song = categorical({'m\_497\_17\_s1','m\_497\_17\_s1\_2fast','m\_497\_17\_s1\_shift\_1', ...
%     'm\_497\_17\_s1\_shift\_6', 'm\_497\_17\_s1\_slow'});
% figure; 
% % hold on; 
% data = zeros(5,5);
% for j=1:5
%     data(:,j) = [cc74(j) cc83(j) cc512s8(j) cc512s16(j) cc1820(j)];
% end
% ax1=subplot(2,1,1);bar(song,data');
% title('correlation coefficient');
% ylabel('correlation coefficient');
% legend('electrode 7\_4 trained with 16 songs', 'electrode 8\_3 trained with 16 songs', ...
%     'electrode 5\_12 trained with 8 songs', 'electrode 5\_12 trained with 16 songs', ...
%     'electrode 18\_20 trained with 16 songs');

%% visualise % data
song = categorical({'m\_497\_17\_s1','m\_497\_17\_s1\_2fast','m\_497\_17\_s1\_shift\_1', ...
    'm\_497\_17\_s1\_shift\_6', 'm\_497\_17\_s1\_slow'});
% figure; 
% hold on; 
percdata = zeros(4,5);
for j=1:5
    percdata(:,j) = [perc74(j) perc83(j) perc512(j) perc1820(j)];
end
ax1=subplot(2,1,1); 
bar(song,percdata');
title('percentage correlation coefficient in 20ms resolution');
ylabel('% correlation coefficient');
legend('electrode 7\_4 trained with 16 songs', 'electrode 8\_3 trained with 16 songs', ...
   'electrode 5\_12 trained with 16 songs', ...
    'electrode 18\_20 trained with 16 songs');

