%%
clc; close all; 

%% plot
figure
subplot(2,2,1)
corre_coef = [corre_coef18_20; corre_coef5_12; corre_coef7_4; corre_coef8_3];
% corre_coef= sort(corre_coef);
hold on
for i=1:4
    scatter(corre_coef(i,2),corre_coef(i,1),100,'o','filled');
end
xlabel('fast-warped','FontSize',20);
ylabel('normal speed','FontSize',20);
axis([-0.2 0.2 -0.2 0.2]);
% refline(1,0);

subplot(2,2,2)
hold on
corre_coef = [corre_coef18_20; corre_coef5_12; corre_coef7_4; corre_coef8_3];
for i=1:4
    scatter(corre_coef(i,3),corre_coef(i,1),100,'o','filled');
end
xlabel('20% shifted up in frequency','FontSize',20);
ylabel('normal speed','FontSize',20);
axis([-0.2 0.2 -0.2 0.2]);

subplot(2,2,3)
hold on
corre_coef = [corre_coef18_20; corre_coef5_12; corre_coef7_4; corre_coef8_3];
for i=1:4
    scatter(corre_coef(i,4),corre_coef(i,1),100,'o','filled');
end
xlabel('20% shifted down in frequency','FontSize',20);
ylabel('normal speed','FontSize',20)
axis([-0.2 0.2 -0.2 0.2]);

subplot(2,2,4)
hold on
corre_coef = [corre_coef18_20; corre_coef5_12; corre_coef7_4; corre_coef8_3];
for i=1:4
    scatter(corre_coef(i,5),corre_coef(i,1),100,'o','filled');
end
xlabel('slow-warped','FontSize',20);
ylabel('normal speed','FontSize',20)
axis([-0.2 0.2 -0.2 0.2]);


