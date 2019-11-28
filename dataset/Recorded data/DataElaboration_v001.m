clc
clear all
Data = readtable('test_16channels_ProbeD_2.xlsx');

figure;
subplot(3,1,1);
plot(Data(:,1),Data(:,2)+1,'.r')

ylim([0.5 16.5]);
xlim([2000 12000]);
set(gca, 'YTick', 1:16)
subplot(3,1,2);
plot(Data(:,1),Data(:,3),'.')
xlim([2000 12000]);
subplot(3,1,3);
plot(Data(:,1),Data(:,4),'.')
xlim([2000 12000]);
figure;
plot(diff(Data(:,1)))