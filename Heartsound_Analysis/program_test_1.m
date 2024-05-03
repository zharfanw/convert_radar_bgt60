clc
clear all
close all
load('Dataset\ppg_1_txt.mat')

% data_csv = readmatrix('Dataset\BGT60TR13C_record_20230726-101820_raw_txt.csv');
data_csv = readmatrix('Dataset\BGT60TR13C_record_20230726-103055_raw_txt.csv');
rx1 = data_csv(2:(64*10000)+1,2);
rx2 = data_csv(2:(64*10000)+1,3);
rx3 = data_csv(2:(64*10000)+1,4);

t = (0:length(rx2)-1) / 81920;

plot_raw_rx = true;
plot_raw_ppg = true;

if(plot_raw_rx)
    figure();
    
    subplot(3,1,1);
    plot(t,rx1)
    xlim([0 2])
    subplot(3,1,2);
    plot(t,rx2)
    subplot(3,1,3);
    plot(t,rx3)
end

if(plot_raw_ppg)
    figure();
    subplot(4,1,1);
    plot(time,ch_1)
    subplot(4,1,2);
    plot(time,ch_2)
    subplot(4,1,3);
    plot(time,ch_3)
    subplot(4,1,4);
    plot(time,ch_4)
end

