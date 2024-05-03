clc;
clear all;

% Add the EMD toolkit to the MATLAB path
%addpath('path/to/emd_toolkit');

%% Load data
data_csv = readmatrix('Dataset\BGT60TR13C_record_20230726-101820_raw_txt.csv');
load('Dataset\ppg_1_txt.mat');
rx1 = data_csv(2:(64*10000)+1,2);
rx2 = data_csv(2:(64*10000)+1,3);
rx3 = data_csv(2:(64*10000)+1,4);

%% Parameters
fs = 1000; % Sampling frequency, adjust as per your data
fs_new = 10000; % New sampling frequency
f_low = 0.5; % Low frequency for bandpass filter
f_high = 4; % High frequency for bandpass filter

%% Bandpass filter
[b, a] = butter(2, [f_low, f_high] / (fs/2), 'bandpass');
rx2_filtered = filter(b, a, rx2);

%% EEMD
% Ensure the length of rx2_filtered is an integer
rx2_filtered = rx2_filtered(1:floor(end));

% Set EEMD parameters
ensemble_size = 100; % Number of ensembles (integer)
noise_strength = 0.2; % Noise strength

% Example input arguments
y = rx2_filtered; % your data
Nstd = 0.2; % noise strength
NR = 100; % number of ensembles
aim = 10; % maximum number of intrinsic mode functions (IMFs)

% Calling the eemd1 function with the example input arguments
IMFs = eemd1(y, Nstd, NR, aim);

%% Time vector
t = (0:length(rx2_filtered)-1) / fs;


% Plot filtered data
figure;
subplot(4,1,1); % Change to 4,1,1 to accommodate EEMD plot
plot(rx2_filtered);
grid on;
title('Filtered Radar Data');

subplot(4,1,2);
plot(time,ch_2);
grid on;
title('PPG_data');

% Plot EEMD IMFs
subplot(4,1,3);
for i = 1:size(IMFs,1)
    plot(t, IMFs(i,:));
    hold on;
end
grid on;
title('EEMD IMFs');

% Plot FFT
subplot(4,1,4); % Change to 4,1,4 to accommodate EEMD plot
plot(f, P1);
title('Single-Sided Amplitude Spectrum');
xlabel('Frequency (Hz)');
ylabel('|P1(f)|');

% ...
