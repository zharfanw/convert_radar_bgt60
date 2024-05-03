clc;
clear all;

%% Load data
% [Your existing data loading code...]
data_csv = readmatrix('Dataset\BGT60TR13C_record_20230726-101820_raw_txt.csv');
load('Dataset\ppg_1_txt.mat');
sample_multiply = 8192;
rx1 = data_csv(2:(64*sample_multiply)+1,2);
rx2 = data_csv(2:(64*sample_multiply)+1,3);
rx3 = data_csv(2:(64*sample_multiply)+1,4);

%% Parameters
fs = 1000; 
f_low = 0.5; 
f_high = 4; 

%% Bandpass filter
[b, a] = butter(2, [f_low, f_high] / (fs/2), 'bandpass'); 
rx2_filtered = filter(b, a, rx2);

%% Time vector
t = (0:length(rx2_filtered)-1) / fs;

%% FFT
% [Your existing FFT code...]

%% Plot filtered data and FFT
% [Your existing plotting code...]

%% PPG DATA PROCESSING
% Assuming ch_2 is your PPG data
ppg_data = ch_2;

% Preprocessing: Bandpass filter
filtered_ppg = filter(b, a, ppg_data);

% Ensure t and filtered_ppg have the same length
t_ppg = (0:length(filtered_ppg)-1) / fs;

% Peak Detection
threshold_ppg = 0.5; % [Adapt] Define a suitable threshold for your data
min_distance_ppg = 100; % [Adapt] Define a suitable minimum peak distance for your data
[pks_ppg, locs_ppg] = findpeaks(filtered_ppg, 'MinPeakHeight', threshold_ppg, 'MinPeakDistance', min_distance_ppg);

% Interval Calculation
pp_intervals_ppg = diff(locs_ppg) / fs; 

% Pulse Variability Analysis: Example with SDNN
sdnn_ppg = std(pp_intervals_ppg);

%% RADAR DATA PROCESSING
% Using rx2_filtered as your radar data
threshold_radar = 0.5; % [Adapt] Define a suitable threshold for your data
min_distance_radar = 100; % [Adapt] Define a suitable minimum peak distance for your data
[pks_radar, locs_radar] = findpeaks(rx2_filtered, 'MinPeakHeight', threshold_radar, 'MinPeakDistance', min_distance_radar);

% Interval Calculation
pp_intervals_radar = diff(locs_radar) / fs; 

% Pulse Variability Analysis: Example with SDNN
sdnn_radar = std(pp_intervals_radar);


%% Plotting PPG and Radar data with detected peaks
figure;
subplot(2,1,1);
plot(t_ppg, filtered_ppg);
hold on;
plot(t_ppg(locs_ppg), filtered_ppg(locs_ppg), 'rx');
title('Filtered PPG Data and Detected Peaks');
xlabel('Time (s)');
ylabel('Amplitude');

subplot(2,1,2);
plot(t, rx2_filtered);
hold on;
plot(t(locs_radar), rx2_filtered(locs_radar), 'rx');
title('Filtered Radar Data and Detected Peaks');
xlabel('Time (s)');
ylabel('Amplitude');

%% Bland-Altman Plot
% [Your existing Bland-Altman Plot code...]

%% Display the detected pulse rate
% [Your existing code for displaying detected pulse rate...]

%% Display Pulse Variability
disp(['PPG SDNN: ', num2str(sdnn_ppg), ' s']);
disp(['Radar SDNN: ', num2str(sdnn_radar), ' s']);
