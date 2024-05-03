clc;
clear all;

%% Load data
data_csv = readmatrix('Dataset\BGT60TR13C_record_20230726-101820_raw_txt.csv');
load('Dataset\ppg_1_txt.mat');
rx1 = data_csv(2:(64*10000)+1,2);
rx2 = data_csv(2:(64*10000)+1,3);
rx3 = data_csv(2:(64*10000)+1,4);

%% Parameters
fs_radar = 1000; % Original sampling frequency for radar data
fs_ppg = 125; % Original sampling frequency for PPG data
fs_new = 10000; % New common sampling frequency

%% Resample data
% If the original sampling frequencies are not equal, resample the data to a common sampling frequency
if fs_radar ~= fs_ppg
    [P, Q] = rat(fs_new / fs_radar);
    rx1 = resample(rx1, P, Q);
    rx2 = resample(rx2, P, Q);
    rx3 = resample(rx3, P, Q);
    
    [P, Q] = rat(fs_new / fs_ppg);
    ch_2 = resample(ch_2, P, Q);
end

%% Parameters
fs = 1000; % Sampling frequency, adjust as per your data
fs_new = 10000; % New sampling frequency
f_low = 0.5; % Low frequency for bandpass filter
f_high = 4; % High frequency for bandpass filter

%% Bandpass filter
[b, a] = butter(2, [f_low, f_high] / (fs/2), 'bandpass');
rx2_filtered = filter(b, a, rx2);

%% Time vector for Radar
% Create a time vector based on the sampling frequency and length of radar data
time_radar = (0:length(rx2)-1) / fs_radar;

%% Time vector for PPG
% Assuming fs_new is the new sampling frequency for PPG data
time = (0:length(ch_2)-1) / fs_new;

%% Threshold based simple peak finding
threshold = 0.5; % Set a threshold
min_distance = fs; % Minimum distance between peaks

%% Find peaks for radar signal
radar_data_above_threshold = rx2_filtered > threshold;
diff_radar_data = diff(radar_data_above_threshold);
radar_peak_starts = find(diff_radar_data == 1);
radar_peak_ends = find(diff_radar_data == -1);
radar_peak_locs = radar_peak_starts + round((radar_peak_ends - radar_peak_starts) / 2);

%% Remove peaks that are too close to each other
radar_peak_locs(diff([0; radar_peak_locs]) < min_distance) = [];

%% Find peaks for PPG signal
ppg_data_above_threshold = ch_2 > threshold;
diff_ppg_data = diff(ppg_data_above_threshold);
ppg_peak_starts = find(diff_ppg_data == 1);
ppg_peak_ends = find(diff_ppg_data == -1);
ppg_peak_locs = ppg_peak_starts + round((ppg_peak_ends - ppg_peak_starts) / 2);

%% Remove peaks that are too close to each other
ppg_peak_locs(diff([0; ppg_peak_locs]) < min_distance) = [];

%% Calculate the PAT for each pair of peaks
num_peaks = min(length(radar_peak_locs), length(ppg_peak_locs));
pat = (radar_peak_locs(1:num_peaks) - ppg_peak_locs(1:num_peaks)) / fs;

%% Plot filtered data
figure;
subplot(4,1,1);
plot(rx2_filtered);
hold on;
plot(radar_peak_locs, rx2_filtered(radar_peak_locs), 'r*');
grid on;
title('Filtered Radar Data');

% Plot PPG_data
subplot(4,1,2);
plot(time, ch_2);
hold on;
plot(ppg_peak_locs / fs_new, ch_2(ppg_peak_locs), 'r*'); % Adjust the time for peak locations
grid on;
title('PPG_data');

% Plot PAT
subplot(4,1,3);
plot(pat, 'o-');
grid on;
title('Pulse Arrival Time (PAT)');
xlabel('Peak number');
ylabel('PAT (seconds)');

%% Plot Radar Data in Time Domain
figure;
plot(time_radar, rx2_filtered);
title('Radar Data in Time Domain');
xlabel('Time (s)');
ylabel('Amplitude');

% Display the calculated average PAT
disp(['Average Calculated Pulse Arrival Time: ', num2str(mean(pat)), ' seconds']);
