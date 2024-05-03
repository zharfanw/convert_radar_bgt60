clc;
clear all;

% Load data
data_csv = readmatrix('Dataset\BGT60TR13C_record_20230726-103055_raw_txt.csv');
load('Dataset\ppg_1_txt.mat');
rx1 = data_csv(2:(64*10000)+1,2);

% Parameters
fs = 1000; % Original sampling frequency
fs_new = 2000; % New sampling frequency
f_low = 0.5; % Low frequency for bandpass filter
f_high = 4; % High frequency for bandpass filter

% Resample the radar signal
rx1_resampled = resample(rx1, fs_new, fs);

% Bandpass filter
[b, a] = butter(2, [f_low, f_high] / (fs_new/2), 'bandpass');
rx1_filtered = filter(b, a, rx1_resampled);

% Time vector for resampled signal
t_resampled = (0:length(rx1_filtered)-1) / fs_new;

% FFT
N = length(rx1_filtered);
f = (0:N/2)*(fs_new/N); % Define f only up to N/2+1 points
Y = fft(rx1_filtered);
P2 = abs(Y/N);
P1 = P2(1:N/2+1);
P1(2:end-1) = 2*P1(2:end-1);

% Find the peak
[~, idx] = max(P1(f >= f_low & f <= f_high));
f_peak = f(idx);

% Plot filtered data
figure;
subplot(3,1,1);
plot(t_resampled, rx1_filtered); % Use t_resampled here
grid on
title('Filtered Radar Data');

subplot(3,1,2);
plot(time,ch_2)
grid on
title('PPG_data');

% Plot FFT
subplot(3,1,3);
plot(f, P1); % Now f and P1 have the same length
title('Single-Sided Amplitude Spectrum');
xlabel('Frequency (Hz)');
ylabel('|P1(f)|');

% Plot filtered data in time domain
figure;
plot(t_resampled, rx1_filtered); % Use t_resampled here
title('Filtered Radar Data in Time Domain');
xlabel('Time (s)');
ylabel('Amplitude');

% Display the detected pulse rate
disp(['Detected Pulse Rate: ', num2str(f_peak), ' Hz']);
