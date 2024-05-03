clc;
clear all;

%% Load data
data_csv = readmatrix('Dataset\BGT60TR13C_record_20230726-102220_raw_txt');
load('Dataset\ppg_1_txt.mat');
sample_multiply = 8192;
rx1 = data_csv(2:(64*sample_multiply)+1,2);
rx2 = data_csv(2:(64*sample_multiply)+1,3);
rx3 = data_csv(2:(64*sample_multiply)+1,4);

%% Parameters
% fs = 1000; % Sampling frequency, adjust as per your data
fs = 1000; % Sampling frequency, adjust as per your data
% fs = 81920;
fs_new = 1000; % New sampling frequency
f_low = 0.5; % Low frequency for bandpass filter
f_high = 4; % High frequency for bandpass filter

%% Bandpass filter
[b, a] = butter(2, [f_low, f_high] / (fs/2), 'bandpass'); % Use 'bandpass' instead of 'band'
rx2_filtered = filter(b, a, rx2);

%% Time vector
% t = (0:length(rx2_filtered)-1) / fs;
t = (0:length(rx2_filtered)-1) / 81920;

% Time vector for resampled signal
t_resampled = (0:length(rx2_filtered)-1) / fs_new;

% FFT
N = length(rx2_filtered);
f = (0:N/2)*(fs/N); % Define f only up to N/2+1 points
Y = fft(rx2_filtered);
P2 = abs(Y/N);
P1 = P2(1:N/2+1);
P1(2:end-1) = 2*P1(2:end-1);


% Find the peak
[~, idx] = max(P1(f >= f_low & f <= f_high));
f_peak = f(idx);

% Plot filtered data
figure;
subplot(3,1,1);
plot(rx2_filtered);
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

% Plot Wrist Pulse
    %figure;
   % subplot(4,1,1);
   % plot(time,ch_1)
   % subplot(4,1,2);
   % plot(time,ch_2)
   % subplot(4,1,3);
   % plot(time,ch_3)
   % subplot(4,1,4);
    %plot(time,ch_4)

    % Plot filtered data in time domain
figure;
subplot (2,1,1)
plot(t, rx2_filtered);
title('Filtered Radar Data in Time Domain');
xlabel('Time (s)');
ylabel('Amplitude');
 
subplot (2,1,2)
plot(time,ch_2)
grid off
title('PPG_data');
xlabel('Time (s)');
ylabel('Amplitude');
xlim([0 25])

%% Bland-Altman Plot
% Assuming rx2_filtered and ch_2 are the measurements you want to compare
% Ensure they have the same length before proceeding
%if length(rx2_filtered) ~= length(ch_2)
  %  error('Measurement vectors must be of equal length for Bland-Altman plot.');
%end

%means = (rx2_filtered + ch_2) / 2;
%differences = rx2_filtered - ch_2;

%figure;
%plot(means, differences, 'o');
%xlabel('Mean of measurements');
%ylabel('Difference between measurements');
%title('Bland-Altman Plot');

%mean_diff = mean(differences);
%std_diff = std(differences);
%hold on;
%plot([min(means), max(means)], [mean_diff, mean_diff], 'k--');
%plot([min(means), max(means)], [mean_diff + 1.96*std_diff, mean_diff + 1.96*std_diff], 'r--');
%plot([min(means), max(means)], [mean_diff - 1.96*std_diff, mean_diff - 1.96*std_diff], 'r--');
%hold off;

% Display the detected pulse rate
%disp(['Detected Pulse Rate: ', num2str(f_peak), ' Hz']);
