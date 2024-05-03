clc;
clear all;


% Load data
data_csv = readmatrix('Dataset\BGT60TR13C_record_20230726-101820_raw_txt.csv');
load('Dataset\ppg_1_txt.mat');
rx1 = data_csv(2:(64*10000)+1,2);

% Parameters
fs = 1000; % Sampling frequency, adjust as per your data
fs_new = 10000; % New sampling frequency
f_low = 0.5; % Low frequency for bandpass filter
f_high = 4; % High frequency for bandpass filter



% Bandpass filter
[b, a] = butter(2, [f_low, f_high] / (fs/2), 'bandpass'); % Use 'bandpass' instead of 'band'
rx1_filtered = filter(b, a, rx1);

% Time vector
t = (0:length(rx1_filtered)-1) / fs;

% Time vector for resampled signal
t_resampled = (0:length(rx1_filtered)-1) / fs_new;

% FFT
N = length(rx1_filtered);
f = (0:N/2)*(fs/N); % Define f only up to N/2+1 points
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
plot(rx1_filtered);
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
plot(t, rx1_filtered);
title('Filtered Radar Data in Time Domain');
xlabel('Time (s)');
ylabel('Amplitude');
 

% Display the detected pulse rate
disp(['Detected Pulse Rate: ', num2str(f_peak), ' Hz']);
