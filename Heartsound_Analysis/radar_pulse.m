clc;
clear all;

% If you have a separate FMCW processing function, call it here.
% e.g., processFMCW();

%% Load data
try
    data_csv = readmatrix('Dataset\BGT60TR13C_record_20230726-101820_raw_txt.csv');
    load('Dataset\ppg_1_txt.mat');
    rx1 = data_csv(2:(64*10000)+1,2);
catch e
    disp('Error loading data:');
    disp(e.message);
    return;
end

%% Parameters
fs = 1000;
f_low = 0.5;
f_high = 4;

%% Bandpass filter
try
    [b, a] = butter(2, [f_low, f_high] / (fs/2), 'bandpass');
    rx1_filtered = filter(b, a, rx1);
catch e
    disp('Error in bandpass filter:');
    disp(e.message);
    return;
end

%% Time vector and FFT
try
    [t, P1, f] = performFFT(rx1_filtered, fs);
catch e
    disp('Error in FFT:');
    disp(e.message);
    return;
end

%% Find the peak
try
    f_peak = findPeakFrequency(P1, f, f_low, f_high);
catch e
    disp('Error in finding peak frequency:');
    disp(e.message);
    return;
end

%% Plotting
try
    plotData(rx1_filtered, t, P1, f, time, ch_2);
catch e
    disp('Error in plotting data:');
    disp(e.message);
    return;
end

%% Display the detected pulse rate
disp(['Detected Pulse Rate: ', num2str(f_peak), ' Hz']);

% Additional Functions
function [t, P1, f] = performFFT(rx1_filtered, fs)
    t = (0:length(rx1_filtered)-1) / fs;
    N = length(rx1_filtered);
    f = (0:N/2)*(fs/N);
    Y = fft(rx1_filtered);
    P2 = abs(Y/N);
    P1 = P2(1:N/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
end

function f_peak = findPeakFrequency(P1, f, f_low, f_high)
    [~, idx] = max(P1(f >= f_low & f <= f_high));
    f_peak = f(idx);
end

function plotData(rx1_filtered, t, P1, f, time, ch_2)
    figure;
    subplot(3,1,1);
    plot(rx1_filtered);
    grid on;
    title('Filtered Radar Data');

    subplot(3,1,2);
    plot(time,ch_2);
    grid on;
    title('PPG_data');

    subplot(3,1,3);
    plot(f, P1);
    title('Single-Sided Amplitude Spectrum');
    xlabel('Frequency (Hz)');
    ylabel('|P1(f)|');

    figure;
    plot(t, rx1_filtered);
    title('Filtered Radar Data in Time Domain');
    xlabel('Time (s)');
    ylabel('Amplitude');
end
