% Load the data from CSV file
data_csv = readmatrix('Dataset\BGT60TR13C_record_20230726-101820_raw_txt.csv');

% Extract received signals from the CSV data
rx1 = data_csv(2:(64*10000)+1,2);
rx2 = data_csv(2:(64*10000)+1,3);
rx3 = data_csv(2:(64*10000)+1,4);

% Parameters
f_start = 59.9e9; % Start frequency
f_stop = 60.1e9; % Stop frequency
T = 1e-3; % Pulse duration
c = 3e8; % Speed of light

% Time vector
t = linspace(0, T, length(rx1));

% Generate Chirp Signal
k = (f_stop - f_start) / T; % Chirp rate
transmitted_signal = exp(1j * 2 * pi * (f_start * t + 0.5 * k * t.^2));

% Process each received signal
for rx = {rx1, rx2, rx3}
    received_signal = rx{1};
    
    % Signal Processing
    mixed_signal = transmitted_signal .* conj(received_signal);
    processed_signal = abs(fft(mixed_signal));
    
    % Find the delay
    [~, idx] = max(processed_signal);
    calculated_delay = idx * T / length(rx1);
    
    % Calculate the distance
    distance = c * calculated_delay / 2;
    
    % Display the result
    disp(['Calculated Distance: ' num2str(distance) ' m']);
end
