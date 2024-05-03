close all
clear all
clc


clear;

% [retVal] = readDCA1000_1('./04060010_1573837157_Raw_0.bin');

% BGT60TR13C_record_20230726-102220.raw.txt x ppg_3_txt.mat
% BGT60TR13C_record_20230726-101820.raw.txt x ppg_2_txt.mat
% BGT60TR13C_record_20230726-102748.raw.txt x ppg_4_txt.mat

% [FileName,PathName] = uigetfile('*.txt','Select the .txt file from radar');
% load('E:\XSDUX\Heartsound\Heartsound_Analysis\Dataset\ppg_2_txt.mat');
% 
% [Radar_Parameter,Frame_Number,NumRXAntenna,Frame]=data_import([PathName FileName]);

% 'E:\XSDUX\Heartsound\Heartsound_Analysis\Dataset\Wrist Pulse Only 2672023\
% % Data Pair 1
load('E:\XSDUX\Heartsound\Heartsound_Analysis\Dataset\ppg_3_txt.mat');
%[Radar_Parameter,Frame_Number,NumRXAntenna,Frame]=data_import("E:\XSDUX\Heartsound\Heartsound_Analysis\Dataset\Wrist Pulse Only 2672023\BGT60TR13C_record_20230726-102748.raw.txt");

[Radar_Parameter,Frame_Number,NumRXAntenna,Frame]=data_import("E:\XSDUX\Radar Data\60 GHz chest\BGT60TR13C_record_20230726-130841.raw.txt");


% % Data Pair 2
% load('E:\XSDUX\Heartsound\Heartsound_Analysis\Dataset\ppg_2_txt.mat');
% [Radar_Parameter,Frame_Number,NumRXAntenna,Frame]=data_import("E:\XSDUX\Heartsound\Heartsound_Analysis\Dataset\Wrist Pulse Only 2672023\BGT60TR13C_record_20230726-101820.raw.txt");


% Data Pair 3
% load('E:\XSDUX\Heartsound\Heartsound_Analysis\Dataset\ppg_4_txt.mat');
% [Radar_Parameter,Frame_Number,NumRXAntenna,Frame]=data_import("E:\XSDUX\Heartsound\Heartsound_Analysis\Dataset\Wrist Pulse Only 2672023\BGT60TR13C_record_20230726-102748.raw.txt");

%% Signal Processing
c = 3e8; % Speed of light (m/s)
fc = (Radar_Parameter.Lower_RF_Frequency_kHz+Radar_Parameter.Upper_RF_Frequency_kHz)/2*1000; % Center frequency (Hz)
CPI = 1/250; % Coherent processing interval (s). Seconds to be processed (CPI must be lower to N*fs)
CRR = 1/Radar_Parameter.Chirp_Time_sec; % Chirp repetition rate (Hz)
FRR=1/Radar_Parameter.Frame_Period_sec;% Frame repetition rate (Hz)
BW = (Radar_Parameter.Upper_RF_Frequency_kHz-Radar_Parameter.Lower_RF_Frequency_kHz)*1000; % Bandwidth (Hz)
gamma = BW*CRR; % Chirp rate (Hz/s)
f0=60000000; % initial frequency
lambda=c/f0; %radar signal wavelength
range_res = c/(2*BW);
max_range = range_res*fix(Radar_Parameter.Sampling_Frequency_kHz*1e3/CRR)/2;
global NL
NL = 1024;


%% FFT
for RXAntenna = 1: NumRXAntenna  % RX 
    
    tau = (0:(Frame_Number-1))/FRR; % Slow time (s)
    raw_data_matrix =zeros(Frame_Number,Radar_Parameter.Samples_per_Chirp);
    for FrameN = 1:Frame_Number
        raw_data_matrix((FrameN),:)= Frame(:,1,RXAntenna,FrameN)';
    end
    
    % subtract DC
    avgDC=mean(raw_data_matrix,2);               
    for jj = 1:size(raw_data_matrix,1)
        raw_data_matrix(jj,:) = raw_data_matrix(jj,:) - avgDC(jj);
    end
    
    [a,b]= size(raw_data_matrix);
    eje_dis = (0:NL-1)/NL*c*Radar_Parameter.Sampling_Frequency_kHz*1000/2/gamma;
    win=rectwin(Radar_Parameter.Samples_per_Chirp);
    win_2=win(:,ones(Frame_Number,1)); %add window
    raw_data_matrix_2 = raw_data_matrix.*win_2';
    per = fft(raw_data_matrix_2,NL,2); % Range profiles
    raw_data_matrix_antenna(:,:,RXAntenna)=raw_data_matrix_2;
    % MTI处理
    per_MTI=diff(per,1,1);
    per_data_matrix_antenna(:,:,RXAntenna)=per_MTI;    
end

global numChirps;
global numADCSamples;%Sampling points
numChirps = Frame_Number*Radar_Parameter.Chirps_per_Frame
numADCSamples = Radar_Parameter.Samples_per_Chirp
% RX1data = reshape(retVal(1,:),numADCSamples,numChirps);   %RX1
% RX2data = reshape(retVal(2,:),numADCSamples,numChirps);   %RX2
% RX3data = reshape(retVal(3,:),numADCSamples,numChirps);   %RX3
% RX4data = reshape(retVal(4,:),numADCSamples,numChirps);   %RX4

RX1data = transpose(raw_data_matrix_antenna(:,:,1));  %RX1
RX2data = transpose(raw_data_matrix_antenna(:,:,2));   %RX2
RX3data = transpose(raw_data_matrix_antenna(:,:,3));   %RX3
% RX4data = reshape(retVal(4,:),numADCSamples,numChirps);   %RX4

% c=3.0e8;  
% slope=60e12;   %FM slope
% Tc=50e-6;      %chirpÖÜÆÚ
% B=slope*Tc;    %FM bandwidth
% Fs=4e6;        %Sampling Rate
% f0=60.36e9;    %initial frequency
% lambda=c/f0;   %radar signal wavelength
% d=lambda/2;    %Antenna Array Spacing
% frame=400;     %number of frames
% Tf=0.05;       %frame period
% N=1024;        %FFT

c=3.0e8;  
slope=(Radar_Parameter.Upper_RF_Frequency_kHz-Radar_Parameter.Lower_RF_Frequency_kHz)*1e3;   %FM slope
Tc=Radar_Parameter.Frame_Period_sec/(Radar_Parameter.Samples_per_Frame/3);      %chirp
B=BW;    %FM bandwidth
Fs=Radar_Parameter.Sampling_Frequency_kHz*1e3;       %Sampling Rate
% f0=f0;    %initial frequency
lambda=c/f0;   %radar signal wavelength
d=lambda/2;    %Antenna Array Spacing
frame=Frame_Number/2;     %number of frames
Tf=Radar_Parameter.Frame_Period_sec;       %frame period
N=1024;        %FFTµãÊý

%hamming window
range_win = hamming(numADCSamples); %generate hamming window
for k=1:1:frame
    din_win(:,k)=RX1data(:,2*k-1).*range_win; %do to the signal Range-fft
    datafft(:,k)=fft(din_win(:,k));
end

%Find Range-bin peaks
rangeBinStartIndex=3;%·Ö±æÂÊ0.1m
rangeBinEndIndex=10;
for k=1:1:frame
    for j=rangeBinStartIndex:1:numADCSamples 
        if(abs(datafft(j,k))==max(abs(datafft((rangeBinStartIndex:rangeBinEndIndex),k)))) 
            data(:,k)=datafft(j,k);
        end
    end
end

%Get the real and imaginary parts of the signal
for k=1:frame
    data_real(:,k)=real(data(:,k));
    data_imag(:,k)=imag(data(:,k));
end

%Calculate signal phase
for k=1:frame
    signal_phase(:,k)=atan(data_imag(:,k)/data_real(:,k));
end

%Phase unwrapping
for k=2:frame
    diff=signal_phase(:,k)-signal_phase(:,k-1);
    if diff>pi/2
        signal_phase(:,(k:end))=signal_phase(:,(k:end))-pi;
    elseif diff<-pi/2
        signal_phase(:,(k:end))=signal_phase(:,(k:end))+pi;
    end
end

%Calculate the phase difference
for k=1:frame-1
    delta_phase(:,k)=signal_phase(:,k+1)-signal_phase(:,k);
end

%Original signal bandpass filtering
thresh=0.8;
for k=1:frame-3
    phaseUsedComputation(:,k)=filter_RemoveImpulseNoise(delta_phase(:,k),delta_phase(:,k+1),delta_phase(:,k+2),thresh);
end
index=1:1:frame-3;
index=index*Tf;

%Original signal bandpass filtering
filter_delta_phase=filter(bpf_vitalsign,phaseUsedComputation);

vital_sign=filter_delta_phase;

%Original Signal Time Domain Plot
figure(1);
plot(index,vital_sign);
xlabel('Time(s)','FontWeight','bold');
ylabel('Amplitude','FontWeight','bold');
title('cardiopulmonary signal','FontWeight','bold');

%do fft on the original signal
vital_sign_fft=fft(vital_sign,N);

%Convert double sideband signal to single sideband
freq=(0:1:N/2)/Tf/N;  %The vital sign signal sampling rate is the number of frames
P2 = abs(vital_sign_fft/(N-1));
P1 = P2(1:N/2+1);   %Select the first half at this time, because after fft is a symmetrical bilateral spectrum
P1(2:end-1) = 2*P1(2:end-1);

%Original signal frequency domain plot
figure(2);
plot(freq,P1);
xlim([0,2]);
xlabel('Frequency(Hz)','FontWeight','bold');
ylabel('Amplitude','FontWeight','bold');
title('Cardiopulmonary signal spectrogram','FontWeight','bold');

% %Breathing signal bandpass filtering
filter_delta_phase_breathe=filter(bpf_breathe,phaseUsedComputation);

breathe=filter_delta_phase_breathe;

%Time Domain Diagram of Respiration Signal
figure(3);
plot(index,breathe);
xlabel('Time(s)','FontWeight','bold');
ylabel('Amplitude','FontWeight','bold');
title('breathing signal','FontWeight','bold');

%Do fft on respiration signal
breathe_fft=fft(breathe,N);

%Convert double sideband signal to single sideband
P2_breathe = abs(breathe_fft/(N-1));
P1_breathe = P2_breathe(1:N/2+1);   %Select the first half at this time, because after fft is a symmetrical bilateral spectrum
P1_breathe(2:end-1) = 2*P1_breathe(2:end-1);

%Respiratory signal frequency domain diagram
figure(4);
plot(freq,P1_breathe);
xlim([0,2]);
xlabel('Frequency(Hz)','FontWeight','bold');
ylabel('Amplitude','FontWeight','bold');
title('Respiratory signal spectrogram','FontWeight','bold');


%Heartbeat signal bandpass filtering
filter_delta_phase_heart=filter(bpf_heart,phaseUsedComputation);

heart=filter_delta_phase_heart;

%Heartbeat signal time domain diagram
% figure(5);
figure("Name","Selected Heart Signal x Time Domain");
plot(index,heart);
% xlim([0 20])
xlabel('Time(s)','FontWeight','bold');
ylabel('Amplitude','FontWeight','bold ');
title('heartbeat signal','FontWeight','bold');

% figure(6);
figure("Name","Selected Heart Signal ABS x Time Domain");
plot(index,abs(heart));
% xlim([0 20])
xlabel('Time(s)','FontWeight','bold');
ylabel('Amplitude','FontWeight','bold ');
title('heartbeat signal','FontWeight','bold');

% Do fft on the heartbeat signal
heart_fft=fft(heart,N);

%Convert double sideband signal to single sideband
P2_heart = abs(heart_fft/(N-1));
P1_heart = P2_heart(1:N/2+1);   %Select the first half at this time, because after fft is a symmetrical bilateral spectrum
P1_heart(2:end-1) = 2*P1_heart(2:end-1);

%Heartbeat Harmonic Detection
[heart_peaks,heart_peaksnum]=findpeaks_1(P1_heart,0.9,2,N,Tf);%0.9Hz-2Hz
[heart_harmonic_peaks,heart_harmonic_peaksnum]=findpeaks_1(P1_heart,1.8,4,N,Tf);%1.8Hz-4Hz
heart_peaks=heart_peaks/N/Tf;
heart_harmonic_peaks=heart_harmonic_peaks/N/Tf;
[heart_peaks_row,heart_peaks_column]=size(heart_peaks);
[heart_harmonic_peaks_row,heart_harmonic_peaks_column]=size(heart_harmonic_peaks);
for i=1:heart_peaks_column
    if max(P1_heart)-P1_heart(round(heart_peaks(i)*N*Tf)+1)<0.3 % When the difference between the two peaks is within a certain range, the signal that finds the harmonic is expanded
        for j=1:heart_harmonic_peaks_column
            if heart_harmonic_peaks(j)/heart_peaks(i)==2
                P1_heart(round(heart_peaks(i)*N*Tf)+1)=2*P1_heart(round(heart_peaks(i)*N*Tf)+1);
            end
        end
    end
end

%Heartbeat signal frequency domain diagram
figure("Name","Heartbeat signal frequency domain diagram");
plot(freq,P1_heart);
xlim([0,4]);
xlabel('Frequency(Hz)','FontWeight','bold');
ylabel('Amplitude','FontWeight','bold');
title('Spectrum diagram of heartbeat signal','FontWeight','bold');

%% EEMD Radar Processing 
% figure("Name","EEMD Data")
fs = index(1);
[imfEnst,imft] = eemd2(transpose(heart),fs,5,100,5,0);

%% Get All IMF
figure("Name","EEMD Signal")
len_imf = length(imft(:,1));
for c = 1:len_imf
    subplot(len_imf,1,c);
    if(c==2)
        plot(index,imft(c,:)) 
    end
    plot(index,imft(c,:)) 
    s_title = "IMF " + num2str(c);
    title(s_title)
end

%% Get Specific IMF
figure("Name","EEMD 2 Signal")
plot(index,imft(3,:)*3) 
% figure,

%% Get All imf Ensemble IMF
figure("Name","EEMD 3 Signal")
len_imfe = length(imfEnst(:,1));
for c = 1:len_imfe
    subplot(len_imfe,1,c);
    plot(index,imfEnst(c,:)) 
    s_title = "IMFe " + num2str(c);
    title(s_title)
    xlim([30 50])
end

%% Get Spesific omf Ensemble IMF
figure("Name","EEMD 4 Signal")
len_imfe = length(imfEnst(:,1));
for c = 1:len_imfe
    subplot(len_imfe,1,c);
    plot(index,imfEnst(c,:)) 
     s_title = "IMFe " + num2str(c);
    title(s_title)
end

%% Try Sum of Signal
% radar_heartbeat = imfEnst(2,:)+imfEnst(3,:);
radar_heartbeat = imfEnst(3,:);
figure("Name","EEMD 5 Signal")
plot(index,radar_heartbeat)
title('radar extracted heartbeat signal')
xlabel('Time(s)')
ylabel('Amplitude')

%% Try Sum of Signal with PPG
% radar_heartbeat = imfEnst(2,:)+imfEnst(3,:);
radar_heartbeat = imfEnst(3,:);
figure("Name","EEMD 5 Signal")
subplot(2,1,1)
plot(index,abs(radar_heartbeat))
xlim([10 30])
title('radar extracted heartbeat signal')
xlabel('Time(s)')
ylabel('Amplitude')
subplot(2,1,2)
plot(time,ch_2)
xlim([10 30])
title('Biopac Ch2');

heart = radar_heartbeat;

figure,
ax = gca;
yyaxis left
radar_heartbeat = imfEnst(2,:)+imfEnst(5,:);
plot(radar_heartbeat, 'Color', 'k','LineWidth',1)
ylabel('Amplitude','fontweight','bold','fontsize',12)
ylim([-0.03 6e-3])
ax.YColor = 'k';
yyaxis right
plot(time,ch_2, 'Color', 'r','LineWidth',1)
% plot(org_taxis+13, ecgSig, 'Color', 'r','LineWidth',1)
xlabel('Time (s)','fontweight','bold','fontsize',12)
ylabel('Amplitude(mV)','fontweight','bold','fontsize',12', 'Color', 'k')
ylim([-0.03 1.12])
ax.YColor = 'k';
ax.FontWeight = 'bold';
ax.LineWidth = 0.8;
% title('Radar(HR) & ECG')
legend('Radar Heartbeat Signal 3', 'PPG signal')
xlim([10 200])

% ------------------------------------- PPG DATA PROCESSING --------------------------------
% load('E:\XSDUX\Heartsound\Heartsound_Analysis\Dataset\ppg_2_txt.mat')
disp("+-----+");
disp("| PPG |");
disp("+-----+");
figure('Name','PPG Ch1')
subplot(4,1,1)
plot(time,ch_1)
title('Biopac Ch1');
subplot(4,1,2)
plot(time,ch_2)
xlim([0 20])
title('Biopac Ch2');
subplot(4,1,3)
plot(time,ch_3)
title('Biopac Ch3');
subplot(4,1,4)
plot(time,ch_4)
title('Biopac Ch3');

% ------------------------------------- PPG DATA PROCESSING --------------------------------
disp("+-----------+");
disp("| Radar PPG |");
disp("+-----------+");
figure('Name','Radar X PPG')
subplot(3,1,1)
[pks,locs] = findpeaks(abs(heart));
plot(index,abs(heart),(locs/10),pks,"o");
% plot();
xlim([10 20])
title('Radar Breath');

[pks_2,locs_2] = findpeaks(ch_2,"MinPeakHeight",0.1);
subplot(3,1,2)
plot(time,ch_2,(locs_2/1000),pks_2,"o")
xlim([10 20])
subplot(3,1,3)
title('PPG and Peaks');
plot(locs_2,pks_2)
xlim([10 20])
title('Biopac Ch2');

figure('Name','Radar to BPM')
subplot(4,1,1)
plot(time,ch_2)
xlim([10 20])
title('Biopac Ch2');
subplot(4,1,2)
[pks,locs] = findpeaks(abs(heart));
plot(index,abs(heart),(locs/10),pks,"o");
% plot();
xlim([10 20])
title('Radar Breath');

%%
disp("+-----------+");
disp("| Radar PPG |");
disp("+-----------+");
figure('Name','Radar and 2 PPG')
subplot(3,1,1)
[pks,locs] = findpeaks(abs(heart));
plot(index,abs(heart),(locs/10),pks,"o");
% plot();
xlim([13 20])
title('Radar Breath');

[pks_2,locs_2] = findpeaks(ch_2,"MinPeakHeight",0.1);
subplot(3,1,2)
plot(time,ch_2,(locs_2/1000),pks_2,"o")
xlim([13 20])
subplot(3,1,3)
plot(locs_2,pks_2)
xlim([13 20])
title('Biopac Ch2');

figure('Name','Radar to BPM')
subplot(4,1,1)
plot(time,ch_2)
xlim([10 20])
title('Biopac Ch2');
subplot(4,1,2)
[pks,locs] = findpeaks(abs(heart));
plot(index,abs(heart),(locs/10),pks,"o");
% plot();
xlim([10 20])
title('Radar Breath');

%% Counting BPM
time_peak = locs/10;
time_diff=[];
time_diff_x=[];
for i=2:length(time_peak)
    time_diff(i-1)=time_peak(i)-time_peak(i-1);
    time_diff_x(i-1)=time_peak(i);
end
time_peak_ppg = locs_2/1000;
time_diff_ppg=[];
time_diff_x_ppg=[];
for i=2:length(time_peak_ppg)
    time_diff_ppg(i-1)=time_peak_ppg(i)-time_peak_ppg(i-1);
    time_diff_x_ppg(i-1)=time_peak_ppg(i);
end
subplot(4,1,3)
plot(time_diff_x,time_diff)

subplot(4,1,4)
plot(time_diff_x_ppg,time_diff_ppg)

figure("Name","PPG x Radar")
subplot(2,1,1)
plot(time_diff_x,time_diff)
subplot(2,1,2)
plot(time_diff_x_ppg,time_diff_ppg)


% %% detect R peak and store R-R interval
% %%% perform double-threshold method to detect R peak
% % threshold 1: magnitude
% indmag=locs;    % find the index with value larger than 1.5
% 
% % threshold 2: define window
% diffind=indmag(2:end)-indmag(1:end-1);
% indgap=find(diffind>1);
% 
% indmax=[];   % the location of index with maximal value in each cycle
% for i=1:length(indgap)+1
%     if i==1
%         period=indmag(1:indgap(1));
%     elseif i==length(indgap)+1    
%         period=indmag(indgap(i-1)+1:end);
%     else
%         period=indmag(indgap(i-1)+1:indgap(i));
%     end
%     [value,ind]=max(heart(period));
%     indmax(i)=period(ind(1));
% end
% fs=1000;
% indmax(end)=[];
% PPinterval=[];  % heart rate variability (HRV)
% PPinterval=diff(indmax)/fs;
% HRV_radar=60./PPinterval;
% 
% avgPinterval_radar = mean(PPinterval);
% avgHR_radar = mean(HRV_radar);
% % plot the radar RR tachogram to show obvious artifacts:
% Ann_radar = cumsum(PPinterval);
% figure,
% subplot(3,1,1)
% plot(Ann_radar,PPinterval)
% title("Radar RR tachogram")
% % % Filter from artifacts and plot the average heart rate:
% % PPinterval = hrv.RRfilter(PPinterval,0.05);
% % subplot(4,1,2)
% % plot(Ann_radar,PPinterval)
% % title("Radar RR tachogram after filtering the artifacts")
% % Plot the average heart rate:
% subplot(3,1,2)
% plot(Ann_radar, hrv.HR(PPinterval,60))
% title("average heart rate from Radar")
% % Corresponding relative RR intervals:
% subplot(3,1,3)
% pp = hrv.rrx(PPinterval);
% plot(Ann_radar,pp)
% title("radar relative RR intervals")
% %% save RR intervals matrix
% save('rr_subject_ecg.mat','RRinterval')   % save A and B variables
% save('rr_subject_radar.mat','PPinterval')
% %% rewrite data set by HRV
% % samplerate=1/mean(RRinterval);
% % N=length(RRinterval);

PPinterval = time_diff;
RRinterval = time_diff_ppg;
avgRRinterval = mean(RRinterval);
avgPinterval_radar = mean(PPinterval);
HRV=60./RRinterval;
HRV_radar=60./PPinterval;
avgHR = mean(HRV)
avgHR_radar = mean(HRV_radar)

disp("+-----------+");
disp("| Interval RR dan INTERVAL PP|");
disp("+-----------+");
figure("Name","Interval RR dan INTERVAL PP")
plot(RRinterval),
hold on
plot(PPinterval),
% title('Before Handgrip Test Challenge')
xlabel('beats'),ylabel('R-R interval (s)')
legend({ ['Average RR interval of PPG = ',num2str(avgRRinterval),'sec.' ];
    [ 'Average Peak interval of Radar = ',num2str(avgPinterval_radar),'sec.' ]})

figure("Name","Histogram")
histogram(HRV,100)
hold on
histogram(HRV_radar,100)
xlim([0 90])
ylim([0 90])
title('After Handgrip Test Challenge')
legend({ ['Mean HRV of PPG is ' num2str(avgHR) ];
    ['Mean HRV of Radar is ' num2str(avgHR_radar) ]})
xlabel('R-R interval (s)','fontweight','bold','fontsize',12')
ylabel('count beats','fontweight','bold','fontsize',12')
%% RR interval spectral analysis
data = RRinterval;
lpad=2*zeros(size(data));
data = [data lpad];
samplerate=1/mean(RRinterval);
[nfft, faxis, PS_RR] = spectral(data, length(data), samplerate);
figure("Name","HRV Spectral Analysis (FFT)")
plot(faxis,20*log10(PS_RR(1:nfft/2+1)))
title('PPG HRV Spectral Analysis (FFT)')
xlabel('Frequency (Hz)')
ylabel('Power Spectrum (dB)')

% 2. Estimate heart rate
[radarBPM, ppgBPM] = estimateHeartRate(heart, ch_2, fs);

% 3. Display results
fprintf('Estimated Heart Rate from Radar: %.2f BPM\n', radarBPM);
fprintf('Estimated Heart Rate from PPG: %.2f BPM\n', ppgBPM);

% 4. Bland-Altman plot
if(length(RRinterval) > length(PPinterval))
    data_rr = RRinterval(1:length(PPinterval));
    data_pp = PPinterval(1:length(PPinterval));
else
    data_rr = RRinterval(1:length(RRinterval));
    data_pp = PPinterval(1:length(RRinterval));
end
blandAltmanPlot(data_pp, data_rr);

% Function to estimate heart rate
function [radarBPM, ppgBPM] = estimateHeartRate(radar_signal, ch_2_signal, fs)
    % Peak detection for radar signal
    [radarPeaks, ~] = findpeaks(radar_signal, 'MinPeakDistance', fs/3); % Assuming HR > 40 BPM
    radarBPM = (length(radarPeaks) - 1) / (length(radar_signal) / fs) * 60;

    % Peak detection for PPG signal
    [ppgPeaks, ~] = findpeaks(ch_2_signal, 'MinPeakDistance', fs/3); % Assuming HR > 40 BPM
    ppgBPM = (length(ppgPeaks) - 1) / (length(ch_2_signal) / fs) * 60;
end

% Function to generate Bland-Altman plot
function blandAltmanPlot(radarBPM, ppgBPM)
    % Calculate mean and difference for the two measurements
    meanBPM = (radarBPM + ppgBPM) / 2;
    diffBPM = radarBPM - ppgBPM;

    % Calculate mean difference and limits of agreement
    meanDiff = mean(diffBPM);
    stdDiff = std(diffBPM);
    upperLimit = meanDiff + 1.96*stdDiff;
    lowerLimit = meanDiff - 1.96*stdDiff;

    % Plotting
    figure("Name","Bland Altman");
    plot(meanBPM, diffBPM, 'b*'); hold on;
    plot(meanBPM, repmat(meanDiff, size(meanBPM)), 'r-');
    plot(meanBPM, repmat(upperLimit, size(meanBPM)), 'r--');
    plot(meanBPM, repmat(lowerLimit, size(meanBPM)), 'r--');
    title('Bland-Altman Plot');
    xlabel('Mean BPM [(Radar + PPG)/2]');
    ylabel('Difference BPM (Radar - PPG)');
    legend('Data', 'Mean difference', 'Upper 95% limit', 'Lower 95% limit');
    grid on;
end