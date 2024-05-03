clear all;
close all;


clear;

% [retVal] = readDCA1000_1('./04060010_1573837157_Raw_0.bin');

[FileName,PathName] = uigetfile('*.txt','Select the .txt file from radar');

[Radar_Parameter,Frame_Number,NumRXAntenna,Frame]=data_import([PathName FileName]);



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
    avgDC=nanmean(raw_data_matrix,2);               
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
% RX1data = reshape(retVal(1,:),numADCSamples,numChirps);   %RX1Êý¾Ý
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
% N=1024;        %FFTµãÊý

c=3.0e8;  
slope=(Radar_Parameter.Upper_RF_Frequency_kHz-Radar_Parameter.Lower_RF_Frequency_kHz)*1e3;   %FM slope
Tc=Radar_Parameter.Frame_Period_sec/(Radar_Parameter.Samples_per_Frame/3);      %chirpÖÜÆÚ
B=BW;    %FM bandwidth
Fs=Radar_Parameter.Sampling_Frequency_kHz*1e3;       %Sampling Rate
% f0=f0;    %initial frequency
lambda=c/f0;   %radar signal wavelength
d=lambda/2;    %Antenna Array Spacing
frame=Frame_Number/2;     %number of frames
Tf=Radar_Parameter.Frame_Period_sec;       %frame period
N=1024;        %FFTµãÊý


% plus Hamming window
range_win = hamming(numADCSamples);                     %Generate Hamming window
for k=1:1:frame
    din_win(:,k)=RX1data(:,2*k-1).*range_win;           %do to the signal Range-fft
    datafft(:,k)=fft(din_win(:,k));
end

%Find Range-bin peaks
rangeBinStartIndex=3;       %resolution0.1m
rangeBinEndIndex=20;
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

%phase unwrapping
for k=2:frame
    diff=signal_phase(:,k)-signal_phase(:,k-1);
    if diff>pi/2
        signal_phase(:,(k:end))=signal_phase(:,(k:end))-pi;
    elseif diff<-pi/2
        signal_phase(:,(k:end))=signal_phase(:,(k:end))+pi;
    end
end

%Calculate phase difference
for k=1:frame-1
    delta_phase(:,k)=signal_phase(:,k+1)-signal_phase(:,k);
end

%Remove Impulse Noise from a Waveform
thresh=0.8;
for k=1:frame-3
    phaseUsedComputation(:,k)=filter_RemoveImpulseNoise(delta_phase(:,k),delta_phase(:,k+1),delta_phase(:,k+2),thresh);
end

%Real chest displacement signal
for i=1:frame-3
   vital_sign(i)= (phaseUsedComputation(1)+phaseUsedComputation(i))*i/2;
end


vital_sign=phaseUsedComputation;
vital_sign=filter(bpf_vitalsign,vital_sign);

%VMD parameter settings
alpha = 5000;      % moderate bandwidth constraint
tau = 0;            % noise-tolerance (no strict fidelity enforcement)0
K = 3;              % modes
DC = 0;             % no DC part imposed
init = 1;           % initialize omegas uniformly
tol = 1e-6;

%VMD
[u, u_hat, omega] = VMD(vital_sign, alpha, tau, K, DC, init, tol);



%% Respiration signal
%Respiratory signal time domain diagram
respiration=filter(bpf_breathe,u(1,:));
index=1:1:frame-3;
index=index*Tf;
figure(3);
plot(index,respiration,'LineWidth',1,'color','b');
xlabel('Time(s)','FontWeight','bold');
ylabel('Amplitude','FontWeight','bold');
title('breathing signal','FontWeight','bold');
%saveas(gcf,'_BreathingSignal.png')


%Respiratory signal frequency domain diagram
freq=(0:1:N/2)/Tf/N;
figure(4);
imf_breathe_fft=fft(respiration,N);
P2_breathe_imf = abs(imf_breathe_fft/(N-1));
P1_breathe_imf = P2_breathe_imf(1:N/2+1);               %Select the first half at this time, because after fft is a symmetrical bilateral spectrum
P1_breathe_imf(2:end-1) = 2*P1_breathe_imf(2:end-1);
plot(freq,P1_breathe_imf,'Linewidth',1,'color','b');
xlim([0,2]);
xlabel('Frequency(Hz)','FontWeight','bold');
ylabel('Amplitude','FontWeight','bold');
title('Respiratory signal spectrogram','FontWeight','bold');
%saveas(gcf,'_Respiratory signal spectrogram.png')


%% breath area ratio
[inhale_to_exhale, exhale_to_inhale, point30Percent, point70Percent] = exhale_inhale_area_ratio(respiration,Tf);

%% Plot VMD Decomposition
figure(99)
subplot(length(u(:,1))+2,1,1);
indexnya = 1:(length(radar_raw_real)/2);
indexnya = indexnya * (Tc/numADCSamples);
plot(indexnya, radar_raw_real(1:length(radar_raw_real)/2))
ylabel('Raw','fontweight','bold','fontsize',18)
subplot(length(u(:,1))+2,1,2);
plot(index,vital_sign,'LineWidth',1,'color','b')
ylabel('Vital Sign','fontweight','bold','fontsize',18)
for ii=1:length(u(:,1))
    subplot(length(u(:,1))+2,1,ii+2);
    title_imf = append('IMF-',int2str(ii));    
    plot(index,u(ii,:));
    if(ii==length(u(:,1)))
        xlabel('Time(s)')
        xlabel('Time (s)','fontweight','bold','fontsize',20)
    end
    ylabel(title_imf,'fontweight','bold','fontsize',18)
   
end

%Connect the 30%-70% points with a line
figure(3);
hold on;
plot(point30Percent(1, :), point30Percent(2, :), 'r-x');
plot(point70Percent(1, :), point70Percent(2, :), 'g-o');
hold off;

%% heartbeat signal
%Heartbeat signal time domain diagram
heart=filter(bpf_heart,u(2,:));
figure(5);
plot(index,heart,'LineWidth',1,'color','b');
xlabel('Time(s)','FontWeight','bold');
ylabel('Amplitude','FontWeight','bold ');
title('heartbeat signal','FontWeight','bold');
saveas(gcf,'_heartbeat signal.png')



%Heartbeat signal frequency domain diagram
figure(6);
imf_heart_fft=fft(heart,N);
P2_heart_imf = abs(imf_heart_fft/(N-1));
P1_heart_imf = P2_heart_imf(1:N/2+1);   %Select the first half at this time, because after fft is a symmetrical bilateral spectrum
P1_heart_imf(2:end-1) = 2*P1_heart_imf(2:end-1);
plot(freq,P1_heart_imf,'LineWidth',1,'color','b');
xlim([0,4]);
xlabel('Frequency(Hz)','FontWeight','bold');
ylabel('Amplitude','FontWeight','bold');
title('Spectrum diagram of heartbeat signal','FontWeight','bold');
saveas(gcf,'_Spectrum diagram of heartbeat signal.png')

%% HRV
[IBI, MEAN, SDNN, r_MSSD] = mmHRV(heart, Tf);
figure(7);
plot(IBI(1,:),'LineWidth',1,'color','b');
axis auto;
xlabel('IBI Index','FontWeight','bold');
ylabel('IBI(s)','FontWeight','bold');
title('IBI estimation of mmHRV','FontWeight','bold');
saveas(gcf,'_IBI estimation of mmHRV.png')

figure(8)
plot(r_MSSD,size(IBI),'LineWidth',1,'color','b');
axis auto;
xlabel('r_MMSD Index','FontWeight','bold');
ylabel('r_MSSD (s)','FontWeight','bold');
title('r_MMSD estimation of mmHRV','FontWeight','bold');
saveas(gcf,'_r_MMSD estimation of mmHRV.png')

figure(9)
plot(mean(IBI, 2),'LineWidth',1,'color','b');
axis auto;
xlabel('MEAN Index','FontWeight','bold');
ylabel('MEAN (s)','FontWeight','bold');
title('MEAN estimation of mmHRV','FontWeight','bold');
saveas(gcf,'_MEAN estimation of mmHRV.png')

%% original signal
%Vital sign signal reconstruction
vital_sign=filter(bpf_heart,u(2,:))+filter(bpf_breathe,u(1,:));

%%Time Domain Diagram of Vital Signs Signal
figure(1);
plot(index,vital_sign,'LineWidth',1,'color','b');
xlabel('Time(s)','FontWeight','bold');
ylabel('Amplitude','FontWeight','bold');
title('cardiopulmonary signal','FontWeight','bold');
% saveas(gcf,'_cardiopulmonary signal.png')

%%Vital sign signal frequency domain diagram
figure(2);
imf_vital_sign_fft=fft(vital_sign,N);
P2_vital_sign_imf = abs(imf_vital_sign_fft/(N-1));
P1_vital_sign_imf = P2_vital_sign_imf(1:N/2+1);   %Select the first half at this time, because after fft is a symmetrical bilateral spectrum
P1_vital_sign_imf(2:end-1) = 2*P1_vital_sign_imf(2:end-1);
plot(freq,P1_vital_sign_imf,'LineWidth',1,'color','b')
xlim([0,2]);
xlabel('Frequency(Hz)','FontWeight','bold');
ylabel('Amplitude','FontWeight','bold');
title('Cardiopulmonary signal spectrogram','FontWeight','bold');
% saveas(gcf,'_Cardiopulmonary signal spectrogram.png')