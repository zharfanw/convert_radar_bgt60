clc
clear all
close all

addpath(".library\npy-matlab-master\npy-matlab\")
addpath("Heartsound_Analysis")

datanya = [
    glob('.dataset/Data 22-25/Measurement 23042024/*/RadarIfxAvian_00/radar.npy');
    glob('.dataset/Data 22-25/Measurement 25042024/*/RadarIfxAvian_00/radar.npy');
    glob('.dataset/Data 22-25/Measurement 25042024/*/*/RadarIfxAvian_00/radar.npy');
    glob('.dataset/Data 22-25/Measurement 22042024/*/*/RadarIfxAvian_00/radar.npy');
];

datanya_conf = [
    glob('.dataset/Data 22-25/Measurement 23042024/*/RadarIfxAvian_00/config.json');
    glob('.dataset/Data 22-25/Measurement 25042024/*/RadarIfxAvian_00/config.json');
    glob('.dataset/Data 22-25/Measurement 25042024/*/*/RadarIfxAvian_00/config.json');
    glob('.dataset/Data 22-25/Measurement 22042024/*/*/RadarIfxAvian_00/config.json');
];

master_save_path = ".dataset/radar/"

% datanya
% datanya_conf

size_datanya = size(datanya);
num_data = size_datanya(1);
for index=1:num_data
    dum_name = convertCharsToStrings(datanya{index});
    dum_name_arr = dum_name.split("\");
    dum_name_real = dum_name_arr(3) + "_" +dum_name_arr(4) + "_" +dum_name_arr(5);
    dum_name_real = dum_name_real.replace(" ","_");
    master_title = dum_name_real;

    path_radar_datanya = convertCharsToStrings(datanya{index});
    path_radar_conf = convertCharsToStrings(datanya_conf{index});
    
    fprintf("%i Namanya        : %s \n",index,dum_name_real);
    fprintf("%i Datanya        : %s \n",index,path_radar_datanya);
    fprintf("%i Datanya conf   : %s \n\n",index,path_radar_conf);
    % datasetnya = readNPY(".dataset\Data 22-25\Measurement 25042024\BGT60TR13C_record_2504202420240425-101640\RadarIfxAvian_00\radar.npy");
    % conf_fname =         '.dataset\Data 22-25\Measurement 25042024\BGT60TR13C_record_2504202420240425-101640\RadarIfxAvian_00\config.json'; 
%% Setup RADAR
    close all

    anal_name = dum_name_real;
    datasetnya = readNPY(path_radar_datanya);
    conf_fname = path_radar_conf;

    conf_fid = fopen(conf_fname); 
    conf_raw = fread(conf_fid,inf); 
    conf_str = char(conf_raw'); 
    fclose(conf_fid); 
    conf_confignya = jsondecode(conf_str);
    
    fprintf("Chirps Reptition Time(s): %d", conf_confignya.device_config.fmcw_single_shape.chirp_repetition_time_s)
    fprintf("num_chirps_per_frame: %i", conf_confignya.device_config.fmcw_single_shape.num_chirps_per_frame)
    fprintf("num_samples_per_chirp: %i", conf_confignya.device_config.fmcw_single_shape.num_samples_per_chirp)
    fprintf("sample_rate_Hz: %i", conf_confignya.device_config.fmcw_single_shape.sample_rate_Hz)
    fprintf("Frame Rate(Hz): %i", 1 / conf_confignya.device_config.fmcw_single_shape.frame_repetition_time_s)

    Radar_Parameter.Num_Tx_Antennas = conf_confignya.device_config.fmcw_single_shape.tx_antennas;
    Radar_Parameter.Num_Rx_Antennas= length(conf_confignya.device_config.fmcw_single_shape.rx_antennas);
    Radar_Parameter.Mask_Tx_Antennas = 1;
    Radar_Parameter.Mask_Rx_Antennas = 7;
    Radar_Parameter.Are_Rx_Antennas_Interleaved = 1;
    Radar_Parameter.Modulation_Type_Enum = 1;
    Radar_Parameter.Chirp_Shape_Enum= 0;
    Radar_Parameter.Lower_RF_Frequency_kHz = conf_confignya.device_config.fmcw_single_shape.start_frequency_Hz;
    Radar_Parameter.Upper_RF_Frequency_kHz = conf_confignya.device_config.fmcw_single_shape.end_frequency_Hz;
    Radar_Parameter.Sampling_Frequency_kHz = conf_confignya.device_config.fmcw_single_shape.sample_rate_Hz/1000;
    Radar_Parameter.ADC_Resolution_Bits=12;
    Radar_Parameter.Are_ADC_Samples_Normalized =1;
    Radar_Parameter.Data_Format_Enum=0;
    Radar_Parameter.Chirps_per_Frame=conf_confignya.device_config.fmcw_single_shape.num_chirps_per_frame;
    Radar_Parameter.Samples_per_Chirp= conf_confignya.device_config.fmcw_single_shape.num_samples_per_chirp;
    Radar_Parameter.Samples_per_Frame=Radar_Parameter.Chirps_per_Frame*Radar_Parameter.Samples_per_Chirp*Radar_Parameter.Num_Rx_Antennas;
    Radar_Parameter.Chirp_Time_sec=conf_confignya.device_config.fmcw_single_shape.frame_repetition_time_s;
    Radar_Parameter.Pulse_Repetition_Time_sec=conf_confignya.device_config.fmcw_single_shape.chirp_repetition_time_s;
    Radar_Parameter.Frame_Period_sec=conf_confignya.device_config.fmcw_single_shape.frame_repetition_time_s;
    
    dummy_size = size(datasetnya);
    Frame_Number = dummy_size(1);
    NumRXAntenna = Radar_Parameter.Num_Rx_Antennas;
    Frame = datasetnya;

%% PROCSESS FMCW    
    used_range_low = 50;
    used_range_up = 60;
    

    % Set device configuration for presence sensing
    f_start = Radar_Parameter.Lower_RF_Frequency_kHz * 1e7; %FMCW chirp start frequency in Hz
    f_end   = Radar_Parameter.Upper_RF_Frequency_kHz * 1e7; %FMCW chirp end frequency in Hz
    
    samples_per_chirp  = Radar_Parameter.Samples_per_Chirp;
    chirps_per_frame   = Radar_Parameter.Chirps_per_Frame;
    
    zero_padding_factor = 2;
    
    % mat = squeeze(Frame(:,:,2,1));
    mat = squeeze(Frame(1,2,:,:));
    % mat = transpose(mat);
    
    
    zero_padded = zeros(size(mat,1),size(mat,2)*2^(zero_padding_factor));
    result = zeros(1,200);
    
    hFig = figure('Name','Range Doppler','numbertitle', 'off');
    % In this example, show frames till figure is closed.
    frame_number = 0;
    
    size_frame = size(Frame);
    
    % mti_history = zeros((self.num_chirps_per_frame, num_samples, num_ant))
    mti_history = zeros(Radar_Parameter.Chirps_per_Frame, Radar_Parameter.Samples_per_Chirp, Radar_Parameter.Num_Rx_Antennas);
    mti_alpha = 0.8;
    
    n_distance = 2;
    n_velocity = 50;
    data_vs = [];
    
    countdown = 0;
    old_status = 0;
    % processed_fmcw = zeros(Radar_Parameter.Chirps_per_Frame,Radar_Parameter.Samples_per_Chirp,size_frame(4));
    % processed_fmcw = zeros(Radar_Parameter.Chirps_per_Frame,Radar_Parameter.Samples_per_Chirp,size_frame(1));
    processed_fmcw = zeros(Radar_Parameter.Chirps_per_Frame,Radar_Parameter.Samples_per_Chirp,size_frame(1),Radar_Parameter.Num_Rx_Antennas);
    filtered_processed_fmcw = zeros(Radar_Parameter.Chirps_per_Frame,Radar_Parameter.Samples_per_Chirp,size_frame(1),Radar_Parameter.Num_Rx_Antennas);
    data_vs = zeros(size_frame(1),Radar_Parameter.Num_Rx_Antennas);
    data_vs_filt = zeros(size_frame(1),Radar_Parameter.Num_Rx_Antennas);
    %##########################################################################
    % STEP 6: Continue fetching Radar data and run desired algorithm on it.
    %##########################################################################
    % while ishandle(hFig)
    figure("Name","Frame Surface")
    % for loop=1:size_frame(4)
    for cur_antenna=1:Radar_Parameter.Num_Rx_Antennas
        for loop=1:size_frame(1)
            range_window = window(@blackmanharris,samples_per_chirp);
            doppler_window = window(@blackmanharris,Radar_Parameter.Chirps_per_Frame);
            % data = Frame(:,:,1,loop);
            data = double(squeeze(Frame(loop,1,:,:)));
            % data = transpose(data);
            data = data - mean(data);
    
            % data_mti = data - mti_history;
            mti_history(:, :, 1) = data * mti_alpha + mti_history(:, :, 1) * (1 - mti_alpha);
    
                % Mean removal and windowing of raw data
            avg_windowed = (data - mean(data,2)).*range_window';
    
            zero_padded(:,1:samples_per_chirp) = avg_windowed;
    
            % Range FFT
            range_fft = (fft(zero_padded'))/samples_per_chirp; % normalized FFT
            range_fft = 2*range_fft.'; % compensating for the negative half the spectrum
            % remove the redundant negative spectrum since FFT input was real
            fft1d = range_fft(:,1:(size(range_fft,2)/2));
    
            % fft1d = transpose(fft1d);
    
            fft1d = fft1d .* doppler_window;
    
            zero_padded_2(:,1:Radar_Parameter.Chirps_per_Frame) = transpose(fft1d);
    
            % Range FFT
            range_fft_2 = (fft(zero_padded_2'))/Radar_Parameter.Chirps_per_Frame; % normalized FFT
            range_fft = 2*range_fft_2.'; % compensating for the negative half the spectrum
            % remove the redundant negative spectrum since FFT input was real
            fft2d = range_fft_2(:,1:(size(range_fft_2,2)/2));
            sum_range = fft2d(1,:) + fft2d(2,:) + fft2d(3,:) + fft2d(4,:);
            sum_all = sum(sum(sum_range));
            processed_fmcw(:,:,loop,cur_antenna) = fft2d;
            % data_vs = [data_vs max(max(fft2d))];
            % data_vs = [data_vs max(max(sum_all))];
            
            % Filtered
            [M,I] = max(abs(processed_fmcw(:,:,loop,1)),[],"all","linear");
            [dim1, dim2] = ind2sub(size(processed_fmcw(:,:,loop,1)),I);
        
            [min_M,min_I] = min(abs(processed_fmcw(:,:,loop)),[],"all","linear");
            [min_dim1, min_dim2] = ind2sub(size(processed_fmcw(:,:,loop,1)),min_I);
            size_frame_proc=size(processed_fmcw);
        
            data_treshold = M*0.9;
            new_frame=zeros(size_frame_proc(1),size_frame_proc(2));
            % half_dist_size = round(size_frame(2)/2);
            half_dist_size = size_frame_proc(2);
            for d1_frame=1:size_frame_proc(1)
                for d2_frame=1:half_dist_size
                    if(abs(processed_fmcw(d1_frame,d2_frame,loop,1))>data_treshold)
                        % if(
                        new_frame(d1_frame,d2_frame) = abs(processed_fmcw(d1_frame,d2_frame,loop,1));
                    end
                end
            end
            % dum_dat = new_frame(120,:)+new_frame(121,:)+new_frame(122,:)+new_frame(123,:)+new_frame(124,:)+new_frame(125,:)+new_frame(126,:)+new_frame(127,:)+new_frame(128,:);
            % data_raw(loop) = max(new_frame)
            % new_frame = new_frame - mean(new_frame);
            % data_raw = max(dum_dat);
            filt_dum_dat = zeros(size_frame_proc(1),1);
            for sel_dat=5:20
                % dum_dat = dum_dat+new_frame(:,sel_dat);
                filt_dum_dat = filt_dum_dat+abs(processed_fmcw(:,sel_dat,loop,cur_antenna));
            end
            data_raw(loop) = sum(filt_dum_dat);
    
            dum_dat = zeros(size_frame_proc(2),1);
            for sel_dat=used_range_low:used_range_up
                % dum_dat = dum_dat+new_frame(:,sel_dat);
                dum_dat = dum_dat + processed_fmcw(:,sel_dat,loop,cur_antenna);
            end
            % data_raw(loop) = sum(dum_dat);
            data_vs(loop,cur_antenna) = sum(dum_dat);
            % data_vs = [data_vs sum(dum_dat)];
    
    
        end
    end
%% Visualize
    size_frame = size(processed_fmcw);

    data_raw = zeros(size_frame(3),1);
    
    
    figure("Name","Animation ")
    for loop=1:size_frame(3)
        [M,I] = max(abs(processed_fmcw(:,:,loop,1)),[],"all","linear");
        [dim1, dim2] = ind2sub(size(processed_fmcw(:,:,loop,1)),I);
    
        [min_M,min_I] = min(abs(processed_fmcw(:,:,loop)),[],"all","linear");
        [min_dim1, min_dim2] = ind2sub(size(processed_fmcw(:,:,loop,1)),min_I);
    
        data_treshold = M*0.9;
        new_frame=zeros(size_frame(1),size_frame(2));
        % half_dist_size = round(size_frame(2)/2);
        half_dist_size = size_frame(2);
        for d1_frame=1:size_frame(1)
            for d2_frame=1:half_dist_size
                if(abs(processed_fmcw(d1_frame,d2_frame,loop,1))>data_treshold)
                    % if(
                    new_frame(d1_frame,d2_frame) = abs(processed_fmcw(d1_frame,d2_frame,loop,1));
                end
            end
        end
        % dum_dat = new_frame(120,:)+new_frame(121,:)+new_frame(122,:)+new_frame(123,:)+new_frame(124,:)+new_frame(125,:)+new_frame(126,:)+new_frame(127,:)+new_frame(128,:);
        % data_raw(loop) = max(new_frame)
        % new_frame = new_frame - mean(new_frame);
        % data_raw = max(dum_dat);
        dum_dat = zeros(size_frame(1),1);
        for sel_dat=5:20
            % dum_dat = dum_dat+new_frame(:,sel_dat);
            dum_dat = dum_dat+abs(processed_fmcw(:,sel_dat,loop,1));
        end
        data_raw(loop) = sum(dum_dat);
    
    
    
        % subplot(1,3,1)
        % surface(abs(processed_fmcw(:,:,loop,1)),"EdgeColor","none")
        % xlabel("Range")
        % ylabel("Velocity")
        % title("Antenna 1")
        % view(2)
        % subplot(1,3,2)
        % surface(abs(processed_fmcw(:,:,loop,2)),"EdgeColor","none")
        % xlabel("Range")
        % ylabel("Velocity")
        % title("Antenna 2")
        % view(2)
        % subplot(1,3,3)
        % surface(abs(processed_fmcw(:,:,loop,3)),"EdgeColor","none")
        % xlabel("Range")
        % ylabel("Velocity")
        % title("Antenna 3")
        % view(2)
    
        % subplot(1,3,2)
        % surface(new_frame,"EdgeColor","none")
        % view(2)
        % time
        % drawnow
        % pause(0.0001)
    end
    figure("Name","Doppler FFT "+master_title)
    subplot(1,3,1)
    surface(abs(processed_fmcw(:,:,round(size_frame(3)/2),1)),"EdgeColor","none")
    xlabel("Range")
    ylabel("Velocity")
    title("Antenna 1")
    view(2)
    subplot(1,3,2)
    surface(abs(processed_fmcw(:,:,round(size_frame(3)/2),2)),"EdgeColor","none")
    xlabel("Range")
    ylabel("Velocity")
    title("Antenna 2")
    view(2)
    subplot(1,3,3)
    surface(abs(processed_fmcw(:,:,round(size_frame(3)/2),3)),"EdgeColor","none")
    xlabel("Range")
    ylabel("Velocity")
    title("Antenna 3")
    view(2)

    savefig(master_save_path+"Figure_dopplerfft_"+master_title);
%% Vital Sign Analysis
    % % figure("Name","Vital Sign")
    % t_vs = 1:length(data_vs)
    % t_vs = t_vs * Radar_Parameter.Frame_Period_sec;
    % 
    % %Get the real and imaginary parts of the signal
    % data = zeros(length(data_vs));
    % data_real=zeros(length(data_vs));
    % data_imag=zeros(length(data_vs));
    % for k=1:length(data_vs)
    %     data(k) = data_vs(k);
    %     data_real(k)=real(data(k));
    %     data_imag(k)=imag(data(k));
    % end
    % 
    % %Calculate signal phase
    % for k=1:length(data_vs)
    %     signal_phase(k)=atan(data_imag(k)/data_real(k));
    % end
    % 
    % %Phase unwrapping
    % for k=2:length(data_vs)
    %     diff=signal_phase(:,k)-signal_phase(:,k-1);
    %     if diff>pi/2
    %         signal_phase(:,(k:end))=signal_phase(:,(k:end))-pi;
    %     elseif diff<-pi/2
    %         signal_phase(:,(k:end))=signal_phase(:,(k:end))+pi;
    %     end
    % end
    % 
    % %Calculate the phase difference
    % for k=1:length(data_vs)-1
    %     % delta_phase(k)=signal_phase(k+1)-signal_phase(k);
    %     phaseUsedComputation(k)=signal_phase(k+1)-signal_phase(k);
    % end
    % 
    % %Original signal bandpass filtering
    % % thresh=0.000000000001;
    % % for k=1:length(data_vs)-3
    % %     phaseUsedComputation(k)=filter_RemoveImpulseNoise(delta_phase(k),delta_phase(k+1),delta_phase(k+2),thresh);
    % % end
    % % index=1:1:length(data_vs)-1;
    % % index=1:1:length(phaseUsedComputation);
    % % index=index*Tf;
    % 
    % N=1024; 
    % Tf=Radar_Parameter.Frame_Period_sec;       %frame period
    % %Original signal bandpass filtering
    % filter_delta_phase=filter(bpf_vitalsign,phaseUsedComputation);
    % 
    % 
    % vital_sign=filter_delta_phase;
    % index = t_vs;
    % index = index(1:length(index)-1);
    % %Original Signal Time Domain Plot
    % % figure("Name","Vital Sign");
    % % plot(index,vital_sign);
    % % xlabel('Time(s)','FontWeight','bold');
    % % ylabel('Amplitude','FontWeight','bold');
    % % title('cardiopulmonary signal','FontWeight','bold');
    % 
    % %do fft on the original signal
    % vital_sign_fft=fft(vital_sign,N);
    % 
    % %Convert double sideband signal to single sideband
    % freq=(0:1:N/2)/Tf/N;  %The vital sign signal sampling rate is the number of frames
    % P2 = abs(vital_sign_fft/(N-1));
    % P1 = P2(1:N/2+1);   %Select the first half at this time, because after fft is a symmetrical bilateral spectrum
    % P1(2:end-1) = 2*P1(2:end-1);
    % 
    % %Original signal frequency domain plot
    % figure("Name","Original Signal Frequency");
    % plot(freq,P1);
    % xlim([0,2]);
    % xlabel('Frequency(Hz)','FontWeight','bold');
    % ylabel('Amplitude','FontWeight','bold');
    % title('Cardiopulmonary signal spectrogram','FontWeight','bold');
    % 
    % %% Show Pre VITAL SIGN Process
    % figure("Name","Non VITAL BAND");
    % plot(t_vs,abs(data_vs));
    % xlabel('Time(s)','FontWeight','bold');
    % ylabel('Amplitude','FontWeight','bold');
    % title('Radar signal','FontWeight','bold');
    % 
    % % %Breathing signal bandpass filtering
    % % filter_delta_phase_breathe=filter(bpf_breathe,data_vs);
    % filter_delta_phase_breathe=filter(bpf_breathe,phaseUsedComputation);
    % 
    % breathe=filter_delta_phase_breathe;
    % 
    % %Time Domain Diagram of Respiration Signal
    % figure("Name","Breathing Data");
    % plot(index,breathe);
    % xlabel('Time(s)','FontWeight','bold');
    % ylabel('Amplitude','FontWeight','bold');
    % title('breathing signal','FontWeight','bold');
    % 
    % %Do fft on respiration signal
    % breathe_fft=fft(breathe,N);
    % 
    % %Convert double sideband signal to single sideband
    % P2_breathe = abs(breathe_fft/(N-1));
    % P1_breathe = P2_breathe(1:N/2+1);   %Select the first half at this time, because after fft is a symmetrical bilateral spectrum
    % P1_breathe(2:end-1) = 2*P1_breathe(2:end-1);
    % 
    % %Respiratory signal frequency domain diagram
    % figure("Name","Respiratory signal spectrogram");
    % plot(freq,P1_breathe);
    % xlim([0,2]);
    % xlabel('Frequency(Hz)','FontWeight','bold');
    % ylabel('Amplitude','FontWeight','bold');
    % title('Respiratory signal spectrogram','FontWeight','bold');
    % 
    % 
    % %Heartbeat signal bandpass filtering
    % filter_delta_phase_heart=filter(bpf_heart,phaseUsedComputation);
    % 
    % heart=filter_delta_phase_heart;
    % 
    % %Heartbeat signal time domain diagram
    % % figure(5);
    % figure("Name","Selected Heart Signal x Time Domain");
    % plot(index,heart);
    % % xlim([0 20])
    % xlabel('Time(s)','FontWeight','bold');
    % ylabel('Amplitude','FontWeight','bold ');
    % title('heartbeat signal','FontWeight','bold');
    % 
    % % figure(6);
    % figure("Name","Selected Heart Signal ABS x Time Domain");
    % plot(index,abs(heart));
    % % xlim([0 20])
    % xlabel('Time(s)','FontWeight','bold');
    % ylabel('Amplitude','FontWeight','bold ');
    % title('heartbeat signal','FontWeight','bold');
    % 
    % % Do fft on the heartbeat signal
    % heart_fft=fft(heart,N);
    % 
    % figure("Name","Radar Heart and Breath");
    % subplot(2,1,1);
    % plot(index,heart);
    % subplot(2,1,2);
    % plot(index,breathe);
    % savefig(master_save_path+"Figure_RadarHeartBreath_"+master_title);

end