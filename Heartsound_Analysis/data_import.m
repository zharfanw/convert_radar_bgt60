function [Radar_Parameter,Frame_Number,NumRXAntenna,Frame]=data_import(FileName)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Infineon XENSIV 60GHz Radar Matlab Interface   %
%                                                 %
%  Li Zhang  & Prateek                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% IFX_radar
%% load data_parameters
    opts = delimitedTextImportOptions("NumVariables", 1);
    opts.DataLines = [9, 29];
    opts.Delimiter = ",";
    opts.VariableNames = "parameters";
    opts.VariableTypes = "string";
    opts.ExtraColumnsRule = "ignore";
    opts.EmptyLineRule = "read";
    %[FileName,PathName] = uigetfile('*.txt','Select the .txt file from radar');
    %IFX_radar_parameters = readtable([PathName,FileName], opts);
    IFX_radar_parameters = readtable([FileName], opts);
    IFX_radar_parameters_index = isstrprop(IFX_radar_parameters.parameters,'digit'); 
    clear opts

    Num_Tx_Antennas = str2num(IFX_radar_parameters.parameters{1}(IFX_radar_parameters_index{1,1}));
    Num_Rx_Antennas = str2num(IFX_radar_parameters.parameters{2}(IFX_radar_parameters_index{2,1}));
    Mask_Tx_Antennas = str2num(IFX_radar_parameters.parameters{3}(IFX_radar_parameters_index{3,1}));
    Mask_Rx_Antennas = str2num(IFX_radar_parameters.parameters{4}(IFX_radar_parameters_index{4,1}));
    Are_Rx_Antennas_Interleaved = str2num(IFX_radar_parameters.parameters{5}(IFX_radar_parameters_index{5,1}));
    Modulation_Type_Enum = str2num(IFX_radar_parameters.parameters{6}(IFX_radar_parameters_index{6,1}));                                  % Modulation_Type_Enum_Def = {DOPPLER = 0, FMCW = 1}
    Chirp_Shape_Enum = str2num(IFX_radar_parameters.parameters{7}(IFX_radar_parameters_index{7,1}));                                      %{UP_CHIRP = 0, DOWN_CHIRP = 1, UP_DOWN_CHIRP = 2, DOWN_UP_CHIRP = 3}
    Lower_RF_Frequency_kHz = str2num(IFX_radar_parameters.parameters{8}(IFX_radar_parameters_index{8,1}))/10^3;
    Upper_RF_Frequency_kHz = str2num(IFX_radar_parameters.parameters{9}(IFX_radar_parameters_index{9,1}))/10^3;
    Sampling_Frequency_kHz = str2num(IFX_radar_parameters.parameters{10}(IFX_radar_parameters_index{10,1}));
    ADC_Resolution_Bits = str2num(IFX_radar_parameters.parameters{11}(IFX_radar_parameters_index{11,1}));
    Are_ADC_Samples_Normalized = str2num(IFX_radar_parameters.parameters{12}(IFX_radar_parameters_index{12,1}));
    Data_Format_Enum = str2num(IFX_radar_parameters.parameters{13}(IFX_radar_parameters_index{13,1}));                                    %Data_Format_Enum_Def = {DATA_REAL = 0, DATA_COMPLEX = 1, DATA_COMPLEX_INTERLEAVED = 2}
    Chirps_per_Frame = str2num(IFX_radar_parameters.parameters{14}(IFX_radar_parameters_index{14,1}));
    Samples_per_Chirp = str2num(IFX_radar_parameters.parameters{15}(IFX_radar_parameters_index{15,1}));
    Samples_per_Frame = str2num(IFX_radar_parameters.parameters{16}(IFX_radar_parameters_index{16,1}));
    Chirp_Time_sec = str2num(IFX_radar_parameters.parameters{17}(IFX_radar_parameters_index{17,1}))/10^(length(IFX_radar_parameters.parameters{17}(IFX_radar_parameters_index{17,1}))-1);
    Pulse_Repetition_Time_sec = str2num(IFX_radar_parameters.parameters{18}(IFX_radar_parameters_index{18,1}))/10^(length(IFX_radar_parameters.parameters{18}(IFX_radar_parameters_index{18,1}))-1);
    Frame_Period_sec = str2num(IFX_radar_parameters.parameters{19}(IFX_radar_parameters_index{19,1}))/10^(length(IFX_radar_parameters.parameters{19}(IFX_radar_parameters_index{19,1}))-1);
    %Frame_Number = str2num(IFX_radar_parameters.parameters{21}(IFX_radar_parameters_index{21,1}));
    Radar_Parameter= struct(  "Num_Tx_Antennas", Num_Tx_Antennas,     ...
                              "Num_Rx_Antennas", Num_Rx_Antennas,     ...
                              "Mask_Tx_Antennas", Mask_Tx_Antennas,     ...
                              "Mask_Rx_Antennas", Mask_Rx_Antennas,     ...
                              "Are_Rx_Antennas_Interleaved",Are_Rx_Antennas_Interleaved,     ...
                              "Modulation_Type_Enum", Modulation_Type_Enum,     ...
                              "Chirp_Shape_Enum", Chirp_Shape_Enum,     ...
                              "Lower_RF_Frequency_kHz", Lower_RF_Frequency_kHz, ...
                              "Upper_RF_Frequency_kHz", Upper_RF_Frequency_kHz, ...
                              "Sampling_Frequency_kHz", Sampling_Frequency_kHz, ...
                              "ADC_Resolution_Bits", ADC_Resolution_Bits, ...
                              "Are_ADC_Samples_Normalized", Are_ADC_Samples_Normalized, ...
                              "Data_Format_Enum", Data_Format_Enum, ...
                              "Chirps_per_Frame", Chirps_per_Frame, ...
                              "Samples_per_Chirp", Samples_per_Chirp, ...
                              "Samples_per_Frame", Samples_per_Frame, ...
                              "Chirp_Time_sec", Chirp_Time_sec, ...
                              "Pulse_Repetition_Time_sec", Pulse_Repetition_Time_sec, ...
                              "Frame_Period_sec", Frame_Period_sec);



% load singal data
    opts = delimitedTextImportOptions("NumVariables", 1);
    opts.DataLines = [30, Inf];
    opts.Delimiter = ",";
    opts.VariableNames = "data";
    opts.VariableTypes = "double";
    opts.ExtraColumnsRule = "ignore";
    opts.EmptyLineRule = "read";
    %IFX_radar_data = readtable([PathName,FileName], opts);
    IFX_radar_data = readtable([FileName], opts);
    clear opts

% Antenna
    switch Mask_Rx_Antennas
        case 1     % 001
            NumRXAntenna = 1;
            RXAntenna_1 = 1;
            RXAntenna_2 = 0;
            RXAntenna_3 = 0;
        case 2     % 010
            NumRXAntenna = 1;
            RXAntenna_1 = 0;
            RXAntenna_2 = 1;
            RXAntenna_3 = 0;
        case 3     % 011
            NumRXAntenna = 2;
            RXAntenna_1 = 1;
            RXAntenna_2 = 1;
            RXAntenna_3 = 0;
        case 4     % 100
            NumRXAntenna = 1;
            RXAntenna_1 = 0;
            RXAntenna_2 = 0;
            RXAntenna_3 = 1;
        case 5     % 101
            NumRXAntenna = 2;
            RXAntenna_1 = 1;
            RXAntenna_2 = 0;
            RXAntenna_3 = 1;
        case 6     % 110
            NumRXAntenna = 2;
            RXAntenna_1 = 0;
            RXAntenna_2 = 1;
            RXAntenna_3 = 1;
        case 7     % 111
            NumRXAntenna = 3;
            RXAntenna_1 = 1;
            RXAntenna_2 = 1;
            RXAntenna_3 = 1;
    end
    IFX_radar_data_noNan=IFX_radar_data.data(~any(isnan(IFX_radar_data.data),2),:);
    Frame_Number= floor(length(IFX_radar_data_noNan)/(Samples_per_Chirp*Chirps_per_Frame*NumRXAntenna));
    sn = 0:Samples_per_Chirp-1; % zero based sample number
    Frame=zeros(Samples_per_Chirp, Chirps_per_Frame, NumRXAntenna,Frame_Number);
% dispatch data
    for nf= 0: Frame_Number-1
        Chirp = zeros(Samples_per_Chirp, Chirps_per_Frame, NumRXAntenna);
        for nc = 0:Chirps_per_Frame-1
            for na = 0:NumRXAntenna-1
%                 IData = IFX_radar_data_noNan(1+ sn*NumRXAntenna + na + Samples_per_Chirp*nc**NumRXAntenna+Samples_per_Frame*nf); % real
                  IData = IFX_radar_data_noNan(1+ sn*NumRXAntenna + na + Samples_per_Chirp*nc*NumRXAntenna+Samples_per_Frame*nf); % real
               try
                Chirp(:,nc+1,na+1) = IData;
               end
            end
        end
        Frame(:,:,:,nf+1) = Chirp;
    end
end
