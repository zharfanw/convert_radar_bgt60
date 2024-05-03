function [SDNN, RMSSD, HTI] = HRV(NN)
    SDNN = std(NN);
    
    RMSSD = 0;
    for i = 1:length(NN)-1
        RMSSD = RMSSD + (NN(i) - NN(i+1))^2;
    end
    RMSSD = RMSSD/length(NN);
    RMSSD = sqrt(RMSSD);
    
    [H, ~] = histcounts(NN, linspace(500,1250,96));
    HTI = length(NN)/max(H);
end

