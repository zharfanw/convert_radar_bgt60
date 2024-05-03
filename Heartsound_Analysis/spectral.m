function [nfft, faxis, PS] = spectral(signal, N, fs)
    nfft = 2^nextpow2(N); % Next power of 2 from length of y
    sig_freq=fft(signal,nfft);
    PS=abs(sig_freq).^2;
    PS=PS/max(PS);  % normalize PS to its maximum
    faxis=fs/2*linspace(0,1,nfft/2+1);
