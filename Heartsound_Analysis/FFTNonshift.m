

function [X,w]=FFTNonshift(x,N)

N1=length(x);
X=fft(x,N)/N1;
ri=2*pi/N; rad=0:ri:2*pi-ri;

w=rad;