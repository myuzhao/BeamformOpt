%%%% 6.2 频域DFT波束形成器
clear all; clc
close all

%-----------产生线性调频信号---------------
% f0 = 1000; fL = f0/2; fU = f0; fs = 5*f0;
fL = 800; fU = 1800; fs = 8192;
T = 512/fs;
t = 0:1/fs:(T-1/fs);
% s = chirp(t,fL,T,fU);
s = sin(2*pi*(fL+(fU-fL)/(2*T)*t).*t);
figure; plot( t*fs, s , 'LineWidth', 1.0 );  axis([0 512 -1.5 1.5])

NFFT = 2^nextpow2(length(t)); % Next power of 2 from length of y
S = abs(fft(s,NFFT)).^2;
flfm = fs/2*linspace(0,1,NFFT/2+1);
% XdB = 10*log10(S(1:NFFT/2+1)/max(S(1:NFFT/2+1)));
XdB = 10*log10(S(1:NFFT/2+1));
figure; plot(flfm,XdB); xlabel('Frequency (Hz)'); ylabel('|S(f)|');grid on

f=flfm;
index1 = find(f>=0.0); index2 = find(f>=2500);
findex1 = index1(1,1); findex2 = index2(1,1);
f1 = f(1,findex1); f2 = f(1,findex2);
FN = findex2 - findex1 +1;

%  ifft 通过频谱求时域信号
SS = zeros(1,NFFT/2+1);
L=length(SS);
SS(1,findex1:findex2) = S(1,findex1:findex2);
sigifft = ifft(SS, NFFT) * ( 2 * L );
timeifft = 1/fs * linspace(0,NFFT/2,NFFT/2+1);
signalifft = real(sigifft(1,1:NFFT/2+1));
figure; plot(timeifft, signalifft);  grid on
Sifft = fft(signalifft,NFFT)/length(signalifft);
fi = fs/2*linspace(0,1,NFFT/2+1);
figure; plot(fi,2*abs(Sifft(1:NFFT/2+1))); grid on
SifftdB = 20*log10(abs(Sifft(1:NFFT/2+1))/max(abs(Sifft(1:NFFT/2+1))));
figure; plot(fi,SifftdB); grid on
