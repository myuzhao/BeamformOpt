%%%% 6.5.1 最小加权误差准则
clear all; clc
close all

L = 15;                   % 滤波器长度
D = (L-1)/2;
tau = 0.12345;
fd = linspace(0, 0.5, 100);
%-----------期望滤波器响应---------------
Hd = exp(-1i*2*pi*fd*(D+tau));

figure(1); hold on; grid on
plot(fd(1,1:80), 20*log10(abs(Hd(1,1:80))),'*'); axis([0 0.5 -1.5 0.2])
figure(2); hold on;
plot(fd(1,1:80), angle(Hd(1,1:80))*180/pi,'*');  axis([0 0.5 -200 200])

%---------- e(f) -------------
e = exp(-1i*2*pi*(0:L-1)'*fd);   

Lambda = [ones(1,80) zeros(1,20)];
%----- 利用yalmip求解 L无穷范数 -----%
y0 = sdpvar(1);
h0 = sdpvar(L,1,'full','real');

constraint0 = set( cone( Lambda(1,1)*(e(:,1).'*h0-Hd(1,1) ),y0 ) );
for ii = 2:length(fd)
    constraint0 = constraint0 + set( cone( Lambda(1,ii)*(e(:,ii).'*h0-Hd(1,ii) ),y0 ) );
end
constraints = constraint0;
obj = y0;
solvesdp(constraints,obj);
h_0 = double(h0);
%----- 利用yalmip求解 -----%
H0 = e.'*h_0;
% figure(1); plot(fd, 20*log10(abs(H0)),'b'); 
% figure(2); plot(fd, angle(H0)*180/pi,'b');  


%----- 利用yalmip求解 L1范数 -----%
h1 = sdpvar(L,1,'full','real');
y1 = sdpvar(length(fd),1);

constraint1 = set( cone( Lambda(1,1)*(e(:,1).'*h1-Hd(1,1) ),y1(1,1) ) );
% obj = y1;
for ii = 2:length(fd)
%     y1 = sdpvar(1);
    constraint1 = constraint1 + set( cone( Lambda(1,ii)*(e(:,ii).'*h1-Hd(1,ii) ),y1(ii,1) ) );
%     obj = obj + y1;
end
constraints = constraint1;
obj = sum(y1);
solvesdp(constraints,obj);
h_1 = double(h1);
%----- 利用yalmip求解 -----%
H1 = e.'*h_1;
% figure(1); plot(fd, 20*log10(abs(H1)),'r'); 
% figure(2); plot(fd, angle(H1)*180/pi,'r'); 


%----- 利用yalmip求解 L2范数 -----%
y2 = sdpvar(1);
h2 = sdpvar(L,1,'full','real');

zz = [];
for ii = 1 : length(fd)
    zz = [zz;sqrt(Lambda(1,ii))*e(:,ii).'*h2 - sqrt(Lambda(1,ii))*Hd(1,ii)];
end
constraint2 = set(cone(zz,y2));

constraints = constraint2;
obj = y2;
solvesdp(constraints,obj);
h_2 = double(h2);
%----- 利用yalmip求解 -----%
H2 = e.'*h_2;
figure(1); plot(fd, 20*log10(abs(H2)),'k'); 
figure(2); plot(fd, angle(H2)*180/pi,'k');  

%-----------产生线性调频信号---------------
f0 = 1000; fL = f0/2; fU = f0; fs = 5*f0;
% fL = 800; fU = 1800; fs = 8192;
T = 512/fs;
t = 0:1/fs:(T-1/fs);
% s = chirp(t,fL,T,fU);
s = sin(2*pi*(fL+(fU-fL)/(2*T)*t).*t);
figure; plot( t*fs, s , 'LineWidth', 1.0 );  axis([0 512 -1.5 1.5])
ylabel('幅度');

NFFT = 2^nextpow2(length(t)); % Next power of 2 from length of y
S = abs(fft(s,NFFT)).^2;
flfm = fs/2*linspace(0,1,NFFT/2+1);
XdB = 10*log10(S(1:NFFT/2+1)/max(S(1:NFFT/2+1)));
figure; plot(flfm,XdB); xlabel('Frequency (Hz)'); ylabel('|S(f)|');grid on
%-----------产生线性调频信号---------------

a = 1;
y = filter(h_2,a,s);
[hh,ww]=freqs(h_2,1);
temp = zeros(1,length(y));
temp(1,1:(length(y)-D)) = y(1,(D+1):length(y));
y = temp;
figure(5);subplot(3,1,2); plot( t*fs, y , 'LineWidth', 1.0 );  axis([0 512 -1.5 1.5])
ylabel('幅度');
title('滤波器延迟波形');

Ts = 1/fs;
tt = t - tau*Ts;
s_delay = sin(2*pi*(fL+(fU-fL)/(2*T)*tt).*tt);
figure(5);subplot(3,1,1); plot( t*fs, s_delay , 'LineWidth', 1.0 );  axis([0 512 -1.5 1.5])
ylabel('幅度');
title('理想延迟波形');
err = y - s_delay;
figure(5);subplot(3,1,3); plot( t*fs, err , 'LineWidth', 1.0 );  axis([0 512 -5e-3 5e-3])
xlabel('采样点数');
ylabel('幅度');
title('波形误差');