 %%传感器阵列波束优化设计与应用
 %%20170815
 %%myuzhao
 %%
clc;
clear ;
close all;

%扫描计算范围
freq = 1000;  %信号频率
fs = 10000; % 采样频率
c0 = 344;
snaps = 500;

% %阵元位置
M = 10;
d_lamda = 1/2;%阵元间距d与波长lamda的关系
d = d_lamda*c0/freq*[0:1:M-1];
%%扫描参数设置
angle0 = 10;
angle = linspace(-90,90,10000);

a0 = exp(-1i*2*pi*freq*d'*sind(angle0)/c0)/M;
w = exp(-1i*2*pi*freq*d'*sind(angle)/c0);   
p = w'*a0;

energy_cbf_P=20*log10(abs(p));
figure
plot(angle,energy_cbf_P)
xlabel('方位/(^o)')
ylabel('波束/dB')
ylim([-70 3])
grid on
title('波束图')

%%纯信号时波束扫描方位谱
snr = 0;
s = 10^(snr/20)*exp(-1i*2*pi*freq*[1/fs:1/fs:snaps/fs]);
as = exp(-1i*2*pi*freq*d'*sind(angle0)/c0)*s;
R = as*as'/snaps;
w = exp(1i*2*pi*freq*d'*sind(angle)/c0)/M;   
p = diag(w'*R*w);
energy_cbf_P = 10*log10(abs(p));
figure
plot(angle,energy_cbf_P)
xlabel('方位/(^o)')
ylabel('波束/dB')
ylim([-70 3])
grid on
title('纯信号时波束扫描方位谱')

%%SNR=30dB信噪比时波束扫描方位谱
snr = 30;
s = 10^(snr/20)*exp(-1i*2*pi*freq*[1/fs:1/fs:snaps/fs]);
as = exp(-1i*2*pi*freq*d'*sind(angle0)/c0)*s;
noise = 1/sqrt(2)*(randn(M,snaps)+1i*randn(M,snaps));

as = as + noise;
R = as*as'/snaps;
w = exp(1i*2*pi*freq*d'*sind(angle)/c0)/M;  
p= diag(w'*R*w);
energy_cbf_P=10*log10(abs(p));
figure
plot(angle,energy_cbf_P)
xlabel('方位/(^o)')
ylabel('波束/dB')
ylim([-10 33])
grid on
title('SNR=30dB时波束扫描方位谱')


%%不同方位时波束阵增益 10*logM
a0 = exp(1i*2*pi*freq*d'*sind(angle0)/c0);
pn = eye(10,10);
w = exp(1i*2*pi*freq*d'*sind(angle)/c0);
pn = diag(w'*pn*w);
ps = diag((w'*a0)*(w'*a0)');
G = abs(ps./pn);
G_dB = 10*log10(G);
figure
plot(angle,G_dB)
xlabel('方位/(^o)')
ylabel('阵增益/dB')
title('不同方位时波束阵增益')
ylim([-50 13])
grid on