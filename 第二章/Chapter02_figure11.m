%%传感器阵列波束优化设计与应用
 %%20170813
 %%myuzhao
 %%
clc;
clear;
close all;

%扫描计算范围
freq =1000;  %信号频率
fs = 10000; % 采样频率
c0 = 344;
snaps = 500;

% %阵元位置
M = 10;
d_lamda=1/2;%阵元间距d与波长lamda的关系
d=d_lamda*c0/freq*[0:1:M-1];
%%扫描参数设置
angle0 = -10;
angle1 = 10;
snr0 = 15;
snr1 = 30;
angle = linspace(-90,90,10000);
 
noise=1/sqrt(2)*(randn(M,snaps)+1i*randn(M,snaps));

s0 = 10^(snr0/20)*exp(1i*2*pi*freq*randn(1,snaps));
s1 = 10^(snr1/20)*exp(1i*2*pi*freq*randn(1,snaps));%不同频不同相位

as0=exp(1i*2*pi*freq*d'*sind(angle0)/c0)*s0;
as1=exp(1i*2*pi*freq*d'*sind(angle1)/c0)*s1;


aps = as0+as1+noise;
R = aps*aps'/snaps;
% R=(ap0*ap0'+ap1*ap1'+noise*noise')/snapshots_N;
ap = exp(1i*2*pi*freq*d'*sind(angle)/c0);
wp = ap/(M);   
p = diag(wp'*R*wp);

energy_cbf_P=10*log10(abs(p));
figure
plot(-10,15,'bo',10,30,'bo')
hold on
plot(angle,energy_cbf_P)
xlabel('方位/(^o)')
ylabel('波束/dB')
ylim([-20 40])
grid on
legend('真实信号','方位谱')


%%MVDR
iR=inv(R);
p_mvdr=abs(diag(1./(ap'*iR*ap)));
energy_mvdr_P=10*log10(abs(p_mvdr));
figure
plot(-10,15,'bo',10,30,'bo')
hold on
plot(angle,energy_mvdr_P)
xlabel('方位/(^o)')
ylabel('波束/dB')
% ylim([-20 40])
grid on
legend('真实信号','方位谱')

