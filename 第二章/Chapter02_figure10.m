%%传感器阵列波束优化设计与应用
 %%20170815
 %%myuzhao
 %%
clc;
clear;
close all;

%扫描计算范围
freq = 1000;  %信号0频率
fs = 10000; % 采样频率
c0 = 344;
snaps = 500;

% %阵元位置
M = 10;
d_lamda = 1/2;%阵元间距d与波长lamda的关系
d = d_lamda * c0/freq * [0:1:M-1];
%%扫描参数设置
angle0 = 10;
snr = 30;
angle = linspace(-90,90,10000);
 
noise = 1/sqrt(2)*(randn(M,snaps)+1i*randn(M,snaps));
s1 = 10^(snr/20)*exp(1i*2*pi*freq*[1/fs:1/fs:snaps/fs]);
as = exp(1i*2*pi*freq*d'*sind(angle0)/c0)*s1; 

as = as+noise;
R = as*as'/snaps;

%%%%%%%%%%%%%%%%%%%%%%%%%期望方向 10度%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
as0 = exp(1i*2*pi*freq*d'*sind(angle0)/c0);%期望波束
iR = eye(M);
w0 = iR*as0/(as0'*iR*as0);
as = exp(1i*2*pi*freq*d'*sind(angle)/c0);   
p = as'*w0;
energy_mvdr_P1 = 20*log10(abs(p));

%%%%%%%%%%%%%%%%%%%%%%%%%期望方向-10度%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
angle1 = -10;
as1 = exp(1i*2*pi*freq*d'*sind(angle1)/c0);%期望波束
% iR=inv(R);
w = iR*as1/(as1'*iR*as1);
as = exp(1i*2*pi*freq*d'*sind(angle)/c0);   
p = as'*w;
energy_mvdr_P2 = 20*log10(abs(p));

figure
plot(angle,energy_mvdr_P1,'b--',angle,energy_mvdr_P2,'k-')
xlabel('方位/(^o)')
ylabel('波束/dB')
ylim([-100 30])
grid on
legend('10','-10')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%方位谱%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%MVDR
iR=inv(R);
as = exp(1i*2*pi*freq*d'*sind(angle)/c0);  
p_mvdr = diag(1./abs(as'*iR*as));
energy_mvdr_P = 10*log10(abs(p_mvdr));
figure
plot(10,30,'bo')
hold on
plot(angle,energy_mvdr_P)
xlabel('方位/(^o)')
ylabel('波束/dB')
ylim([-20 40])
grid on
legend('真实信号','方位谱')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%加权范数%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:length(angle)
    as = exp(1i*2*pi*freq*d'*sind(angle(i))/c0);%期望波束
    if angle(i) ==angle0
         iR = eye(M);
    else
         iR = inv(R);
    end
    w = iR*as/(as'*iR*as);
    w_2(i)=10*log10(norm(w,2));
end
figure
plot(angle,w_2)