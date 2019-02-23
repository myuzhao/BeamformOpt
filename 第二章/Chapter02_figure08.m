 %%传感器阵列波束优化设计与应用
 %%20170813
 %%myuzhao
clc;
clear;
close all;

%扫描计算范围
freq = 1000;  %信号频率
fs = 10000; % 采样频率
c0 = 344;

% %阵元位置
M = 10;
d_lamda=1/2;%阵元间距d与波长lamda的关系

%%扫描参数设置
angle0 = 0;
angle = linspace(-90,90,10000);


a0 = exp(1i*2*pi*d_lamda*[0:M-1]'*sind(angle0))/M;
w = exp(1i*2*pi*d_lamda*[0:M-1]'*sind(angle));   
p = w'*a0;

energy_cbf_P = 20*log10(abs(p));
figure
plot(angle,energy_cbf_P)
xlabel('方位/(^o)')
ylabel('波束/dB')
ylim([-60 3])
grid on

figure
plot(angle,energy_cbf_P)
xlabel('方位/(^o)')
ylabel('波束/dB')
ylim([-50 0])
xlim([-20 20])
grid on
hold on
plot(angle,-3*ones(length(angle),1),'k--')
plot(angle,-13*ones(length(angle),1),'k--')
y=[-50:2:-20];
x1=11.53*ones(length(y),1);
x2=-11.53*ones(length(y),1);
plot(x1,y,'k--',x2,y,'k--')
xx=[-9.316 9.316];yy=[-13 -13];
annotation('doublearrow',[0.294270833333333 0.741666666666667],...
    [0.599207684319834 0.595015576323987])
annotation('doublearrow',[0.336458333333333 0.6984375],...
    [0.711357217030114 0.712357217030114]);
annotation('doublearrow',[0.419270833333333 0.616666666666667],...
    [0.875389408099688 0.876427829698858]);
