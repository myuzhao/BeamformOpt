%%传感器阵列波束优化设计与应用
 %%20170818
 %%myuzhao@163.com
 %%存在单干扰时的MVDR波束图 凹槽噪声法的基础
clc;
clear;
close all;

%扫描计算范围
snaps = 5000;

% %阵元位置
M = 10;
%%扫描参数设置
angle0 = 0;%主轴方向
angle1 = -30;%干扰方向
angle = linspace(-90,90,18000);
format = ['k.';'r-';'b-'];
inr = [-10 0 10]; 
noise = 1/sqrt(2)*(randn(M,snaps)+1i*randn(M,snaps));
fs = 10000;
f = 1000;
for i=1:3
    inr1 = inr(i);%干燥比
    s1 = 10^(inr1/20)*randn(1,snaps);
    ap = exp(-1i*pi*(0:1:M-1)'*sind(angle1)); 
    apis = ap*s1;  %%%阵列接收干扰信号

    apisn = apis + noise;
    Rin = apisn*apisn'/snaps;
    iRin = inv(Rin);
    
    as0 = exp(-1i*pi*(0:1:M-1)'*sind(angle0))
    w0 = iRin*as0/(as0'*iRin*as0);  
    as = exp(-1i*pi*(0:1:M-1)'*sind(angle));
    p = as'*w0;
    enegry_mvdr = 20*log10(abs(p));
    figure(1)
    hold on
    plot(angle,enegry_mvdr,format(i))
end
figure(1)
xlabel('方位/(^o)')
ylabel('波束/dB')
ylim([-60 3])
grid on
title('存在单干扰时的MVDR波束图')
legend('INR=-10dB','INR=0dB','INR=10dB')
% Create arrow
annotation('arrow',[0.401785714285714 0.401785714285714],...
    [0.854761904761905 0.65]);
