%%传感器阵列波束优化设计与应用
 %%20170818
 %%myuzhao@163.com
 %%DC加权在标准线列阵中的应用
 %%指定旁瓣级,非0方向
clc;
clear ;
close all;

%扫描计算范围
freq =1000;  %信号频率
c0 = 344;
lamd=c0/freq;

% %阵元位置
element_num=10;
d_lamda=1/2;%阵元间距d与波长lamda的关系

%%扫描参数设置
theta0=-10*pi/180;
theta1=linspace(-90,90,10000);
theta=theta1*pi/180;
 
p=zeros(length(theta),1);
%%DC
M=element_num;
theta11=asin(2/M);
sll=-20;
lamd=c0/freq;
M=10;
d=1/2*lamd;
type=2;
w_dc=DC_win(theta11,sll,d,M,lamd,type)
w_dc=w_dc/sum(w_dc);
wp=w_dc.*exp(1i*2*pi*d_lamda*[0:element_num-1]'*sin(theta0))
ap=exp(1i*2*pi*d_lamda*[0:element_num-1]'*sin(theta));   
p=wp'*ap;
energy_cbf_P=20*log10(abs(p));
figure(1)
hold on
plot(theta1,energy_cbf_P,'.')
xlabel('方位/(^o)')
ylabel('DC加权波束/dB')
ylim([-60 3])
grid on

sll=-30;
w_dc=DC_win(theta11,sll,d,M,lamd,type)
w_dc=w_dc/sum(w_dc);
wp=w_dc.*exp(1i*2*pi*d_lamda*[0:element_num-1]'*sin(theta0))
ap=exp(1i*2*pi*d_lamda*[0:element_num-1]'*sin(theta));   
p=wp'*ap;
energy_cbf_P=20*log10(abs(p));
figure(1)
hold on
plot(theta1,energy_cbf_P,'--')
xlabel('方位/(^o)')
ylabel('DC加权波束/dB')
ylim([-60 3])
grid on


sll=-40;
w_dc=DC_win(theta11,sll,d,M,lamd,type)
w_dc=w_dc/sum(w_dc);
wp=w_dc.*exp(1i*2*pi*d_lamda*[0:element_num-1]'*sin(theta0))
ap=exp(1i*2*pi*d_lamda*[0:element_num-1]'*sin(theta));   
p=wp'*ap;
energy_cbf_P=20*log10(abs(p));
figure(1)
hold on
plot(theta1,energy_cbf_P,'-')
xlabel('方位/(^o)')
ylabel('DC加权波束/dB')
ylim([-60 3])
grid on

legend('sll=-20dB','sll=-30dB','sll=-40dB')
title('非零方向不同旁瓣DC加权对比')