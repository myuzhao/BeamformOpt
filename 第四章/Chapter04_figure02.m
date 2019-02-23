%%传感器阵列波束优化设计与应用
 %%20170818
 %%myuzhao@163.com
 %%DC加权在标准线列阵中的应用
 %%指定主瓣宽度
clc;
clear ;
close all;
% %阵元位置
M = 10;
d_lamda = 1/2;%阵元间距d与波长lamda的关系
%%参数设置
angle0 = 0;
angle = linspace(-90,90,10000);
%%DC
sll=[];
type=1;
  as = exp(-1i*pi*[0:M-1]'*sind(angle0));
for angle11 = [asin(2/M) asin(3/M) asin(4/M)]
    w_dc = DC_win(angle11,sll,d_lamda,M,type);
    w_dc = w_dc/sum(w_dc);
    w = w_dc.*exp(-1i*pi*[0:M-1]'*sind(angle));   
    p = w'*as;
    enegry_cbf = 20*log10(abs(p));
    figure(1)
    hold on
    plot(angle,enegry_cbf,'.')
end
legend('sin=2/M','sin=3/M','sin=4/M')
ylabel('DC加权波束/dB')
xlabel('方位/(^o)')
set(gca,'fontsize',16)
ylim([-80 3])
grid on
