%%传感器阵列波束优化设计与应用
 %%20170818
 %%myuzhao@163.com
 %%存在单干扰时的MVDR波束图 凹槽噪声法的基础
clc;
clear;
close all;

%扫描计算范围
freq1 =1000;  %信号频率
freq0 =1000;  %信号频率
f=freq1;
fs = 10000; % 采样频率
c0 = 344;
snapshots_N=5000;

% %阵元位置
element_num=10;
d_lamda=1/2;%阵元间距d与波长lamda的关系
d=d_lamda*c0/f*[0:1:element_num-1];
%%扫描参数设置
theta0=0*pi/180;%主轴方向
theta1=-30*pi/180;%干扰方向
type1=['k.';'r-';'b-'];
snr=[-10 0 10];
for i=1:3
snr1=snr(i);%干燥比
theta=linspace(-90,90,10000);
theta=theta*pi/180;
 
noise=1/sqrt(2)*(randn(element_num,snapshots_N)+1i*randn(element_num,snapshots_N));
s1=10^(snr1/20)*exp(1i*2*pi*f*[1/fs:1/fs:snapshots_N/fs]);
ap=exp(1i*2*pi*f*d'*sin(theta1)/c0); 
aps=ap*s1;

aps=aps+noise;
R=aps*aps'/snapshots_N;
% R=(ap0*ap0'+ap1*ap1'+noise*noise')/snapshots_N;

as=exp(1i*2*pi*f*d'*sin(theta0)/c0);%期望波束
iR=inv(R);
w=iR*as/(as'*iR*as);
ap=exp(1i*2*pi*f*d'*sin(theta)/c0);   
p=ap'*w;
energy_mvdr_P1=20*log10(abs(p));


figure1=figure(1)
hold on
plot(theta*180/pi,energy_mvdr_P1,type1(i))
xlabel('方位/(^o)')
ylabel('波束/dB')
ylim([-60 3])
grid on

end


title('存在单干扰时的MVDR波束图')
legend('INR=-10dB','INR=0dB','INR=10dB')
% Create arrow
annotation('arrow',[0.401785714285714 0.401785714285714],...
    [0.854761904761905 0.65]);
