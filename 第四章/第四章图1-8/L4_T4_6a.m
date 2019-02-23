%%传感器阵列波束优化设计与应用
 %%20170818
 %%myuzhao@163.com
 %%扇面内存在多个干扰时的 MVDR 波束图
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
theta1=linspace(-50,-20,16)*pi/180;%干扰方向
type1=['k.';'r-';'b-'];
snr=[-10 10];
for i=1:2
snr1=snr(i);%干燥比
theta=linspace(-90,90,10000);
theta=theta*pi/180;
 
noise=1/sqrt(2)*(randn(element_num,snapshots_N)+1i*randn(element_num,snapshots_N));
% s1=10^(snr1/20)*exp(1i*2*pi*f*[1/fs:1/fs:snapshots_N/fs]);
s1=10^(snr1/20)*exp(1i*2*pi*f*randn(16,snapshots_N));%%假设干扰信号不相干
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
ylim([-100 3])
grid on

end
title('扇面内存在多个干扰时的 MVDR 波束图')
legend('INR=-10dB','INR=10dB')


annotation(figure1,'arrow',[0.388534783406755 0.388534783406756],...
    [0.898546405360363 0.763832119646074]);

% Create arrow
annotation(figure1,'arrow',[0.408195668135096 0.408195668135097],...
    [0.897852124045927 0.763137838331637]);

% Create arrow
annotation(figure1,'arrow',[0.398045154185022 0.398045154185023],...
    [0.897424070904696 0.762709785190406]);

% Create arrow
annotation(figure1,'arrow',[0.378384269456682 0.378384269456683],...
    [0.899402511642826 0.764688225928536]);

% Create arrow
annotation(figure1,'arrow',[0.368694933920705 0.368694933920706],...
    [0.898280177187158 0.763565891472868]);

% Create arrow
annotation(figure1,'arrow',[0.358083241556534 0.358083241556535],...
    [0.900686671066521 0.765972385352231]);

% Create arrow
annotation(figure1,'arrow',[0.3490937041116 0.349093704111601],...
    [0.901380952380957 0.766666666666667]);

% Create arrow
annotation(figure1,'arrow',[0.324999999999999 0.325],...
    [0.901380952380957 0.766666666666667]);

% Create arrow
annotation(figure1,'arrow',[0.339583333333333 0.339583333333333],...
    [0.898013949013954 0.763299663299664]);