close all;clear all;
%题目：例6.1 P169
%      12元均匀线列阵，阵元间隔f0对应的半波长，-30°方向LFM信号，fl=0.5*f0;fu=f0;fs=5*f0;宽带波束形成，加权系数采用常规波束形成
%还有什么没做：1、验证分段不重复时的效果；2、功率谱归一化
%注1：Figure2与书中有细微差别（红线第二、第三点较明显）的原因是12个阵元接收到的信号构造不同，书中存在一个假想阵元0，12个信号依次延时，本仿真将阵元1定成参考点
%1、分段有什么要求2、哪些点是无效的，怎么拼在一起3、分段部分程序的改进4、做多少点的dft
%% 基本参数
f0=50;fl=0.5*f0;fu=f0;fs=5*f0;
N=512;%LFM信号一个周期的长度
N1=256;%缓存长度
T=1/fs*(N-1);
t0=0:1/fs:T;%处理LFM一个周期的所有数据
M=12;%阵元数
theta=-30*pi/180;%信号方向
%% 造各路LFM信号
m=0:M-1;%见注1
tau=-sin(theta)/2/fu*m;   %d应该取最大频率的半波长，如果阵元间距不够下，会有什么结果？                       
c=1500;%声速
d=c/2/fu;%阵元间距为半波长 %c=1500;%声速
                          %d=c/2/fu;%阵元间距为半波长
                          %theta=10*pi/180;%信号方向
                          %tau=-d*sin(theta)/c*0:M-1;                       
[T0 Tau0]=meshgrid(t0,tau);
tt=T0-Tau0;
t1=tt;               %各路延时处理
t1(tt>T)=0;t1(tt<0)=0;%将延迟前后多余部分置零

figure(1);%各阵元接收的波形
st_all=sin(2*pi*(fl+(fu-fl)/(2*T)*t1).*t1);%每一个行向量表示每一个阵元采样的结果
for i=1:12
    plot(1:512,st_all(i,1:512)/2+i)
    hold on
end
xlabel('i');ylabel('阵元号m');title('各阵元接收的波形');
axis([0 512 0 13]);

t(1:M,1:N1,1)=t1(1:M,1:N1);t(1:M,1:N1,2)=t1(1:M,0.5*N1+1:1.5*N1);t(1:M,1:N1,3)=t1(1:M,N1+1:N);%信号的重叠率为0.5
st=sin(2*pi*(fl+(fu-fl)/(2*T)*t).*t);%每一个行向量表示每一个阵元采样的结果
%% 利用dft公式求傅立叶变换，对每个频率进行波束形成，IDFT后进行拼接
L=256;%DFT点数
l=0:L-1;%l=0:size(t,2)-1;
%k1=floor(L*fl/fs);k2=ceil(L*fu/fs);%k1小于fl的最大频率对应的整数，k2大于fu的最小频率对应的整数
k1=21;k2=56;
k=k1:k2;%不同的k对应不同的频率
xk=zeros(M,size(k,2),3);
for i=1:3
xk(1:M,1:size(k,2),i)=st(1:M,1:size(l,2),i)*exp(-2i*pi*l.'*k/L);%公式6.5，每个行向量表示频谱系数
end
%对每一列进行波束形成
jiao0=-30;%期望信号方向
theta0=pi/2-jiao0*pi/180;%转换角度，并用弧度表示
v0=-[cos(theta0);sin(theta0)];
p0=[0:M-1;zeros(1,M)]';
fk=fs*k/L;
wc=1/M*exp(1i*2*pi/c*d*p0*v0*fk);%列向量，式（2.69）带入式（2.11），实际加权值取共轭
                                 %每一个列向量表示不同频率（与fk对应）的加权向量
                                 
k3=0:L-1;
yk1=zeros(1,size(k,2),3);yk2=zeros(1,L/2,3);yk=zeros(1,L,3);yt=zeros(1,N1,3);%声明变量大小
for i=1:3
yk1(1,1:size(k,2),i)=sum(xk(1:M,1:size(k,2),i).*wc);%对每一个频率进行波束形成
yk2(1,1:L/2,i)=[zeros(1,k1),yk1(1,1:size(k,2),i),zeros(1,L/2-1-k2)];%对没有计算的频率进行补零
yk(1,1:L,i)=[yk2(1,1:L/2,i),0,conj(fliplr(yk2(1,2:L/2,i)))];%为了使输出为实数，负频率带取正频率的共轭
%可以将上两句换成注视掉的两句                %  yk(1,(k1+1):(k2+1),i)=yk1(1,1:size(k,2),i);
                                           %  yk(1,(L-k2+1):(L-k1+1),i)=conj(fliplr(yk1(1,1:size(k,2),i)));
yt(1,1:N1,i)=yk(1,1:L,i)*exp(2i*pi*k3.'*l/L)/L;%IDFT
end

figure(2);%重叠率0.5时的频域数据
plot(k1:k2,20*log10(abs(yk1(1,1:size(k,2),1))),':dk',...
                                    'LineWidth',2,...
                                 	'MarkerEdgeColor','k',...
                                    'MarkerFaceColor','k',...
                                    'MarkerSize',5);
hold on;
plot(k1:k2,20*log10(abs(yk1(1,1:size(k,2),2))),'--ob',...
                                 	'LineWidth',2,...
                                   	'MarkerEdgeColor','b',...
                                   	'MarkerFaceColor','b',...
                                 	'MarkerSize',5);
plot(k1:k2,20*log10(abs(yk1(1,1:size(k,2),3))),'-dr',...
                                  	'LineWidth',2,...
                                   	'MarkerEdgeColor','r',...
                                   	'MarkerFaceColor','r',...
                                  	'MarkerSize',5);
legend('{\itn}=1','{\itn}=2','{\itn}=3','Location','South');
xlabel('\itk');
ylabel('|{{\itY}^{({\itn})}}({\itk})|/dB');
axis([k1 k2 -40 40]);
yout=[yt(1,1:192,1),yt(1,65:192,2),yt(1,65:256,3)];

figure(3);%输入LFM功率谱
st0=sin(2*pi*(fl+(fu-fl)/(2*T)*t0).*t0);%第一个阵元接收到的LFM的一个周期T内的信号，st_all(1,:)
xk_st0=fft(st0(1,:));
P_st0=xk_st0(1:N/2).*conj(xk_st0(1:N/2))/(N/2);
plot((0:N/2-1).*fs/N/f0,10*log10(P_st0),'LineWidth',2);
xlabel('频率({\itf}_0)');
ylabel('功率谱/(dB/Hz)')
grid on;

figure(4);plot(t0,real(yout)-st0);axis([0 2 -1.5 1.5 ]);
xlabel('\iti');
ylabel('{\ity}({\iti})-{\its}({\iti})');
title('波束输出序列与信号源波形失真大小');
%% 验证
ykk=fft(yout,5000);
f = fftshift(ykk);
w = linspace(-fs/2,fs/2,5000);%频率坐标，单位Hz 
figure(5),plot(w,abs(f)); title('宽带波束形成输出信号的频谱'); xlabel('频率(Hz)');