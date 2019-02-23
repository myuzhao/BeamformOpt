%%传感器阵列波束优化设计与应用
 %%20171020
 %%myuzhao@163.com
 %%直线阵列L2准则优化抗阵列流型误差
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
theta0=0*pi/180;
theta1=linspace(-90,90,180);
theta=theta1*pi/180;

M=element_num;
d=1/2*lamd;


index_temp_r=find(theta1<12);
index_temp_l=find(theta1>-12);

theta2=[theta1(1:index_temp_l(1)) theta1(index_temp_r(end):end)];
theta3=theta1(index_temp_l(1):index_temp_r(end));
%导向
ap0=exp(1i*2*pi*d_lamda*[0:element_num-1]'*sind(theta0));
ap0_error=ap0+abs(ap0).*rand(size(ap0))*0.05;
%%扫描
ap=exp(1i*2*pi*d_lamda*[0:element_num-1]'*sind(theta1));
ap_error=ap+abs(ap).*rand(size(ap))*0.1;

aps=exp(1i*2*pi*d_lamda*[0:element_num-1]'*sind(theta2));  
aps_error=aps+abs(aps).*rand(size(aps))*0.05;  %%旁瓣

apm=exp(1i*2*pi*d_lamda*[0:element_num-1]'*sind(theta3));
apm_error=apm+abs(apm).*rand(size(apm))*0.05;  %%主瓣

figure
plot(theta1,20*log10(abs(ap0'/M*ap)),'k-.')
ylim([-60 10])
hold on


% % sedumi()
% % b=[1 0 zeros(1,2*M)].';
% % aa=[real(ap).' imag(ap).'];
% % A=[0 1 zeros(1,2*M);0 0 aa;1 -1 zeros(1,2*M);zeros(2*M,1) zeros(2*M,1) 0.05*diag(ones(1,2*M));0 0 aa;zeros(2*M,1) zeros(2*M,1) 0.05*diag(ones(1,2*M))].';
% % % A=[0 1 zeros(1,2*M);zeros(72,1) zeros(72,1) aa;1 -1 zeros(1,2*M);zeros(2*M,1) zeros(2*M,1) 0.05*diag(ones(1,2*M));zeros(72,1) zeros(72,1) aa;zeros(2*M,1) zeros(2*M,1) 0.05*diag(ones(1,2*M))].';
% % c=[0;zeros(1,1);0;zeros(2*M,1);-1*ones(1,1);zeros(2*M,1)];
% % K.q=[2 2*M+1 2*M+1];
% % [x,y]=sedumi(A,b,c,K);

b=[-1 0 zeros(1,2*M)].';
aa=[real(aps).' imag(aps).'];
aa1=[real(ap0).' imag(ap0).'];
temp_A1=[];
temp_A2=[];
SS=0.05;
for i=1:size(aa,1)
    temp_A1=[temp_A1;0 -1 zeros(1,2*M);0 0 -aa(i,:)];
end
for i=1:size(aa1,1)
    temp_A2=[temp_A2;0 0 -aa1(i,:);zeros(2*M,1) zeros(2*M,1) -SS*diag(ones(1,2*M))];
end

A=[temp_A1;-1 1 zeros(1,2*M);zeros(2*M,1) zeros(2*M,1) -SS*diag(ones(1,2*M));temp_A2].';
c=[zeros(size(aa,1)*2,1);0;zeros(2*M,1);repmat([-1 zeros(1,2*M)].',size(aa1,1),1)];
K.q=[2*ones(1,size(aa,1)) 2*M+1 (2*M+1)*ones(1,size(aa1,1))];
[x,y]=sedumi(A,b,c,K);
ww=y(3:12)+1i*y(13:22);
