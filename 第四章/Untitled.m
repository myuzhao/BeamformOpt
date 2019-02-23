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


%%最小化旁瓣
cvx_begin %
        variable w_msl(M) complex
        variable s_msl(1)
    minimize(s_msl)
    subject to
          abs(w_msl'*aps)<=s_msl;
          w_msl'*ap0==1;
%            s_msl<=-20;
cvx_end

hold on
energy_cbf_P_msl=20*log10(abs(w_msl'*ap));
plot(theta1,energy_cbf_P_msl,'b')


%%l2准则稳健
ss=0.05*(M^0.5);
cvx_begin %%
        variable w_l2(M) complex
        variable s_l2(1)
    minimize(s_l2)
    subject to
          for i=1:size(aps,2)
            abs(aps(:,i)'*w_l2)+ss*norm(w_l2,2)<=s_l2;
          end
          for i=1:size(ap0,2)
          {ss*w_l2,ap0(:,i)'*w_l2-1} <In> complex_lorentz(M);      
          end
%           s<=-20;
cvx_end

%%l1准则稳健
ss=0.05;
sum_w_l1=0;
cvx_begin %%
        variable w_l1(M) complex
        variable s_l1(1)
    minimize(s_l1)
    subject to
%           for ii=1:M
%                sum_w_l1=abs(w_l1(ii))+sum_w_l1;
%           end
          sum_w_l1=sum(abs(w_l1));
          for i=1:size(aps,2)
              abs(aps(:,i)'*w_l1)+ss*sum_w_l1<=s_l1;
          end
          for i=1:size(ap0,2)
          {ss*w_l1,ap0(:,i)'*w_l1-1} <In> complex_lorentz(M);      
          end
%           s<=-20;
cvx_end

hold on
energy_cbf_P_l2=20*log10(abs(w_l2'*ap));
plot(theta1,energy_cbf_P_l2,'g')

energy_cbf_P_l1=20*log10(abs(w_l1'*ap));
plot(theta1,energy_cbf_P_l1,'r--')


figure%% error
plot(theta1,20*log10(abs(w_l2'*ap_error)),'g-')
hold on
plot(theta1,20*log10(abs(w_msl'*ap_error)),'b-.')
plot(theta1,20*log10(abs(ap0'/M*ap_error)),'k-.')
plot(theta1,20*log10(abs(w_l1'*ap_error)),'r--')
