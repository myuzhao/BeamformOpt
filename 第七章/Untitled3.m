clc;
clear all;
close all;
%%
array_num=12;%%传感器数目
f0=1000;
f1=f0/2;%%最小f
fu=f0;%%最大f
f_jg=0.01;%归一化间隔
fs=3.125*f0;%采样f
c=1500;%声速
d=c/f0/2*[0:array_num-1];%%半波长距离
Fpb=(0.16:0.01:0.32);%%通带频率
Fsb=[(0:0.01:0.13) (0.35:0.01:0.5)];%%阻带频率
Ftb=[(0.14:0.01:0.15) (0.33:0.01:0.34)];%%过渡带频率
fpb=Fpb*fs;%%通带频率
fsb=Fsb*fs;%阻带频率
ftb=Ftb*fs;%%过渡带频率
theta=(-90:2:90);%%角度
thetaML=(-8:2:28);%%主瓣角度
thetaSL=[(-90:2:-12) (32:2:90)];%%旁瓣角度
theta0=10;%%导向角度
%%%求期望波束%%%
SL=-25;
wd=exp(1i*2*pi*(f0/2)*d'*sind(theta0)/c)/array_num;
a=exp(1i*2*pi*(f0/2)*d'*sind(theta)/c);
cbf_p=wd'*a;
energy_cbf_P=20*log10(abs(cbf_p));
energy_cbf_PML=energy_cbf_P(:,42:60);
energy_cbf_PSL=[energy_cbf_P(:,1:40) energy_cbf_P(:,62:end)];
energy_cbf_PSL(:,1:end)=SL;
figure(1)
hold on
plot(theta,energy_cbf_P,'k-');
hold on
scatter(thetaML,energy_cbf_PML,'*');
plot(thetaSL,energy_cbf_PSL,'r');
legend('常规波束','期望主瓣','期望旁瓣')
title('期望波束')
xlabel('方位/(^o)')
ylabel('波束/dB')
ylim([-60 3])
grid on
%%
f=(0:0.01:0.5)*fs;%%%所有频率
p_d_ML=cbf_p(:,42:60);%%%主瓣
L=25;%%滤波器长度
Ts=1/fs;%%采样间隔
taus=-d'*sind(theta0)/c;
Tm=-(taus/Ts+(L-1)/2)*(Ts);


sum=0;
cvx_begin quiet
    variable hh(N*L) 
    variable s_s(1)
    minimize(s_s)
    subject to
        for k=1:length(Fpb)     
            fk=Fpb(k)*fs;
            ak_ML=exp(1i*2*pi*fk*d'*sind(thetaML)/c); %N*Nml
            ak_SL=exp(1i*2*pi*fk*d'*sind(thetaSL)/c); 
            ek=exp(-1i*(0:L-1)'*2*pi*(fk/fs));     %L*1
            ka=exp(-1i*2*pi*fk*Tm);                %N*1
            temp_ML=cheng(ak_ML,ka);                     %N*Nml
            u_ML=kron(ek,temp_ML);                       %(L*N)*(Nml) 
            pk_ML(k,:)=hh'*u_ML;                              %1*Nml
            temp2(k,:)=pk_ML(k,:)-p_d_ML  ;                    %1*Nml
            temp_SL=cheng(ak_SL,ka);  
            u_SL=kron(ek,temp_SL);
            pk_SL(k,:)=hh'*u_SL;
        end
        max(abs(temp2))<=s_s;
        abs(pk_SL)<=10^(SL/20);
%         norm(hh,2)<=0.25;
cvx_end

for k=1:length(fpb)
    fk=fpb(k);
    ak=exp(1i*2*pi*fk*d'*sind(theta)/c);
    ek=exp(-1i*(0:L-1)'*2*pi*(fk/fs));
    ka=exp(-1i*2*pi*fk*Tm);
    temp=cheng(ak,ka);
    u=kron(ek,temp);
    pk=hh'*u;
    energy_p(k,:)=20*log10(abs(pk));
    figure(5)
    plot(theta,20*log10(abs(pk)));
    hold on
  end

 figure(7)
 [ff,tt]=meshgrid(theta,fpb/fs);
 surf(ff,tt,energy_p);
 zlim([-80 0])
 xlabel('方位/(^o)')
 ylabel('f/fs')
 zlabel('波束/dB')
 title('宽带波束图');










 
 
 
 
   
   
   
   
   
   
   
   
   
   