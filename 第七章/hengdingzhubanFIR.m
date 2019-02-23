clc;
clear all;
close all;
%%
N=12;
f0=1000;
f1=f0/2;
fu=f0;

fs=3.125*f0;
c=1500;
d=c/f0/2*[0:N-1];
Fpb=(0.16:0.01:0.32);
Fsb=[(0:0.01:0.13) (0.35:0.01:0.5)];
Ftb=[(0.14:0.01:0.15) (0.33:0.01:0.34)];
fpb=Fpb*fs;
fsb=Fsb*fs;
ftb=Ftb*fs;
theta=(-90:2:90);
thetaML=(-8:2:28);
thetaSL=[(-90:2:-12) (32:2:90)];
thetas=10;
%%
%%%求期望波束%%%
    SL=-25;
    wd=exp(1i*2*pi*(f0/2)*d'*sind(thetas)/c)/N;
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
      hold on
    plot(thetaSL,energy_cbf_PSL,'r');
      hold on
    legend('常规波束','期望主瓣','期望旁瓣')
    title('期望波束')
    xlabel('方位/(^o)')
    ylabel('波束/dB')
    ylim([-60 3])
    grid on
 %%
f=(0:0.01:0.5)*fs;
p_d_ML=cbf_p(:,42:60);
L=25;
Ts=1/fs;
taus=-d'*sind(thetas)/c;
Tm=-round(taus/Ts)*(Ts);


sum=0;
cvx_begin quiet
    variable hh(N*L) complex
    variable s_s(1)
    minimize(s_s)
    subject to
        for k=1:length(Fpb)     
            fk=Fpb(k)*fs;
            ak_ML=exp(1i*2*pi*fk*d'*sind(thetaML)/c); %N*Nml
            ak_SL=exp(1i*2*pi*fk*d'*sind(thetaSL)/c); 
            ek=exp(-1i*(0:L-1)'*2*pi*(fk/fs));     %L*1
%             ka=1;%exp(-1i*2*pi*fk*Tm);                %N*1
            temp_ML=ak_ML;%cheng(ak_ML,ka);                     %N*Nml
            u_ML=kron(ek,temp_ML);                       %(L*N)*(Nml) 
            pk_ML(k,:)=hh'*u_ML;                              %1*Nml
            temp2(k,:)=pk_ML(k,:)-p_d_ML  ;                    %1*Nml
            temp_SL=ak_SL;%cheng(ak_SL,ka);  
            u_SL=kron(ek,temp_SL);
            pk_SL(k,:)=hh'*u_SL;
        end
        max(abs(temp2))<=s_s;
        abs(pk_SL)<=10^(SL/20);
        norm(hh,2)<=0.25;
cvx_end

for k=1:length(fpb)
    fk=fpb(k);
    ak=exp(1i*2*pi*fk*d'*sind(theta)/c);
    ek=exp(-1i*(0:L-1)'*2*pi*(fk/fs));
    ka=1;%exp(-1i*2*pi*fk*Tm);
    temp=ak;%cheng(ak,ka);
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
 zlim([-200 0])
 xlabel('方位/(^o)')
 ylabel('f/fs')
 zlabel('波束/dB')
 title('宽带波束图');










 
 
 
 
   
   
   
   
   
   
   
   
   
   