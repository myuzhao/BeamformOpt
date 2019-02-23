clc;
clear all;
% close all;
%%
%%%%参数设置%%%
N=12;
fu=1000;
f1=fu/2;
fs=3.125*fu;
f=linspace(f1,fu,33);
K=length(f);
f0=(f1+fu)/2;  
Ts=1/fs;
c=1500;
d=c/fu/2*[0:N-1];
theta=(-90:2:90);
thetaML=(-8:2:28);
thetaSL=[(-90:2:-12) (32:2:90)];
thetas=10;
SL=-30;
%%
L=15;
Ts=1/fs;
taus=-d'*sind(thetas)/c;
Tm=-round(taus/Ts)*(Ts);
%%
%%%已知量求解
a0_s=exp(1i*2*pi*f0*d'*sind(thetas)/c); 
a0_ML=exp(1i*2*pi*f0*d'*sind(thetaML)/c); 
e0=exp(-1i*(0:L-1)'*2*pi*(f0/fs));  
ka0=exp(-1i*2*pi*f0*Tm); 
temps=cheng(a0_s,ka0);                    
u_s=kron(e0,temps);
temp_ML=cheng(a0_ML,ka0);                    
u_ML=kron(e0,temp_ML);   
for k=1:length(f)     
    fk=f(k);
    ak_ML=exp(1i*2*pi*fk*d'*sind(thetaML)/c); 
    ak_SL=exp(1i*2*pi*fk*d'*sind(thetaSL)/c); 
    ek=exp(-1i*(0:L-1)'*2*pi*(fk/fs));     
    ka=exp(-1i*2*pi*fk*Tm);  
    temp_ML=cheng(ak_ML,ka);                     
    u_MLk=kron(ek,temp_ML);   
    temp_SL=cheng(ak_SL,ka);
    u_SL=kron(ek,temp_SL);
    U_SL(:,:,k)=u_SL;
    delta_ML(:,:,k)=u_MLk-u_ML;
end

%%%%旁瓣峰值约束Minimax主瓣差异波束设计（式7.45）
cvx_begin quiet
     variable hh(N*L) complex
     variable s_s(1)
  minimize(s_s)
  subject to
     u_s'*hh==1;
for k=1:K
    delta_p_ML(:,k)=delta_ML(:,:,k)'*hh;
%     abs(delta_p_ML(:,k))<=s_s;
    pk_SL(:,k)=U_SL(:,:,k)'*hh;
%      abs(pk_SL(:,k))<=10^(SL/20);
end
max(abs(delta_p_ML))<=s_s;
abs(pk_SL(:,1:K))<=10^(SL/20);
  
% norm(hh,2)<=3;
cvx_end
%%
for k=1:length(f)
     fk=f(k);
     ak=exp(1i*2*pi*fk*d'*sind(theta)/c); 
     ek=exp(-1i*(0:L-1)'*2*pi*(fk/fs));
     ka=exp(-1i*2*pi*fk*Tm);
     temp=cheng(ak,ka);
     u=kron(ek,temp);
     pk=hh'*u;
     energy_p(k,:)=20*log10(abs(pk));
     figure(1)
     plot(theta,20*log10(abs(pk)));
     title('旁瓣峰值约束Minimax主瓣差异波束设计');
     hold on 
end
 



