%%paper
%%Robust Adaptive Beamforming Using Worst-Case Performance Optimization A Solution to the Signal Mismatch Problem
%%Sergiy A. Vorobyov
%%code 
%%20171205
%%myuzhao
clc;
clear;
% close all;

%扫描计算范围
freq1 =1000;  %信号频率
freq0 =1000;  %信号频率
f=freq1;
fs = 10000; % 采样频率
c0 = 344;
ii=1;

snapshots_N=1:1:100;
SINR_OUT_SMI=zeros(size(snapshots_N));
SINR_OUT_LSMI=zeros(size(snapshots_N));
SINR_OUT_OPT=zeros(size(snapshots_N));
SINR_OUT_PRJ=zeros(size(snapshots_N));
SINR_OUT_PRJ=zeros(size(snapshots_N));
SINR_OUT_PRJ1=zeros(size(snapshots_N));
SINR_OUT_MPRJ=zeros(size(snapshots_N));
SNR=-30:1:25;
SINR_OUT_SMI=zeros(size(SNR));
SINR_OUT_LSMI=zeros(size(SNR));
SINR_OUT_OPT=zeros(size(SNR));
SINR_OUT_PRJ=zeros(size(SNR));
SINR_OUT_PRJ1=zeros(size(SNR));
SINR_OUT_MPRJ=zeros(size(SNR));
SINR_OUT_GLC=zeros(size(SNR));
SINR_OUT_WCO=zeros(size(SNR));
SINR_OUT_RECON=zeros(size(SNR));
for i_snap=100%1:1:length(snapshots_N)
    % %阵元位置
element_num=10;
M=element_num;
d_lamda=1/2;%阵元间距d与波长lamda的关系
d=d_lamda*c0/f*[0:1:element_num-1];
%%扫描参数设置
theta0=0*pi/180;%信号方向
% SNR=20;


for SNR=-30:1:25
%     SNR=-10;
for i_m=1:100  %%%
    theta1=20*pi/180;%干扰方向
    INR1=20;%干噪比
    
    theta2=60*pi/180;%干扰方向
    INR2=20;%干噪比

    theta=linspace(-90,90,360);
    theta=theta*pi/180;
    
    randn('state',30);
    noise=1/sqrt(2)*(randn(element_num,snapshots_N(i_snap))+1i*randn(element_num,snapshots_N(i_snap)));%
    %%signal
    s0=exp(1i*2*pi*f*[1/fs:1/fs:snapshots_N(i_snap)/fs]);%1/sqrt(2)*(randn(1,snapshots_N(i_snap))+i*randn(1,snapshots_N(i_snap)));%
    ap=10^(SNR/20)*exp(1i*2*pi*f*d'*sin(theta0)/c0); 
    aps=ap*s0;
%     s0=sqrt((10^(SNR/10)))*(randn(1,snapshots_N(i_snap))+1i*randn(1,snapshots_N(i_snap)))/sqrt(2);
%     ap=exp(1i*(0:M-1)'*pi*sin(theta0));
%     aps=ap*s0;
    %%IS
    s1=exp(1i*2*pi*(f+200)*[1/fs:1/fs:snapshots_N(i_snap)/fs]+1i*2*pi*rand(1));%1/sqrt(2)*(randn(1,snapshots_N(i_snap))+1i*randn(1,snapshots_N(i_snap)));%
    ap1=10^(INR1/20)*exp(1i*2*pi*f*d'*sin(theta1)/c0); 
    apn1=ap1*s1;
%     s1=sqrt((10^(INR1/10)))*(randn(1,snapshots_N(i_snap))+1i*randn(1,snapshots_N(i_snap)))/sqrt(2);
%     ap1=exp(1i*(0:M-1)'*pi*sin(theta1));
%     apn1=ap1*s1;

    %%
    s2=exp(1i*2*pi*(f+100)*[1/fs:1/fs:snapshots_N(i_snap)/fs]+1i*2*pi*rand(1));%1/sqrt(2)*(randn(1,snapshots_N(i_snap))+1i*randn(1,snapshots_N(i_snap)));%
    ap2=10^(INR2/20)*exp(1i*2*pi*f*d'*sin(theta2)/c0); 
    apn2=ap2*s2;
%     s2=sqrt((10^(INR2/10)))*(randn(1,snapshots_N(i_snap))+1i*randn(1,snapshots_N(i_snap)))/sqrt(2);
%     ap2=exp(1i*(0:M-1)'*pi*sin(theta2));
%     apn2=ap2*s2;
    aps_all=aps+apn1+apn2+noise; %%接收的数据
    X=aps_all;
    R = (aps_all*aps_all')/snapshots_N(i_snap);%%
%     aps_all=aps+apn1+apn2; %%接收的数据
%     R=aps_all*aps_all'/snapshots_N(i_snap)+eye(10);%%

%    R=(aps*aps'+apn1*apn1'+apn2*apn2')/snapshots_N(i_snap)+eye(10);
   

    
    Rin=((apn1+apn2+noise)*(apn1+apn2+noise)')/snapshots_N(i_snap);%%干扰和噪声协方差矩阵
  %  Rin=((apn1*apn1')+(apn2*apn2'))/snapshots_N(i_snap)+eye(10);%%干扰和噪声协方差矩阵
    iRin=inv(Rin);

    [v,ds]=eig(R);
    [ds_new,index]=sort(diag(ds),'descend');
    
%   as=exp(1i*2*pi*f*d'*sin((theta0+0)*pi/180)/c0);%期望波束+5*pi/180
    as=exp(1i*(0:M-1)'*pi*sin(theta0+rand(1)*5*pi/180));
    %%%MVDR CAPON
    w_opt=iRin*as/(as'*iRin*as);%;%
    %SMI 
    w_smi=inv(R)*as/(as'*inv(R)*as);% ;% 
    %LSMI DL 
    iR=inv(R+10*ds_new(end)*eye(size(R)));
    w_lsmi=iR*as/(as'*iR*as);%;%
     
    %%ESB PRJ 
    v=v(:,index);
    A_s_i=v(:,1:3);
    w_prj=inv(R)*A_s_i*A_s_i'*as;
%   w_prj1=A_s_i*diag(1./ds_new(1:3))*A_s_i'*as;

    %%GLC 线形组和R和I
    %%Fully Automatic Computation of Diagonal Loading Levels for Robust Adaptive Beamforming
    v_mprj=trace(R)/M;
    p_glc=0;
    p_sum=0;
    for i_glc=1:snapshots_N(i_snap)
        p_sum=p_sum+norm(aps_all(:,i_glc),'fro')^4;
    end
    p_glc=p_sum/(snapshots_N(i_snap)^2)-norm(R,'fro')^2/snapshots_N(i_snap);
    alfa=min(v_mprj*p_glc/(norm(R-v_mprj*eye(size(R)))^2),v_mprj);
    blfa=1-alfa/v_mprj;
    R_glc=alfa*eye(size(R))+blfa*R;
    w_glc=inv(R_glc)*as/(as'*inv(R_glc)*as);
    
    %%%%M_ESB M_PRJ  modify prj 
    %%%Modified projection approach for robust adaptive array beamforming
    %%%在小误差，低SNR时，与PRJ方法相比，有比较大的改善
    [v_mprj,ds_mprj]=eig(R_glc);
    [ds_mprj_new,index_mprj]=sort(diag(ds_mprj),'descend');
    v_mprj=v_mprj(:,index_mprj);
    p_mprj=0.8;
    A_s_i_mprj=[];
    sum_mprj=0;
    i_mprj=1;
    sum_mprj=abs(v_mprj(:,index_mprj(i_mprj))'*as)/i_mprj;
    while  (sum_mprj/i_mprj<p_mprj)&(i_mprj<10)
         i_mprj=i_mprj+1;
         sum_mprj=sum_mprj+abs(v_mprj(:,index_mprj(i_mprj))'*as);
    end
    A_s_i_mprj=v_mprj(:,1:(i_mprj-1));
    w_mprj=inv(R_glc)*A_s_i_mprj*A_s_i_mprj'*as;
    
%     w_prj1=A_s_i*diag(1./ds_new(1:3))*A_s_i'*as;
    %%%%%Worst-Case Performance Optimization
    %%%%%将误差限定在一个范围内
    %%%%%Robust adaptive bamforming using worst-case performance optimization: a solution to the signal mismatch problem’
    ss_wco=sqrt(3);%%%阵列流型误差范数上届
    U_chol=chol(R); %U'*U=R
    cvx_begin quiet
        variable w_wco(M) complex
        minimize(norm(U_chol * w_wco))
        subject to
                % w_wco'*as>=ss_wco*+1;
                {ss_wco*w_wco,w_wco'*as-1} <In> complex_lorentz(M)
                % imag(w_wco'*as)==0;
    cvx_end
   
    %%%Robust Adaptive Beamforming Based on Interference Covariance Matrix Reconstruction and Steering Vector Estimation   Yujie Gu and Amir Leshem   
    %%%Recon-based
    P_Recon_based=zeros(size(theta));
    iR_Recon=inv(R);
    %%%估计干扰噪声区域
%     for i_theat=1:length(theat)
%         ap_Recon=exp(1i*2*pi*f*d'*sin(theta(i_theat))/c0);   
%         P_Recon_based(i_theat)=1./(ap_Recon'*iR_Recon*ap_Recon);
%     end
        
%        P_Recon_based=10*log10(P_Recon_based/max(P_Recon_based));
     %  theta_index=find(P_Recon_based<-3);%%%
     
    %%%%%% %%%%%% %%%%%% %%%%%% %%%%%% %%%%%% %%%%%%
       theta_index=find(theta<theta0-10*pi/180|theta>theta0+10*pi/180); %%%%%%干扰噪声区域
       RIN_Recon=0;
       for i_index=1:length(theta_index)
            ap_Recon=exp(1i*2*pi*f*d'*sin(theta(i_index))/c0); 
            P_Recon_based(i_index)=1./(ap_Recon'*iR_Recon*ap_Recon)
            RIN_Recon=RIN_Recon+P_Recon_based(i_index)*ap_Recon*ap_Recon';
       end
       
       w_recon=inv(RIN_Recon)*as/(as'*inv(RIN_Recon)*as);
       
    
    
%     ap=exp(1i*2*pi*f*d'*sin(theta)/c0);   
%     p1=ap'*w_opt;
%     energy_mvdr_P1=20*log10(abs(p1));
%     figure
% %     figure1=figure(1);
%     hold on
%     plot(theta*180/pi,energy_mvdr_P1)
%     xlabel('方位/(^o)')
%     ylabel('波束/dB')
%     title('MVDR')
% %     ylim([-50 20])
%     grid on
%    
%     %     %%%%
%     ap=exp(1i*2*pi*f*d'*sin(theta)/c0);   
%     p2=ap'*w_smi;
%     energy_mvdr_P2=20*log10(abs(p2));
% %     figure2=figure(2);
%     hold on
%     plot(theta*180/pi,energy_mvdr_P2)
%     xlabel('方位/(^o)')
%     ylabel('波束/dB')
% %     ylim([-50 20])
%     title('SMI')
%     grid on
%     
%     ap=exp(1i*2*pi*f*d'*sin(theta)/c0);   
%     p3=ap'*w_lsmi;
%     energy_lsim_P3=20*log10(abs(p3));
% %     figure3=figure(3);
%     hold on
%     plot(theta*180/pi,energy_lsim_P3)
%     xlabel('方位/(^o)')
%     ylabel('波束/dB')
% %     ylim([-50 20])
%     title('LSMI')
%     grid on
% 
% %     %%%%
%     ap=exp(1i*2*pi*f*d'*sin(theta)/c0);   
%     p4=ap'*w_prj;
%     energy_prj_P4=20*log10(abs(p4));
% %     figure4=figure(4);
%     hold on
%     plot(theta*180/pi,energy_prj_P4)
%     xlabel('方位/(^o)')
%     ylabel('波束/dB')
% %     ylim([-50 20])
%     title('PRJ')
%     grid on
%     title(num2str(SNR))
%     legend('MVDR','SMI','LSMI','PRJ')
    SINR_OUT_SMI(ii)=SINR_OUT_SMI(ii)+10*log10(abs(10^(SNR/10)*(w_smi'*(as*as')*w_smi)/(w_smi'*Rin*w_smi)));
    SINR_OUT_LSMI(ii)=SINR_OUT_LSMI(ii)+10*log10(abs(10^(SNR/10)*(w_lsmi'*(as*as')*w_lsmi)/(w_lsmi'*Rin*w_lsmi)));
    SINR_OUT_OPT(ii)=SINR_OUT_OPT(ii)+10*log10(abs(10^(SNR/10)*(w_opt'*(as*as')*w_opt)/(w_opt'*Rin*w_opt)));%(w_opt'*(as*as')*w_opt)/(w_opt'*R*w_opt);
   
    SINR_OUT_PRJ(ii)=SINR_OUT_PRJ(ii)+10*log10(abs(10^(SNR/10)*(w_prj'*(as*as')*w_prj)/(w_prj'*Rin*w_prj)));
%     SINR_OUT_PRJ1(ii)=SINR_OUT_PRJ1(ii)+10*log10(abs(10^(SNR/10)*(w_prj1'*(as*as')*w_prj1)/(w_prj1'*Rin*w_prj1)));
    SINR_OUT_GLC(ii)=SINR_OUT_GLC(ii)+10*log10(abs(10^(SNR/10)*(w_glc'*(as*as')*w_glc)/(w_glc'*Rin*w_glc)));
    SINR_OUT_MPRJ(ii)=SINR_OUT_MPRJ(ii)+10*log10(abs(10^(SNR/10)*(w_mprj'*(as*as')*w_mprj)/(w_mprj'*Rin*w_mprj)));
    SINR_OUT_WCO(ii)=SINR_OUT_WCO(ii)+10*log10(abs(10^(SNR/10)*(w_wco'*(as*as')*w_wco)/(w_wco'*Rin*w_wco)));
   
    SINR_OUT_RECON(ii)=SINR_OUT_RECON(ii)+10*log10(abs(10^(SNR/10)*(w_recon'*(as*as')*w_recon)/(w_recon'*Rin*w_recon)));
    
end
    SINR_OUT_SMI(ii)=SINR_OUT_SMI(ii)/i_m;
    SINR_OUT_LSMI(ii)=SINR_OUT_LSMI(ii)/i_m;
    SINR_OUT_OPT(ii)=SINR_OUT_OPT(ii)/i_m;
    SINR_OUT_PRJ(ii)=SINR_OUT_PRJ(ii)/i_m;
    SINR_OUT_PRJ1(ii)=SINR_OUT_PRJ1(ii)/i_m;
    SINR_OUT_GLC(ii)=SINR_OUT_GLC(ii)/i_m;
    SINR_OUT_MPRJ(ii)=SINR_OUT_MPRJ(ii)/i_m;
    SINR_OUT_WCO(ii)=SINR_OUT_WCO(ii)/i_m;
    SINR_OUT_RECON(ii)=SINR_OUT_RECON(ii)/i_m;
    ii=ii+1;
end
end
% figure
% plot(snapshots_N,SINR_OUT_SMI,'-Ko')
% hold on
% plot(snapshots_N,SINR_OUT_LSMI,'-b*')
% plot(snapshots_N,SINR_OUT_OPT,'--r')
% plot(snapshots_N,SINR_OUT_PRJ,'-v')
% plot(snapshots_N,SINR_OUT_PRJ1,'->')
% 
SNR=-30:1:25;
figure
plot(SNR,SINR_OUT_SMI,'-ro')
hold on
plot(SNR,SINR_OUT_LSMI,'-*y')
plot(SNR,SINR_OUT_OPT,'--g')
plot(SNR,SINR_OUT_PRJ,'-vb')
plot(SNR,SINR_OUT_MPRJ,'->')
plot(SNR,SINR_OUT_GLC,'--square')
plot(SNR,SINR_OUT_WCO,'-.')
plot(SNR,SINR_OUT_RECON,'--o')

grid on
legend('SMI','LSMI','OPT','PRJ','MPRJ','GLC','WCO')
%figure1=figure(1)
% hold on
% plot(theta*180/pi,energy_mvdr_P1,type1(1))
% xlabel('方位/(^o)')
% ylabel('波束/dB')
% ylim([-60 3])
% grid on

% end

% 
% title('存在单干扰时的MVDR波束图')
% legend('INR=-10dB','INR=0dB','INR=10dB')
% % Create arrow
% annotation('arrow',[0.401785714285714 0.401785714285714],...
%     [0.854761904761905 0.65]);
