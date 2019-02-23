%%传感器阵列波束优化设计与应用
%%20181215
%%myuzhao
clc;
clear;
close all;
M = 10;
snaps = 2*M;
lnr_all = -20:5:70;
Monte_Carlo = 200;
snr_all = [-20 -10 0 30];
SINR_LSMI = zeros(length(snr_all),length(lnr_all));
SINR_OPT = zeros(length(snr_all),length(lnr_all));
w2 = zeros(length(snr_all),length(lnr_all));
rand_state = 0;
rng('default');
rng(rand_state);
for i_snr = 1 :length(snr_all)
    snr0 = snr_all(i_snr);
    for i_lnr = 1:1:length(lnr_all)
        LNR = lnr_all(i_lnr);
        for i_M = 1:1:Monte_Carlo
            noise = 1/sqrt(2)*(randn(M,snaps)+1i*randn(M,snaps));%高斯复噪声
           
            angle0 = 0;
            s0 = randn(1, snaps) ;%;exp(1i*2*pi*f*(1/fs:1/fs:snaps/fs));%
            as0 = 10^(snr0/20)*exp(-1i*pi*(0:1:M-1)'*sind(angle0)) * s0; 

            inr1 = 30;
            angle1 = 20;
            s1 = randn(1, snaps) ;%;exp(1i*2*pi*(f+100)*(1/fs:1/fs:snaps/fs)+1i*2*pi*rand(1));%randn(1, snaps);
            an1 = 10^(inr1/20)*exp(-1i*pi*(0:1:M-1)'*sind(angle1)) * s1; 

            inr2 = 35;
            angle2 = 35;
            s2 = randn(1, snaps) ;%;exp(1i*2*pi*(f+200)*(1/fs:1/fs:snaps/fs)+1i*2*pi*rand(1));%randn(1, snaps);
            an2 = 10^(inr2/20)*exp(-1i*pi*(0:1:M-1)'*sind(angle2)) * s2; 

            as_all = as0 + an1 + an2 + noise; %%接收的数据
            R = (as_all*as_all')/(snaps);     %%数据协方差矩阵
            
%           Rin = ((an1+an2+noise)*(an1+an2+noise)')/snaps;%%干扰和噪声协方差矩阵
            Rin = (an1*an1' + an2*an2'+ noise*noise')/snaps;%%干扰和噪声协方差矩阵
            as = exp(-1i*pi*(0:M-1)'*sind(angle0));
            a_int1 = exp(-1i*pi*(0:M-1)'*sind(angle1));
            a_int2 = exp(-1i*pi*(0:M-1)'*sind(angle2));
            %%%MVDR CAPON
            iRin = inv(Rin);
            w_opt = iRin*as/(as'*iRin*as);%;%
            %%%LSMI 
            iR = inv(R + 10^(LNR/10) * eye(M) );
            w_lsmi = iR*as/(as'*iR*as);% ;%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %               不能使用这种方法，会存在R的计算误差。
% %             SINR_SMI(i_snr, i_snaps) = SINR_SMI(i_snr, i_snaps) +...
% %                    abs(w_smi'*R*w_smi)/abs(w_smi'*Rin*w_smi);%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            SINR_LSMI(i_snr, i_lnr) = SINR_LSMI(i_snr, i_lnr) +...
                  10^(snr0/10)*abs(w_lsmi'*as)^2/...
                  (10^(inr1/10)*abs(w_lsmi'*a_int1)^2+10^(inr2/10)*abs(w_lsmi'*a_int2)^2+w_lsmi'*w_lsmi);  
%             w2(i_lnr, i_snr)  = norm(w_lsmi,2)^2;    
        end
        SINR_LSMI(i_snr, i_lnr) = SINR_LSMI(i_snr, i_lnr)/Monte_Carlo;
    end
    figure(1)
    plot(lnr_all ,10*log10(abs(SINR_LSMI(i_snr,:))))
    hold on
   
end
figure(1)
ylim([-40 20])
xlim([min(lnr_all) max(lnr_all)])
xlabel('LNR/dB');ylabel('SINR/dB');
legend('SNR = -20dB','SNR = -10dB','SNR = 0dB','SNR = 30dB');
set(gca,'fontsize',16)
grid on