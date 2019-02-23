%%传感器阵列波束优化设计与应用
%%20181224
%%myuzhao
clc;
clear;
close all;
M = 10;
snaps = 6*M;
LNR = 10;
snr_all = -20:5:40;
Monte_Carlo = 200;

rand_state = 0;
rng('default');
rng(rand_state);
for error = [1 0]
    SINR_LSMI = zeros(1,length(snr_all));
    SINR_SMI = zeros(1,length(snr_all));
    for i_snr = 1:1:length(snr_all)
     
        for i_M = 1:1:Monte_Carlo
            noise = 1/sqrt(2)*(randn(M,snaps)+1i*randn(M,snaps));%高斯复噪声
            snr0 = snr_all(i_snr);
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
            if error
               as = exp(-1i*pi*(0:M-1)'*sind(angle0+1));
            else
               as = exp(-1i*pi*(0:M-1)'*sind(angle0)); 
            end
            a_int1 = exp(-1i*pi*(0:M-1)'*sind(angle1));
            a_int2 = exp(-1i*pi*(0:M-1)'*sind(angle2));
            %%%SMI
            iR = inv(R);
            w_smi = iR*as/(as'*iR*as);%;%
            %%%LSMI 
            iLR = inv(R + 10^(LNR/10) * eye(M) );
            w_lsmi = iLR*as/(as'*iLR*as);% ;% 
            
            as = exp(-1i*pi*(0:M-1)'*sind(angle0)); %%%%重要
            
            SINR_LSMI(1, i_snr) = SINR_LSMI(1, i_snr) +...
                  10^(snr0/10)*abs(w_lsmi'*as)^2/...
                  (10^(inr1/10)*abs(w_lsmi'*a_int1)^2+10^(inr2/10)*abs(w_lsmi'*a_int2)^2+w_lsmi'*w_lsmi);  
            
            SINR_SMI(1, i_snr) = SINR_SMI(1, i_snr) +...
                  10^(snr0/10)*abs(w_smi'*as)^2/...
                  (10^(inr1/10)*abs(w_smi'*a_int1)^2+10^(inr2/10)*abs(w_smi'*a_int2)^2+w_smi'*w_smi);  
        end
        SINR_LSMI(1, i_snr) = SINR_LSMI(1, i_snr)/Monte_Carlo;
        
        SINR_SMI(1, i_snr) = SINR_SMI(1, i_snr)/Monte_Carlo;
    end
    figure(1)
    plot(snr_all ,10*log10(abs(SINR_SMI)))
    hold on
    plot(snr_all ,10*log10(abs(SINR_LSMI)))
end
figure(1)
ylim([-20 50])
xlim([min(snr_all) max(snr_all)])
xlabel('SNR/dB');ylabel('SINR/dB');
legend('SMI,\theta_s =1','LSMI,\theta_s =1','SMI,\theta_s =0','LSMI,\theta_s =0');
set(gca,'fontsize',16)
