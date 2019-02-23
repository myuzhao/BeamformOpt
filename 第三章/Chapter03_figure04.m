clc;
clear;
close all;
M = 10;
snaps_all = [2*M 6*M ];
LNR = 10;
angle = -90:1:90;
format = ['k.-';'b--'];
rand_state = 0;
rng('default');
rng(rand_state);
for i_snaps = 1:length(snaps_all)
    snaps = snaps_all(i_snaps);
    noise = 1/sqrt(2)*(randn(M,snaps)+1i*randn(M,snaps));%高斯复噪声
    snr0 = 0;
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
    Rin = (an1*an1' + an2*an2'+ noise*noise')/snaps;%%干扰和噪声协方差矩阵
    Rn = (noise*noise')/snaps;%%噪声协方差矩阵
    
    p_smi = zeros(1,length(angle));
    p_lsmi = zeros(1,length(angle));
    p_cbf = zeros(1,length(angle));
    
    LR = Rn + 10^(LNR/10) * eye(M);
    
    as_0 = exp(-1i*pi*(0:1:M-1)'*sind(angle0)) ;
    for i_angle = 1:length(angle)
        as = exp(-1i*pi*(0:1:M-1)'*sind(angle(i_angle)));
        w_smi = Rn\as/(as'*(Rn\as));
        p_smi(1, i_angle) =  w_smi' * as_0;
        
        w_lsmi = LR\as/(as'*(LR\as));
        p_lsmi(1, i_angle) =  w_lsmi' * as_0;
        
        w_cbf = as/M;
        p_cbf(1, i_angle) =  w_cbf' * as_0;
    end
    figure(1)
    hold on
    plot(angle ,20*log10(abs(p_smi)), format(i_snaps,:));
  
    
    figure(2)
    plot(angle ,20*log10(abs(p_lsmi)),format(i_snaps,:));
    hold on
end

for i_figure =1:2
    figure(i_figure)
    plot(angle ,20*log10(abs(p_cbf)), 'r-');
    xlabel('方位（^o）')
    ylabel('波束（dB）')
    legend({'N=2M','N=6M','CBF'})
    set(gca,'fontsize',16)
    ylim([-60 0])
    grid on
end
