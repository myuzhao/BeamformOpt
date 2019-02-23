%%%% 例7.1 全局优化法
clear all; clc
% close all

M = 12;                   % 阵元数
f0 = 2000; c = 1500;
fL = f0/2; fU = f0;              %% 期望频率
fs = 3.125*f0; Ts = 1/fs;
lambda = c/f0;
d = lambda/2;
d_lambda = 1/2;
angle_radian = 1/180*pi;            %角度转换为弧度
theta_s = 10*angle_radian;
delta0 = 2*angle_radian;

%----------array manifold-------------
a_s  = exp(-1i*2*pi*fL*(0:M-1)'*d*sin(theta_s)/c);   

%-----------期望波束图---------------
theta = -pi/2:delta0:pi/2;
w_cbf = a_s/M;    
ad = exp(-1i*2*pi*fL*(0:M-1)'*d*sin(theta)/c); 
Pd = w_cbf'*ad;              %% 期望波束图对应于 f=f0/2;

figure; plot(theta/angle_radian,20*log10(abs(Pd)),'b','linewidth',1.5)
axis([-90 90 -60 0]); grid on

%---- 拟合波束 对应于频率 f0/2~f0 -------%
theta_wl = -pi/2:delta0:pi/2;
theta_ml = [(-8*angle_radian) : delta0 : (28*angle_radian)];
theta_sl = [-pi/2:delta0:(-12*angle_radian), (32*angle_radian):delta0:pi/2];
Pml = Pd(1,find(abs(theta - theta_ml(1,1)) <= 1e-10));
for ii = 2 : length(theta_ml)
    Pml = [Pml Pd(1,find(abs(theta - theta_ml(1,ii)) <= 1e-10))];
end

tau_s = (0:M-1)'*d*sin(theta_s)/c;
L = 64;               % 滤波器长度
T = -round(tau_s/Ts + (L-1)/2) * Ts;
tau_s  = (0:M-1)'*d*sin(theta_s)/c;
tau_ml = (0:M-1)'*d*sin(theta_ml)/c;
tau_sl = (0:M-1)'*d*sin(theta_sl)/c;
tau_wl = (0:M-1)'*d*sin(theta_wl)/c;

%----- 式(7.9) -----%
deltaf = 0.01;                     % 离散化频率间隔
PB = 0.16 : deltaf : 0.32;      % 通带带宽 其中 0.16 = fL/fs；0.32 = fU/fs；
SB = [    0:deltaf:0.13   0.35:deltaf:0.50 ];      % 阻带带宽
TB = [ 0.14:deltaf:0.15   0.33:deltaf:0.34 ];      % 过渡带宽
WB =  0 : deltaf : 0.5;                            % 全部带宽
e_PB = exp(-1i*2*pi*(0:L-1)'*PB);  
e_SB = exp(-1i*2*pi*(0:L-1)'*SB);  
e_TB = exp(-1i*2*pi*(0:L-1)'*TB); 
e_WB = exp(-1i*2*pi*(0:L-1)'*WB);

K_PB = exp(-1i*2*pi*T*PB*fs);
K_SB = exp(-1i*2*pi*T*SB*fs);
K_TB = exp(-1i*2*pi*T*TB*fs);
K_WB = exp(-1i*2*pi*T*WB*fs);

%----- 利用yalmip求解 L2范数 -----%
y = sdpvar(1);
h = sdpvar(M*L,1,'full','real');

%----- PB + ML -------  %   
for pb = 1 : length(PB)
    for ml = 1 : length(theta_ml)
        a_PBML = exp(-1i*2*pi*PB(1,pb)*fs*(0:M-1)'*d*sin(theta_ml(1,ml))/c);
        u_PBML = kron(e_PB(:,pb),conj(a_PBML).*K_PB(:,pb));
        if pb ==1 && ml ==1
            constraint_PBML = set( cone( (u_PBML.'*h-Pml(1,1)),y ) );
        else
            constraint_PBML = constraint_PBML + set( cone( (u_PBML.'*h-Pml(1,ml)),y ) );
        end
    end
end

%----- PB + SL -------  %  
delta_PBSL = 10^(-25/20);
for pb = 1 : length(PB)
    for sl = 1: length(theta_sl)
        a_PBSL = exp(-1i*2*pi*PB(1,pb)*fs*(0:M-1)'*d*sin(theta_sl(1,sl))/c);
        u_PBSL = kron(e_PB(:,pb),conj(a_PBSL).*K_PB(:,pb));
        if pb ==1 && sl ==1
            constraint_PBSL = set(cone(u_PBSL.'*h,delta_PBSL));
        else 
            constraint_PBSL = constraint_PBSL + set(cone(u_PBSL.'*h,delta_PBSL));
        end
    end
end

%----- TB + SL -------  %  
delta_TBSL = 10^(-25/20);
for tb = 1 : length(TB)
    for sl = 1: length(theta_sl)
        a_TBSL = exp(-1i*2*pi*TB(1,tb)*fs*(0:M-1)'*d*sin(theta_sl(1,sl))/c);
        u_TBSL = kron(e_TB(:,tb),conj(a_TBSL).*K_TB(:,tb));
        if tb ==1 && sl ==1
            constraint_TBSL = set(cone(u_TBSL.'*h,delta_TBSL));
        else
            constraint_TBSL = constraint_TBSL + set(cone(u_TBSL.'*h,delta_TBSL));
        end
    end
end

%----- SB + WL -------  %  
delta_SBWL = 10^(-25/20);
for sb = 1 : length(SB)
    for wl = 1: length(theta_wl)
        a_SBWL = exp(-1i*2*pi*SB(1,sb)*fs*(0:M-1)'*d*sin(theta_wl(1,wl))/c);
        u_SBWL = kron(e_WB(:,sb),conj(a_SBWL).*K_SB(:,sb));
        if sb ==1 && wl ==1
            constraint_SBWL = set(cone(u_SBWL.'*h,delta_SBWL));
        else
            constraint_SBWL = constraint_SBWL + set(cone(u_SBWL.'*h,delta_SBWL));
        end
    end
end

% constraint_h = set(cone(h,0.25));

constraints = constraint_PBML + constraint_PBSL;% + constraint_TBSL + constraint_SBWL + constraint_h;
obj = y; 
solvesdp(constraints,obj);
h_opt = double(h);
%----- 利用yalmip求解 -----%


%----- 式(7.10) -----%
p_fir = zeros(length(PB),length(theta_wl));
pdB_fir = zeros(length(PB),length(theta_wl));
figure; hold on
for wb = 1 : length(PB)
    for wl = 1: length(theta_wl)
        U = zeros(M,L);
        for m = 1 : M
            for LL = 1 : L
                U(m,LL) = exp(-1i*2*pi*PB(1,wb)*fs*(T(m,1)+tau_wl(m,wl)+(LL-1)*Ts));
            end
        end
        u = reshape(U,numel(U),1);  
        p_fir(wb,wl) = u.' * h_opt;
        pdB_fir(wb,wl) = 20*log10(abs(p_fir(wb,wl)));
    end
    plot(theta/angle_radian,pdB_fir(wb,:),'b','linewidth',1.0)
%     axis([-90 90 -60 0])
    grid on        
end
figure; 
mesh(theta/angle_radian,PB,pdB_fir)
% zlim([-60 0])           
 
