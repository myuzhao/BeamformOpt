%%%% 例7.9 多约束恒定主瓣响应波束设计
clear all; clc
close all

M = 12;                   % 阵元数
f0 = 2000; c = 1500;
fL = f0/2; fU = f0;              % 期望频率
fM = (fL + fU)/2;                % 参考频率
fs = 3.125*f0; Ts = 1/fs;
lambda = c/f0;
d = lambda/2;
d_lambda = 1/2;
angle_radian = 1/180*pi;            % 角度转换为弧度
delta0 = 2*angle_radian;

%---- 拟合波束 对应于频率 f0/2~f0 -------%
theta_wl = -pi/2:delta0:pi/2;
theta_ml = [(-8*angle_radian) : delta0 : (28*angle_radian)];
theta_sl = [-pi/2:delta0:(-12*angle_radian), (32*angle_radian):delta0:pi/2];

L = 32;               % 滤波器长度
theta_s = 10*angle_radian;
tau_s  = (0:M-1)'*d*sin(theta_s)/c;
tau_ml = (0:M-1)'*d*sin(theta_ml)/c;
tau_sl = (0:M-1)'*d*sin(theta_sl)/c;
tau_wl = (0:M-1)'*d*sin(theta_wl)/c;
T = -round(tau_s/Ts + (L-1)/2) * Ts;

%----- 式(7.10) -----%
deltaf = 0.005;                     % 离散化频率间隔
PB = 0.16 : deltaf : 0.32;      % 通带带宽 其中 0.16 = fL/fs；0.32 = fU/fs；

%--------- 生成宽带协方差矩阵 ----------%
B = fU - fL; fc = (fU + fL)/2;

%--------- 信号 -----------------------%
SNRs = 0;
Rs = zeros(M*L,M*L);
for ma = 1 : M
    for la = 1 : L
        tau_s_mala = T(ma,1)+tau_s(ma,1)+(la-1)*Ts;
        for mb = 1 : M
            for lb = 1 : L
                tau_s_mblb = T(mb,1)+tau_s(mb,1)+(lb-1)*Ts;
                row = ma + (la - 1) * M; col = mb + (lb - 1) * M;
                Rs(row,col) = 1 / B * sinc(pi*B*(tau_s_mala-tau_s_mblb)) * cos(2*pi*fc*(tau_s_mala-tau_s_mblb));
            end
        end
    end
end

%--------- 干扰 1 -----------------------%
SNRd1 = 30;
theta_d1 = -50*angle_radian;
tau_d1  = (0:M-1)'*d*sin(theta_d1)/c;

Rd1 = zeros(M*L,M*L);
for ma = 1 : M
    for la = 1 : L
        tau_d1_mala = T(ma,1)+tau_d1(ma,1)+(la-1)*Ts;
        for mb = 1 : M
            for lb = 1 : L
                tau_d1_mblb = T(mb,1)+tau_d1(mb,1)+(lb-1)*Ts;
                row = ma + (la - 1) * M; col = mb + (lb - 1) * M;
                Rd1(row,col) = 1 / B * sinc(pi*B*(tau_d1_mala-tau_d1_mblb)) * cos(2*pi*fc*(tau_d1_mala-tau_d1_mblb));
            end
        end
    end
end

%--------- 干扰 2 -----------------------%
SNRd2 = 30;
theta_d2 = -30*angle_radian;
tau_d2  = (0:M-1)'*d*sin(theta_d2)/c;

Rd2 = zeros(M*L,M*L);
for ma = 1 : M
    for la = 1 : L
        tau_d2_mala = T(ma,1)+tau_d2(ma,1)+(la-1)*Ts;
        for mb = 1 : M
            for lb = 1 : L
                tau_d2_mblb = T(mb,1)+tau_d2(mb,1)+(lb-1)*Ts;
                row = ma + (la - 1) * M; col = mb + (lb - 1) * M;
                Rd2(row,col) = 1 / B * sinc(pi*B*(tau_d2_mala-tau_d2_mblb)) * cos(2*pi*fc*(tau_d2_mala-tau_d2_mblb));
            end
        end
    end
end

%--------- 噪声 -----------------------%
SNRn = 0;
Rn = zeros(M*L,M*L);
for ma = 1 : M
    for la = 1 : L
        for mb = 1 : M
            for lb = 1 : L
                row = ma + (la - 1) * M; col = mb + (lb - 1) * M;
                Rn(row,col) = 1 / B * (ma == mb) * sinc(pi*B*(la-lb)*Ts) * cos(2*pi*fc*(la-lb)*Ts);
            end
        end
    end
end
R = 10^(SNRs/10) * Rs + 10^(SNRd1/10) * Rd1 + 10^(SNRd2/10) * Rd2 + 10^(SNRn/10) * Rn;


%----- 利用yalmip求解 L2范数 -----%
y = sdpvar(1);
h = sdpvar(M*L,1,'full','real');

%----- 最大增益约束 -------  %
constraint_R = [( h' * R * h <= y )];

%----- 无失真约束 -------  %
U_S = zeros(M,L);
for m = 1 : M
    for LL = 1 : L
        U_S(m,LL)    = exp(-1i*2*pi*fM*(T(m,1)+tau_s(m,1)+(LL-1)*Ts));
    end
end
u_S = reshape(U_S,numel(U_S),1);
constraint_S = [(  u_S.'*h == 1 )];

% % %------- 对每个频率都施加 无失真约束 -------%
% for pb = 1 : length(PB)
%     for s = 1: length(theta_s)
%         U_S = zeros(M,L);
%         for m = 1 : M
%             for LL = 1 : L
%                 U_S(m,LL) = exp(-1i*2*pi*PB(1,pb)*fs*(T(m,1)+tau_s(m,s)+(LL-1)*Ts));
%             end
%         end
%         u_S = reshape(U_S,numel(U_S),1);
%         if pb ==1 && s ==1
%             constraint_S = set( u_S.'*h == 1 );
%         else
%             constraint_S = constraint_S + set( u_S.'*h == 1 );
%         end
%     end
% end

%----- PB + ML 最小差异恒定主瓣：均方主瓣约束 -------  %
delta_PBML = 10^(-50.0/10);
UU_PBML = zeros(M*L,M*L);
for pb = 1 : length(PB)
    for ml = 1: length(theta_ml)
        U_PBML = zeros(M,L);
        U_D    = zeros(M,L);
        for m = 1 : M
            for LL = 1 : L
                U_PBML(m,LL) = exp(-1i*2*pi*PB(1,pb)*fs*(T(m,1)+tau_ml(m,ml)+(LL-1)*Ts));
                U_D(m,LL)    = exp(-1i*2*pi*fM*(T(m,1)+tau_ml(m,ml)+(LL-1)*Ts));
            end
        end
        u_PBML = reshape(U_PBML,numel(U_PBML),1);
        u_D    = reshape(U_D,numel(U_D),1);
        UU_PBML = UU_PBML + real( (u_PBML - u_D) * (u_PBML - u_D)' );
    end
end
UU_PBML = UU_PBML / (length(PB)*length(theta_ml));
constraint_PBML =[( h' * UU_PBML * h <= delta_PBML )];    

% %----- PB + SL 均方旁瓣约束-------  %  
% delta_PBSL = 10^(-33.0/10);
% UU_PBSL = zeros(M*L,M*L);
% for pb = 1 : length(PB)
%     for sl = 1: length(theta_sl)
%         U_PBSL = zeros(M,L);
%         for m = 1 : M
%             for LL = 1 : L
%                 U_PBSL(m,LL) = exp(-1i*2*pi*PB(1,pb)*fs*(T(m,1)+tau_sl(m,sl)+(LL-1)*Ts));
%             end
%         end
%         u_PBSL = reshape(U_PBSL,numel(U_PBSL),1);
%         UU_PBSL = UU_PBSL + real( u_PBSL * u_PBSL' );
%     end
% end
% UU_PBSL = UU_PBSL / (length(PB)*length(theta_sl));
% constraint_PBSL = set( h' * UU_PBSL * h <= delta_PBSL );  

%----- PB + SL 旁瓣峰值约束 -------  %  
delta_PBSL = 10^(-25.0/20);
UU_PBSL = zeros(M*L,M*L);
for pb = 1 : length(PB)
    for sl = 1: length(theta_sl)
        U_PBSL = zeros(M,L);
        for m = 1 : M
            for LL = 1 : L
                U_PBSL(m,LL) = exp(-1i*2*pi*PB(1,pb)*fs*(T(m,1)+tau_sl(m,sl)+(LL-1)*Ts));
            end
        end
        u_PBSL = reshape(U_PBSL,numel(U_PBSL),1);
        if pb ==1 && sl ==1
            constraint_PBSL = [( cone( u_PBSL.'*h, delta_PBSL ) )];
        else
            constraint_PBSL = constraint_PBSL + ( cone( u_PBSL.'*h, delta_PBSL ) );
        end
    end
end

%------- 稳健性约束 ------%
X = 10^(-10/20);
constraint_h = [(cone(h,X))];

constraints = constraint_R + constraint_S + constraint_PBML + constraint_PBSL + constraint_h;% 
obj = y; 
solvesdp(constraints,obj);
h_opt = double(h);
y_opt = 10*log10(double(y));
%----- 利用yalmip求解 -----%

%----- 式(7.10) -----%
p_PB = zeros(length(PB),length(theta_wl));
pdB_PB = zeros(length(PB),length(theta_wl));
figure; hold on
for pb = 1 : length(PB)
    for wl = 1: length(theta_wl)
        U = zeros(M,L);
        for m = 1 : M
            for LL = 1 : L
                U(m,LL) = exp(-1i*2*pi*PB(1,pb)*fs*(T(m,1)+tau_wl(m,wl)+(LL-1)*Ts));
            end
        end
        u = reshape(U,numel(U),1);  
        p_PB(pb,wl) = u.' * h_opt;
        pdB_PB(pb,wl) = 20*log10(abs(p_PB(pb,wl)));
    end
    plot(theta_wl/angle_radian,pdB_PB(pb,:),'b','linewidth',1.0)
    axis([-90 90 -60 0])
    grid on        
end
% figure; 
% mesh(theta_wl/angle_radian,PB,pdB_PB)
% zlim([-60 0])           
  
