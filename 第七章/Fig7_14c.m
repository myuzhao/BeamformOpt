%%%% 例7.5 多干扰情况FIR波束形成器
clear all;
clc
close all

M = 12;                   % 阵元数
f0 = 2000; c = 1500;
fL = f0/2; fU = f0;              %% 期望频率
fs = 3.125*f0; Ts = 1/fs;
lambda = c/f0;
d = lambda/2;
d_lambda = 1/2;
angle_radian = 1/180*pi;            %角度转换为弧度
delta0 = 2*angle_radian;
deltaf = 0.005;                     % 离散化频率间隔
PB = 0.16 : deltaf : 0.32;      % 通带带宽 其中 0.16 = fL/fs；0.32 = fU/fs；

L = 31;               % 滤波器长度
theta_s = 10*angle_radian;
tau_s  = (0:M-1)'*d*sin(theta_s)/c;
T_real = -round(tau_s/Ts + (L-1)/2) * Ts;   % 滤波器系数为实数
T_comp = -round(tau_s/Ts) * Ts;             % 滤波器系数为复数
% theta_sl = [-pi/2:delta0:(-12*angle_radian), (32*angle_radian):delta0:pi/2];
% tau_sl = (0:M-1)'*d*sin(theta_sl)/c;
angle_sl = [ linspace(-15,-5,length(PB)).' linspace(35,25,length(PB)).' ];
angle_sl = angle_sl * angle_radian;

%--------- 生成宽带协方差矩阵 ----------%
B = fU - fL; fc = (fU + fL)/2;

%--------- 信号 -----------------------%
SNRs = 0;
Rs_real = zeros(M*L,M*L);
Rs_comp = zeros(M*L,M*L);
for ma = 1 : M
    for la = 1 : L
        tau_s_real_mala = T_real(ma,1)+tau_s(ma,1)+(la-1)*Ts;
        tau_s_comp_mala = T_comp(ma,1)+tau_s(ma,1)+(la-1)*Ts;
        for mb = 1 : M
            for lb = 1 : L
                tau_s_real_mblb = T_real(mb,1)+tau_s(mb,1)+(lb-1)*Ts;
                tau_s_comp_mblb = T_comp(mb,1)+tau_s(mb,1)+(lb-1)*Ts;
                row = ma + (la - 1) * M; col = mb + (lb - 1) * M;
                Rs_real(row,col) = 1 / B * sinc(pi*B*(tau_s_real_mala-tau_s_real_mblb)) * cos(2*pi*fc*(tau_s_real_mala-tau_s_real_mblb));
                Rs_comp(row,col) = 1 / B * sinc(pi*B*(tau_s_comp_mala-tau_s_comp_mblb)) * exp(-1i*2*pi*fc*(tau_s_comp_mala-tau_s_comp_mblb));
            end
        end
    end
end

%--------------- 干扰 1 -----------------%
SNRd1 = 30;
theta_d1 = -50*angle_radian;
tau_d1  = (0:M-1)'*d*sin(theta_d1)/c;

Rd1_real = zeros(M*L,M*L);
Rd1_comp = zeros(M*L,M*L);
for ma = 1 : M
    for la = 1 : L
        tau_d1_real_mala = T_real(ma,1)+tau_d1(ma,1)+(la-1)*Ts;
        tau_d1_comp_mala = T_comp(ma,1)+tau_d1(ma,1)+(la-1)*Ts;
        for mb = 1 : M
            for lb = 1 : L
                tau_d1_real_mblb = T_real(mb,1)+tau_d1(mb,1)+(lb-1)*Ts;
                tau_d1_comp_mblb = T_comp(mb,1)+tau_d1(mb,1)+(lb-1)*Ts;
                row = ma + (la - 1) * M; col = mb + (lb - 1) * M;
                Rd1_real(row,col) = 1 / B * sinc(pi*B*(tau_d1_real_mala-tau_d1_real_mblb)) * cos(2*pi*fc*(tau_d1_real_mala-tau_d1_real_mblb));
                Rd1_comp(row,col) = 1 / B * sinc(pi*B*(tau_d1_comp_mala-tau_d1_comp_mblb)) * exp(-1i*2*pi*fc*(tau_d1_comp_mala-tau_d1_comp_mblb));
            end
        end
    end
end

%--------- 干扰 2 -----------------------%
SNRd2 = 30;
theta_d2 = -30*angle_radian;
tau_d2  = (0:M-1)'*d*sin(theta_d2)/c;

Rd2_real = zeros(M*L,M*L);
Rd2_comp = zeros(M*L,M*L);
for ma = 1 : M
    for la = 1 : L
        tau_d2_real_mala = T_real(ma,1)+tau_d2(ma,1)+(la-1)*Ts;
        tau_d2_comp_mala = T_comp(ma,1)+tau_d2(ma,1)+(la-1)*Ts;
        for mb = 1 : M
            for lb = 1 : L
                tau_d2_real_mblb = T_real(mb,1)+tau_d2(mb,1)+(lb-1)*Ts;
                tau_d2_comp_mblb = T_comp(mb,1)+tau_d2(mb,1)+(lb-1)*Ts;
                row = ma + (la - 1) * M; col = mb + (lb - 1) * M;
                Rd2_real(row,col) = 1 / B * sinc(pi*B*(tau_d2_real_mala-tau_d2_real_mblb)) * cos(2*pi*fc*(tau_d2_real_mala-tau_d2_real_mblb));
                Rd2_comp(row,col) = 1 / B * sinc(pi*B*(tau_d2_comp_mala-tau_d2_comp_mblb)) * exp(-1i*2*pi*fc*(tau_d2_comp_mala-tau_d2_comp_mblb));
            end
        end
    end
end

%--------- 噪声 -----------------------%
SNRn = 0;
Rn_real = zeros(M*L,M*L);
Rn_comp = zeros(M*L,M*L);
for ma = 1 : M
    for la = 1 : L
        for mb = 1 : M
            for lb = 1 : L
                row = ma + (la - 1) * M; col = mb + (lb - 1) * M;
                Rn_real(row,col) = 1 / B * (ma == mb) * sinc(pi*B*(la-lb)*Ts) * cos(2*pi*fc*(la-lb)*Ts);
                Rn_comp(row,col) = 1 / B * (ma == mb) * sinc(pi*B*(la-lb)*Ts) * exp(-1i*2*pi*fc*(la-lb)*Ts);
            end
        end
    end
end
R_real = 10^(SNRs/10) * Rs_real + 10^(SNRd1/10) * Rd1_real + 10^(SNRd2/10) * Rd2_real + 10^(SNRn/10) * Rn_real;
R_comp = 10^(SNRs/10) * Rs_comp + 10^(SNRd1/10) * Rd1_comp + 10^(SNRd2/10) * Rd2_comp + 10^(SNRn/10) * Rn_comp;

%------滤波器系数为实数-----------%
%----- 利用yalmip求解 L2范数 -----%
y_real = sdpvar(1);
h_real = sdpvar(M*L,1,'full','real');

constraint_real_R = set( h_real'*R_real*h_real <= y_real );

%----- PB + theta_s 无失真约束 -------  %   
for pb = 1 : length(PB)
	U_real_PBS = exp(-1i*2*pi*PB(1,pb)*fs*((T_real+tau_s)*ones(1,L)+ones(M,1)*(0:L-1)*Ts));
	u_real_PBS(:,pb) = vec(U_real_PBS);
end
cvx_begin 
	variable h_real(M*L);
    minimize(h_real'*R_real*h_real)
    subject to
        u_real_PBS'*h_real==1;
        norm(h_real)<=0.25;
cvx_end
% %----- PB + theta_sl 旁瓣控制 -------  %   
% delta_PBSL = 10^(-30/20);
% for pb = 1 : length(PB)
%     theta_sl = [ -pi/2:delta0:angle_sl(pb,1), angle_sl(pb,2):delta0:pi/2 ];
%     tau_sl = (0:M-1)'*d*sin(theta_sl)/c;
%     for sl = 1: length(theta_sl)
%         U_real_PBSL = zeros(M,L);
%         for m = 1 : M
%             for LL = 1 : L
%                 U_real_PBSL(m,LL) = exp(-1i*2*pi*PB(1,pb)*fs*(T_real(m,1)+tau_sl(m,sl)+(LL-1)*Ts));
%             end
%         end
%         u_real_PBSL = reshape(U_real_PBSL,numel(U_real_PBSL),1);
%         if pb ==1 && sl ==1
%             constraint_real_PBSL = set( cone( u_real_PBSL.'*h_real, delta_PBSL ) );
%         else
%             constraint_real_PBSL = constraint_real_PBSL + set( cone( u_real_PBSL.'*h_real, delta_PBSL ) );
%         end
%     end
% end


%----- 式(7.10) -----%
theta_wl = -pi/2:delta0:pi/2;
tau_wl = (0:M-1)'*d*sin(theta_wl)/c;
p_real_PB = zeros(length(PB),length(theta_wl));
pdB_real_PB = zeros(length(PB),length(theta_wl));
figure; hold on
for pb = 1 : length(PB)
    for wl = 1: length(theta_wl)
        U_real = exp(-1i*2*pi*PB(1,pb)*fs*((T_real+tau_wl(:,wl))*ones(1,L)+ones(M,1)*(0:L-1)*Ts));
        u_real = vec(U_real);  
        p_real_PB(pb,wl) = u_real.' * h_real;
        pdB_real_PB(pb,wl) = 20*log10(abs(p_real_PB(pb,wl)));
    end
    plot(theta_wl/angle_radian,pdB_real_PB(pb,:),'b','linewidth',1.0)
    axis([-90 90 -100 5])
    grid on        
end
figure; 
mesh(theta_wl/angle_radian,PB*fs/10,pdB_real_PB)
zlim([-100 5])           
 
% %------滤波器系数为复数-----------%
% %----- 利用yalmip求解 L2范数 -----%
% y_comp = sdpvar(1);
% h_comp = sdpvar(M*L,1,'full','complex');
% 
% constraint_comp_R = set( h_comp.'*R_comp*conj(h_comp) <= y_comp );
% 
% %----- PB + theta_s 无失真约束 -------  %   
% for pb = 1 : length(PB)
%     for s = 1: length(theta_s)
%         U_comp_PBS = zeros(M,L);
%         for m = 1 : M
%             for LL = 1 : L
%                 U_comp_PBS(m,LL) = exp(-1i*2*pi*PB(1,pb)*fs*(T_comp(m,1)+tau_s(m,s)+(LL-1)*Ts));
%             end
%         end
%         u_comp_PBS = reshape(U_comp_PBS,numel(U_comp_PBS),1);
%         if pb ==1 && s ==1
%             constraint_comp_PBS = set( u_comp_PBS.'*h_comp == 1  );
%         else
%             constraint_comp_PBS = constraint_comp_PBS + set( u_comp_PBS.'*h_comp == 1 );
%         end
%     end
% end
% 
% %----- PB + theta_sl 旁瓣控制 -------  %   
% % angle_sl = [ linspace(-15,-5,length(PB)).' linspace(35,25,length(PB)).' ];
% % angle_sl = angle_sl * angle_radian;
% delta_PBSL = 10^(-30/20);
% for pb = 1 : length(PB)
%     theta_sl = [ -pi/2:delta0:angle_sl(pb,1), angle_sl(pb,2):delta0:pi/2 ];
%     tau_sl = (0:M-1)'*d*sin(theta_sl)/c;
%     for sl = 1: length(theta_sl)
%         U_comp_PBSL = zeros(M,L);
%         for m = 1 : M
%             for LL = 1 : L
%                 U_comp_PBSL(m,LL) = exp(-1i*2*pi*PB(1,pb)*fs*(T_comp(m,1)+tau_sl(m,sl)+(LL-1)*Ts));
%             end
%         end
%         u_comp_PBSL = reshape(U_comp_PBSL,numel(U_comp_PBSL),1);
%         if pb ==1 && sl ==1
%             constraint_comp_PBSL = set( cone( u_comp_PBSL.'*h_comp, delta_PBSL ) );
%         else
%             constraint_comp_PBSL = constraint_comp_PBSL + set( cone( u_comp_PBSL.'*h_comp, delta_PBSL ) );
%         end
%     end
% end
% 
% % -------- 宽带白噪声增益约束 ---------%
% constraint_comp_h = set( cone(h_comp,0.25) );
% 
% constraints = constraint_comp_R + constraint_comp_PBS + constraint_comp_PBSL + constraint_comp_h; 
% obj = y_comp; 
% solvesdp(constraints,obj);
% h_comp_opt = double(h_comp);
% %----- 利用yalmip求解 -----%
% 
% 
% %----- 式(7.10) -----%
% theta_wl = -pi/2:delta0:pi/2;
% tau_wl = (0:M-1)'*d*sin(theta_wl)/c;
% p_comp_PB = zeros(length(PB),length(theta_wl));
% pdB_comp_PB = zeros(length(PB),length(theta_wl));
% figure; hold on
% for pb = 1 : length(PB)
%     for wl = 1: length(theta_wl)
%         U_comp = zeros(M,L);
%         for m = 1 : M
%             for LL = 1 : L
%                 U_comp(m,LL) = exp(-1i*2*pi*PB(1,pb)*fs*(T_comp(m,1)+tau_wl(m,wl)+(LL-1)*Ts));
%             end
%         end
%         u_comp = reshape(U_comp,numel(U_comp),1);  
%         p_comp_PB(pb,wl) = u_comp.' * h_comp_opt;
%         pdB_comp_PB(pb,wl) = 20*log10(abs(p_comp_PB(pb,wl)));
%     end
%     plot(theta_wl/angle_radian,pdB_comp_PB(pb,:),'b','linewidth',1.0)
% %     axis([-90 90 -100 5])
%     grid on        
% end
% % figure; 
% % mesh(theta_wl/angle_radian,PB,pdB_comp_PB)
% % zlim([-100 5])  