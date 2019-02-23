%%%% 6.6.3 恒定主瓣响应FIR波束形成器 + 凹槽； 与 图5.13对应
clear all; clc
close all

M = 12;                   % 阵元数
f0 = 2000; c = 1500;
lambda = c/f0;
d = lambda/2;
fU = f0;
fL = f0/2;                 %% 期望频率
d_lambda = 1/2;
angle_radian = 1/180*pi;            %角度转换为弧度
theta_s = 10*angle_radian;
delta0 = 1*angle_radian;

%----------array manifold-------------
a_s  = exp(-1i*2*pi*fL*(0:M-1)'*d*sin(theta_s)/c);   

%-----------期望波束图---------------
theta = -pi/2:delta0:pi/2;
w_cbf = a_s/M;    
ad = exp(-1i*2*pi*fL*(0:M-1)'*d*sin(theta)/c); 
Pd1 = w_cbf'*ad;              %% 期望波束图对应于 f=f0/2;
Pd = 20*log10(abs(Pd1));

figure; hold on
plot(theta/angle_radian,Pd,'b','linewidth',1.5)
axis([-90 90 -60 0])
grid on

%---- 拟合波束 对应于频率 f0/2~f0 -------%
theta_ml = [(-8*angle_radian) : delta0 : (28*angle_radian)];
theta_sl = [-pi/2 : delta0:(-50*angle_radian),(-40*angle_radian) : delta0 : (-12*angle_radian), (32*angle_radian) : delta0 : pi/2];
theta_ac = [(-50*angle_radian) : delta0 : (-40*angle_radian)];
Pml = Pd1(1,find(abs(theta - theta_ml(1,1)) <= 1e-10));
for ii = 2 : length(theta_ml)
    Pml = [Pml Pd1(1,find(abs(theta - theta_ml(1,ii)) <= 1e-10))];
end

FN = 33;
f = linspace(fL,fU,FN);
w_opt = zeros(M,FN);
a = zeros(M,length(theta),FN);
for jj = 1 : FN
    a_ml = exp(-1i*2*pi*f(jj)*(0:M-1)'*d*sin(theta_ml)/c);
    a_sl = exp(-1i*2*pi*f(jj)*(0:M-1)'*d*sin(theta_sl)/c);
    a_ac = exp(-1i*2*pi*f(jj)*(0:M-1)'*d*sin(theta_ac)/c);
    %----- 利用yalmip求解 L2范数 -----%
    y2 = sdpvar(1);
    w2 = sdpvar(M,1,'full','complex');
    Lml = ones(1,length(theta_ml));
    
    % z = 0;
    % for ii = 1 : length(theta_ml)
    %     z = z + (abs( sqrt(Lml(1,ii))*a_ml(:,ii).'*w2 - sqrt(Lml(1,ii))*Pml(1,ii) ))^2;
    % end
    % constraint2_ml = set(z <= y2);
    
    zz = [];
    for ii = 1 : length(theta_ml)
        zz = [zz;sqrt(Lml(1,ii))*a_ml(:,ii).'*w2 - sqrt(Lml(1,ii))*Pml(1,ii)];
    end
    constraint_ml = set(cone(zz,y2));
    
    delta_sl = 10^(-25/20);
    constraint_sl = set(cone(a_sl(:,1).'*w2,delta_sl));
    for ii = 2:length(theta_sl)
        constraint_sl = constraint_sl + set(cone(a_sl(:,ii).'*w2,delta_sl));
    end
    
    delta_ac = 10^(-50/20);
    constraint_ac = set(cone(a_ac(:,1).'*w2,delta_ac));
    for ii = 2:length(theta_ac)
        constraint_ac = constraint_ac + set(cone(a_ac(:,ii).'*w2,delta_ac));
    end
    
    Xi = 10^(-7.5/20);
    constraint_w = set(cone(w2,Xi));
    % Gwd0 = 5;
    % X0 = sqrt((1/M)*10^(Gwd0/20));
    % constraint_w = set(cone(w2,X0));
    
    constraints = constraint_ml + constraint_sl + constraint_ac + constraint_w;
    obj = y2;
    solvesdp(constraints,obj);
    w_opt(:,jj) = double(w2);
    %----- 利用yalmip求解 -----%
    
    a(:,:,jj) = exp(-1i*2*pi*f(jj)*(0:M-1)'*d*sin(theta)/c);
end

P = zeros(FN,pi/delta0+1);
PdB = zeros(FN,pi/delta0+1);
figure; hold on
for jj = 1 : FN
    P(jj,:) = w_opt(:,jj).'*a(:,:,jj);
    PdB(jj,:) = 20*log10(abs(P(jj,:)));
    plot(theta/angle_radian,PdB(jj,:),'b','linewidth',1.0)
    axis([-90 90 -60 0])
    grid on
end
% figure; 
% mesh(theta/angle_radian,f,PdB)
% zlim([-60 0])

%--------- 期望滤波器响应  -------%
tau = (0:M-1)'*d*sin(theta_s)/c;
fs = 3.125*f0; Ts = 1/fs;
L = 64;               % 滤波器长度
T = -round(tau/Ts + (L-1)/2) * Ts;
Hd = conj(w_opt) .* exp(1i*2*pi*T*f);

% fnorm = linspace(fL/fs,fU/fs,FN);
% figure(3);hold on; plot(fnorm,20*log10(abs(Hd(2,:))),'*k'); axis([0 0.5 -60 0])
% figure(4);hold on; plot(fnorm, angle(Hd(2,:))*180/pi,'*k'); axis([0 0.5 -200 200])
%----------- 期望滤波器响应  -------%

%----------- 设计滤波器  -------%
deltaf = (fU-fL)/(FN-1)/fs;                     % 离散化频率间隔
PB = 0.16 : deltaf : 0.32;    % 通带带宽 其中 0.16 = fL/fs；0.32 = fU/fs；
SB = [    0:deltaf:0.13   0.35:deltaf:0.50 ];      % 阻带带宽
TB = [ 0.14:deltaf:0.15   0.33:deltaf:0.34 ];      % 过渡带宽
WB =  0 : deltaf : 0.5;                            % 全部带宽
e_PB = exp(-1i*2*pi*(0:L-1)'*PB);  
e_SB = exp(-1i*2*pi*(0:L-1)'*SB);  
e_TB = exp(-1i*2*pi*(0:L-1)'*TB); 
e_WB = exp(-1i*2*pi*(0:L-1)'*WB);

H = zeros(M,FN);
for m = 1 : M
    %----- 阻带峰值约束最小均方通带方法 -----%
    y2 = sdpvar(1);
    h2 = sdpvar(L,1,'full','real');
    LPB = ones(1,length(PB));
    
    zz = [];
    for ii = 1 : length(PB)
        zz = [zz;sqrt(LPB(1,ii))*e_PB(:,ii).'*h2 - sqrt(LPB(1,ii))*Hd(m,ii)];
    end
    constraint2 = set(cone(zz,y2));
    
    delta_SB = 10^(-40/20);
    constraint_SB = set(cone(e_SB(:,1).'*h2,delta_SB));
    for ii = 2:length(SB)
        constraint_SB = constraint_SB + set(cone(e_SB(:,ii).'*h2,delta_SB));
    end
    
    constraints = constraint2 + constraint_SB;
    obj = y2;
    solvesdp(constraints,obj);
    h_2 = double(h2);
    %----- 利用yalmip求解 -----%
    % fB = [ 0:deltaf/fs:0.5 ];
    % e = exp(-1i*2*pi*(0:L-1)'*fB);
    % H2 = e.'*h_2;
    % figure(3); plot(fB, 20*log10(abs(H2)),'k');
    % figure(4); plot(fB, angle(H2)*180/pi,'k');
    %----- 阻带峰值约束最小均方通带方法 -----%
    
%     %----- 通带均方误差约束最低最带方法 -----%
%     %----- 利用yalmip求解 L2范数 -----%
%     y3 = sdpvar(1);
%     h3 = sdpvar(L,1,'full','real');
%     LPB = ones(1,length(PB));
%     
%     zz = [];
%     for ii = 1 : length(PB)
%         zz = [zz;sqrt(LPB(1,ii))*e_PB(:,ii).'*h3 - sqrt(LPB(1,ii))*Hd(1,ii)];
%     end
%     delta_PB = 10^(-60/20);
%     constraint2 = set(cone(zz,delta_PB));
%     
%     constraint_SB = set(cone(e_SB(:,1).'*h3,y3));
%     for ii = 2:length(SB)
%         constraint_SB = constraint_SB + set(cone(e_SB(:,ii).'*h3,y3));
%     end
%     
%     constraints = constraint2 + constraint_SB;
%     obj = y3;
%     solvesdp(constraints,obj);
%     h_3 = double(h3);
%     %----- 利用yalmip求解 -----%
%     fB = [ 0:deltaf/fs:0.5 ];
%     e = exp(-1i*2*pi*(0:L-1)'*fB);
%     H3 = e.'*h_3;
%     figure(1); plot(fB, 20*log10(abs(H3)),'k');
%     figure(2); plot(fB, angle(H3)*180/pi,'k');
%     %----- 通带均方误差约束最低最带方法 -----%
    
    H(m,:) = h_2.'*e_PB;
    % figure(3); plot(PB, 20*log10(abs(H(m,:))),'r');
    % figure(4); plot(PB, angle(H(m,:))*180/pi,'r');
end

w_fir = H .* exp(-1i*2*pi*T*f);
w_fir = conj(w_fir);

P_fir = zeros(FN,pi/delta0+1);
PdB_fir = zeros(FN,pi/delta0+1);
figure; hold on
for jj = 1 : FN
    P_fir(jj,:) = w_fir(:,jj).'*a(:,:,jj);
    PdB_fir(jj,:) = 20*log10(abs(P_fir(jj,:)));
    plot(theta/angle_radian,PdB_fir(jj,:),'b','linewidth',1.0)
    axis([-90 90 -60 0])
    grid on
end
% figure; 
% mesh(theta/angle_radian,f,PdB_fir)
% zlim([-60 0])
