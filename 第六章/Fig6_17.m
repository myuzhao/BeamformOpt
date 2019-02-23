%%%% 6.6.3 恒定主瓣响应FIR波束形成器
clear all; clc
close all
jay=sqrt(-1);

M = 12;                   % 阵元数
f0 = 2000; c = 1500;
fL = f0/2; fU = f0;              %% 期望频率
lambda = c/f0;
d = lambda/2;
d_lambda = 1/2;
angle_radian = 1/180*pi;            %角度转换为弧度
theta_s = 10*angle_radian;
delta0 = 2*angle_radian;

%----------array manifold-------------
a_s  = exp(-jay*2*pi*fL*(0:M-1)'*d*sin(theta_s)/c);   

%-----------期望波束图---------------
theta = -pi/2:delta0:pi/2;
w_cbf = a_s/M;    
ad = exp(-jay*2*pi*fL*(0:M-1)'*d*sin(theta)/c); 
Pd1 = w_cbf'*ad;              %% 期望波束图对应于 f=f0/2;
Pd = 20*log10(abs(Pd1));

figure; hold on
plot(theta/angle_radian,Pd,'b','linewidth',1.5)
axis([-90 90 -60 0])
grid on

%---- 拟合波束 对应于频率 f0/2~f0 -------%
theta_ml = [(-8*angle_radian) : delta0 : (28*angle_radian)];
theta_sl = [-pi/2:delta0:(-12*angle_radian), (32*angle_radian):delta0:pi/2];
Pml = Pd1(1,find(abs(theta - theta_ml(1,1)) <= 1e-10));
for ii = 2 : length(theta_ml)
    Pml = [Pml Pd1(1,find(abs(theta - theta_ml(1,ii)) <= 1e-10))];
end

FN = 33;
f = linspace(fL,fU,FN);
w_opt = zeros(M,FN);
a = zeros(M,length(theta),FN);
for jj = 1 : FN
    a_ml = exp(-jay*2*pi*f(jj)*(0:M-1)'*d*sin(theta_ml)/c);
    a_sl = exp(-jay*2*pi*f(jj)*(0:M-1)'*d*sin(theta_sl)/c);
    
	cvx_begin quiet
        variable w2(M) complex;
            minimize(norm(w2'*a_ml-Pml))
        subject to
            max(abs(w2'*a_sl))<=10^(-25/20);
            norm(w2)<=10^(-7.5/20);
    cvx_end
    w_opt(:,jj) = w2;
 
    a(:,:,jj) = exp(-jay*2*pi*f(jj)*(0:M-1)'*d*sin(theta)/c);
end

P = zeros(FN,pi/delta0+1);
PdB = zeros(FN,pi/delta0+1);
figure; hold on
for jj = 1 : FN
    P(jj,:) = w_opt(:,jj)'*a(:,:,jj);
    PdB(jj,:) = 20*log10(abs(P(jj,:)));
    plot(theta/angle_radian,PdB(jj,:),'b','linewidth',1.0)
    axis([-90 90 -60 0])
    grid on
end
figure; 
mesh(theta/angle_radian,f,PdB)
zlim([-60 0])

%--------- 期望滤波器响应  -------%
tau = (0:M-1)'*d*sin(theta_s)/c;
fs = 3.125*f0; Ts = 1/fs;
L = 65;               % 滤波器长度
D = (L-1)/2;
T = -round(tau/Ts + D) * Ts;
Hd = conj(w_opt) .* exp(jay*2*pi*T*f);

%----------- 设计滤波器  -------%
deltaf = (fU-fL)/(FN-1)/fs;                     % 离散化频率间隔
PB = 0.16 : deltaf : 0.32;      % 通带带宽 其中 0.16 = fL/fs；0.32 = fU/fs；
SB = [    0:deltaf:0.13   0.35:deltaf:0.50 ];      % 阻带带宽
WB =  0 : deltaf : 0.5;                            % 全部带宽
e_PB = exp(-jay*2*pi*(0:L-1)'*PB);  
e_SB = exp(-jay*2*pi*(0:L-1)'*SB);  
e_WB = exp(-jay*2*pi*(0:L-1)'*WB);

NN = length(WB);
H = zeros(M,NN);
Tsig = 512/fs;
t = 0:1/fs:(Tsig-1/fs);
s = sin(2*pi*(fL+(fU-fL)/(2*Tsig)*t).*t);
for m = 1 : M
%----- 阻带峰值约束最小均方通带方法 -----%
    cvx_begin 
        variable h_2(L);
            minimize(norm(h_2.'*e_PB-Hd(m,:)))   
        subject to
            max(abs(h_2.'*e_SB))<=10^(-40/20);
    cvx_end

    H(m,:) = h_2.'*e_WB;
    figure(3); plot(WB, 20*log10(abs(H(m,:))), PB, 20*log10(abs(Hd(m,:))), '.');
    figure(4); plot(WB, angle(H(m,:))*180/pi, PB, angle(Hd(m,:))*180/pi, '.');

    h(:,m) = h_2;
    tt = t-tau(m);
    s_in(m,:) = sin(2*pi*(fL+(fU-fL)/(2*Tsig)*tt).*tt);
end

%--------------------波束形成输出序列---------------------%
N=length(s);
for m = 1:M,
    Tm = -(round(tau(m,1)/Ts)+D);
    if Tm>=0,
        s_Tm(m,:)=[zeros(1,Tm) s_in(m,1:N-Tm)];
    else
        s_Tm(m,:)=[s_in(m,-Tm+1:N) zeros(1,-Tm)];
    end
    s_fil(m,:)=filter(h(:,m),1,s_Tm(m,:));
end
y=sum(s_fil);
err = y - s;
figure; plot( t*fs, err , 'LineWidth', 1.0 );  %axis([0 512 -0.02 0.02])

%-------------------波束图-------------------%
hh = vec(h.');
for jj = 1 : length(PB),
    for ii = 1:length(theta),
        U = exp(-jay*2*pi*f(jj)*((T+(0:M-1)'*d*sin(theta(ii))/c)*ones(1,L)+ones(M,1)*(0:L-1)*Ts));
        u = vec(U);
        P(jj,ii) = hh.'*u;
    end
    figure(6); hold on
    plot(theta/angle_radian,20*log10(abs(P(jj,:))),'b','linewidth',1.0);
end