clear all;
clc;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  变量初始化 %%%%%%%%%%%%%%%%%%%%%%%%%%
M = 12;       	% 12元均匀线列阵

% 线性调频信号LPM参数
f0 = 1000;    	% f0
f1 = 0.5*f0;    
fu = f0;       	% f1和fu为调频信号源的上下边界频率 
fs = 5*f0;     	% 采样频率fs
T = 512/fs;    	% 信号持续时间，满足T*fs = 512个采样点
N = ceil(T*fs);	% 总的采样点数512；此处如果设定N=T*fs，若不是整数倍的采样频率，则易出错

% λ0=c/f0
c = 3*10^8;
lambda0 = c/f0;

% 取出频域数据中k = 21,...56的子带数据，对应频率范围为[0.4102f0,1.0938f0]
% fk = fs*k/L      , k = 1,...,L/2-1
%    = fs*(k-L)/L  , k = L/2,...L-1
L = N/2; % 缓存数据块长度256点
k = 21:1:56;
fk = zeros(1,length(k));
for i=1:length(k)
    if k(i)< L/2
        fk(i) = fs*k(i)/L;   
    else
        fk(i) = fs*(k(i)-L)/L;
    end
end

% 信号入射方向为θ0
phi = 90*pi/180; % xOy平面，phi = 90°
theta0 = -30*pi/180; % θ0 = -30°

% 阵元坐标向量p_m
p_element_coordinate = ones(2,1,M);
p_element_coordinate(1,1,:) = 0; % x = 0
p_element_coordinate(2,1,:) = (0:M-1).*lambda0/2; % y_m = (m-1)/2*λ

% 信号传播方向的单位向量v(θ0)
v_theta0 = -[sin(phi)*cos(theta0),sin(phi)*sin(theta0)].'; % “.'”非共轭转置，“'”共轭转置

%%%%%%%%%%%%%%%%%%%% part2: 常规波束形成器加权向量w_c(fk) %%%%%%%%%%%%%%%%%%
% 期望方向为θ0时的阵列流形向量a(θ0)
omega = 2*pi*fk; % w = 2πfk
a_theta0 = zeros(M,length(k));
k_theta0 = zeros(2,length(k));
for i=1:length(k)
	k_theta0(:,i) = omega(i)/c*v_theta0; 
    for m=1:M
        a_theta0(m,i) = exp(-1i*k_theta0(:,i).'*p_element_coordinate(:,:,m));
    end
end

% 窄带常规波束形成器的加权向量w_c(fk)
w_c = a_theta0/M;

%%%%%%%%%%%%%%%%%%% part3: 阵元接收到有效信号s_receive_LFM %%%%%%%%%%%%%%%%%
% 阵元间延时，阵元1接收到数据时阵元M还未接收到数据，τm = v(θ0).'*pm/c
tau = zeros(1,M);
for m=1:M
    tau(m) = v_theta0.'*p_element_coordinate(:,:,m)/c;
end
delay = floor(tau*fs);  % 延时采样间隔取整数，τm/Ts = τm*fs

% 经过采样后的线性调频信号s_LFM(i) = sin(2*pi*(f1+(fu-f1)/(2*T)*t(i))*t(i));
te = (0:1:N-1);
t = (0:1:N+delay(M)-1)./fs;
s_receive_LFM = zeros(M,N); %有效信号持续长度512点
s_receive_LFM_tmp = zeros(M,length(t));  % 阵元接收数据
for m=1:M	
	for i=1:length(t)
        s_receive_LFM_tmp(m,i) = sin(2*pi*(f1+(fu-f1)/(2*T)*(t(i)-tau(m)))*(t(i)-tau(m)));  % 阵元m接收到的信号sm(t)=s(t-τm)
        
        % 采样时间 0<=t<T
        if t(i)-tau(m)<0
            s_receive_LFM_tmp(m,i) = 0;
        elseif t(i)-tau(m)>=T
            s_receive_LFM_tmp(m,i) = 0;
        end
        
        % 取前512个点作为采样数据
        if i<=N
            s_receive_LFM(m,i)= s_receive_LFM_tmp(m,i);
        end
    end
end

%%%%%%%%%%%%%%%%% part4: 阵元有效数据分段进行DFT：无重叠分为2段%%%%%%%%%%%%%%
%      |     ,     , ... ,    |          |     ,     , ... ,     | → k0
%      |     ,     , ... ,    |          |     ,     , ... ,     | → k1
% x(n)=|     ,     , ... ,    |   X(k) = |     ,     , ... ,     | → k2
%      |     ,     , ... ,    |          |     ,     , ... ,     | ...
%      |     ,     , ... ,    |          |     ,     , ... ,     | → kL-1
%        ↑    ↑          ↑               ↑    ↑          ↑
%       x1(n) x2(n)   ... xM(n)           X1(k) X2(k)  ...  XM(k)

xn_2_seg = zeros(L,M,2);
Xk_2_seg = zeros(L,M,2);
for m=1:M
    for i=1:L
        xn_2_seg(i,m,1) = s_receive_LFM(m,i);
        xn_2_seg(i,m,2) = s_receive_LFM(m,i+L);
    end
%     Xk_2_seg(:,m,1) = dft(xn_2_seg(:,m,1),L);
%     Xk_2_seg(:,m,2) = dft(xn_2_seg(:,m,2),L);
      Xk_2_seg(:,m,1) = fft(xn_2_seg(:,m,1));
      Xk_2_seg(:,m,2) = fft(xn_2_seg(:,m,2));
end

Yk_2_seg = zeros(2,length(k));
for i=1:length(k)
    Yk_2_seg(1,i) = w_c(:,i)'*Xk_2_seg(k(i)+1,:,1).';   % k(i)+1 才对应X(21)→X(56)
    Yk_2_seg(2,i) = w_c(:,i)'*Xk_2_seg(k(i)+1,:,2).';
end

%%%%%%%%%%%%%%%%% part5: 2段波束输出频域数据进行IDFT，得到时域输出 %%%%%%%%%
Yk_2_seg_full = zeros(L,2);
for i=1:length(k)
    Yk_2_seg_full(k(i)+1,1) = Yk_2_seg(1,i);            % k(i)+1 才对应Y(21)→Y(56)
    Yk_2_seg_full(L-k(i)+1,1) = conj(Yk_2_seg(1,i));
    
    Yk_2_seg_full(k(i)+1,2) = Yk_2_seg(2,i);            
    Yk_2_seg_full(L-k(i)+1,2) = conj(Yk_2_seg(2,i));    % L-k(i)+1对应的负频率
end

yn_2_seg_full = zeros(L,2);
% yn_2_seg_full(:,1) = idft(Yk_2_seg_full(:,1),L);
% yn_2_seg_full(:,2) = idft(Yk_2_seg_full(:,2),L);
yn_2_seg_full(:,1) = ifft(Yk_2_seg_full(:,1));
yn_2_seg_full(:,2) = ifft(Yk_2_seg_full(:,2));

% 数据无重叠进行拼接
yn_2_seg = zeros(1,N);
for i=1:L
    yn_2_seg(i) = yn_2_seg_full(i,1);
    yn_2_seg(i+L) = yn_2_seg_full(i,2);
end

%%%%%%%%%%%%%%%%% part6: 阵元有效数据分段进行DFT：50%重叠分为3段%%%%%%%%%%%%%%
xn_3_seg = zeros(L,M,3);
Xk_3_seg = zeros(L,M,3);
for m=1:M
    for i=1:L
        xn_3_seg(i,m,1) = s_receive_LFM(m,i);
        xn_3_seg(i,m,2) = s_receive_LFM(m,i+L/2);
        xn_3_seg(i,m,3) = s_receive_LFM(m,i+L);
    end
%     Xk_3_seg(:,m,1) = dft(xn_3_seg(:,m,1),L);
%     Xk_3_seg(:,m,2) = dft(xn_3_seg(:,m,2),L);
%     Xk_3_seg(:,m,3) = dft(xn_3_seg(:,m,3),L);
    Xk_3_seg(:,m,1) = fft(xn_3_seg(:,m,1));
    Xk_3_seg(:,m,2) = fft(xn_3_seg(:,m,2));
    Xk_3_seg(:,m,3) = fft(xn_3_seg(:,m,3));
end

Yk_3_seg = zeros(3,length(k));
for i=1:length(k)
    Yk_3_seg(1,i) = w_c(:,i)'*Xk_3_seg(k(i)+1,:,1).';
    Yk_3_seg(2,i) = w_c(:,i)'*Xk_3_seg(k(i)+1,:,2).';
    Yk_3_seg(3,i) = w_c(:,i)'*Xk_3_seg(k(i)+1,:,3).';
end

%%%%%%%%%%%%%%%%% part7: 3段波束输出频域数据进行IDFT，得到时域输出 %%%%%%%%%
Yk_3_seg_full = zeros(L,3);
for i=1:length(k)
    Yk_3_seg_full(k(i)+1,1) = Yk_3_seg(1,i);
    Yk_3_seg_full(L-k(i)+1,1) = conj(Yk_3_seg(1,i));
    
    Yk_3_seg_full(k(i)+1,2) = Yk_3_seg(2,i);
    Yk_3_seg_full(L-k(i)+1,2) = conj(Yk_3_seg(2,i));
    
    Yk_3_seg_full(k(i)+1,3) = Yk_3_seg(3,i);
    Yk_3_seg_full(L-k(i)+1,3) = conj(Yk_3_seg(3,i));
end

yn_3_seg_full = zeros(L,3);
% yn_3_seg_full(:,1) = idft(Yk_3_seg_full(:,1),L);
% yn_3_seg_full(:,2) = idft(Yk_3_seg_full(:,2),L);
% yn_3_seg_full(:,3) = idft(Yk_3_seg_full(:,3),L);
yn_3_seg_full(:,1) = ifft(Yk_3_seg_full(:,1));
yn_3_seg_full(:,2) = ifft(Yk_3_seg_full(:,2));
yn_3_seg_full(:,3) = ifft(Yk_3_seg_full(:,3));

% 数据按照50%重叠进行拼接，顺序不能变动，而且不能在一个循环里完成
yn_3_seg = zeros(1,N);
for i=1:L
    yn_3_seg(i) = yn_3_seg_full(i,1);  	   
end

for i=1:L
    yn_3_seg(i+L/2) = yn_3_seg_full(i,2);
end

for i=1:L
	yn_3_seg(i+L) = yn_3_seg_full(i,3);  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%% part8: plot all figures %%%%%%%%%%%%%%%%%%%%%%%
% Figure1: 采样后信号波形
plot(te,s_receive_LFM(1,:))
xlabel('\iti');
ylabel('{\its}({\iti})');
axis([0 N-1 -1.25 1.25]);

% Figure2/3: 信号的功率谱密度
% algorithm 1:直接法求功率谱：|Xk|^2/N估计功率谱,matlab给出的例程
figure;
Xk = fft(s_receive_LFM(1,:));
Px = abs(Xk).^2/N/(fs/f0)*2;        % 注意:*2 是关键，因为此处仅取[0,fs/2)频带，并且归一化
Hpsd = dspdata.psd(Px(1:length(Px)/2),'Fs',fs/f0); 
F_1 = Hpsd.Frequencies;
PSD_1 = Hpsd.Data;
plot(F_1,10*log10(PSD_1),'LineWidth',2);
xlabel('频率({\itf}_0)');
ylabel('功率谱/(dB/Hz)')
title('|{\itX}({\itk})|^2估计功率谱密度');
axis([0 fs/f0/2 -50 10]);
grid on;

% algorithm 2: 利用periodogram周期图法求,其实和algorithm 1无本质区别
figure;
Hs = spectrum.periodogram;
Hpsd = psd(Hs,s_receive_LFM(1,:),'Fs',fs/f0);   % fs/f0归一化处理
F_2 = Hpsd.Frequencies;
PSD_2 = Hpsd.Data;
plot(F_2,10*log10(PSD_2),'LineWidth',2);
xlabel('频率({\itf}_0)');
ylabel('功率谱/(dB/Hz)')
title('Periodogram周期图法估计功率谱密度');
axis([0 fs/f0/2 -50 10]);
grid on;

% Figure4: 各个阵元接收到的信号波形
figure;
for m=1:M
    t = 1:1:(N+delay(M));
	plot(t,s_receive_LFM_tmp(m,:)/2+m); 
    hold all;
end
xlabel('\iti');
ylabel('阵元号\itm');
axis([0 N+delay(M)-1 0 13]);

% Figure5: 阵元数据无重叠分为2段，波束输出频域数据幅度|Y(k)|/dB
figure;
plot(k,20*log10(abs(Yk_2_seg(1,:))),':dk',...
                                  	'LineWidth',2,...
                                   	'MarkerEdgeColor','k',...
                                   	'MarkerFaceColor','k',...
                                   	'MarkerSize',5);
hold on;
plot(k,20*log10(abs(Yk_2_seg(2,:))),'-ob',...
                                  	'LineWidth',2,...
                                	'MarkerEdgeColor','b',...
                                  	'MarkerFaceColor','b',...
                                  	'MarkerSize',5);
legend('{\itn}=1','{\itn}=2','Location','South');
xlabel('\itk');
ylabel('|{{\itY}^{({\itn})}}({\itk})|/dB');
axis([k(1)-3 k(length(k))+3 -40 40]);

% Figure6: 2段IDFT时域输出数据拼接
figure;
plot(te,real(yn_2_seg));
xlabel('\iti');
ylabel('{\ity}({\iti})');
axis([0 N-1 -1.5 1.5]);

% Figure7: 阵元数据无重叠分为3段，波束输出频域数据幅度|Y(k)|/dB
figure;
plot(k,20*log10(abs(Yk_3_seg(1,:))),':dk',...
                                    'LineWidth',2,...
                                 	'MarkerEdgeColor','k',...
                                    'MarkerFaceColor','k',...
                                    'MarkerSize',5);
hold on;
plot(k,20*log10(abs(Yk_3_seg(2,:))),'--ob',...
                                 	'LineWidth',2,...
                                   	'MarkerEdgeColor','b',...
                                   	'MarkerFaceColor','b',...
                                 	'MarkerSize',5);
hold on;
plot(k,20*log10(abs(Yk_3_seg(3,:))),'-dr',...
                                  	'LineWidth',2,...
                                   	'MarkerEdgeColor','r',...
                                   	'MarkerFaceColor','r',...
                                  	'MarkerSize',5);
legend('{\itn}=1','{\itn}=2','{\itn}=3','Location','South');
xlabel('\itk');
ylabel('|{{\itY}^{({\itn})}}({\itk})|/dB');
axis([k(1)-3 k(length(k))+3 -40 40]);

% Figure8: 3段IDFT时域输出数据拼接
figure;
plot(te,real(yn_3_seg));
xlabel('\iti');
ylabel('{\ity}({\iti})');
axis([0 N-1 -1.5 1.5]);

% Figure9: DFT波束输出序列与信号源波形失真大小
figure;
plot(te,real(yn_2_seg)-s_receive_LFM(1,:),'--k','LineWidth',2);
hold on;
plot(te,real(yn_3_seg)-s_receive_LFM(1,:),'-b','LineWidth',2);
xlabel('\iti');
ylabel('{\ity}({\iti})-{\its}({\iti})');
legend('不重叠','重叠50%','Location','SouthWest');
axis([0 N-1 -1.5 1.5]);
