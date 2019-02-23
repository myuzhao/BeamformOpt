%%%% 例6.5 时域宽带常规波束形成
clear; clc
close all;
jay=sqrt(-1);

M = 20;                   % 阵元数
f0 = 100; c = 340;
lambda = c/f0;
d = lambda/2;
angle_radian = 1/180*pi;            %角度转换为弧度
theta_s = 10*angle_radian;

%-----------产生线性调频信号---------------
fL = f0/2; fU = 200; fs = 5*f0; Ts = 1/fs;
% fL = 800; fU = 1800; fs = 8192;
Tsig = 512/fs;
t = 0:1/fs:(Tsig-1/fs);
% s = chirp(t,fL,T,fU);
s = sin(2*pi*(fL+(fU-fL)/(2*Tsig)*t).*t);
% figure; plot( t*fs, s , 'LineWidth', 1.0 );  axis([0 512 -1.5 1.5])

% NFFT = 2^nextpow2(length(t)); % Next power of 2 from length of y
% S = fft(s,NFFT);
% flfm = fs/2*linspace(0,1,NFFT/2+1);
% XdB = 20*log10(abs(S(1:NFFT/2+1)));
% figure; plot(flfm,XdB); xlabel('Frequency (Hz)'); ylabel('|S(f)|');grid on
%-----------产生线性调频信号---------------

L = 65;                   % 滤波器长度
D = (L-1)/2;
% posi=[0 5.5 10.5 15 19 23.5 28.5 34];
posi=(0:(M-1));
tau = posi'*d*sin(theta_s)/c; 
deltaf = (fU - fL)/40;
fd = fL/fs : deltaf/fs : fU/fs;
ww=chebwin(M,30);
for m = 1:M,
%-----------期望滤波器响应---------------
    Hd(m,:) = (ww(m)/sum(ww))*exp(-jay*2*pi*fd*(D+round(tau(m,1)/Ts)-tau(m,1)/Ts));

%---------- e(f) -------------
    PB = fL/fs : deltaf/fs : fU/fs;
    SB = [ 0:deltaf/fs:(fL-8*deltaf)/fs  (fU+8*deltaf)/fs:deltaf/fs:0.5 ];
    e_PB = exp(-jay*2*pi*(0:L-1)'*PB);  
    e_SB = exp(-jay*2*pi*(0:L-1)'*SB);  
%----- 阻带峰值约束最小均方通带方法 -----% 
    cvx_begin
        variable h_2(L);
        minimize(norm(e_PB.'*h_2-Hd(m,:).'))
        subject to
            max(abs(e_SB.'*h_2))<=10^(-40/20);
    cvx_end

    h(:,m) = h_2;
    tt = t - tau(m,1);
    s_in(m,:) = sin(2*pi*(fL+(fU-fL)/(2*Tsig)*tt).*tt);
end

figure(1); hold on; grid on
plot(fd, 20*log10(abs(Hd(2,:))),'*'); axis([0 0.5 -60 -20])
figure(2); hold on;
plot(fd, angle(Hd(2,:))*180/pi,'*');  axis([0 0.5 -200 200])

fB = [ 0:deltaf/fs:0.5 ];
e = exp(-1i*2*pi*(0:L-1)'*fB);  
H2 = e.'*h(:,2);
figure(1); plot(fB, 20*log10(abs(H2)),'k'); 
figure(2); plot(fB, angle(H2)*180/pi,'k');

%------------------波束输出
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

%----------------波束图
theta = -90:2:90;
T = -round(tau/Ts + D) * Ts;
hh = vec(h.');
for jj = 1 : length(PB),
    for ii = 1:length(theta),
        U = exp(-jay*2*pi*PB(jj)*fs*((T+posi'*d*sind(theta(ii))/c)*ones(1,L)+ones(M,1)*(0:L-1)*Ts));
        u = vec(U);
        P(jj,ii) = hh.'*u;
    end
    figure(5); hold on
    plot(theta,20*log10(abs(P(jj,:))),'b','linewidth',1.0);
    set(gca,'YLim',[-100 5])
end