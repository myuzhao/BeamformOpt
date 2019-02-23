 %%传感器阵列波束优化设计与应用
 %%20170818  20181225重写
 %%myuzhao@163.com
 %%几种窗函数幅度加权对比
 %%均匀加权 余弦加权 hanning hamming 
 %%%书上的hanning 和hamming结论写反了
clc;
clear ;
close all;
%扫描计算范围
c0 = 340;
M = 10;
angle0 = 0*pi/180;
angle = linspace(-90,90,180);
m = 1:1:M;

as0 = exp(-1i*pi*(0:1:M-1)'*sind(angle0));
w = exp(-1i*pi*(0:1:M-1)'*sind(angle));   
w_uniform = 1/M *ones(M,1);

p_uniform = (w_uniform.*w)'*as0;
enegry_p_uniform = 20*log10(abs(p_uniform));
figure(1)
subplot(2,2,1)
stem(m, w_uniform)
ylim([0 1/M*2])
xlim([0  M+1])
hold on

figure(2)
plot(angle,enegry_p_uniform,'k-.')
hold on

%%cossin加权
w_cossin = cos(pi*(m-(M+1)/2)/M);
w_cossin = w_cossin.'/sum(w_cossin); %%%归一化
p_cossin = (w_cossin.*w)'*as0;
enegry_p_cossin=20*log10(abs(p_cossin));
figure(1)
subplot(2,2,2)
stem(m, w_cossin)
ylim([0 1/M*2])
xlim([0  M+1])
figure(2)
plot(angle,enegry_p_cossin,'b.')

%%hanning加权
w_hanning = hanning(M);
w_hanning = w_hanning/sum(w_hanning);
p_hanning = (w_hanning.*w)'*as0;
enegry_p_hanning=20*log10(abs(p_hanning));
figure(1)
subplot(2,2,3)
stem(m, w_hanning)
ylim([0 1/M*2])
xlim([0  M+1])

figure(2)
plot(angle,enegry_p_hanning,'r--')

%%hamming加权
w_hamming = hamming(M);
w_hamming = w_hamming/sum(w_hamming);
p_hamming = (w_hamming.*w)'*as0;
enegry_p_hamming = 20*log10(abs(p_hamming));
figure(1)
subplot(2,2,4)
stem(m, w_hamming)
ylim([0 1/M*2])
xlim([0  M+1])

figure(2)
plot(angle,enegry_p_hamming,'g-')
xlabel('方位/(^o)')
ylabel('波束/dB')
ylim([-80 3])
grid on

legend('uniform','cossin','hanning','hamming')