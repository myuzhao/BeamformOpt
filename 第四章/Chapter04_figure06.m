%%传感器阵列波束优化设计与应用
 %%20170818
 %%myuzhao@163.com
 %%扇面内存在多个干扰时的 MVDR 波束图
clc;
clear;
close all;

snaps=5000;

% %阵元位置
M = 10;
%%扫描参数设置
angle0 = 0;%主轴方向
angle1 = linspace(-50,-20,16);%干扰方向
angle =linspace(-90,90,10000);
type1 = ['k.';'r-';'b-'];
inr=[-10 10];
noise=1/sqrt(2)*(randn(M,snaps)+1i*randn(M,snaps));
for i=1:2
    inr1=inr(i);%干燥比
    s1 = 10^(inr1/20)*randn(16,snaps);
    ap=exp(-1i*pi*(0:1:M-1)'*sind(angle1)); 
    apis=ap*s1;
    
    apisn = apis+noise;
    R = apisn*apisn'/snaps;

    as0 = exp(-1i*pi*(0:1:M-1)'*sind(angle0));%期望波束
    iR = inv(R);
    w = iR*as0/(as0'*iR*as0);
    ap = exp(-1i*pi*(0:1:M-1)'*sind(angle));   
    p = ap'*w;
    enegry_mvdr=20*log10(abs(p));

    figure1 = figure(1)
    hold on
    plot(angle,enegry_mvdr,type1(i))
    xlabel('方位/(^o)')
    ylabel('波束/dB')
    ylim([-100 3])
    grid on

end

title('扇面内存在多个干扰时的 MVDR 波束图')
legend('INR=-10dB','INR=10dB')
annotation(figure1,'arrow',[0.388534783406755 0.388534783406756],...
    [0.898546405360363 0.763832119646074]);
annotation(figure1,'arrow',[0.408195668135096 0.408195668135097],...
    [0.897852124045927 0.763137838331637]);
annotation(figure1,'arrow',[0.398045154185022 0.398045154185023],...
    [0.897424070904696 0.762709785190406]);
annotation(figure1,'arrow',[0.378384269456682 0.378384269456683],...
    [0.899402511642826 0.764688225928536]);
annotation(figure1,'arrow',[0.368694933920705 0.368694933920706],...
    [0.898280177187158 0.763565891472868]);
annotation(figure1,'arrow',[0.358083241556534 0.358083241556535],...
    [0.900686671066521 0.765972385352231]);
annotation(figure1,'arrow',[0.3490937041116 0.349093704111601],...
    [0.901380952380957 0.766666666666667]);
annotation(figure1,'arrow',[0.324999999999999 0.325],...
    [0.901380952380957 0.766666666666667]);
annotation(figure1,'arrow',[0.339583333333333 0.339583333333333],...
    [0.898013949013954 0.763299663299664]);