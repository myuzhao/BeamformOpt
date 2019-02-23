clear
clc
close all

format = ['r-.';'b-+';' k.'];
%%%%%%%%%%%%%%%%%%%%%%信号参数%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C0 = 340;
Snapshots = 15;
Fs = 10000;
Angle = [0  90;    %%%水平角  俯仰角
%          -5 90
         ];
SNR =   [0];
Freq = [1000;
         ];
Freqs_ref = 1000;
rand_state = 1;
rng('default');
rng(rand_state);
%%%%%%%%%%%%%%%%%%%%%%阵列参数%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Sensor = 2; %%Sound_Pressur%%Acoustic_Vector_Sensor
if Sensor == 2
    Sensor_Channel_Num = 4;
else
    Sensor_Channel_Num = 1;
end

M = 10;
Pos_array = zeros(M,3);
Pos_array(:,1) = (0:1:M-1) * (C0/Freqs_ref/2);
Array_Loc = Pos_array;

Angle0 = Angle;
uu =[sind(Angle0(2)) * sind(Angle0(1));
     sind(Angle0(2)) * cosd(Angle0(1));
     cosd(Angle0(2))];
     as = exp(-1i * 2 * pi * Freq /C0 * (Array_Loc * uu) );
     if( Sensor == 2)
         avxs= as * uu(1);
         avys= as * uu(2);
         avzs= as * uu(3);
         steer_vector =[as; avxs; avys; avzs];%
     else
         steer_vector = as;
     end
a0 = steer_vector;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
resolution = 1;
scan_angle1 = [-90:resolution:90]; %%%划分时不存在网格误差
scan_angle2 = Angle(1,2);          %%%固定俯仰角
a_all =[]; 
for i_1 = 1:length(scan_angle1) %%%水平角
        Angle0 = [scan_angle1(i_1), scan_angle2];
        uu =[sind(Angle0(2)) * sind(Angle0(1));
             sind(Angle0(2)) * cosd(Angle0(1));
             cosd(Angle0(2))];
        as = exp(-1i * 2 * pi * Freq /C0 * (Array_Loc * uu) );
        if( Sensor == 2)
            avxs= as * uu(1);
            avys= as * uu(2);
            avzs= as * uu(3);
            steer_vector =[as; avxs; avys; avzs];% 
        else
            steer_vector = as;
        end
        a_all= [a_all steer_vector];
end
p = 20*log10(abs(a_all'*a0)/max(abs(a_all'*a0)));

angle_ml = -12:1:12;
angle_sl1 = [-90:1:-12];
angle_sl2 = [12:1:90];
angle_sl =[angle_sl1 angle_sl2];

% a_head = zeros(size(a_all));
index_ml = find(scan_angle1>=angle_ml(1) & scan_angle1<=angle_ml(end))
       
a_ml = a_all(:,index_ml);%%%主瓣导向向量
p_ml = p(index_ml);


index_s1 = find(scan_angle1>=angle_sl1(1) & scan_angle1<=angle_sl1(end));     
a_sl1  = a_all(:,index_s1);

index_s2 = find(scan_angle1>=angle_sl2(1) & scan_angle1<=angle_sl2(end));    
a_sl2  = a_all(:,index_s2);

index_sl = [index_s1 index_s2];
a_sl = [a_sl1 a_sl2];%%%旁瓣导向向量
p_sl = p(index_sl);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
plot(scan_angle1 , p,format(1,:))
hold on
plot(angle_ml , p_ml,format(2,:))
plot(angle_sl , p_sl,format(3,:))
ylim([-90 0])
xlabel('角度/^o','fontsize',14)
ylabel('波束响应/dB','fontsize',14)
legend({'波束响应','主瓣区域','旁瓣区域'},'fontsize',14,'Location','best')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%主瓣无约束，最小旁瓣波束响应%%%%%%%
cvx_begin quiet%
   variable w_msl(M*4) complex
   variable s_msl(1)
   minimize(s_msl)
   subject to
        abs(w_msl'*a_sl)<=s_msl;
        w_msl'* a0 == 1;
cvx_end 

figure
plot(scan_angle1 , p,format(1,:))
p_msl =20*log10(abs(w_msl' * a_all)/max(abs(w_msl' * a_all)));
hold on
plot(scan_angle1,p_msl ,format(2,:))
ylim([-90 0])
xlabel('角度/^o','fontsize',14)
ylabel('波束响应/dB','fontsize',14)
legend({'DAS波束响应','最低旁瓣波束响应'},'fontsize',14,'Location','best')
norm(w_msl)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%主瓣无约束，稳健最小旁瓣波束响应%%%%%%%
cvx_begin quiet%
   variable w_msl(M*4) complex
   variable s_msl(1)
   minimize(s_msl)
   subject to
        abs(w_msl'*a_sl)<=s_msl;
        w_msl'* a0 == 1;
        norm(w_msl) <= 1/(M*4)*10;
cvx_end 

figure
plot(scan_angle1 , p,format(1,:))
p_msl =20*log10(abs(w_msl' * a_all)/max(abs(w_msl' * a_all)));
hold on
plot(scan_angle1,p_msl ,format(2,:))
ylim([-90 0])
xlabel('角度/^o','fontsize',14)
ylabel('波束响应/dB','fontsize',14)
legend({'DAS波束响应','稳健最低旁瓣波束响应'},'fontsize',14,'Location','best')
norm(w_msl)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%旁瓣控制主瓣最小误差%%%%%%%
max_sl = -20; %%%
max_sl=10.^(max_sl/20);
p_ml = 10.^(p_ml/20); %%%期望主瓣响应
cvx_begin quiet%
    variable w_msl(M*4) complex
    p_temp=(w_msl' * a_all);
    p_temp_ml=p_temp(index_ml);
    p_temp_ml=p_temp_ml(:); %%%优化主瓣响应
    p_ml=p_ml(:);%%%期望主瓣响应
    p_temp_sl=p_temp(index_sl);%%%优化旁瓣响应
    
    minimize(norm(p_temp_ml-p_ml,2))  %ML 2范数
subject to
    abs(p_temp_sl) <= max_sl;%%sl L_inf
    w_msl'*a0 == 1;
%     norm(w_msl,2) <=55;% 2.5;%$1/(M*4)*12500;
cvx_end

figure
plot(scan_angle1 , p,format(1,:))
p_msl =20*log10(abs(w_msl' * a_all)/max(abs(w_msl' * a_all)));
hold on
plot(scan_angle1,p_msl ,format(2,:))
ylim([-90 0])
xlabel('角度/^o','fontsize',14)
ylabel('波束响应/dB','fontsize',14)
legend({'DAS波束响应','旁瓣控制主瓣最小误差'},'fontsize',14,'Location','best')
norm(w_msl,2)