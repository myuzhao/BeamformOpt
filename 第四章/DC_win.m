function [ w_dc] = DC_win(angle11,sll,d_lamda,M,type)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% angle11:指定主瓣宽度
% sll:指定旁瓣级
% d_lamda:阵元间距d与波长lamda的关系
% M:阵元数
% type:1:指定主瓣宽度 else:指定旁瓣级
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%code : myuzhao@163.com  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
w=zeros(M,1);
if(type==1)
    z = cos(pi/(2*(M-1)))/cos(pi*d_lamda*sin(angle11));%指定主瓣宽度
else
    z=cosh(acosh(sqrt(10^(-sll/10)))/(M-1)); %指定旁瓣级
end

w(1)=(z^(M-1))/2;
temp=M/2+1;
m=2;
while m<temp
    for i=1:m-1
        tenpm=0.5*(M-1)*factorial(m-2)*factorial(M-i-1)/(factorial(m-i)*factorial(i-1)*factorial(m-i-1)*factorial(M-m))...
            *(z^(M-2*m+1))*((z^2-1)^(m-i));
        w(m)=w(m)+tenpm;
    end
    m=m+1;
end

for m=m:1:M
    w(m)=w(M+1-m);
end
w_dc=w;


