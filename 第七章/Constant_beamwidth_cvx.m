function [ output_args ] = Constant_beamwidth_cvx(array_x,array_y,array_z,f_ml,f_sl,f_jg,fs,theta_ml,theta_sl,theta_jg,sl,c0)
    if (nargin==0)
       fs=10000;
       f_ml=[500 1000];
       f_sl=[490 1010];
       f_jg=10;
       fml=f_ml(1):f_jg:f_ml(2);
       fsl=[0:f_jg:f_ml(1) f_ml(2):f_jg:fs/2;
       
       theta_jg=1;  
       theta=-90:theta_jg:90;
       
       thetaml= theta_ml(1):theta_jg:theta_ml(2);
       thetasl=[-90:theta_jg:theta_ml(1)  theta_ml(2):theta_jg:90];
     
    end


end

