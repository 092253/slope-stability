clc
H_s=input('please input a number H_s:')%H_f的一半
L=input('please input a number L:')
c=5;phi=18*pi/180;%湿润锋处的抗剪强度
co=3;phio=18*pi/180;%有效抗剪强度饱和区
ZL=10;
gam_d=13.5;gam_s=19.5;
cta_s=0.45;cta_r=0.015;cta_i=0.15;alpha=3.5;m=0.333;n=1.5;
q=0.036*3600;%rainfall intensive 单位m/h
u_a=0;
beta=30*pi/180;%slope angle
Ks=1e-6;
hs=0.5;%单位m
%cta_zf=cta_i
%cta_zf=cta_r+(cta_s-cta_r)/(1+(alpha*hf)^n)^m;
%cta_z=cta_zf+(cta_s-cta_zf)*sqrt(1-(z/zf)^2);%the relationship between cta_z and zf 
H_c=hs/(q*cos(beta)/Ks-cos(beta));%临界入渗深度 饱和层的厚度，湿润锋的厚度应再乘以2
I_c=(cta_s-cta_i)*(1+pi/4)*H_c;%临界累积入渗量
t_c=I_c/(q*cos(beta));
M=-L*cos(beta)-sqrt(L^2*(cos(beta))^2+4*L*hs*sin(beta));
N=2*(1+pi/4)*sin(beta);
if H_s<H_c
    
    t=(1+pi/4)*(cta_s-cta_i)*H_s/(q*cos(beta));
else
    %t=L*(1+pi/4)(cta_s-cta_i)/(Ks*sin(beta))*(N/(M-N)*log(((1+pi/4)*H_s+N))/(H_c+N))-M/(M-N)*log(((1+pi/4)*H_s+M)/(H_c+M))+t_c
    t=L*(1+pi/4)*(cta_s-cta_i)/(Ks*sin(beta))*(N/(M-N)*log(((1+pi/4)*H_s+N)/(H_c+N))-M/(M-N)*log(((1+pi/4)*H_s+M)/(H_c+M)))+t_c
    %t=((cta_zf-cta_i)+pi*(cta_s-cta_zf)/4))*(zf-zp-hf/cos(beta)*log((zp*cos(beta)+hf)/(zp*cos(beta)+hf))/3600+t_p;
end 
%饱和区的安全系数
Wo=L*H_s*gam_s;
Wno=Wo*sin(beta);
sigmao=H_s*gam_s*cos(beta);
J=10*L*H_s*sin(beta);
SF2=(co+sigmao*tan(phio))*L/(Wno*sin(beta)+J)
SFS=(co+gam_s*H_s*cos(beta)*tan(phio))/(H_s*sin(beta)*(10+gam_s*sin(beta)))
%湿润锋处的安全系数
W=gam_s*H_s*L+H_s*L*gam_d*(1+cta_i+pi/4*(cta_s-cta_i)^2);
sigmn=W*cos(beta)/L;
R=(c+sigmn*tan(phi)+ZL)*L;
Wn=W*sin(beta);
FSW=R/(Wn+J)


%sigmas=-(cta_zf-cta_r)*hf/(cta_s-cta_r);
%sigmap=(4*(1+cta_zf)+pi*(cta_s-cta_zf))*gam_d*zf*cos(beta)/4-sigmas;
%if 
%FS=tan(phi)/tan(beta)+c/((4*(1+cta_zf)+pi*(cta_s-cta_zf))*gam_d*z_f*sin(beta)/4)-...
    %(u_a+simas)*tan(phi)/((4*(1+cta_zf)+pi*(cta_s-cta_zf))*gam_d*z_f*sin(beta/4))
    %FS=4*(c+sigmap*tan(phi))/((4*(1+cta_zf)+pi*(cta_s-cta_zf))*gam_d*zf*sin(beta)/4)

    
