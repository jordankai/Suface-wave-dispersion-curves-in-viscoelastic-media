function [vs, vp] = visco_model(vs1,vp1,freq)
%速度模型
freq=reshape(freq,1,length(freq));
Tao11=0.0344113;Tao12=0.0028676;Tao21=0.0414902;Tao22=0.0034575;%%分子松弛时间
tao11=0.0287482;tao12=0.0023957;tao21=0.0287482;tao22=0.0023957;%%分母松弛时间

M1=2*(vp1.^2-vs1.^2);
M2=2*vs1.^2;
fa1=(1+1i*2*pi*freq*Tao11)./(1+1i*2*pi*freq*tao11);
fa2=(1+1i*2*pi*freq*Tao12)./(1+1i*2*pi*freq*tao12);
fb1=(1+1i*2*pi*freq*Tao21)./(1+1i*2*pi*freq*tao21);
fb2=(1+1i*2*pi*freq*Tao22)./(1+1i*2*pi*freq*tao22);
MC1=M1.*(0.5*(fa1+fa2));  %此处速度为列向量，频率为数组，两者的乘积即为速度矩阵。
MC2=M2.*(0.5*(fb1+fb2));
vp=((MC1+MC2)/2).^0.5;
vs=(MC2/2).^0.5;

end

