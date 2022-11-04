function [vs, vp] = model_KK_constant(cvs,cvp,Qs,Qp,freq)
%速度模型
% 参考频率为1Hz，参考lai论文以及aki的书 second order approximation
vp=zeros(length(cvp),length(freq));
vs=zeros(length(cvs),length(freq));
Qs = reshape(Qs,length(Qs),1);
Qp = reshape(Qp,length(Qp),1);
om = 2*pi*freq;
for jj=1:length(freq)
    p1 = cvp;
    s1 = cvs;
    vp(:,jj)=p1.*(1+1i./(2*Qp))./(1+(0.5./Qp).^2);
    vs(:,jj)=s1.*(1+1i./(2*Qs))./(1+(0.5./Qs).^2);

end