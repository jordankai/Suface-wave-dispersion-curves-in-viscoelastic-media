function [vs, vp] = model_KD_constant(cvs,cvp,Qs,Qp,freq)
%速度模型
% 参考频率为1Hz，参考lai论文以及aki的书
vp=zeros(length(cvp),length(freq));
vs=zeros(length(cvs),length(freq));
Qs = reshape(Qs,length(Qs),1);
Qp = reshape(Qp,length(Qp),1);
om = 2*pi*freq;
for jj=1:length(freq)
    p1 = cvp;
    s1 = cvs;
    hp = sqrt(1+1./(Qp.^2));
    hs = sqrt(1+1./(Qs.^2));
    vp(:,jj)=p1.*((1+hp)/2+1i./(2*Qp))./hp;
    vs(:,jj)=s1.*((1+hs)/2+1i./(2*Qs))./hs;

end