function [vs ,vp] = model_simple(cvs,cvp,Qs,Qp,freq)
%速度模型
 
Qs = reshape(Qs,length(Qs),1);
Qp = reshape(Qp,length(Qp),1);
for jj=1:length(freq)
    vp(:,jj)=cvp.*(1+1i./(2*Qp));
    vs(:,jj)=cvs.*(1+1i./(2*Qs));

end
