function y=fast_ht_kai(k,om,thk,dns,vp,vs)

%% Thomson–Haskell Method
%**************************************************************************
%    Written and Copyleft by Zhang Kai,IGGE, zhangkai@igge.cn 
%    Version 1, 2020/4/11
%    You may use and modify this code, provided you acknowledge the source
%**************************************************************************

n=length(vs);
epsilon = 0.0001;
 while any(abs(om/k-vs)<epsilon) || any(abs(om/k-vp)<epsilon)
    k = k * (1+epsilon);
 end
x2=(om/k)^2;
vs2=vs.^2;


  
%   for i=1:n
%     if(x<vp(i))
%         rp(i)=sqrt(1-x^2/vp(i)^2);
%     else
%         rp(i)=1i*(sqrt(x^2/vp(i)^2-1));
%     end
%     if(x<vs(i))
%         rs(i)=sqrt(1-x^2/vs(i)^2);
%     else
%         rs(i)=1i*(sqrt(x^2/vs(i)^2-1));
%     end
%    end
%     index = find(imag(-1i*nus) > 0);
%      nus(index) = -nus(index);
 %  gammas = nus/k;
   
   rs = sqrt(1-x2./vs2);
    rp = sqrt(1-x2./vp.^2);
%       index = find(imag(-1i*nup)> 0);
%       nup(index) = -nup(index);
 %  gammap= nup/k;


    t=2-x2./(vs2);
    u=dns.*vs2/1000000;

U=[0 0 1 0;0 0 0 1];
V=[1 rs(n);rp(n) 1;2*u(n)*rp(n) u(n)*t(n);u(n)*t(n) 2*u(n)*rs(n)];
E=zeros(4,4);

for i=1:n-1
 %% according to "A New Misfit Function for Multimode Dispersion Curve Inversion of Rayleigh Wave"by Cai Wei et.al
 %% "新的瑞雷波多模式频散曲线反演目标函数" 蔡伟等    
    tmp_2=1;%abs(exp(0.5*k*(rp(i)+rs(i))*thk(i)));  
    E(1,1)=exp(k*(rp(i))*thk(i))/tmp_2;
    E(2,2)=exp(-k*rp(i)*thk(i))/tmp_2;
    E(3,3)=exp(k*rs(i)*thk(i))/tmp_2;
    E(4,4)=exp(-k*rs(i)*thk(i))/tmp_2;
    Q=[1 1 rs(i) -rs(i);
        rp(i) -rp(i) 1 1;
        2*rp(i)*u(i) -2*rp(i)*u(i) t(i)*u(i) t(i)*u(i);
        t(i)*u(i) t(i)*u(i) 2*rs(i)*u(i) -2*rs(i)*u(i)];
    tmp1=0.5/(2*rp(i)-rp(i)*t(i));
    tmp2=1/(t(i)-2);
    tmp3=0.5/(2*rs(i)-rs(i)*t(i));
    Q_inv=[ -tmp2 -t(i)*tmp1 tmp1/u(i) 0.5*tmp2/u(i);
            -tmp2 t(i)*tmp1 -tmp1/u(i) 0.5*tmp2/u(i);
            -t(i)*tmp3 -tmp2 0.5*tmp2/u(i) tmp3/u(i);
            t(i)*tmp3 -tmp2 0.5*tmp2/u(i) -tmp3/u(i)];
  
    T=Q*E*Q_inv;
%    T=Q*E/Q;
    U=U*T;
end
y=((det(U*V)));