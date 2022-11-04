function dd = Fast_Delta(k,om,thk,dns,cvp,cvs)

%  2011.4 by Kai
% 参开1993年S.Ivansson 的文章
% Check to see if the trial phase velocity is equal to the shear wave velocity
% or compression wave velocity of one of the layers
epsilon = 0.0001;
while any(abs(om/k-cvs)<epsilon) || any(abs(om/k-cvp)<epsilon)
   k = k * (1+epsilon);
end   
%     k = om/vr;
    
   cvs2 = cvs.^2; cvp2 = cvp.^2;
   mu = dns.*cvs2;

   k2 = k^2; om2 = om^2;
   p = k/om;
   p2 = p^2;
   r1 = cvs2./((om/k)^2);
   for j=1:length(thk)
       r_si(j) = dns(j+1)/dns(j);
       r_en(j) = 2*(r1(j)-r_si(j)*r1(j+1));
   end
   r_a = r_si + r_en;
   r_a1 = r_a - 1;
   r_b = 1-r_en;
   r_b1 = r_b - 1;
   
   
   rk = sqrt(1-(om/k)^2./cvp2);
   sk = sqrt(1-(om/k)^2./cvs2);
   tk = 2-(om/k)^2./cvs2;

   Ca = cosh(k.*rk(1:length(thk)).*thk);
   Sa = sinh(k.*rk(1:length(thk)).*thk);
   Cb = cosh(k.*sk(1:length(thk)).*thk);
   Sb = sinh(k.*sk(1:length(thk)).*thk);
   
   X = mu(1)^2*[2*tk(1),-tk(1)^2,0,0,-4,2*tk(1)];
   
   for i=1:length(thk)
       pp1 = Cb(i)*X(2)+sk(i)*Sb(i)*X(3);
       pp2 = Cb(i)*X(4)+sk(i)*Sb(i)*X(5);
       pp3 = (1/sk(i))*Sb(i)*X(2)+Cb(i)*X(3);
       pp4 = (1/sk(i))*Sb(i)*X(4)+Cb(i)*X(5);
       
       qq1 = Ca(i)*pp1-rk(i)*Sa(i)*pp2;
       qq2 = -(1/rk(i))*Sa(i)*pp3+Ca(i)*pp4;
       qq3 = Ca(i)*pp3-rk(i)*Sa(i)*pp4;
       qq4 = -(1/rk(i))*Sa(i)*pp1+Ca(i)*pp2;
       
       yy1 = r_a1(i)*X(1)+r_a(i)*qq1;
       yy2 = r_a(i)*X(1)+r_a1(i)*qq2;
       
       zz1 = r_b(i)*X(1)+r_b1(i)*qq1;
       zz2 = r_b1(i)*X(1)+r_b(i)*qq2;
       
       X_new(1) = r_b1(i)*yy1+r_b(i)*yy2;
       X_new(2) = r_a(i)*yy1+r_a1(i)*yy2;
       X_new(3) = r_si(i)*qq3;
       X_new(4) = r_si(i)*qq4;
       X_new(5) = r_b1(i)*zz1+r_b(i)*zz2;
       X_new(6) = X_new(1);
       
       X = X_new;
   end
%        dd = real(X(2) + sk(end)*X(3) - rk(end)*(X(4)+sk(end)*X(5)));  
       dd = (X(2) + sk(end)*X(3) - rk(end)*(X(4)+sk(end)*X(5)));  
  