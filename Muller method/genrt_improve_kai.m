function [Td,Rd,e11,e12,e21,e22,du,mu,nus,nup] = genrt_improve_kai(thk,dns,cvp,cvs,om,k)

% This function calculates the E and Lambda matrices (up-going and 
% down-going matrices) for the P-SV case. Note that a separate function,
% updown, is provided for calculating the Lambda matrices for use in
% determining the displacement-stress vectors.

% Copyright 1999 by Glenn J. Rix and Carlo G. Lai

%  根据快速RT方法改进，2008，Donghong Pei，Improvement on Computation of RT method
%  2011.4 by  kai
% 根据文章：袁腊梅, 凡友华, 孙书荣. 2009. 广义反射2透射系数算法的无量纲化. 
% 地震学报, 31 (4) :37722384.将E修改后的结果

cvs2 = cvs.^2; cvp2 = cvp.^2;
mu = dns.*cvs2;
mu1 = sum(mu)/length(mu);
mu2 = mu/mu1;

e11 = zeros(2,2,length(cvs));
e12 = zeros(2,2,length(cvs));
e21 = zeros(2,2,length(cvs));
e22 = zeros(2,2,length(cvs));
du = zeros(2,2,length(thk));

if om == 0

   kappa = (1.0 + cvs2./cvp2)./(1.0 - cvs2./cvp2);
   kmu = k*mu;
   
   e11(1,1,:) = ones(1,length(cvs));
   e11(1,2,:) = e11(1,1,:);
   e12(1,1,:) = e11(1,1,:);
   e12(1,2,:) = e11(1,1,:);
   e11(2,1,:) = -(kappa - 1.0);
   e11(2,2,:) = e11(1,1,:);
   e12(2,1,:) = -e11(2,1,:);
   e12(2,2,:) = -e11(1,1,:);
   e21(1,1,:) = (kappa - 3.0).*kmu;
   e21(1,2,:) = -2*kmu;
   e22(1,1,:) = -e21(1,1,:);
   e22(1,2,:) = -e21(1,2,:);
   e21(2,1,:) = (kappa - 1.0).*kmu;
   e21(2,2,:) = -2*kmu;
   e22(2,1,:) = e21(2,1,:);
   e22(2,2,:) = e21(2,2,:);
   
   du(1,1,:) = exp(-k*thk);
   du(2,2,:) = exp(-k*thk);
   du(2,1,:) = -k*thk.*exp(-k*thk);

else
   
   k2 = k^2; om2 = om^2;
   tk = 1-(om2./k2)./(2*cvs2);
   
   ks2 = om2./cvs2;
   nus = sqrt(k2-ks2);
   index = find(imag(-1i*nus) > 0);
   nus(index) = -nus(index);
   gammas = nus/k;
   
   kp2 = om2./cvp2;
   nup = sqrt(k2-kp2);
   index = find(imag(-1i*nup) > 0);
   nup(index) = -nup(index);
   gammap= nup/k;

   chi = 2.0*k - ks2/k;
   

   e11(1,1,:) = ones(1,length(cvs));
   e11(1,2,:) = -gammas;
   e12(1,1,:) = e11(1,1,:);
   e12(1,2,:) = gammas;
   e11(2,1,:) = -gammap;
   e11(2,2,:) = e11(1,1,:);
   e12(2,1,:) = gammap;
   e12(2,2,:) = e11(1,1,:);
   e21(1,1,:) = -mu2.*gammap;
   e21(1,2,:) = mu2.*tk;
   e22(1,1,:) = -e21(1,1,:);
   e22(1,2,:) = e21(1,2,:);
   e21(2,1,:) = e21(1,2,:);
   e21(2,2,:) = -mu2.*gammas;
   e22(2,1,:) = e21(1,2,:);
   e22(2,2,:) = -e21(2,2,:);
   
   du(1,1,:) = exp(-nup(1:length(thk)).*thk);
   du(2,2,:) = exp(-nus(1:length(thk)).*thk);
   
   
   % 矩阵逆
   mofa = 1./(2*(tk-1));
   
   e11_inv(1,1,:) = -mofa;
   e11_inv(1,2,:) = -mofa.*tk./gammap;
   e12_inv(1,1,:) = mofa./(gammap.*mu2);
   e12_inv(1,2,:) = mofa./mu2;
   
   e11_inv(2,1,:) = -mofa.*tk./gammas;
   e11_inv(2,2,:) = e11_inv(1,1,:);
   e12_inv(2,1,:) =  e12_inv(1,2,:);
   e12_inv(2,2,:) = mofa./(gammas.*mu2);
   
   e21_inv(1,1,:) = e11_inv(1,1,:);
   e21_inv(1,2,:) = - e11_inv(1,2,:);
   e22_inv(1,1,:) = -e12_inv(1,1,:);
   e22_inv(1,2,:) = e12_inv(1,2,:);
   
   e21_inv(2,1,:) = -e11_inv(2,1,:);
   e21_inv(2,2,:) = e11_inv(1,1,:);
   e22_inv(2,1,:) = e12_inv(2,1,:);
   e22_inv(2,2,:) = -e12_inv(2,2,:);
   
end

 N = length(thk);
% Initialize 2x2xN matrices
% L11 = zeros(2,2);
% L12 = zeros(2,2);
% L21 = zeros(2,2);
% L22 = zeros(2,2);
Td = zeros(2,2,N);
Rd = zeros(2,2,N+1);

% Calculate the Td and Rd matrices for the Nth layer
du(1:2,1:2,N+1) = 1; % 任意值都行
Rd(1:2,1:2,N+1) = 0;

% Loop through the first N-1 layers in reverse order
for j = N:-1:1
   A = [e11_inv(:,:,j) e12_inv(:,:,j); e21_inv(:,:,j) e22_inv(:,:,j)]; 
   B = [e11(:,:,j+1) e12(:,:,j+1); e21(:,:,j+1) e22(:,:,j+1)];
   L = A*B;
   L11 = L(1:2,1:2); 
   L12 = L(1:2,3:4);
   L21 = L(3:4,1:2);  L22 = L(3:4,3:4);
   
   Td(:,:,j) = (L11+L12*du(:,:,j+1)*Rd(:,:,j+1))\du(:,:,j);
   Rd(:,:,j) = L21*Td(:,:,j) + L22*du(:,:,j+1)*Rd(:,:,j+1)*Td(:,:,j);
end



