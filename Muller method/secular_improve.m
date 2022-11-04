function dd = secular_improve(k,om,thk,dns,cvp,cvs)

% This function calculates the absolute value of the secular function for
% a particular frequency and wavenumber.

% Copyright 1999 by Glenn J. Rix and Carlo G. Lai

% Check to see if the trial phase velocity is equal to the shear wave velocity
% or compression wave velocity of one of the layers
epsilon = 0.0001;
while any(abs(om/k-cvs)<epsilon) || any(abs(om/k-cvp)<epsilon)
   k = k * (1+epsilon);
end   

[Td,Rd,e11,e12,e21,e22,du,mu,nus,nup] = genrt_improve_kai(thk,dns,cvp,cvs,om,k);

% Note that the absolute value of the secular function is calculated

%  dd = (det(e21(:,:,1) + e22(:,:,1)*du(:,:,1)*Rd(:,:,1))/(nus(1)*nup(1)*mu(1)^2));
 
 dd = det(e21(:,:,1) + e22(:,:,1)*du(:,:,1)*Rd(:,:,1));
