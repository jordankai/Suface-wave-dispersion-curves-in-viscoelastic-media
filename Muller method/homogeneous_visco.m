function cvr1 = homogeneous_visco(cvp,cvs)

% This function calculates the Rayleigh phase velocity in a homogeneous
% viscoelastic half space

% nu = 0.5*((cvp*cvp-2*cvs*cvs)/(cvp*cvp-cvs*cvs));

% Define Coefficients of Rayleigh's Equation
a =  1;
b = -8;
c =  8*(3-2*(cvs.*cvs)./(cvp.*cvp));
d = 16*((cvs.*cvs)./(cvp.*cvp)-1);
cvr=zeros(3,length(cvs));
% Solution of Rayleigh Equation
for ii=1:length(cvs)
    p   = [a b c(ii) d(ii)];
    x   = roots(p);
    cr  = cvs(ii)*sqrt(x);
    cvr(1:length(cr),ii)=1./real(1./cr);
    cvr1(ii)=cr(3);
end
    %cvr1=cr(3,:);
    

% % Determine which of the roots is correct using the estimated velocity (Achenbach, 1973)
% crest = cvs*((0.862+1.14*nu)/(1+nu));
% index = find(abs(cr-crest) == min(abs(cr-crest)));
% cvr = cr(index);
% if isempty(cvr)
%    error('No root found for homogeneous half space')
% end