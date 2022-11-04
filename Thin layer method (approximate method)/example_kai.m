%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% PROGRAM:
% numerical_tests.m
%
% PROGRAMMERS:
% Matt Haney and Victor Tsai
%
% Last revision date:
% 26 April 2017
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This code is distributed as part of the source-code package 
%                   raylee_inversion_codes 
% that accompanies Haney and Tsai (2017). The package can be downloaded 
% from the Geophysics source-code archive at 
%                   http://software.seg.org/2017/0003/index.html
% Use of this code is subject to acceptance of the terms and conditions
% that can be found at http://software.seg.org/disclaimer.txt 
% Copyright 2017 by The Society of Exploration Geophysicists (SEG)
% Reference:
% Haney, M. M., Tsai, V. C. (2017) Perturbational and nonperturbational 
% inversion of Rayleigh-wave velocities, Geophysics, 82(3), F15-F28.
% doi: 10.1190/geo2016-0397.1
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Program numerical_tests is a Matlab script that compares the Jacobian 
% matrix computed with the Raylee codes to Jacobian matrices shown in 
% Table 2 of Cercato (2007, GJI). It also computes the derivative of 
% phase velocity with respect to the thickness of a crustal layer using 
% perturbation theory and compares the result to the brute force method 
% of making the crust slightly thinner and slightly thicker and taking  
% finite differences.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% this is a script
% modified by Kai Zhang, 2022/11/4
%  This code is used to compare with the muller method.

clear;

vp=[800 1200];
vs=[200 400];
dns=[2 2];
Qs=[10 10];
Qp=[10 10];
vp=vp.*(1+1i./(2*Qp));
vs=vs.*(1+1i./(2*Qs));

% grid spacings (non-uniform) in meters
h = [0.1*ones(1,100) 2*ones(1,100)];
% number of nodes
Nn = length(h);
% frequencies (Hz)
fks =1:50; 
% mode number for each frequency
% 1=fundamental, 2=first higher mode, etc
modnv = ones(1,10);

% the Novotny crust/mantle model in Cercato
vpv = zeros(1,Nn);
vsv = zeros(1,Nn);
rhov = zeros(1,Nn);

% layer 1
vpv(1:100) = vp(1);
vsv(1:100) = vs(1);
rhov(1:100) = dns(1);

% layer 2
vpv(101:end) = vp(2);
vsv(101:end) = vs(2);
rhov(101:end) = dns(2);


% no fluid layer on top of this model
Nnf = 0; vpfv = 1; rhofv =1; hfv = 1;

% compute phase and group velocities
countr = 0;
modn = 1:5;
tic

for f=fks
    
    kk = raylee_lysmer_kai(Nn,vsv,vpv,rhov,f,h,modn,Nnf,vpfv,rhofv,hfv);
    countr = countr + 1;
    vpp2(countr,:) = kk;
end

cr_real_thin=fks'*2*pi./real(vpp2);

cr_imag_thin=imag(vpp2);

cr_imag_thin(cr_real_thin>max(real(vs)))=0;
cr_real_thin(cr_real_thin>max(real(vs)))=0;


cr_real_thin(cr_real_thin==0)=NaN;
cr_imag_thin(cr_imag_thin==0)=NaN;

toc
