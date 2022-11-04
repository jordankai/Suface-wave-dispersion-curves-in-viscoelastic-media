%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% PROGRAM:
% raylee_lysmer.m
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
% Program raylee_lysmer is a Matlab function that computes Rayleigh/Scholte 
% wave phase and group velocities and mode shapes for a given model based 
% on the finite element method of Lysmer (BSSA, 1970). It performs the 
% forward problem at a single frequency.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Input:
% Nn            number of elements in solid part of model
% Nnf           number of elements in fluid part of model
% hv            vector of grid spacings for solid (meters)
% hfv           vector of grid spacings for fluid (meters)
% f             frequency (Hz)
% modn          which mode (1=fundamental, 2=first overtone, etc)
% vsv           shear velocity model in solid, a vector (m/s)
% vpv           compressional velocity model in solid, a vector (m/s)
% rhov          density model in solid, a vector (kg/m^3)
% vpfv          compressional velocity model in fluid, a vector (m/s)
% rhofv         density model in fluid, a vector (kg/m^3)
%
% Output:
% kk            wavenumber for the Rayleigh wave at this frequency
% vpk           phase velocity for the Rayleigh wave at this frequency
% vgk           group velocity for the Rayleigh wave at this frequency
% ev            vertical and horizontal displacement eigenfunctions 
%               (mode shapes), note these are scrambled
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  kk = ... 
    raylee_lysmer_kai(Nn,vsv,vpv,rhov,f,hv,modn,Nnf,vpfv,rhofv,hfv)



% the number of nodes in the fluid, based on the number of elements
if (Nnf > 0)
    Nnfo = Nnf + 1;
else
    Nnfo = 0;
end

% make fluid portion of model

% make kappaf, the modulus
kappafv = rhofv.*vpfv.*vpfv;

% make angular frequency
omga = 2*pi*f;

% initialize some local matrices
L1 = sparse(2,2);
L3 = sparse(2,2);
M1 = sparse(2,2);

% initialize the global matrix
Ka1 = sparse(Nnfo+(2*Nn),Nnfo+(2*Nn));
Ka3 = sparse(Nnfo+(2*Nn),Nnfo+(2*Nn));
M = sparse(Nnfo+(2*Nn),Nnfo+(2*Nn));

% for all elements
for ii=1:Nnf
    
    % grab grid interval of current element
    h = hfv(ii);
    
    % grab material properties of current element
    rhof = rhofv(ii);
    kappaf = kappafv(ii);
    
    % make elemental mass matrix
    M1 = sparse(2,2);
    M1(1,1) = h/(2*kappaf);
    M1(2,2) = h/(2*kappaf);
    
    % make elemental stiffness matrices
    L1 = sparse(2,2);
    L3 = sparse(2,2);

    % some alternate variables from Lysmer
    alph = 1/(6*rhof);
    bet = 1/(6*rhof);
    
    % the 4 entries of the 2x2 elemental stiffness matrices of Lysmer

    L1(1,1) = 2*alph*h; 
    L3(1,1) = (6*bet/h);
    
    L1(1,2) = alph*h;
    L3(1,2) = -(6*bet/h);
    
    L1(2,1) = L1(1,2);
    L3(2,1) = L3(1,2);
    
    L1(2,2) = L1(1,1);
    L3(2,2) = L3(1,1);
    
    % assemble mass and stiffness matrices from elemental matrices
    M((1*(ii-1)+1):(1*(ii+1)),(1*(ii-1)+1):(1*(ii+1))) = ...
        M((1*(ii-1)+1):(1*(ii+1)),(1*(ii-1)+1):(1*(ii+1))) + M1;
    Ka1((1*(ii-1)+1):(1*(ii+1)),(1*(ii-1)+1):(1*(ii+1))) = ...
        Ka1((1*(ii-1)+1):(1*(ii+1)),(1*(ii-1)+1):(1*(ii+1))) + L1;
    Ka3((1*(ii-1)+1):(1*(ii+1)),(1*(ii-1)+1):(1*(ii+1))) = ...
        Ka3((1*(ii-1)+1):(1*(ii+1)),(1*(ii-1)+1):(1*(ii+1))) + L3;
    
end

M(1,1) = M(1,1)*2;
Ka1(1,1) = Ka1(1,1)*2;
Ka3(1,1) = Ka3(1,1)*2;

% make solid portion of model

% make mu and lambda
muv = rhov.*vsv.*vsv;
lamdav = rhov.*vpv.*vpv - 2*muv;

%initialize some matrices
Ka2 = sparse(Nnfo+(2*Nn),Nnfo+(2*Nn));

L1 = sparse(4,4);
L2 = sparse(4,4);
L3 = sparse(4,4);
M1 = sparse(4,4);

% for all elements
for ii=1:Nn
    
    % grab grid interval of current element
    h = hv(ii);
    
    % grab material properties of current element
    mu = muv(ii);
    lamda = lamdav(ii);
    
    % make elemental mass matrix
    M1 = sparse(4,4);
    M1(1,1) = h*rhov(ii)/2;
    M1(2,2) = h*rhov(ii)/2;
    M1(3,3) = h*rhov(ii)/2;
    M1(4,4) = h*rhov(ii)/2;
    
    % make elemental stiffness matrices
    L1 = sparse(4,4);
    L2 = sparse(4,4);
    L3 = sparse(4,4);
    
    % some alternate variables from Lysmer
    alph = ((2*mu)+lamda)/6;
    bet = mu/6;
    theta = (mu+lamda)/4;
    psi = (mu-lamda)/4;
    
    % the 16 entries of the 4x4 elemental stiffness matrices of Lysmer
    
    L1(1,1) = 2*alph*h; 
    L3(1,1) = (6*bet/h);
    
    L2(1,2) = 2*psi;
    
    L1(1,3) = alph*h;
    L3(1,3) = -(6*bet/h);
    
    L2(1,4) = 2*theta;
    
    L2(2,1) = L2(1,2);
    
    L1(2,2) = 2*bet*h;
    L3(2,2) = (6*alph/h);
    
    L2(2,3) = -2*theta;
    
    L1(2,4) = bet*h;
    L3(2,4) = -(6*alph/h);
    
    L1(3,1) = L1(1,3);
    L3(3,1) = L3(1,3);
    
    L2(3,2) = L2(2,3);
    
    L1(3,3) = L1(1,1);
    L3(3,3) = L3(1,1);
    
    L2(3,4) = -2*psi;
    
    L2(4,1) = L2(1,4);
    
    L1(4,2) = L1(2,4);
    L3(4,2) = L3(2,4);
    
    L2(4,3) = L2(3,4);
    
    L1(4,4) = L1(2,2);
    L3(4,4) = L3(2,2);
    
    % assemble mass and stiffness matrices from elemental matrices
    if (ii == Nn)
    M((Nnfo+(2*(ii-1)+1)):(Nnfo+(2*ii)),(Nnfo+(2*(ii-1)+1)):(Nnfo+(2*ii))) = ...
        M((Nnfo+(2*(ii-1)+1)):(Nnfo+(2*ii)),(Nnfo+(2*(ii-1)+1)):(Nnfo+(2*ii))) + M1(1:2,1:2);
    Ka1((Nnfo+(2*(ii-1)+1)):(Nnfo+(2*ii)),(Nnfo+(2*(ii-1)+1)):(Nnfo+(2*ii))) = ...
        Ka1((Nnfo+(2*(ii-1)+1)):(Nnfo+(2*ii)),(Nnfo+(2*(ii-1)+1)):(Nnfo+(2*ii))) + L1(1:2,1:2);
    Ka2((Nnfo+(2*(ii-1)+1)):(Nnfo+(2*ii)),(Nnfo+(2*(ii-1)+1)):(Nnfo+(2*ii))) = ...
        Ka2((Nnfo+(2*(ii-1)+1)):(Nnfo+(2*ii)),(Nnfo+(2*(ii-1)+1)):(Nnfo+(2*ii))) + L2(1:2,1:2);
    Ka3((Nnfo+(2*(ii-1)+1)):(Nnfo+(2*ii)),(Nnfo+(2*(ii-1)+1)):(Nnfo+(2*ii))) = ...
        Ka3((Nnfo+(2*(ii-1)+1)):(Nnfo+(2*ii)),(Nnfo+(2*(ii-1)+1)):(Nnfo+(2*ii))) + L3(1:2,1:2);
    else
    M((Nnfo+(2*(ii-1)+1)):(Nnfo+(2*(ii+1))),(Nnfo+(2*(ii-1)+1)):(Nnfo+(2*(ii+1)))) = ...
        M((Nnfo+(2*(ii-1)+1)):(Nnfo+(2*(ii+1))),(Nnfo+(2*(ii-1)+1)):(Nnfo+(2*(ii+1)))) + M1;
    Ka1((Nnfo+(2*(ii-1)+1)):(Nnfo+(2*(ii+1))),(Nnfo+(2*(ii-1)+1)):(Nnfo+(2*(ii+1)))) = ...
        Ka1((Nnfo+(2*(ii-1)+1)):(Nnfo+(2*(ii+1))),(Nnfo+(2*(ii-1)+1)):(Nnfo+(2*(ii+1)))) + L1;
    Ka2((Nnfo+(2*(ii-1)+1)):(Nnfo+(2*(ii+1))),(Nnfo+(2*(ii-1)+1)):(Nnfo+(2*(ii+1)))) = ...
        Ka2((Nnfo+(2*(ii-1)+1)):(Nnfo+(2*(ii+1))),(Nnfo+(2*(ii-1)+1)):(Nnfo+(2*(ii+1)))) + L2;
    Ka3((Nnfo+(2*(ii-1)+1)):(Nnfo+(2*(ii+1))),(Nnfo+(2*(ii-1)+1)):(Nnfo+(2*(ii+1)))) = ...
        Ka3((Nnfo+(2*(ii-1)+1)):(Nnfo+(2*(ii+1))),(Nnfo+(2*(ii-1)+1)):(Nnfo+(2*(ii+1)))) + L3;
    end
    
end

% construct the coupling matrix
if (Nnf > 0)
Cm = sparse(Nnfo+(2*Nn),Nnfo+(2*Nn));
Cm(Nnfo,Nnfo+2) = 1;
Cm(Nnfo+2,Nnfo) = 1;
else
    Cm = sparse(Nnfo+(2*Nn),Nnfo+(2*Nn));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       
% find the rayleigh/scholte wave speed which would exist if the solid model 
% were a halfspace with the minimum model velocity
%
% this is a lower bound on the velocity that can be passed to EIGS, based 
% on ARPACK
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (Nnf > 0)
    [msval msloc] = min(vsv);
    [mfval mfloc] = min(vpfv);
    vsmay = msval;
    vpmay = vpv(msloc);
    vpfmay = mfval;
    rhofmay = rhofv(mfloc);
    rhomay = rhov(msloc);
    rspd = stoneley_vel(vpmay,vsmay,vpfmay,rhofmay,rhomay);
else
    [msval msloc] = min(vsv);
    vsmay = msval;
    vpmay = vpv(msloc);
    % coefficients of rayleigh's polynomial
    t1 = 1/(vsmay^6);
    t2 = -8/(vsmay^4);
    t3 = ((24/(vsmay^2))-(16/(vpmay^2)));
    t4 = -16*(1-((vsmay/vpmay)^2));
    % rayleigh wave speed
    rspd = sqrt(min(roots([ t1 t2 t3 t4 ]))); 
end
    

% find the eigenvalue closest to the upper-bound eigenvalue

for ii=1:length(modn)
    mn=modn(ii);
    opts.disp = 0;
    [xp,dp]=eigs([sparse((Nnfo+(2*Nn)),(Nnfo+(2*Nn))) speye((Nnfo+(2*Nn)),(Nnfo+(2*Nn))); ((omga*omga*M)-Ka3-(omga*Cm)) Ka2],...
        [speye((Nnfo+(2*Nn)),(Nnfo+(2*Nn))) sparse((Nnfo+(2*Nn)),(Nnfo+(2*Nn))); sparse((Nnfo+(2*Nn)),(Nnfo+(2*Nn))) Ka1],...
        mn,omga/rspd,opts);
    % pick the mode of interest
    x = xp(:,mn);
    d = dp(mn,mn);
    kk(1,ii) = d;
end









