
 thk=10;
 dns=[2 2];
 vs=[200 400];
 vp=[800 1200];
 Qs=[10 10];
 Qp=[10 10];

 freq =1:50; 

[cr,cr_real,cr_imag] = Rayleigh_DC(freq,thk,dns,vs,vp,Qs,Qp);

figure
plot(freq,cr_real)
xlabel('Frequency (Hz)')
ylabel('Phase velocity (m/s)')

figure
plot(freq,cr_imag)
xlabel('Frequency (Hz)')
ylabel('Attenuation coefficient')








