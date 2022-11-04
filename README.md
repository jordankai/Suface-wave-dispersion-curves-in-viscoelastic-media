# Suface-wave-dispersion-curves-in-viscoelastic-media
I recode the muller method for calculating the Rayleigh wave dispersion curves in viscoelastic media.
It is simple and can work well in most models. 
I also give an example from the thin layer method written by Matt Haney and Victor Tsai.
Generally, both the two methods show good agreements.
Since the thin layer method is a finite element method, its results are approximate. 
When we consider higher modes, the thin layer method will probably be not very suitable while the muller
method provide exact surface wave eignvalues in viscoealstic layered half space.


Here are some other methods or articles on this topic:

1. M. J. S. Lowe, "Matrix techniques for modeling ultrasonic waves in multilayered media," in IEEE Transactions on Ultrasonics, Ferroelectrics, and Frequency Control, 
   vol. 42, no. 4, pp. 525-542, July 1995.    http://www.disperse.software

2. Carlo G. Lai, Glenn J. Rix; Solution of the Rayleigh Eigenproblem in Viscoelastic Media. Bulletin of the Seismological Society of America 2002;; 92 (6): 2297–2309.
    Lai et.al. solved the Rayleigh wave eignvalue problem in viscoelastic media based on the argument principle from the calculus of residues.

3. Michael Wayne Morrison, IN-SITU VISCOELASTIC SOIL PARAMETER ESTIMATION USING LOVE WAVE INVERSION. In the end of his dissertation, there are some matlab codes used for
   solve the Love wave eignvalue problem in viscoelastic media which is based on the argument principle from the calculus of residues. 

4. Cao, Z., and H. Dong (2010). Attenuation dispersion of Love waves in a viscoelastic multilayered half-space, 2010 SEG Annual Meeting,OnePetro.  

5. Orta, Adil Han, et al. "A comparative study for calculating dispersion curves in viscoelastic multi-layered plates." Composite Strutures (2022): 115779.
   matlab app available:  https://github.com/adilorta/The-Dispersion-Box 

6. Lei Pan, Shichuan Yuan, Xiaofei Chen; Modified Generalized R/T Coefficient Method for Surface‐Wave Dispersion‐Curve Calculation in Elastic and Viscoelastic Media. 
   Bulletin of the Seismological Society of America 2022;; 112 (5): 2280–2296. doi: https://doi.org/10.1785/0120210294


If you can share more other methods (codes, maybe more efficient and more rigorous mathematical processing) 
for solving the surface wave eignvalues in viscoealstic layered half space, I would be very grateful for your contact and email.

Zhang Kai  naturekai@126.com  2022/11/5
