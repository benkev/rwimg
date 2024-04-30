//
//              Kuniji Saito's Plasma Density
//
// Plasma density routine for the beam_path.c 
// The plasma density and its gradient are calculated as:
// - exponential for the chromosphere, i.e. below 9000 km over the photosphere;
// - Kuniji Saito's distribution for the corona above 11000 km height;
// - a "stitch function", smoothly merging the chromosphere and corona
//     in the transitional region between 9000 and 11000 km altitude
//

#include <stdio.h>
#include <math.h>
#include "raytrace.h"

void plasma_density(int nRay, double Pos_DI[3][nRay], 
                    double Rho_I[nRay], 
                    double GradRho_DI[3][nRay],
		    double DS_I[nRay], 
		    short Flags_I[nRay]) {
  //
  // Coronal density model by Kuniji Saito
  //
  // 
  //
  double extern const DeltaS; 
            
  double extern const cPi;               // = 3.1415926535897931;
  double extern const ProtonChargeSGSe;  // = 4.8e-10, StatCoulombs, SGSe
  double extern const cProtonMass;       // = 1.672621636E-24; //g
  double extern const cElectronMass;     // = 9.10938215E-28; //g
  double const static R0 = 6.950e5;   // km, Solar radius

  double static Pos_D[3], r2;
  int iRay, i;
  //double const static a1 = 3.09E8, a2 = 1.58E8, a3 = 0.0251E6;
  double static t1, t2, t3;
  double static x, y, z, r, rhoCyl2, rhoCyl, rhoCylm1;
  double static rm1, rm2, rm2d5, rm3d5, rm6, rm7, rm16, rm17;
  double static cosTh, sinTh, cosPh, sinPh, absCosTh;
  double static Ne, Ner, Neth, Nerho, Nex, Ney, Nez;
  double static expo, expo1, dNdr, b0, b2;
  double static a0, a1, a2;
  double static saito_rcor, pr, pt, prt, gr, gt, grt;
  int static cnt0 = 100, cnt1 = 500, cnt2 = 500, cnt3 = 500;

  // The Saito's function coefficients
  double const static  g1 = 3.09e8, g2 = 1.58e8, g3 = 2.5e4;

  //double const static r_chr_cor = 1.014388489208633;  // *R0 = R0+1e4, km

  // The merging points in upper chromosphere, r_chromo and
  // lower corona, r_corona, between which the el. density
  // Ne is smoothly replaced by a polynomial
  // The chromosphere height is supposed to be ~ 10000 km,
  // h_chrom = 10000, and the merging points are 1000km
  // away from the chromospheric boundary on both sides:
  // r_chromo = 1 + (h_chrom - 1000)/R0,
  // r_corona = 1 + (h_chrom + 1000)/R0,
  double const static r_chromo = 1.01294964029;  // *R0 = R0 + 9000, km
  double const static r_corona = 1.01582733813;  // *R0 = R0 + 11000, km
  double static delr = 0.00287769784; //= r_corona - r_chromo; *R0 = 2000 km
  double static delrm1 = 347.5; //= 1/(r_corona - r_chromo)
  double static delrm2 = 120756.2501449; //= 1/(r_corona - r_chromo)^2

  // Bits for Flags_I
  short const static Penetr = 0x0001;  // Eps < 0; Critical surface search
  short const static WasRefl = 0x0002; // First step after reflection

  for (iRay = 0; iRay < nRay; iRay++) {
    x = Pos_D[0] = Pos_DI[0][iRay];
    y = Pos_D[1] = Pos_DI[1][iRay];
    z = Pos_D[2] = Pos_DI[2][iRay];
    //r2 = dot_product(Pos_D, Pos_D);
    r2 = x*x + y*y + z*z;
    r = sqrt(r2);
    if ((cnt0--) > 0) printf("x, y, z = %g %g %g, r2 = %10.8f\n", x, y, z, r2);
    
    rm1 = 1.0/r;                 // r^(-1)
    rm2 = rm1*rm1;               // r^(-2)
    rm2d5 = pow(r,-2.5);         // r^(-2.5)
    rm3d5 = rm2d5*rm1;           // r^(-3.5)
    rm6 = pow(r,-6);             // r^(-6)
    rm7 = rm6*rm1;               // r^(-7)
    rm16 = pow(r,-16);           // r^(-16)
    rm17 = rm16*rm1;             // r^(-17)
    cosTh = fabs(z*rm1);         // cos(theta) = z/r
    sinTh = rhoCyl*rm1;          // sin(theta) = rho/r
    cosPh = x*rhoCylm1;           // cos(phi) = x/rho
    sinPh = y*rhoCylm1;           // sin(phi) = y/rho
    absCosTh = fabs(cosTh);        // abs(cos(theta)) = abs(z/r)
    
    t1 = 1.0 - 0.5*absCosTh;          // (1 - 0.5cos(theta))
    t2 = 1.0 - 0.95*absCosTh;         // (1 - 0.95cos(theta))
    t3 = 1.0 - sqrt(absCosTh); 
    
    if (r < r_chromo) { 
      
      // Chromospheric exponential density distribution
      expo = exp(-7.7e-4*(R0*(r-1) - 500));
      Ne = 5.7e11*expo;  // Ne, electron # density
      dNdr = -438.9e6*R0*expo*rm1;
      Nex = dNdr*x;                   // Ne derivative by x
      Ney = dNdr*y;                   // Ne derivative by y
      Nez = dNdr*z;                   // Ne derivative by z

      printf("Nex, Ney, Nez = %g %g %g, r = %10.8f\n", Nex, Ney, Nez);

    }

    else if (r < r_corona) { // r >= r_chromo, but still r < r_corona:
      
      // "Stitching" function 
      //    p(r,th) = a0 + a1*r + a2*r^2 + (r-r1)/(r2-r1)*
      // Find its coefficients a0, a1, and a2 solving the system
      // M*a = b
      expo1 = exp(-7.7e-4*(R0*(r_chromo-1) - 500));
      b0 = 5.7e11*expo1;  // Ne, electron # density
      b2 = -438.9e6*R0*expo1;
      a2 = -delrm2*b0 - delrm1*b2;
      a1 = -2.0*r_chromo*a2 + b2;
      a0 = b0 - r_chromo*a1 - r_chromo*r_chromo*a2;
      // Value of stitching function as number density
      saito_rcor = saito(r_corona, cosTh);
      Ne = a0 + r*(a1 + r*a2) + (r - r_chromo)*delrm1*saito_rcor;
      // Gradient components
      pr = (a1 + 2.0*a2*r + delrm1*saito_rcor)*rm1; // dp/dr*dr/dx = pr*x etc.
      pt = (r - r_chromo)*delrm1*rm2*
	((0.5*g1*pow(r_corona,-16) + 0.95*g2*pow(r_corona,-6))*cosTh 
	 + 0.5*g3*pow(r_corona,-2.5)*sqrt(cosTh)); // dp/dth*dth/dx = pt*y, x
      prt = pr + pt;
      Nex = prt*x;
      Ney = prt*y;
      Nez = prt*z - (r - r_chromo)*delrm1*rm1*
	(0.5*g1*pow(r_corona,-16) + 0.95*g2*pow(r_corona,-6)
	 + 0.5*g3*pow(r_corona,-2.5)/sqrt(cosTh));

      if (cnt2>0) {printf("r_chromo <= r < r_corona: r = %10.8f\n", r); cnt2--;}

    }

    else { // r > r_corona:
      // The Saito's density distribution
      Ne = saito(r, cosTh);
      gr = -16.0*g1*pow(r,-18)*(1.0 - 0.5*cosTh) 
	- 6.0*g2*pow(r,-8)*(1.0 - 0.95*cosTh)
	- 2.5*g3*pow(r,-4.5)*(1.0 - sqrt(cosTh));// dg/dr*dr/dx = gr*x etc.
      gt = 0.5*g1*pow(r,-18)*cosTh 
	+ 0.95*g2*pow(r,-8)*cosTh
	+ 0.5*g3*pow(r,-4.5)*sqrt(cosTh);    // dg/dth*dth/dx = gt*x etc.
      grt = gr + gt;
      Nex = grt*x;
      Ney = grt*y;
      Nez = grt*z - (0.5*g1*pow(r,-17)
		     + 0.95*g2*pow(r,-7)
		     + 0.5*g3*pow(r,-3.5)/sqrt(cosTh));
      if (cnt3>0) {printf("r >= r_corona: r = %10.8f \n", r); cnt3--;}

    }


    //
    // Play with the Ne: 
    //

    if (0) { // Ten times Ne in gramms cm^-3
      Rho_I[iRay] = 10.0*Ne*cProtonMass;
      GradRho_DI[0][iRay] = 10.0*Nex*cProtonMass;
      GradRho_DI[1][iRay] = 10.0*Ney*cProtonMass;
      GradRho_DI[2][iRay] = 10.0*Nez*cProtonMass;
    }
    else if (1) { // Just Ne in gramms cm^-3
      Rho_I[iRay] = Ne*cProtonMass;
      GradRho_DI[0][iRay] = Nex*cProtonMass;
      GradRho_DI[1][iRay] = Ney*cProtonMass;
      GradRho_DI[2][iRay] = Nez*cProtonMass;
    }
    else if (0) {  // Ne in cm^-3
      Rho_I[iRay] = Ne;
      GradRho_DI[0][iRay] = Nex;
      GradRho_DI[1][iRay] = Ney;
      GradRho_DI[2][iRay] = Nez;
    }
   

    // Revert the step to original if it is not in
    // the process of finding the critical surface
    if (!(Flags_I[nRay]&Penetr)) {
      if (r > 2.0) 
        DS_I[iRay] = DeltaS;
      else if (r > (r_corona+DeltaS)) 
        DS_I[iRay] = 0.1*DeltaS;
      else
        DS_I[iRay] = 0.01*DeltaS;
    }
  }
}  // void plasma_density()

double inline saito(double r, double cosTh) {
  // The Saito's function coefficients
  double const static  g1 = 3.09e8, g2 = 1.58e8, g3 = 2.5e4;
  double g;

  //cosTh = cos(th);
  g = g1*pow(r,-16)*(1.0 - 0.5*cosTh) + g2*pow(r,-6)*(1.0 - 0.95*cosTh)
    + g3*pow(r,-2.5)*(1.0 - sqrt(cosTh));
  return g;
}
