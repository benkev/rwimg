#include <stdio.h>
#include <math.h>
#include "raytrace.h"

//
// Allen-Baumbach's plasma density, stitched with the 
// chromospheric exponential distribution
//

void plasma_density(int nRay, double Pos_DI[3][nRay], 
		    double Density_I[nRay], 
		    double GradDensity_DI[3][nRay],
		    double DS_I[nRay], 
		    short Flags_I[nRay]) {
  //
  // After S. F. Smerd, 1945
  // Chromospheric density model N = 10^8(1.55rho^-6 + 2.99rho^-16) cm^-3
  //                        rho = R/R0, R0 = 6.95e5 km is the solar radius
  // Coronal density model by Baumbach
  // 
  //double extern const cPi;              // = 3.1415926535897931;
  //double extern const ProtonChargeSGSe; // = 4.8e-10; // StatCoulombs, SGSe
  double extern const cProtonMass;      // = 1.672621636E-24; //g
  //double extern const cElectronMass;    // = 9.10938215E-28; //g
  double extern const DeltaS;             

  int iRay, i;
  double const static R0 = 6.950e5;   // km, Solar radius
  // The merging points in upper chromosphere, r_chromo and
  // lower corona, r_corona, between which the el. density
  // Ne is smoothly replaced by a polynomial
  // The chromosphere height is supposed to be ~ 10000 km,
  // h_chrom = 10000, and the merging points are 1000km
  // away from the chromospheric boundary on both sides:
  // r_chromo = 1 + (h_chrom - 1000)/R0,
  // r_corona = 1 + (h_chrom + 1000)/R0,
  double const static r_chromo = 1.01294964029;  // *R0, km
  double const static r_corona = 1.01582733813; // *R0, km
  // The "sewing" polynomial coefficients a0..a3, where
  //    p(r) = a0 + a1*r + a2*r^2 + a3*r^3, are
  //double const static a0 = 1.68684609e+16, a1 = -4.98108520e+16;
  //double const static a2 = 4.90288005e+16, a3 = -1.60863438e+16;
  double const static 
    a0 = 1.686846091695834800e+16, 
    a1 = -4.981085198383761600e+16,
    a2 =  4.902880054853674400e+16, 
    a3 = -1.608634376506040800e+16;

  double static x, y, z, r, r2;
  double static rm1, rm2, rm6, rm8, rm16, rm18;
  double static Ne, dNdr, expo;
  short deb = 0;
  // Bits for Flags_I
  short const static Penetr = 0x0001;  // Eps < 0; Critical surface search
  short const static WasRefl = 0x0002; // First step after reflection

  for (iRay = 0; iRay < nRay; iRay++) {
    x = Pos_DI[0][iRay];
    y = Pos_DI[1][iRay];
    z = Pos_DI[2][iRay];
    r2 = x*x + y*y + z*z;
    r = sqrt(r2); // Distance from the centre of Sun       
    
    if (r < r_chromo) {
      expo = exp(-7.7e-4*(R0*(r-1) - 500));
      Ne = 5.7e11*expo;
      dNdr = -438.9e6*R0*expo/r;
    }
    
    else if (r < r_corona) { // but r >= r_chromo
      Ne = a0 + r*(a1 + r*(a2 + r*a3));
      // dNdr = (1/r)*dNe/dr
      dNdr = a1/r + 2.0*a2 + 3*a3*r;
    }
    
    else { // r >= r_corona:
      // rm1 = 1/r;                   // r^(-1)
      rm2 = 1/r2;                  // r^(-2)
      rm6 = pow(r,-6);             // r^(-6)
      rm8 = rm6*rm2;               // r^(-8)
      rm16 = pow(r,-16);           // r^(-16)
      rm18 = rm16*rm2;             // r^(-18)
      Ne = 1e8*(1.55*rm6 + 2.99*rm16);
      // dNdr = (1/r)*dNe/dr
      dNdr = -1e8*(9.3*rm8 + 47.84*rm18);
    }
    
    Density_I[iRay] = Ne*cProtonMass;
    GradDensity_DI[0][iRay] = dNdr*x*cProtonMass;
    GradDensity_DI[1][iRay] = dNdr*y*cProtonMass;
    GradDensity_DI[2][iRay] = dNdr*z*cProtonMass;
    
    // Revert the step to original if it is not in
    // the process of finding the critical surface
    if (!(Flags_I[nRay]&Penetr)) {
      //DS_I[iRay] = DeltaS;
      if (r > 2.0) 
	DS_I[iRay] = DeltaS;
      else if (r > (r_corona+DeltaS)) 
	DS_I[iRay] = 0.1*DeltaS;
      else
	DS_I[iRay] = 0.01*DeltaS;
    }
  } 
  
}  // void plasma_density()


