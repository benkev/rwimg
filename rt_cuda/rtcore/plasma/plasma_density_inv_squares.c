#include <stdio.h>
#include <math.h>
#include "raytrace.h"

void plasma_density_inv_squares(int nRay, 
				double PhysConst_I[],
				double Pos_ID[nRay][3], 
				double Rho_I[nRay], 
				double GradRho_ID[nRay][3], 
				double DS_I[nRay], 
				short Flags_I[nRay]) {
  //
  // DensityAtSolarSurface = Ne(@SolarSurface)*ProtonMass 
  //                       = 2x10^8(cm^-3)*1.6726x10^-24(g) = 3.3452e-16(g/cm^3)
  //
  double DeltaS_Rsun =      PhysConst_I[0];
  /* double ProtonChargeCGSe = PhysConst_I[1]; */
  /* double ProtonMass_g =     PhysConst_I[2]; */
  /* double ElectronMass_g =   PhysConst_I[3]; */
  /* double Rsun_cm =          PhysConst_I[4]; */
  /* double Rsun_km =          PhysConst_I[5];  */
  double RhoSurf =          PhysConst_I[6]; /* = 2.5*3.3452E-16 g/cm^3 */
  /* double h_chromo_km =      PhysConst_I[7]; */
  /* double Te_corona_K =      PhysConst_I[8]; */
  /* double Te_chromo_K =      PhysConst_I[9]; */
  /* double AU_m =             PhysConst_I[10]; */
  /* double AU_Rsun =          PhysConst_I[11]; */

  double Pos_D[3], r2, r4;
  int iRay;

  for (iRay = 0; iRay < nRay; iRay++) {
    Pos_D[0] = Pos_ID[iRay][0];
    Pos_D[1] = Pos_ID[iRay][1];
    Pos_D[2] = Pos_ID[iRay][2];
    r2 = dot_product(Pos_D, Pos_D); /* Solar distance squared */
    Rho_I[iRay] = RhoSurf/r2;
    r4 = pow(r2,2);    /* Solar distance to the power of four */
    GradRho_ID[iRay][0] = -2.*RhoSurf*Pos_D[0]/r4;
    GradRho_ID[iRay][1] = -2.*RhoSurf*Pos_D[1]/r4;
    GradRho_ID[iRay][2] = -2.*RhoSurf*Pos_D[2]/r4;
    DS_I[iRay] = DeltaS_Rsun;       // Revert the step to original
  }

} 

/* // */
/* // Scalar Product of 3-Dimensional Vectors a and b */
/* // */
/* inline double dot_product(double a[3], double b[3]) { */
/*   return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]; */
/* }  */
