#include <stdio.h>
#include <math.h>
#include "raytrace.h"



void plasma_density(int nRay, 
		    double Pos_DI[3][nRay], 
		    double Density_I[nRay], 
		    double GradDensity_DI[3][nRay],
		    double DS_I[nRay], 
		    short Flags_I[nRay]) {
  //
  // Plasma density and its gradient at specified locations, Pos_DI
  //
  // DensityAtSolarSurface = Ne(@SolarSurface)*ProtonMass 
  //                       = 2x10^8(cm^-3)*1.6726x10^-24(g) = 3.3452e-16(g/cm^3)
  //
  double extern const cPi;              // = 3.1415926535897931;
  double extern const ProtonChargeSGSe; // = 4.8e-10; // StatCoulombs, SGSe
  double extern const cProtonMass;      // = 1.672621636E-24; //g
  double extern const cElectronMass;    // = 9.10938215E-28; //g
  double extern const DeltaS;

  //double const static OneAU = 215.0; // solar radii
  double const static DensityAtSolarSurface = 2.5*3.3452E-16;    // g/cm^3
  double Pos_D[3], SolarDistSqr, SolarDistQuad;
  int iRay, i;
  short deb = 0;

  for (iRay = 0; iRay < nRay; iRay++) {
    Pos_D[0] = Pos_DI[0][iRay];
    Pos_D[1] = Pos_DI[1][iRay];
    Pos_D[2] = Pos_DI[2][iRay];
    SolarDistSqr = dot_product(Pos_D, Pos_D);
    Density_I[iRay] = DensityAtSolarSurface/SolarDistSqr;
    SolarDistQuad = SolarDistSqr*SolarDistSqr;
    GradDensity_DI[0][iRay] = -2.*DensityAtSolarSurface*Pos_D[0]/SolarDistQuad;
    GradDensity_DI[1][iRay] = -2.*DensityAtSolarSurface*Pos_D[1]/SolarDistQuad;
    GradDensity_DI[2][iRay] = -2.*DensityAtSolarSurface*Pos_D[2]/SolarDistQuad;
    DS_I[iRay] = DeltaS;       // Revert the step to original
  }

}  // void plasma_density()


