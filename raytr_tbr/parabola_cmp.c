#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "raytrace.h"

// parabola_cmp
//
// Trace rays ising beam_path() in 2-dimensional medium with linearly distributed dielectric permittivity
// and compare the prajectories with exact parabolic solutions
//

int  main(int argc, char* argv[])
{
  double extern const cPi;                         // = 3.1415926535897931;
  short static const FALSE = 0, TRUE = 1;
  double const dtor = cPi/180.0;                     // Degrees-to-Radians conversion constant

  int const nIter = 1000;
  int const nRay = 5;
  double DensityCr;
  double const al[] = {0.1*dtor, 1*dtor, 10*dtor, 30*dtor, 60*dtor};
  short ExcludeRay_I[nRay];   // A ray is excluded from processing if it is .true.
  double Pos_DI[3][nRay], Slope_DI[3][nRay];
  double Intensity_I[nRay], DeltaS_I[nRay];
  double Slope_D[3], SlopeCrossREarth_D[3], Pos_D[3];
  short RayFlag_I[nRay];         // 1 if a ray is OK, 0 otherwise
  short NewEntry;           // Must be set to .true. before a call with new value of nRay
  double Tolerance = 0.01, DeltaS = 1.0;
  int iIter, i, j, k, iRay, iTraj, pTraj;
  FILE *fh0, *fh1, *fh2, *fh3, *fh4, *fh5;

  fh0 = fopen("ray0.txt", "w");
  fh1 = fopen("ray1.txt", "w");
  fh2 = fopen("ray2.txt", "w");
  fh3 = fopen("ray3.txt", "w");
  fh4 = fopen("ray4.txt", "w");
  fh5 = fopen("ray5.txt", "w");

  //
  // Calculate the critical density from the frequency
  //
  //DensityCr = cPi*cProtonMass*cElectronMass*pow(RadioFrequency/ProtonChargeSGSe,2);
  DensityCr = 1.0;
  printf("rho_cr = %g\n", DensityCr);

  for (i = 0; i < nRay; i++) {
    Intensity_I[i] = 0.0;
    RayFlag_I[i] = FALSE;
    ExcludeRay_I[i] = FALSE;
    DeltaS_I[i] = DeltaS;
    Pos_DI[0][i] = .0;
    Pos_DI[1][i] = .0;
    Pos_DI[2][i] = .0;
    Slope_DI[0][i] = cos(al[i]);
    Slope_DI[1][i] = sin(al[i]);
    Slope_DI[2][i] = 0.0;
    fprintf(fh5, "%g %g %g\n", cos(al[i]), sin(al[i]));
  } // for

  fprintf(fh0, "%g %g %g\n", Pos_DI[0][0], Pos_DI[1][0], Pos_DI[2][0]);
  fprintf(fh1, "%g %g %g\n", Pos_DI[0][1], Pos_DI[1][1], Pos_DI[2][1]);
  fprintf(fh2, "%g %g %g\n", Pos_DI[0][2], Pos_DI[1][2], Pos_DI[2][2]);
  fprintf(fh3, "%g %g %g\n", Pos_DI[0][3], Pos_DI[1][3], Pos_DI[2][3]);
  fprintf(fh4, "%g %g %g\n", Pos_DI[0][4], Pos_DI[1][4], Pos_DI[2][4]);
 
  Tolerance = 0.1;
  NewEntry = 1;
  pTraj = 0;

 
 
  for (iIter = 0; iIter < nIter; iIter++) {
    if (!(iIter % 100)) {
      printf("Iter %i\n", iIter);
      printf("X, Y, Z: %g %g %g \n", Pos_DI[0][0], Pos_DI[1][0], Pos_DI[2][0]);
      printf("Vx, Vy, Vz: %g %g %g \n", Slope_DI[0][0], Slope_DI[1][0], Slope_DI[2][0]);
    }

    //beam_path(void (* Get_Plasma_Density)(), int nRay, short ExcludeRay_I[nRay], double Pos_DI[3][nRay],
    //          double Slope_DI[3][nRay], double DeltaS_I[nRay], double ToleranceInit, double DensityCr, 
    //          double Intensity_I[nRay], short RayFlag_I[nRay], short *NewEntry)
    beam_path(plasma_density, nRay, ExcludeRay_I, Pos_DI, Slope_DI, DeltaS_I,
	      Tolerance, DensityCr, Intensity_I, RayFlag_I, &NewEntry);

    for (iRay = 0; iRay < nRay; iRay++) {
      if (Pos_DI[0][iRay] < -20.0) ExcludeRay_I[iRay] = TRUE;
      if (RayFlag_I[iRay]) printf("+++ Bad ray %i", iRay);
    }

    if (!ExcludeRay_I[0]) fprintf(fh0, "%g %g %g\n", Pos_DI[0][0], Pos_DI[1][0], Pos_DI[2][0]);
    if (!ExcludeRay_I[1]) fprintf(fh1, "%g %g %g\n", Pos_DI[0][1], Pos_DI[1][1], Pos_DI[2][1]);
    if (!ExcludeRay_I[2]) fprintf(fh2, "%g %g %g\n", Pos_DI[0][2], Pos_DI[1][2], Pos_DI[2][2]);
    if (!ExcludeRay_I[3]) fprintf(fh3, "%g %g %g\n", Pos_DI[0][3], Pos_DI[1][3], Pos_DI[2][3]);
    if (!ExcludeRay_I[3]) fprintf(fh4, "%g %g %g\n", Pos_DI[0][4], Pos_DI[1][4], Pos_DI[2][4]);

  }

  fclose(fh0);
  fclose(fh1);
  fclose(fh2);
  fclose(fh3);
  fclose(fh4);
  fclose(fh5);

  //
  // Free internal memory
  //
  NewEntry = -1;
  beam_path(plasma_density, nRay, ExcludeRay_I, Pos_DI, Slope_DI, DeltaS_I,
	    Tolerance, DensityCr, Intensity_I, RayFlag_I, &NewEntry);

} // int main()

// ========================================================== 

//
// This program returns linearly increasing plasma density along the X axis
// and its gradient, which is constant
// 
void plasma_density(int nRay, double Pos_DI[3][nRay], double Density_I[nRay], double GradDensity_DI[3][nRay],
		    double DeltaS_I[nRay], short RayFlag_I[nRay]) {
  int iRay;
  
  double Ldist = 100.0;               // // L = 100; Distance from x=0 to the critical surface
  //double DensityCr = 2.07756e-18;
  double DensityCr = 1.0;
  double coef = DensityCr/Ldist;

  for (iRay = 0; iRay < nRay; iRay++) {
    Density_I[iRay] = coef*Pos_DI[0][iRay];         // rho = (rho_cr/L)*x
    GradDensity_DI[0][iRay] = DensityCr/Ldist;      // nabla rho = |rho_cr/L; 0; 0|
    GradDensity_DI[1][iRay] = 0.0;
    GradDensity_DI[2][iRay] = 0.0;
    DeltaS_I[iRay] = 1.0;                           // Revert the step to original
  }
}  // void plasma_density()

