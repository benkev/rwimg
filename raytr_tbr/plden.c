#include <stdio.h>
#include <math.h>
#include "raytrace.h"

double extern const cProtonMass;      // = 1.672621636E-24; //g

int main() {
 
  int i, j, nRay = 1;
  double pos[3][1] = {1.1, -1.0, 0.0};
  //double pos[3][1] = {0.57735026919, 0.57735026919, 0.57735026919};
  double p1[3][1], p2[3][1], grdn[3], grd[3], den1[1], den2[1];
  double dx = 1e-6, dy = 1e-6, dz = 1e-6;
  double den[1], gradd[3][1], ds[1], p[3], r;
  short rflg[1];
  FILE *fh, *fh1;

  int static cnt0 = 1000;

//void plasma_density(int nRay, double Position_DI[3][nRay], 
//		    double Density_I[nRay], 
//		    double GradDensity_DI[3][nRay],
//		    double DeltaS_I[nRay], 
//		    short RayFlag_I[nRay])

  fh = fopen("den.txt", "w");
  fh1 = fopen("dengrad.txt", "w");
  for (i = 0; i < 50001; i++) {

    plasma_density(nRay, pos, den, gradd, ds, rflg);

    for (j = 0; j < 3; j++) { grd[j] = gradd[j][0]; } // grd=grdd

    //    printf("pos[x]= %g, N= %g, gradd[x]= %g\n", 
    //           pos[0][0], den[0]/cProtonMass, gradd[0][0]);

    // Numerically calculate gradient
    for (j = 0; j < 3; j++) { p1[j][0] = p2[j][0] = pos[j][0]; } // p1=p2=pos
    p1[0][0] = pos[0][0] - dx;
    p2[0][0] = pos[0][0] + dx;
    plasma_density(nRay, p1, den1, gradd, ds, rflg);
    plasma_density(nRay, p2, den2, gradd, ds, rflg);
    grdn[0] = (den2[0] - den1[0])/(2.0*dx);

    for (j = 0; j < 3; j++) { p1[j][0] = p2[j][0] = pos[j][0]; } // p1=p2=pos
    p1[1][0] = pos[1][0] - dy;
    p2[1][0] = pos[1][0] + dy;
    plasma_density(nRay, p1, den1, gradd, ds, rflg);
    plasma_density(nRay, p2, den2, gradd, ds, rflg);
    grdn[1] = (den2[0] - den1[0])/(2.0*dy);

    for (j = 0; j < 3; j++) { p1[j][0] = p2[j][0] = pos[j][0]; } // p1=p2=pos
    p1[2][0] = pos[2][0] - dz;
    p2[2][0] = pos[2][0] + dz;
    plasma_density(nRay, p1, den1, gradd, ds, rflg);
    plasma_density(nRay, p2, den2, gradd, ds, rflg);
    grdn[2] = (den2[0] - den1[0])/(2.0*dz);

    for (j = 0; j < 3; j++) { p[j] = pos[j][0]; } // p=pos
    r = sqrt(dot_product(p, p));

    //fprintf(fh, "%g %g\n", r,  den[0]);
    fprintf(fh, "%g %g\n", pos[1][0],  den[0]);
    
    fprintf(fh1, "%6f (%g, %g, %g) [%g %g %g] [%g %g %g]\n", r,
	    p[0], p[1], p[2],
	    grdn[0], grdn[1], grdn[2],
	    grd[0],  grd[1],  grd[2]);

    //pos[0][0] += 0.00005;
    pos[1][0] += 0.00005;
    //pos[2][0] -= 0.00005;
  }

  fclose(fh);
  fclose(fh1);

  return 0;
}

