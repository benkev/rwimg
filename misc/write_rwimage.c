#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//#include <string.h>
#include "raytrace.h"
#include <sys/time.h>
#include <sys/resource.h>



double const cPi = 3.1415926535897931;
double const ProtonChargeSGSe = 4.8e-10; // StatCoulombs, SGSe
double const cProtonMass = 1.672621636E-24; //g
double const cElectronMass = 9.10938215E-28; //g


double sum(double a[], int n) {
  double s = 0.0;
  int i;
  for (i = 0; i < n; i++) {
    s += a[i];
  }
  return s;
}

//
// Sum of squared components of vector vec[n]
//
double sum_squares(double vec[], int n) {
  int i;
  double sumsq;
  sumsq = 0.0;
  for (i = 0; i < n; i++) {
    sumsq += vec[i]*vec[i];
  }
  return sumsq;
}

//
// Print 1D array
//
void print1d(int n, double a[n]) {
  int i;
  for (i = 0; i < n; i++) {
    printf("%g ", a[i]);
  }
  printf("\n");
}

//
// Print 2D array
//
void print2d(int m, int n, double a[m][n]) {
  int i, j;
  for (i = 0; i < m; i++) {
    for (j = 0; j < n; j++) {
      printf("%g ", a[i][j]);
    }
    printf("\n");
  }
}

//
// Scalar Product of 3-Dimensional Vectors a and b
//
double inline dot_product(double a[3], double b[3]) {
  //double dot_product(double a[3], double b[3]) {
  return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
} 

//
// Vector cross product
// c = a x b
//
void inline cross_product(double a[3], double b[3],double c[3]) {
  c[0] = a[1]*b[2] - a[2]*b[1];
  c[1] = -a[0]*b[2] + a[2]*b[0];
  c[2] = a[0]*b[1] - a[1]*b[0];
}
//
// Minimum value in array a
//
double minval(double a[], int n) {
  double mina = a[0];
  int i;
  if (n <= 1) return mina;
  for (i = 1; i < n; i++) {
    if (a[i] < mina) mina = a[i];
  }
  return mina;
}
//
// Maximum value in array a
//
double maxval(double a[], int n) {
  double maxa = a[0];
  int i;
  if (n <= 1) return maxa;
  for (i = 1; i < n; i++) {
    if (a[i] > maxa) maxa = a[i];
  }
  return maxa;
}

int main(int argc, char *argv[]) 
{
  //double REarth_D[3] = {150.0, 145.7, 50.0};           // Earth position in the SGI, in solar radii
  double REarth_D[3] = {150.0, 0.0, 0.0};           // Earth position in the SGI, in solar radii
  double RadioFrequency =  42.3e6;                       // MHz
  //double RadioFrequency =  84.5E6;                       // MHz
  double ImageRange_I[4] = {-5.0, -5.0, 5.0, 5.0}; //(x0, y0, x1, y1), i.e. (XLower, YLower, XUpper, YUpper)
  double rIntegration = 25.0;                          // Radius of "integration sphere" in solar radii
  int const nXPixel = 200, nYPixel = 200;                  // Image plane dimensions in pixels
  double Intensity_II[nYPixel][nXPixel];               // Resultant sun image in radiowaves at RadioFrequency
  int i, j, k, nTraj;
  int const nTrajMax = nYPixel/10+1;                         // Maximum number of selected trajectories
  int const lenTrajMax = 1000;                            // Maximum lengths of selected trajectories
  int SelTraj_II[2][nTrajMax]; 
  int LenTraj_I[nTrajMax];                                // Lengths of calculated trajectories 
  double Traj_DII[3][nTrajMax][lenTrajMax];               // The returned trajectories
  int iTraj, lenTraj;
  char fname[100], fnpart[100];
  FILE *fh;

  //SelTraj_II[0][0] = 3; SelTraj_II[1][0] = 4;      // i1, j1
  //SelTraj_II[0][1] = 10; SelTraj_II[1][1] = 10;    // i2, j2
  printf("\nSelected pixels for trajectories:\n");
  nTraj = 0;
  for (i = 0; i < nYPixel; i += 10) {
    SelTraj_II[0][nTraj] = i;             // Select i each 10th
    SelTraj_II[1][nTraj] = nXPixel/2;     // Select j at the middle vertical line
    printf("%3i %3i\n", SelTraj_II[0][nTraj], SelTraj_II[1][nTraj]);
    if (++nTraj == nTrajMax) break;
  }
  printf("nTraj = %i\n", nTraj);
  
  if (argc > 1) {
    RadioFrequency = atof(argv[1]);
    printf("\nFreq = %g\n", RadioFrequency);
  }

  
  beam_intens_traj(REarth_D, RadioFrequency, ImageRange_I, rIntegration, nXPixel, nYPixel, Intensity_II, 
		   nTrajMax, lenTrajMax, nTraj, SelTraj_II, LenTraj_I, Traj_DII);

  fh = fopen("image.txt", "w");
  
  //printf("Image:\n");
  
  for (i = 0; i < nYPixel; i++) {
    for (j = 0; j < nXPixel; j++) {
      //printf("%g ", Intensity_II[i][j]);
      fprintf(fh, "%g ", Intensity_II[i][j]);
    }
    //printf("\n");
    fprintf(fh, "\n");
  }

  fclose(fh);

  //return;
  
  //
  // Save the selected ray trajectories as a table
  // nTraj*SUM(LenTraj_I) lines by 3 columns of X Y Z coordinates
  // 
  //  lenTraj = maxval(LenTraj_I, nTraj);
  //		      int nTraj, int lenTrajMax, int SelTraj_II[2][nTraj], 
  //		      int LenTraj_I[nTraj], double Traj_DII[3][nTraj][lenTrajMax]) 
  printf("#\n# Ray trajectories generated via raytracing\n");
  printf("#\n# Radio wave frequency, Hz:\n");
  printf("%12.4e\n", RadioFrequency);
  printf("#\n# Lengths of rays stored:\n");
  for (j = 0; j < nTraj; j++) printf("%i ", LenTraj_I[j]);
  printf("\n#\n");
  printf("# Data:\n");

  fh = fopen("traj.txt", "w");

  fprintf(fh, "#\n# Ray trajectories generated via raytracing\n");
  fprintf(fh, "#\n# Radio wave frequency, Hz:\n");
  fprintf(fh, "%12.4e\n", RadioFrequency);
  fprintf(fh, "#\n# Lengths of rays stored:\n");
  for (j = 0; j < nTraj; j++) fprintf(fh, "%i ", LenTraj_I[j]);
  fprintf(fh, "\n#\n");
  fprintf(fh, "# Data:\n");

  for (iTraj = 0; iTraj < nTraj; iTraj++) {
    lenTraj = LenTraj_I[iTraj];
    for (i = 0; i < lenTraj; i++) { 
      for (j = 0; j < 3; j++) fprintf(fh, "%12.4e ", Traj_DII[j][iTraj][i]);
      fprintf(fh, "\n");
    }
  }

  fclose(fh);

}
