//
// Save 3d image of the trajectories
//


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "raytrace.h"
#include <sys/time.h>
#include <sys/resource.h>
#include <time.h>


double const cPi = 3.1415926535897931;
double const ProtonChargeSGSe = 4.8e-10;      // StatCoulombs, SGSe
double const cProtonMass = 1.672621636E-24;   //g
double const cElectronMass = 9.10938215E-28;  //g
double const cSolarRadius = 6.95E10;         //cm
double const h_chromo = 10000;         //cm


double const DeltaS = 0.1;         // in R_sun units



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

//
// Fast linear system solver
//
void linsolve(const int N, double a[N][N], double b[N], double x[N]) {
  // 
  // Solving a system of linear equations a*x = b,
  // where a is NxN matrix, and b is N-vector.
  //
  // This is a pure Gauss algorithm. The program does not check if the 
  // matrix is well-conditioned. It does not check if a diagonal element 
  // is zero or not before division. The purpose of this program is 
  // "SPEED at the expence of RELIABILITY", whether you like it or not :)
  //
  // Note that the matrix a and RHS vector b are damaged in the course of
  // computations, so one cannot assume a and b remain unchanged. If you
  // need the values 
  // 
  int i, j, k, kp1, Nm1;
  double c, akk;

  //
  // Reduce system to upper-triangle form
  //
  Nm1 = N - 1;
  for(k = 0; k < Nm1; k++) {
    kp1 = k + 1;
    akk = a[k][k]; // Just to save time not accessing the a array
    for(i = kp1; i < N; i++) {
      c = a[i][k]/akk;
      for(j = kp1; j < N; j++) {
	a[i][j] -= c*a[k][j];
      }
      b[i] -= c*b[k];
    }
  }

  //
  // Back substitution run
  //
  x[Nm1] = b[Nm1]/a[Nm1][Nm1]; // Find the last root
  
  for(i = Nm1-1; i >= 0; i--) {
    c = 0.0;
    for(j = i+1; j < N; j++) {
      c = c + a[i][j]*x[j];
    }
    x[i] = (b[i] - c)/a[i][i];
  }
}

//
// Matrix inversion
//
void minv(const int N, double a[N][N]) {
  // 
  // Matrix inversion by solving N systems of linear equations a*ai = I 
  // for ai, where a is NxN matrix, and I is the identity matrix 
  // (all zeroes except the diagonal elements, which are ones)
  // Input:
  //   N: system size;
  //   a: matrix N by N
  // Output:
  //   a: inverse of a.
  //
  // This is a pure Gauss algorithm. The program does not check if the 
  // matrix a is well-conditioned. It does not check if a diagonal element 
  // is zero or not before division. The purpose of this program is 
  // "SPEED at the expence of RELIABILITY", whether you like it or not :)
  //
  int i, j, k, l, kp1, Nm1;
  double c, akk;
  double b[N][N]; // Identity matrix 
  double x[N][N]; // Inverse of a

  //
  // Prepare the identity matrix
  //
  for(i = 0; i < N; i++) 
    for(j = 0; j < N; j++) 
      if (i == j) b[i][j] = 1.0; else b[i][j] = 0.0;

  //
  // Reduce system to upper-triangle form
  //
  Nm1 = N - 1;
  for(k = 0; k < Nm1; k++) {
    kp1 = k + 1;
    akk = a[k][k]; // Just to save time not accessing the a array
    for(i = kp1; i < N; i++) {
      c = a[i][k]/akk;
      for(j = kp1; j < N; j++) {
	a[i][j] -= c*a[k][j];
      }
      for(l = 0; l < N; l++) b[i][l] -= c*b[k][l];
    }
  }

  //
  // Back substitution run
  //
  for(l = 0; l < N; l++) 
    x[Nm1][l] = b[Nm1][l]/a[Nm1][Nm1]; // Find the last roots of each system
  
  for(i = Nm1-1; i >= 0; i--) {
    for(l = 0; l < N; l++) {
      c = 0.0;
      for(j = i+1; j < N; j++) {
	c = c + a[i][j]*x[j][l];
      }
      x[i][l] = (b[i][l] - c)/a[i][i];
    }
  }

  //
  // Override the a with its inverse from x
  //
  for(i = 0; i < N; i++) 
    for(j = 0; j < N; j++) 
      a[i][j] = x[i][j];
}

//
// Matrix by vector multiplication
//
void mvmul(int N, double a[N][N], double x[N], double b[N]){
  //
  // Multiply matrix a by vector x, return result in vector b.
  //  b = a*x
  //
  int i, j;

  for(i = 0; i < N; i++) {
    b[i] = 0.0;
    for(j = 0; j < N; j++) {
      b[i] += a[i][j]*x[j];
    }
  }
}

int main(int argc, char *argv[]) 
{
  //double REarth_D[3] = {150.0, 145.7, 50.0}; // SGI earth position in R_sol
  double REarth_D[3] = {150.0, 0.0, 0.0};   // SGI earth position in R_sol
  //  double Freq =  231.68e6;       // Hz
  //double Freq =  100e6;            // Hz
  //double Freq =  84.5E6;           // Hz
  //double Freq =  10.0E6;           // Hz
  double Freq = 160.0E6;            // Hz
  // ImageRange is (x0, y0, x1, y1), i.e. (XLower, YLower, XUpper, YUpper)
  //double ImageRange_I[4] = {-1.5, -1.5, 1.5, 1.5}; 
  double ImageRange_I[4] = {-1., -1., 1., 1.}; 
  //double ImageRange_I[4] = {-5.0, -5.0, 5.0, 5.0}; 
  double rIntegration = 5.0;  // Radius of "integration sphere" in solar radii
  int const nXPixel = 7, nYPixel = 7; // Image plane dimensions in pixels
  //int const nXPixel = 640, nYPixel = 640; // Image plane dimensions in pixels
  double Intensity_II[nYPixel][nXPixel];  // Sun radioimage at Freq, Hz
  int i, j, k, nTraj;
  //int const nTrajMax = nYPixel/10+21;// Maximum number of selected trajects
  int const nTrajMax = 100;// Maximum number of selected trajectories
  int const lenTrajMax = 3000;      // Maximum lengths of selected trajectories
  int SelTraj_II[2][nTrajMax];      // (y,x) indices of selected trajectories
  int LenTraj_I[nTrajMax];          // Lengths of calculated trajectories 
  double Traj_DII[3][nTrajMax][lenTrajMax];    // The returned trajectories
  double Tb_II[nTrajMax][lenTrajMax]; // Brightness temperatures 
  int iTraj, lenTraj;
  double dy;
  char fname[100], fnpart[100], cmd[100];
  FILE *fh;

  clock_t cstart, cend;
  double elapsed;
  time_t time1, time2;


  cstart = clock();
  time(&time1);

  //SelTraj_II[0][0] = 3; SelTraj_II[1][0] = 4;   // i1, j1
  //SelTraj_II[0][1] = 10; SelTraj_II[1][1] = 10; // i2, j2
  printf("\nSelected pixels for trajectories:\n");
  nTraj = 0;
  //for (i = 0; i < nYPixel; i += 10) {
  //  SelTraj_II[0][nTraj] = i;         // Select i each 10th
  //  SelTraj_II[1][nTraj] = nXPixel/2; // Select j at the middle vertical line
  //  printf("%3i %3i\n", SelTraj_II[0][nTraj], SelTraj_II[1][nTraj]);
  //  if (++nTraj == nTrajMax) break;
  //}
  //for (i = 0; i < 120; i++) { // Denser near the center
  //  SelTraj_II[0][nTraj] = i + 40;  //nYPixel/3;        
  //  SelTraj_II[1][nTraj] = nXPixel/2; // Select j at the middle vertical line
  //  printf("%3i %3i\n", SelTraj_II[0][nTraj], SelTraj_II[1][nTraj]);
  //  if (++nTraj == nTrajMax) break;
  //}
  for (i = 0; i < nYPixel; i++) { 
    for (j = 0; j < nXPixel; j++) { 
      SelTraj_II[0][nTraj] = i;  //         
      SelTraj_II[1][nTraj] = j;  //
      printf("%3i %3i\n", SelTraj_II[0][nTraj], SelTraj_II[1][nTraj]);
      if (++nTraj == nTrajMax) break;
    }
    if (nTraj == nTrajMax) break;
  }

  printf("nTraj = %i\n", nTraj);
  
  if (argc > 1) {
    Freq = 1e6*atof(argv[1]);
  }
  printf("\nFreq = %g\n", Freq);
  printf("nTrajMax, lenTrajMax, nTraj\n = %i, %i, %i\n", 
	 nTrajMax, lenTrajMax, nTraj);

  
  beam_intens_traj(REarth_D, Freq, ImageRange_I, 
		   rIntegration, nXPixel, nYPixel, Intensity_II, 
		   nTrajMax, lenTrajMax, nTraj, SelTraj_II, 
		   LenTraj_I, Traj_DII, Tb_II);

  fh = fopen("image.txt", "w");
  fprintf(fh, "#\n# Image generated via raytracing\n");
  fprintf(fh, "#\n# Radio wave frequency, Hz:\n");
  fprintf(fh, "%12.4e\n", Freq);
  fprintf(fh, "#\n# Image lower-left and upper-right points,\n");
  fprintf(fh, "# in the SGI coordinates, x0, y0, x1, y1,\n");
  fprintf(fh, "# solar radius units:\n");
  for (j = 0; j < 4; j++) fprintf(fh, "%g ", ImageRange_I[j]);
  fprintf(fh, "\n#\n");
  fprintf(fh, "# Image Y and X sizes in pixels:\n");
  fprintf(fh, "%i %i\n", nYPixel, nXPixel);
  fprintf(fh, "#\n# Data:\n");

  
  for (i = 0; i < nYPixel; i++) {
    for (j = 0; j < nXPixel; j++) {
      fprintf(fh, "%g ", Intensity_II[i][j]);
    }
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
  //		      int LenTraj_I[nTraj], 
  //                  double Traj_DII[3][nTraj][lenTrajMax]) 
  printf("#\n# Ray trajectories generated via raytracing\n");
  printf("#\n# Radio wave frequency, Hz:\n");
  printf("%12.4e\n", Freq);
  printf("#\n# Lengths of rays stored:\n");
  for (j = 0; j < nTraj; j++) printf("%i ", LenTraj_I[j]);
  printf("\n#\n");
  printf("# Data:\n");

  fh = fopen("traj.txt", "w");

  fprintf(fh, "#\n# Ray trajectories generated via raytracing\n");
  fprintf(fh, "#\n# Radio wave frequency, Hz:\n");
  fprintf(fh, "%12.4e\n", Freq);
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

  fh = fopen("tbr.txt", "w");

  fprintf(fh, "#\n# Brightness temperatures generated via raytracing\n");
  fprintf(fh, "#\n# Radio wave frequency, Hz:\n");
  fprintf(fh, "%12.4e\n", Freq);
  fprintf(fh, "#\n# Lengths of data stored:\n");
  for (j = 0; j < nTraj; j++) fprintf(fh, "%i ", LenTraj_I[j]);
  fprintf(fh, "\n#\n");
  fprintf(fh, "# Data:\n");

  for (iTraj = 0; iTraj < nTraj; iTraj++) {
    lenTraj = LenTraj_I[iTraj];
    for (i = 0; i < lenTraj; i++) { 
      fprintf(fh, "%12.4e ", Tb_II[iTraj][i]);
      fprintf(fh, "\n");
    }
  }

  fclose(fh);

  //=============================================
  //char fname[100], fnpart[100];
  //double dy;

  dy = (ImageRange_I[3] - ImageRange_I[1])/nYPixel;
  
  
  //strcpy(fname, "tb_");                   // fname = "tb_";
  //sprintf(fnpart, "%i", (int)(Freq/1e6)); // fnpart = (str)Freq/1e6; (in MHz)
  sprintf(fname, "tb_%i.txt", (int)(Freq/1e6)); // fname = Freq/1e6; (in MHz)
  //sprintf(fnpart, "_%i.txt", (int)(Freq/1e6)); // fnpart = Freq/1e6; (in MHz)
  //strcat(fname, fnpart);                  // fname = "tb_" + "200";
  //strcpy(fnpart, ".txt");                 // fnpart = ".txt";
  //strcat(fname, fnpart);                  // fname = "tb_200" + ".txt";

  fh = fopen(fname, "w");
  
  for (i = nYPixel/2; i < nYPixel; i++) {
    fprintf(fh, "%g %g\n", dy/2.0 + (float)(i-nYPixel/2)*dy, 
	    Intensity_II[i][nXPixel/2]);
  }
  //  fprintf(fh, "\n");
  
  fclose(fh);
  
 
  cend = clock();
  elapsed = ((double)(cend - cstart)) / CLOCKS_PER_SEC;
  printf("\n\nCLOCKS_PER_SEC = %i;\n", CLOCKS_PER_SEC);
  printf("Elapsed time: %3g sec = %02g:%02g\n", 
	 elapsed, floor(elapsed/60.), fmod(elapsed,60.));

  //time(&time2);
  //elapsed = difftime (time2, time1);
  //printf("Elapsed time: %3g sec = %02g:%02g\n\n\n\n", 
  //	 elapsed, floor(elapsed/60.), fmod(elapsed,60.));


  sprintf(cmd, "cp image.txt image_%i.txt", (int)(Freq/1e6));
  //sprintf(cmd, "cp image.txt image_%ia.txt", (int)(Freq/1e6));
  system(cmd);
  sprintf(cmd, "cp traj.txt traj_%i.txt", (int)(Freq/1e6));
  system(cmd);


}
