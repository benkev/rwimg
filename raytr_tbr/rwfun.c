#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <stdio.h>


double const cPi = 3.1415926535897931;
double const ProtonChargeSGSe = 4.8e-10; // StatCoulombs, SGSe
double const cProtonMass = 1.672621636E-24; //g
double const cElectronMass = 9.10938215E-28; //g
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
    printf("%g\t", a[i]);
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
      printf("%g\t", a[i][j]);
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
// Square matrix-by-vector multiplication
//
void mvmul(int N, double a[N][N], double x[N], double b[N]){
  //
  // Multiply square matrix a[NxN] by vector x, return result in vector b.
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

//
// Linear system solving
//
void linsolve(const int N, double a[N][N], double b[N], double x[N]) {
  // 
  // Solving a system of linear equations a*x = b,
  // where a is NxN matrix, and b is N-vector.
  //
  // This is a pure Gauss algorithm. The program does not check if the 
  // matrix is well-conditioned. It does not check if a diagonal element 
  // is zero or not before division. The purpose of this program is 
  // "SPEED at the expence of REAIABILITY", whether you like it or not :)
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
// Matrix inversion "in place"
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
  // "SPEED at the expence of REAIABILITY", whether you like it or not :)
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

