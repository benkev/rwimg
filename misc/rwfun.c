#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <stdio.h>


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

