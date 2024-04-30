#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define min(a, b) ((a) < (b) ? (a) : (b))
#define max(a, b) ((a) > (b) ? (a) : (b))

double sum(double a[], int n) {
  double s = 0.0;
  int i;
  for (i = 0; i < n; i++) {
    s += a[i];
  }
  return s;
}

//
// Vector cross product
// c = a x b
//
void cross_product(double a[3], double b[3],double c[3]) {
  c[0] = a[1]*b[2] - a[2]*b[1];
  c[1] = -a[0]*b[2] + a[2]*b[0];
  c[2] = a[0]*b[1] - a[1]*b[0];
}

//
// Scalar Product of 3-Dimensional Vectors a and b
//
double inline dot_product(double a[3], double b[3]) {
  return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
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
    printf("%f ", a[i]);
  }
  printf("\n");
}

//
// Print 2D array
//
void print2d(int m, int n, double a[m][n]) {
  int i, j;
  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) {
      printf("%f ", a[i][j]);
    }
    printf("\n");
  }
}


int main(int argc, char *argv[]) {
  double p[] = {10, 1, 2};
  double q[] = {6, 2, 9, -3, 0, 1, 1, -7, 9};
  printf("\np = ");
  print1d(3, p);
  printf("\np^2 = %g\n\n", dot_product(p, p));
  printf("sum(q) = %g\n\n", sum(q, 9));
  // min(a, b) ((a) < (b) ? (a) : (b)) 
  printf("min(0.015,0.01) = %g\n\n", min(0.015, 0.01));
  printf("min(0.01,0.015) = %g\n\n", min(0.01, 0.015));
  printf("min(1.,2.) = %g\n\n", min(1.,2.));
  printf("min(5.,3.) = %g\n\n", min(5.,3.));
  double ten = 10.0;
  printf("10^(-3/2) = %g\n\n", pow(ten,-1.5));

  

}
