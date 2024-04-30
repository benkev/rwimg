#include <stdio.h>
//#include "raytrace.h"

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

void testsol(int N, double a[N][N], double x[N], double b[N]){
  //
  // Multiply matrix a by vector x, return result in vector b.
  //  b = a*x
  //
  // Used to test that a*x = b, where x is the solution of linear system 
  //  a*x = b
  //
  int i, j;

  for(i = 0; i < N; i++){
    b[i] = 0.0;
    for(j = 0; j < N; j++){
      b[i] += a[i][j]*x[j];
    }
  }
}

void linsolve(const int N, double a[N][N], double b[N], double x[N]) {
  // 
  // Solving a system of linear equations a*x = b for x,
  // where a is NxN matrix, and b is N-vector.
  // Input:
  //   N: system size;
  //   a: system matrix N by N
  //   b: right-hand-side vector of N elements.
  // Output:
  //   x: solution vector of N elements.
  //
  // This is a pure Gauss algorithm. The program does not check if the 
  // matrix a is well-conditioned. It does not check if a diagonal element 
  // is zero or not before division. The purpose of this program is 
  // "SPEED at the expence of REAIABILITY", whether you like it or not :)
  //
  // Note that the matrix a and RHS vector b are damaged in the course of
  // computations, so one cannot assume a and b remain unchanged. If you
  // need the values in a and b, make their copies before the call. 
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


void linsolve2(const int N, const int M, 
	       double a[N][N], double b[N][M], double x[N][M]) {
  // 
  // Solving N systems of linear equations a*x = b for x,
  // where a is NxN matrix, and b is a set of M columns, each of which
  // is N-vector. 
  // This program can be used to find the inverse matrix Ainv for A by 
  // solving the N systems
  //    A*Ainv = I,
  // where I is the unit NxN matrix (all zeroes except the diagonal
  // elements, which are ones)
  // Input:
  //   N: system size;
  //   M: number of systems to solve
  //   a: system matrix N by N
  //   b: right-hand-side vector of N elements.
  // Output:
  //   x: solutions as M-set of column vectors of N elements.
  //
  // This is a pure Gauss algorithm. The program does not check if the 
  // matrix a is well-conditioned. It does not check if a diagonal element 
  // is zero or not before division. The purpose of this program is 
  // "SPEED at the expence of REAIABILITY", whether you like it or not :)
  //
  // Note that the matrix a and RHS vectors b are damaged in the course of
  // computations, so one cannot assume a and b remain unchanged. If you
  // need the values in a and b, make their copies before the call. 
  // 
  int i, j, k, l, kp1, Nm1;
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
      for(l = 0; l < M; l++) b[i][l] -= c*b[k][l];
    }
  }

  //
  // Back substitution run
  //
  for(l = 0; l < M; l++) 
    x[Nm1][l] = b[Nm1][l]/a[Nm1][Nm1]; // Find the last roots of each system
  
  for(i = Nm1-1; i >= 0; i--) {
    for(l = 0; l < M; l++) {
      c = 0.0;
      for(j = i+1; j < N; j++) {
	c = c + a[i][j]*x[j][l];
      }
      x[i][l] = (b[i][l] - c)/a[i][i];
    }
  }
}



int main(){

  int const N = 100; // Specify here the system size
  int i, j;
  //
  // If you do not have the files with matrix a and the RHS vector b
  // in the files a.txt and b.txt, you can uncomment some of the following
  // lines where the a and b are initialized 
  // (and comment out the file operations): 
  //
  //double a[3][3] = {{2., -3., 1.}, {-1., -1., 2.}, {1., 2., -1.}};
  //double a1[3][3];
  //double a[3][3] = {{ 0.78122748,  0.83681232,  0.03065747},
  //                { 0.11431929,  0.87048356,  0.5410083 },
  //		    { 0.28966147,  0.19748816,  0.2879289 }};
  //double b[3] = {5., 0., 3.};
  //double b[3] = {-0.56016294,  0.9195913 ,  0.45032625};
  //double x[3], b1[3];
  double b[100], b1[100], x[100], a[100][100], a1[100][100];
  //double b[50], b1[50], x[50], a[50][50], a1[50][50];
  FILE *fh;

  //
  // The matrix a and the RHS vector b must be in the text files 
  // a.txt and b.txt, respectively.
  //

  printf("privet1\n");

  fh = fopen("a.txt", "r");
  for(i = 0; i < N; i++) {
    for(j = 0; j < N; j++)
      fscanf(fh, "%lf", &a[i][j]);
  }
  fclose(fh);

  fh = fopen("b.txt", "r");
  for(i = 0; i < N; i++) {
    fscanf(fh, "%lf", &b[i]);
  }
  fclose(fh);

  // Save the values in a and b:
  // Copy a1 := a; b1 := b; 
  for(i = 0; i < N; i++) {
    b1[i] = b[i]; 
    for(j = 0; j < N; j++) 
      a1[i][j] = a[i][j];
  }

  printf("Initial Matrix a:\n");
  print2d(N, N, a);
  printf("Right-hand vector b:\n");
  print1d(N, b);
  printf("\n");

  linsolve(N, a, b, x);

  printf("Result x:\n");
  print1d(N, x);


  printf("Right-hand vector b:\n");
  print1d(N, b1);
  printf("\n");
  
  testsol(N, a1, x, b1);

  printf("Test a*x:\n");
  print1d(N, b1);

  return 0;
}
