#include <stdio.h>
#include <stdlib.h>
#include <math.h>

//
// Scalar Product of 3-Dimensional vectors a and b
//
//double inline dot_product(double a[3], double b[3]) {
//double dot_product(double a[3], double b[3]) {
double dot_product(double a[], double b[]) {
  return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
} 
 

double dsum(double a[], int n) {
  int i;
  double sum = 0.0;
  for (i = 0; i < n; i++) {
    sum += a[i];
  }
  return sum;
}

void fun1(int param) {
  printf("param = %i\n", param);
}

void funar(int m, int n) {
  int a[m][n];
  int i, j;
  printf("m = %i, n = %i\n", m, n);
  for (i = 0; i < m; i++) {
    for (j = 0; j < n; j++) {
      a[i][j] = n*i + j;
    }
  }
  printf("Variable-Length 2D Array a[%i][%i]:\n", m, n);
  for (i = 0; i < m; i++) {
    for (j = 0; j < n; j++) {
      printf("%i ", a[i][j]);
    }
    printf("\n");
  }
  printf("\n");
  //printf("Sum(a) = \n", dsum(a));

}

void func ( void (*f)(int) ) {
  f(10);
  f(-10);
}

void afn(int n, int a[3][n]) {
  int i;
  double static x;
  double* y;
  const static short FALSE = 0, TRUE = 1;
  y = (double *)calloc(3*n, sizeof(double));
  for (i = 0; i < n; i++) {
    a[0][i] = i;
    a[1][i] = i*i;
    a[2][i] = i*i*i;
    printf("a[*][%i] = %i, %i, %i\n", i, a[0][i], a[1][i], a[2][i]);
  }
  printf("static x = %g\n", x); 
  x += 100.0;

  for (i = 0; i < n; i++) {
    y[i] = a[0][i];
    y[n+i] = a[1][i];
    y[2*n+i] = a[2][i];
  }
  for (i = 0; i < n; i++) {
    printf("y[*][%i] = %g, %g, %g\n", i, y[i], y[n+i], y[2*n+i]);
  }
  //printf("dsum(a[0,:] = %g\n", dsum(a, n));
  printf("dsum(y[0,:] = %g\n", dsum(y, n));
  printf("dsum(y[1,:] = %g\n", dsum(y+n, n));
}

int main(int argc, char* argv[]) {
  int m = 10;
  double p, q;
  double a[] = {3, 2, 1};
  double b[] = {3, 2, 1};
  //  int b[3][10];
  //  printf("Hello, World!\n");
  //  func(fun1);
  // afn(m, b);
  //  afn(m, b);
  //  afn(m, b);
  //funar(2,3);
  //funar(5,5);
  //funar(10,8);
  p = dot_product(a, b);
  q = dot_product(a, a);
  printf("p = %g\n", p);
  printf("q = %g\n", q);
  printf("a*a = %g\n", pow(a[0],2) + pow(a[1],2) + pow(a[2],2));
}
