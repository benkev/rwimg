#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void cross_product(double a[3], double b[3],double c[]) {
  c[0] = a[0]*b[1] - a[1]*b[0];
  c[1] = -a[0]*b[2] + a[2]*b[0];
  c[2] = a[1]*b[2] - a[2]*b[1];
}

int main(int argc, char* argv[]) {
  double pi = 3.14159265, al = 30.*pi/180.;
  double r[3];
  double p[] = {cos(al), sin(al), 0.0};
  double q[] = {-sin(al), cos(al), 0.0};
  cross_product(p, q, r);
  printf("p     = [%g, %g, %g]\n", p[0], p[1], p[2]);
  printf("q     = [%g, %g, %g]\n", q[0], q[1], q[2]);
  printf("p x q = [%g, %g, %g]\n", r[0], r[1], r[2]);

}
