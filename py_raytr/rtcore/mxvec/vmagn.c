#include <math.h>
/*
 * Magnitude of vector a[n]
 */
inline double vmagn(double a[], int n) {
  int i;
  double magn = 0.0;
  if (n <= 1) return a[0];

  for (i = 1; i < n; i++) {
    magn += pow(a[i],2);
  }

  return sqrt(magn);
}
