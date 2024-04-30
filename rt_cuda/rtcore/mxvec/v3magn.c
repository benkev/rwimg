#include <math.h>
/*
 * Magnitude of vector a[3]
 */
inline double v3magn(double a[]) {
  return sqrt(pow(a[0],2) + pow(a[1],2) +pow(a[2],2));
}
