/*
 * Return np numbers spaced evenly on a log scale in array sp.
 * The sequence starts at base^start and ends with base^stop.
 */
#include <math.h>

void logspace(int base, double start, double stop, long np, double sp[]) {
  double delx = (stop - start)/(np - 1); 
  double x = start;
  long i;

  for (i = 0; i < np; i++) {
    sp[i] = pow(base, x);
    x = x + delx;
  }
} 
