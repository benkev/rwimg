/*
 * Return np numbers spaced evenly on a log scale in array sp.
 * The sequence starts at 10^start and ends with 10^stop.
 */
#include <math.h>

void log10space(double start, double stop, long np, double sp[]) {
  double delx = (stop - start)/(np - 1); 
  double x = start;
  long i;

  for (i = 0; i < np; i++) {
    sp[i] = pow(10.0, x);
    x = x + delx;
  }
} 
 
