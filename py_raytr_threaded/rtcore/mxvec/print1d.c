#include <stdio.h>
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
