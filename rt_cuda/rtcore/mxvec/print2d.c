#include <stdio.h>
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
