/*
 * Sum of array a elememts
 */
double sum(double a[], int n) {
  double s = 0.0;
  int i;
  for (i = 0; i < n; i++) {
    s += a[i];
  }
  return s;
}
