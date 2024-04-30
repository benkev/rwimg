//
// Minimum value in array a
//
double minval(double a[], int n) {
  double mina = a[0];
  int i;
  if (n <= 1) return mina;
  for (i = 1; i < n; i++) {
    if (a[i] < mina) mina = a[i];
  }
  return mina;
}
