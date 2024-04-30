//
// Maximum value in array a
//
double maxval(double a[], int n) {
  double maxa = a[0];
  int i;
  if (n <= 1) return maxa;
  for (i = 1; i < n; i++) {
    if (a[i] > maxa) maxa = a[i];
  }
  return maxa;
}
