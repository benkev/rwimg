//
// Matrix by vector multiplication
//
void mvmul(int N, double a[N][N], double x[N], double b[N]){
  //
  // Multiply matrix a by vector x, return result in vector b.
  //  b = a*x
  //
  int i, j;

  for(i = 0; i < N; i++) {
    b[i] = 0.0;
    for(j = 0; j < N; j++) {
      b[i] += a[i][j]*x[j];
    }
  }
}
