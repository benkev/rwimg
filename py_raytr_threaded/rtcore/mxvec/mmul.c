//
// Matrix by matrix multiplication
//
void mmul(int const N, int const M, int const L, 
	  double a[N][M], double b[M][L], double r[N][L]){
  //
  // Multiply matrix a[NxM] by matrix b[MxL], 
  // return the result in matrix r[NxL].
  //  r = a*x
  //
  // Used to test that a*x = b, where x is the solution of linear system 
  //  a*x = b
  //
  int i, j, k;

  for(k = 0; k < L; k++) {
    for(i = 0; i < N; i++) {
	r[i][k] = 0.0;
	for(j = 0; j < M; j++) {
	  r[i][k] += a[i][j]*b[j][k];
	}
      }
    }
  }
