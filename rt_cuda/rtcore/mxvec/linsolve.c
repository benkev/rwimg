//
// Fast linear system solver
//
void linsolve(const int N, double a[N][N], double b[N], double x[N]) {
  // 
  // Solving a system of linear equations a*x = b,
  // where a is NxN matrix, and b is N-vector.
  //
  // This is a pure Gauss algorithm. The program does not check if the 
  // matrix is well-conditioned. It does not check if a diagonal element 
  // is zero or not before division. The purpose of this program is 
  // "SPEED at the expence of RELIABILITY", whether you like it or not :)
  //
  // Note that the matrix a and RHS vector b are damaged in the course of
  // computations, so one cannot assume a and b remain unchanged. If you
  // need the values 
  // 
  int i, j, k, kp1, Nm1;
  double c, akk;

  //
  // Reduce system to upper-triangle form
  //
  Nm1 = N - 1;
  for(k = 0; k < Nm1; k++) {
    kp1 = k + 1;
    akk = a[k][k]; // Just to save time not accessing the a array
    for(i = kp1; i < N; i++) {
      c = a[i][k]/akk;
      for(j = kp1; j < N; j++) {
	a[i][j] -= c*a[k][j];
      }
      b[i] -= c*b[k];
    }
  }

  //
  // Back substitution run
  //
  x[Nm1] = b[Nm1]/a[Nm1][Nm1]; // Find the last root
  
  for(i = Nm1-1; i >= 0; i--) {
    c = 0.0;
    for(j = i+1; j < N; j++) {
      c = c + a[i][j]*x[j];
    }
    x[i] = (b[i] - c)/a[i][i];
  }
}
