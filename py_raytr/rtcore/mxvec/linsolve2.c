void linsolve2(const int N, const int M, 
	       double a[N][N], double b[N][M], double x[N][M]) {
  // 
  // Solving N systems of linear equations a*x = b for x,
  // where a is NxN matrix, and b is a set of M columns, each of which
  // is N-vector. 
  // This program can be used to find the inverse matrix Ainv for A by 
  // solving the N systems
  //    A*Ainv = I,
  // where I is the unit NxN matrix (all zeroes except the diagonal
  // elements, which are ones)
  // Input:
  //   N: system size;
  //   M: number of systems to solve 
  //   a: system matrix N by N
  //   b: right-hand-side vector of N elements.
  // Output:
  //   x: solutions as M-set of column vectors of N elements.
  //
  // This is a pure Gauss algorithm. The program does not check if the 
  // matrix a is well-conditioned. It does not check if a diagonal element 
  // is zero or not before division. The purpose of this program is 
  // "SPEED at the expence of REAIABILITY", whether you like it or not :)
  //
  // Note that the matrix a and RHS vectors b are damaged in the course of
  // computations, so one cannot assume a and b remain unchanged. If you
  // need the values in a and b, make their copies before the call. 
  // 
  int i, j, k, l, kp1, Nm1;
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
      for(l = 0; l < M; l++) b[i][l] -= c*b[k][l];
    }
  }

  //
  // Back substitution run
  //
  for(l = 0; l < M; l++) 
    x[Nm1][l] = b[Nm1][l]/a[Nm1][Nm1]; // Find the last roots of each system
  
  for(i = Nm1-1; i >= 0; i--) {
    for(l = 0; l < M; l++) {
      c = 0.0;
      for(j = i+1; j < N; j++) {
	c = c + a[i][j]*x[j][l];
      }
      x[i][l] = (b[i][l] - c)/a[i][i];
    }
  }
}

 
