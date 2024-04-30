//
// Matrix in-place inversion
//
void minv(const int N, double a[N][N]) {
  // 
  // Matrix inversion by solving N systems of linear equations a*ai = I 
  // for ai, where a is NxN matrix, and I is the identity matrix 
  // (all zeroes except the diagonal elements, which are ones)
  // Input:
  //   N: system size;
  //   a: matrix N by N
  // Output:
  //   a: inverse of a.
  //
  // This is a pure Gauss algorithm. The program does not check if the 
  // matrix a is well-conditioned. It does not check if a diagonal element 
  // is zero or not before division. The purpose of this program is 
  // "SPEED at the expence of RELIABILITY", whether you like it or not :)
  //
  int i, j, k, l, kp1, Nm1;
  double c, akk;
  double b[N][N]; // Identity matrix 
  double x[N][N]; // Inverse of a

  //
  // Prepare the identity matrix
  //
  for(i = 0; i < N; i++) 
    for(j = 0; j < N; j++) 
      if (i == j) b[i][j] = 1.0; else b[i][j] = 0.0;

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
      for(l = 0; l < N; l++) b[i][l] -= c*b[k][l];
    }
  }

  //
  // Back substitution run
  //
  for(l = 0; l < N; l++) 
    x[Nm1][l] = b[Nm1][l]/a[Nm1][Nm1]; // Find the last roots of each system
  
  for(i = Nm1-1; i >= 0; i--) {
    for(l = 0; l < N; l++) {
      c = 0.0;
      for(j = i+1; j < N; j++) {
	c = c + a[i][j]*x[j][l];
      }
      x[i][l] = (b[i][l] - c)/a[i][i];
    }
  }

  //
  // Override the a with its inverse from x
  //
  for(i = 0; i < N; i++) 
    for(j = 0; j < N; j++) 
      a[i][j] = x[i][j];
}
