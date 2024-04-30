//
//              Kuniji Saito's Plasma Density
//              plus the Newkirk's coronal streamer
//
// The gradient is calculated in catresian coordinates NUMERICALLY!
//
// Plasma density routine for the beam_path.c 
// The plasma density and its gradient are calculated as:
// - exponential for the chromosphere, i.e. below 9000 km over the photosphere;
// - Kuniji Saito's distribution for the corona above 11000 km height;
// - a "stitch function", smoothly merging the chromosphere and corona
//     in the transitional region between 9000 and 11000 km altitude
//

#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <cuda.h>
#include "simulation.h"
#include "plasma_parameters.h"

__device__ __constant__ double A[4*4], rhsa[4], a[4];

void mvmul_host(int N, double *a, double x[], double b[]){
    //
    // Multiply matrix a by vector x, return result in vector b.
    //  b = a*x
    //
    int i, j;

    for(i = 0; i < N; i++) {
        b[i] = 0.0;
        for(j = 0; j < N; j++) {
            b[i] += a[i*N+j]*x[j];
        }
    }
}

void minv(const int N, double *a) {
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
  double *b, *x;

  b = (double *)malloc(sizeof(double)*N*N);
  x = (double *)malloc(sizeof(double)*N*N);

  //
  // Prepare the identity matrix
  //
  for(i = 0; i < N; i++) 
    for(j = 0; j < N; j++) 
      if (i == j) b[i*N+j] = 1.0; else b[i*N+j] = 0.0;

  //
  // Reduce system to upper-triangle form 
  //
  Nm1 = N - 1;
  for(k = 0; k < Nm1; k++) {
    kp1 = k + 1;
    akk = a[k*N+k]; // Just to save time not accessing the a array
    for(i = kp1; i < N; i++) {
      c = a[i*N+k]/akk;
      for(j = kp1; j < N; j++) {
	a[i*N+j] -= c*a[k*N+j];
      }
      for(l = 0; l < N; l++) b[i*N+l] -= c*b[k*N+l];
    }
  }

  //
  // Back substitution run
  //
  for(l = 0; l < N; l++) 
    x[Nm1*N+l] = b[Nm1*N+l]/a[Nm1*N+Nm1]; // Find the last roots of each system

  for(i = Nm1-1; i >= 0; i--) {
    for(l = 0; l < N; l++) {
      c = 0.0;
      for(j = i+1; j < N; j++) {
	c = c + a[i*N+j]*x[j*N+l];
      }
      x[i*N+l] = (b[i*N+l] - c)/a[i*N+i];
    }
  }

  //
  // Override the a with its inverse from x
  //
  for(i = 0; i < N; i++) 
    for(j = 0; j < N; j++) 
      a[i*N+j] = x[i*N+j];

  free(b);
  free(x);
}

extern "C"{
  void intro(int blocks, int threads, Simulation *sim, Simulation
	     *sim_h) {

    static int first = 0;
    /* static struct plasma_constants p_const_d, p_const_h; */


    if (first==0) {
      first = 1;
      /* double r_corm16_h, r_corm6_h, r_corm2d5_h;
       * double r_corm17_h, r_corm7_h, r_corm3d5_h; */
      double A_h[4*4], rhsa_h[4], a_h[4];
      double const R0 = 6.955e5;


      A_h[0*4+0] = 1.; A_h[0*4+1] = sim_h->prm_h->r_chromo; 
      A_h[0*4+2] = pow(sim_h->prm_h->r_chromo,2); 
      A_h[0*4+3] = pow(sim_h->prm_h->r_chromo,3);
      A_h[1*4+0] = 1.; A_h[1*4+1] = sim_h->prm_h->r_corona; 
      A_h[1*4+2] = pow(sim_h->prm_h->r_corona,2); 
      A_h[1*4+3] = pow(sim_h->prm_h->r_corona,3);
      A_h[2*4+0] = 0.; A_h[2*4+1] = 1.; 
      A_h[2*4+2] = 2.*sim_h->prm_h->r_chromo; 
      A_h[2*4+3] = 3.*pow(sim_h->prm_h->r_chromo,2);
      A_h[3*4+0] = 0.; A_h[3*4+1] = 1.;
      A_h[3*4+2] = 2.*sim_h->prm_h->r_corona; 
      A_h[3*4+3] = 3.*pow(sim_h->prm_h->r_corona,2);


      minv(4, A_h);  // A_h = inv(A_h), repla_hce A_h with its inverse
      rhsa_h[0] = 5.7e11*exp(-7.7e-4*(R0*(sim_h->prm_h->r_chromo - 1) 
				      - 500.0));
      rhsa_h[1] = 0.; rhsa_h[2] = 0.; rhsa_h[3] = 0.; 
      // a_h = A_hinv*rhsa_h: solution to system A_h*a_h = rhsa_h
      mvmul_host(4, A_h, rhsa_h, a_h);  


      cudaMemcpyToSymbol(A,A_h,sizeof(double)*16);
      cudaMemcpyToSymbol(a,a_h,sizeof(double)*4);
      cudaMemcpyToSymbol(rhsa,rhsa_h,sizeof(double)*4);

      printf("cudaMemcpyToSymbol(A,A_h,sizeof(double)*16);\n");
      printf("blocks=%d, threads=%d\n", blocks, threads);

    }

    plasma_parameters<<<blocks,threads>>>(sim);

  }
  
  void cleanup(){

  }

} /* extern "C" { */





__device__ void mvmul(int N, double *a, double x[], double b[]) {
  //
  // Multiply matrix a by vector x, return result in vector b.
  //  b = a*x
  //
  int i, j;

  for(i = 0; i < N; i++) {
    b[i] = 0.0;
    for(j = 0; j < N; j++) {
      b[i] += a[i*N+j]*x[j];
    }
  }
}




__global__ void plasma_parameters(Simulation *sim){
    //
    // Coronal density model by Kuniji Saito
    //
    // 
    //


    double const dx = 1e-6, dy = 1e-6, dz = 1e-6;
    double const ddx = 2e-6, ddy = 2e-6, ddz = 2e-6;

    double den1, den2;
    double r2;
    double x, y, z, r;
    double Ne, Nex, Ney, Nez;

    double mr_3r2; /* 3 times dot product of dipole m and radius r over r^2 */ 
    double Msun_r3; /* Msun/r^3 */


    int iRay = blockDim.x*blockIdx.x+threadIdx.x;


    if(iRay>=sim->nRay)         return; 
 

    printf("plasma_parameters: iRay=%d\n", iRay);
    printf("plasma_parameters: sim->nRay=%d\n", sim->nRay);
    printf("plasma_parameters: BIT_IS_ON(INACTIVE, iRay)=%d\n", 
	   BIT_IS_ON(INACTIVE, iRay));

 
    if (BIT_IS_ON(INACTIVE, iRay)) return;     //------------------------>>

    x = sim->Pos_ID[iRay*3+0];
    y = sim->Pos_ID[iRay*3+1];
    z = sim->Pos_ID[iRay*3+2];


    r2 = x*x + y*y + z*z;
    r = sqrt(r2);

    Ne = density_only(sim->prm, sim->oStreamers, x, y, z);


     /* Nex, Ney, Nez are partial derivatives of Ne, 
     * calculated numerically */
    den1 = density_only(sim->prm, sim->oStreamers, x-dx, y, z);
    den2 = density_only(sim->prm, sim->oStreamers, x+dx, y, z);
    Nex = (den2 - den1)/ddx;

    den1 = density_only(sim->prm, sim->oStreamers, x, y-dy, z);
    den2 = density_only(sim->prm, sim->oStreamers, x, y+dy, z);
    Ney = (den2 - den1)/ddy;

    den1 = density_only(sim->prm, sim->oStreamers, x, y, z-dz);
    den2 = density_only(sim->prm, sim->oStreamers, x, y, z+dz);
    Nez = (den2 - den1)/ddz;

    /* Turn number density into g/cm^3*/

    sim->Rho_I[iRay] = Ne*sim->prm->ProtonMass_g;


    sim->GradRho_ID[iRay*3+0] = Nex*sim->prm->ProtonMass_g;
    sim->GradRho_ID[iRay*3+1] = Ney*sim->prm->ProtonMass_g;
    sim->GradRho_ID[iRay*3+2] = Nez*sim->prm->ProtonMass_g;

    if(sim->rtmode == 3){

        double B[3];

        /* Calculate the magnetic field vector */
        mr_3r2 = 3.0*(sim->prm->dirMx*x + sim->prm->dirMy*y + sim->prm->dirMz*z)/r2;
        Msun_r3 = sim->prm->Msun_G/pow(r,3);

        //Streamers_calculateBField(sim->oStreamers,x,y,z,B);

        B[0]=0;
        B[1]=0;
        B[2]=0;

        sim->Bfield_ID[iRay*3+0] = Msun_r3*(mr_3r2*x - sim->prm->dirMx)+B[0];
        sim->Bfield_ID[iRay*3+1] = Msun_r3*(mr_3r2*y - sim->prm->dirMy)+B[1];
        sim->Bfield_ID[iRay*3+2] = Msun_r3*(mr_3r2*z - sim->prm->dirMz)+B[2];

    }

    /* Revert the step to original if it is not in
     * the process of finding the critical surface */
    if (!(sim->Flags_I[iRay] & PENETR)) {
        if (r > 2.0) 
            sim->DS_New_I[iRay] = sim->prm->DeltaS;
        else if (r > (sim->prm->r_corona + sim->prm->DeltaS)) 
            sim->DS_New_I[iRay] = 0.1*sim->prm->DeltaS;
        else
            sim->DS_New_I[iRay] = 0.01*sim->prm->DeltaS;
    }

}  /* void plasma_density() */


/*
 * This function finds the electron density of the plasma when given
 * the ray, plasma parameters, and cartesian coordinates
 */
__device__ double density_only(struct param *prm, Streamers_T oStreamers,
        double x, double y, double z) {




    double const R0 = 6.955e5;   // km, Solar radius
    // The Saito's function coefficients
    double const  g1 = 3.09e8, g2 = 1.58e8, g3 = 2.51e6;
    // The merging points in upper chromosphere, r_chromo and
    // lower corona, r_corona, besim->oStreamers, tween which the el. density
    // Ne is smoothly replaced by a polynomial
    // The chromosphere height is supposed to be ~ 10000 km,
    // h_chrom = 10000, and the merging points are 1000km
    // away from the chromospheric boundary on both sides:
    // r_chromo = 1 + (h_chrom - 1000)/R0,
    // r_corona = 1 + (h_chrom + 1000)/R0,
    /* double const r_chromo = 1.01294964029;  // *R0 = R0 + 9000, km */
    /* double const r_corona = 1.01582733813;  // *R0 = R0 + 11000, km */
    /* double r_chromo =     sim->prm[10]; /\*  "top" of chromosphere  *\/ */
    /* double r_corona =     sim->prm[12]; /\*  "bottom" of corona *\/ */

    /* Access parameters by structure field names */

    double r, r2;
    double t1, t2, t3;
    double rm2d5, rm6, rm16;
    double cosTh, sqrtCosTh;
    double Ne;



    double saito_rcor;
    // Linear systems A*a = rhsa and A*b = rhsb

    double rhsb[4], b[4]; 


    double r_corm16, r_corm6, r_corm2d5, r_corm17, r_corm7, r_corm3d5;


    r2 = x*x + y*y + z*z;
    r = sqrt(r2);



    //if ((cnt0--) > 0) printf("x, y, z = %g %g %g, r2 = %10.8f\n", x, y, z,r2);


    if (r < prm->r_chromo) { 

      // Chromospheric exponential density distribution
        Ne = 5.7e11*exp(-7.7e-4*(R0*(r-1) - 500));  // Ne, electron # density

    }

    else { // r >= prm->r_chromo
        rm2d5 = pow(r,-2.5);         // r^(-2.5)
        rm6 = pow(r,-6);             // r^(-6)
        rm16 = pow(r,-16);           // r^(-16)
        cosTh = fabs(z/r);         // cos(theta) = z/r: only abs value used
        sqrtCosTh = sqrt(cosTh);

        t1 = 1.0 - 0.5*cosTh;          // (1 - 0.5cos(theta))
        t2 = 1.0 - 0.95*cosTh;         // (1 - 0.95cos(theta))
        t3 = 1.0 - sqrtCosTh;          // (1 - sqrt(cos(theta)))

        if (r < prm->r_corona) { // r >= r_chromo, but still r < r_corona:
	  //\\//\\//////////////////////////////////////////////////////
          //  return 0;
	  //\\//\\//////////////////////////////////////////////////////
            r_corm16 = pow(prm->r_corona,-16);
            r_corm6 = pow(prm->r_corona,-6);
            r_corm2d5 = pow(prm->r_corona,-2.5);
            r_corm17 = pow(prm->r_corona,-16);
            r_corm7 = pow(prm->r_corona,-6);
            r_corm3d5 = pow(prm->r_corona,-2.5);

            // int flag = 0;
            // if(isnan(Ne)) flag = 1;
	    //
            // "Stitching" function 
            //    p(r,th) = a0 + a1*r + a2*r^2 + 
	    //              (r-r1)/(r2-r1)*saito(r_corona,th)
            // The saito function at the point r_corona
            saito_rcor = g1*r_corm16*t1 + g2*r_corm6*t2 + g3*r_corm2d5*t3;

            rhsb[0] = 0.; rhsb[1] = 1.; 
            rhsb[2] = -438.9e6*R0
                *exp(-7.7e-4*(R0*(prm->r_chromo - 0.0) - 500.0))/saito_rcor; 
            rhsb[3] = (- 16.*g1*r_corm17*t1
                    -  6.*g2*r_corm7*t2
                    - 2.5*g3*r_corm3d5*t3)/saito_rcor; 
	    // b = Ainv*rhsb: solution to system A*b = rhsb
            mvmul(4, A, rhsb, b);  

            /* if (iRay == 22281) { */
            /* 	printf("b[0..3] = "); */
            /* 	print1d(4, b); */
            /* 	printf("a[0..3] = "); */
            /* 	print1d(4, a); */
            /* } */

            Ne = a[0] + r*(a[1] + r*(a[2] + r*a[3]))
                + (b[0] + r*(b[1] + r*(b[2] + r*b[3])))*saito_rcor;

            /* if (iRay == 22281) {  */
            /* 	printf("Ne = %20.12e\n", Ne); */
            /* 	printf("pos = %20.12e, %20.12e, %20.12e\n", x, y, z); */
            /* } */

        }

        else { // r >= r_corona:
            // The Saito's density distribution
            Ne = g1*rm16*t1 + g2*rm6*t2 + g3*rm2d5*t3;

        }
    }


    //    Streamers_calculateDensity(oStreamers, x, y, z, r, &Ne);



    return Ne;
}
