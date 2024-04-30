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


void intro(int blocks, int threads, Simulation *simulation, plasma_constants
*p_consts_d){

   plasma_parameters<<<blocks,threads>>>(simulation,p_consts_d);

}


__device__ void mvmul(int N, double *a, double x[], double b[]){
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




__global__ void plasma_parameters(Simulation *simulation, struct
plasma_constants *p_consts_d){
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


    if(iRay>=simulation->nRay)
        return; 

    if (BIT_IS_ON(INACTIVE, iRay)) return;     //------------------------>>

    x = simulation->Pos_ID[iRay*3+0];
    y = simulation->Pos_ID[iRay*3+1];
    z = simulation->Pos_ID[iRay*3+2];


    r2 = x*x + y*y + z*z;
    r = sqrt(r2);

    Ne = density_only(simulation->prm, x, y, z, p_consts_d);

    /* Nex, Ney, Nez are partial derivatives of Ne, 
     * calculated numerically */
    den1 = density_only(simulation->prm, x-dx, y, z, p_consts_d);
    den2 = density_only(simulation->prm, x+dx, y, z, p_consts_d);
    Nex = (den2 - den1)/ddx;

    den1 = density_only(simulation->prm, x, y-dy, z, p_consts_d);
    den2 = density_only(simulation->prm, x, y+dy, z, p_consts_d);
    Ney = (den2 - den1)/ddy;

    den1 = density_only(simulation->prm, x, y, z-dz, p_consts_d);
    den2 = density_only(simulation->prm, x, y, z+dz, p_consts_d);
    Nez = (den2 - den1)/ddz;

    /* Turn number density into g/cm^3*/

    simulation->Rho_I[iRay] = Ne*simulation->prm->ProtonMass_g;


    simulation->GradRho_ID[iRay*3+0] = Nex*simulation->prm->ProtonMass_g;
    simulation->GradRho_ID[iRay*3+1] = Ney*simulation->prm->ProtonMass_g;
    simulation->GradRho_ID[iRay*3+2] = Nez*simulation->prm->ProtonMass_g;

    if(simulation->rtmode == 3){

        double B[3];

        /* Calculate the magnetic field vector */
        mr_3r2 = 3.0*(simulation->prm->dirMx*x + simulation->prm->dirMy*y + simulation->prm->dirMz*z)/r2;
        Msun_r3 = simulation->prm->Msun_G/pow(r,3);

        //Streamers_calculateBField(simulation->oStreamers,x,y,z,B);

        B[0]=0;
        B[1]=0;
        B[2]=0;

        simulation->Bfield_ID[iRay*3+0] = Msun_r3*(mr_3r2*x - simulation->prm->dirMx)+B[0];
        simulation->Bfield_ID[iRay*3+1] = Msun_r3*(mr_3r2*y - simulation->prm->dirMy)+B[1];
        simulation->Bfield_ID[iRay*3+2] = Msun_r3*(mr_3r2*z - simulation->prm->dirMz)+B[2];


    }

    /* Revert the step to original if it is not in
     * the process of finding the critical surface */
    if (!(simulation->Flags_I[iRay] & PENETR)) {
        if (r > 2.0) 
            simulation->DS_New_I[iRay] = simulation->prm->DeltaS;
        else if (r > (simulation->prm->r_corona + simulation->prm->DeltaS)) 
            simulation->DS_New_I[iRay] = 0.1*simulation->prm->DeltaS;
        else
            simulation->DS_New_I[iRay] = 0.01*simulation->prm->DeltaS;
    }
    
}  /* void plasma_density() */


/*
 * This function finds the electron density of the plasma when given
 * the ray, plasma parameters, and cartesian coordinates
 */
__device__ double density_only(struct param *prm, 
        double x, double y, double z, struct plasma_constants *p_consts_d) {




    double const R0 = 6.955e5;   // km, Solar radius
    // The Saito's function coefficients
    double const  g1 = 3.09e8, g2 = 1.58e8, g3 = 2.51e6;
    // The merging points in upper chromosphere, r_chromo and
    // lower corona, r_corona, between which the el. density
    // Ne is smoothly replaced by a polynomial
    // The chromosphere height is supposed to be ~ 10000 km,
    // h_chrom = 10000, and the merging points are 1000km
    // away from the chromospheric boundary on both sides:
    // r_chromo = 1 + (h_chrom - 1000)/R0,
    // r_corona = 1 + (h_chrom + 1000)/R0,
    /* double const r_chromo = 1.01294964029;  // *R0 = R0 + 9000, km */
    /* double const r_corona = 1.01582733813;  // *R0 = R0 + 11000, km */
    /* double r_chromo =     simulation->prm[10]; /\*  "top" of chromosphere  *\/ */
    /* double r_corona =     simulation->prm[12]; /\*  "bottom" of corona *\/ */

    /* Access parameters by structure field names */

    double r, r2;
    double t1, t2, t3;
    double rm2d5, rm6, rm16;
    double cosTh, sqrtCosTh;
    double Ne;

    double saito_rcor;
    // Linear systems A*a = rhsa and A*b = rhsb

    double rhsb[4], b[4]; 


    r2 = x*x + y*y + z*z;
    r = sqrt(r2);



    //if ((cnt0--) > 0) printf("x, y, z = %g %g %g, r2 = %10.8f\n", x, y, z,r2);


    if (r > prm->r_corona) { 
        rm2d5 = pow(r,-2.5);         // r^(-2.5)
        rm6 = pow(r,-6);             // r^(-6)
        rm16 = pow(r,-16);           // r^(-16)
        cosTh = fabs(z/r);         // cos(theta) = z/r: only abs value used
        sqrtCosTh = sqrt(cosTh);

        t1 = 1.0 - 0.5*cosTh;          // (1 - 0.5cos(theta))
        t2 = 1.0 - 0.95*cosTh;         // (1 - 0.95cos(theta))
        t3 = 1.0 - sqrtCosTh;          // (1 - sqrt(cos(theta)))

        // The Saito's density distribution
        Ne = g1*rm16*t1 + g2*rm6*t2 + g3*rm2d5*t3;
    }

    else { // r >= prm->r_chromo



        if (r > prm->r_chromo) { // r >= r_chromo, but still r < r_corona:
            rm2d5 = pow(r,-2.5);         // r^(-2.5)
            rm6 = pow(r,-6);             // r^(-6)
            rm16 = pow(r,-16);           // r^(-16)
            cosTh = fabs(z/r);         // cos(theta) = z/r: only abs value used
            sqrtCosTh = sqrt(cosTh);

            t1 = 1.0 - 0.5*cosTh;          // (1 - 0.5cos(theta))
            t2 = 1.0 - 0.95*cosTh;         // (1 - 0.95cos(theta))
            t3 = 1.0 - sqrtCosTh;          // (1 - sqrt(cos(theta)))
            int flag = 0;
            if(isnan(Ne))
                flag = 1;
            // "Stitching" function 
            //    p(r,th) = a0 + a1*r + a2*r^2 + (r-r1)/(r2-r1)*saito(r_corona,th)
            // The saito function at the point r_corona
            saito_rcor = g1*p_consts_d->r_corm16*t1 + g2*p_consts_d->r_corm6*t2 + g3*p_consts_d->r_corm2d5*t3;

            rhsb[0] = 0.; rhsb[1] = 1.; 
            rhsb[2] = -438.9e6*R0
                *exp(-7.7e-4*(R0*(prm->r_chromo - 0.0) - 500.0))/saito_rcor; 
            rhsb[3] = (- 16.*g1*p_consts_d->r_corm17*t1
                    -  6.*g2*p_consts_d->r_corm7*t2
                    - 2.5*g3*p_consts_d->r_corm3d5*t3)/saito_rcor; 
            mvmul(4, p_consts_d->A, rhsb, b);  // b = Ainv*rhsb: solution to system A*b = rhsb

            /* if (iRay == 22281) { */
            /* 	printf("b[0..3] = "); */
            /* 	print1d(4, b); */
            /* 	printf("p_consts_d->a[0..3] = "); */
            /* 	print1d(4, a); */
            /* } */

            Ne = p_consts_d->a[0] + r*(p_consts_d->a[1] + r*(p_consts_d->a[2] + r*p_consts_d->a[3]))
                + (b[0] + r*(b[1] + r*(b[2] + r*b[3])))*saito_rcor;

            /* if (iRay == 22281) {  */
            /* 	printf("Ne = %20.12e\n", Ne); */
            /* 	printf("pos = %20.12e, %20.12e, %20.12e\n", x, y, z); */
            /* } */

        }

        else { // r >= r_corona:

            // Chromospheric exponential density distribution
            Ne = 5.7e11*exp(-7.7e-4*(R0*(r-1) - 500));  // Ne, electron # density
        }
    }


    //Streamers_calculateDensity(oStreamers, x, y, z, r, &Ne);



    return Ne;
}
