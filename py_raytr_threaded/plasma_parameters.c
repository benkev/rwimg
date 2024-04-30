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

#include <stdio.h>
#include <math.h>
#include "raytrace.h"
#include "streamer.h"


Streamers_T oStreamers = NULL;

/* Function declaration for the streamer bfield calculator */
inline void calculate_streamer_field(double,double,double,double[]);

/*
 * Adds a stremer to the streamer list which will be applied to all future
 * simulations until remove_streamers is called.
 */
void add_streamer(double dRadius, double dTheta, double dPhi,
		  double dOrientation, double dDensity, double dBaseStrength, 
		  double dStalkStrength,double dScaleX,
		  double dScaleY, double dScaleZ){

  if(oStreamers==NULL)
    oStreamers = Streamers_new();

  Streamers_makeStreamer(oStreamers, dRadius, dTheta, dPhi,dOrientation,
			 dDensity,dBaseStrength,dStalkStrength,dScaleX,
			 dScaleY,dScaleZ);
}

/*
 * Reomves all of the streamers from the streamer list, so that none will
 * be applied to future simulations.
 */
void remove_streamers(){
  Streamers_free(oStreamers);
  oStreamers = NULL;
}

void plasma_parameters(int start,
		       int stop,
		       int nRay, 
		       double Param_I[],
		       double Pos_ID[nRay][3], 
		       double Rho_I[nRay], 
		       double GradRho_ID[nRay][3],
		       double Bfield_ID[nRay][3],
		       double DS_I[nRay], 
		       short Flags_I[nRay],
		       int rtmode) {
  //
  // Coronal density model by Kuniji Saito
  //
  // 
  //
  /* Access parameters by structure field names */
  struct param *prm = (void *) Param_I; 
  
  double static const dx = 1e-6, dy = 1e-6, dz = 1e-6;
  double static const ddx = 2e-6, ddy = 2e-6, ddz = 2e-6;

  double den1, den2;
  double r2;
  int iRay;
  double x, y, z, r;
  double Ne, Nex, Ney, Nez;

  double mr_3r2; /* 3 times dot product of dipole m and radius r over r^2 */ 
  double Msun_r3; /* Msun/r^3 */


  for (iRay = start; iRay < stop; iRay++) {
    
    

    x = Pos_ID[iRay][0];
    y = Pos_ID[iRay][1];
    z = Pos_ID[iRay][2];

    r2 = x*x + y*y + z*z;
    r = sqrt(r2);
    
    Ne = density_only(iRay, Param_I, x, y, z);

    /* printf("plasma_parameters: Ne = %g\n", Ne); */

    /* Nex, Ney, Nez are partial derivatives of Ne, 
     * calculated numerically */
    den1 = density_only(iRay, Param_I, x-dx, y, z);
    den2 = density_only(iRay, Param_I, x+dx, y, z);
    Nex = (den2 - den1)/ddx;

    den1 = density_only(iRay, Param_I, x, y-dy, z);
    den2 = density_only(iRay, Param_I, x, y+dy, z);
    Ney = (den2 - den1)/ddy;

    den1 = density_only(iRay, Param_I, x, y, z-dz);
    den2 = density_only(iRay, Param_I, x, y, z+dz);
    Nez = (den2 - den1)/ddz;

    /* Turn number density into g/cm^3 */ 

    Rho_I[iRay] = Ne*prm->ProtonMass_g;


    GradRho_ID[iRay][0] = Nex*prm->ProtonMass_g;
    GradRho_ID[iRay][1] = Ney*prm->ProtonMass_g;
    GradRho_ID[iRay][2] = Nez*prm->ProtonMass_g;

    
    
    if(rtmode == 3){

      /* double B[3] = {0., 0., 0.}; */
      double B[3];

      /* Calculate the magnetic field vector */

      /* mr_3r2 = 3.0*(prm->dirMx*x + prm->dirMy*y + prm->dirMz*z)/r2; */
      /* Msun_r3 = prm->Msun_G/pow(r,3); */
      
      Streamers_calculateBField(oStreamers, x, y, z, B);

      /* Bfield_ID[iRay][0] = Msun_r3*(mr_3r2*x - prm->dirMx) + B[0]; */
      /* Bfield_ID[iRay][1] = Msun_r3*(mr_3r2*y - prm->dirMy) + B[1]; */
      /* Bfield_ID[iRay][2] = Msun_r3*(mr_3r2*z - prm->dirMz) + B[2]; */
      Bfield_ID[iRay][0] = prm->Msun_G;
      Bfield_ID[iRay][1] = 0.;
      Bfield_ID[iRay][2] = 0.;


    }

    /* Revert the step to original if it is not in
     * the process of finding the critical surface */
    if (!(Flags_I[iRay] & PENETR)) {
      if (r > 2.0) 
        DS_I[iRay] = prm->DeltaS;
      else if (r > (prm->r_corona + prm->DeltaS)) 
        DS_I[iRay] = 0.1*prm->DeltaS;
      else
        DS_I[iRay] = 0.01*prm->DeltaS;
    }
  } /* for (iRay = 0; iRay < nRay; iRay++)  */
}  /* void plasma_density() */


/*
 * This function finds the electron density of the plasma when given
 * the ray, plasma parameters, and cartesian coordinates
 */
inline double density_only(int iRay, double Param_I[], 
			   double x, double y, double z) {


  short static FirstEntry = 1;


  double static radeg, dtor;
  double const static R0 = 6.955e5;   // km, Solar radius
  // The Saito's function coefficients
  double const static  g1 = 3.09e8, g2 = 1.58e8, g3 = 2.51e6;
  // The merging points in upper chromosphere, r_chromo and
  // lower corona, r_corona, between which the el. density
  // Ne is smoothly replaced by a polynomial
  // The chromosphere height is supposed to be ~ 10000 km,
  // h_chrom = 10000, and the merging points are 1000km
  // away from the chromospheric boundary on both sides:
  // r_chromo = 1 + (h_chrom - 1000)/R0,
  // r_corona = 1 + (h_chrom + 1000)/R0,
  /* double const static r_chromo = 1.01294964029;  // *R0 = R0 + 9000, km */
  /* double const static r_corona = 1.01582733813;  // *R0 = R0 + 11000, km */
  /* double r_chromo =     Param_I[10]; /\*  "top" of chromosphere  *\/ */
  /* double r_corona =     Param_I[12]; /\*  "bottom" of corona *\/ */

  /* Access parameters by structure field names */
  struct param *prm = (void *) Param_I; 

  double static th0; // = 30.0*dtor;
  double static ph0; // = 90.0*dtor;
  double static s0x; // = sin(th0)*cos(th0);
  double static s0y; // = sin(th0)*sin(th0);
  double static s0z; // = cos(th0);
  double static th1, th2, th3;
  double static s1x, s1y, s1z;
  double static s2x, s2y, s2z;
  double static s3x, s3y, s3z;
  double static ph1, ph2, ph3;
  double sig, dsig2; // sigma, 2*sigma^2
  double static Cn0, Cn1, Cn2, Cn3;
  double bet;  // Arc distance b/w streamer axis and the point (x,y,z)

  double r, r2;
  double t1, t2, t3;
  double rm2d5, rm6, rm16;
  double cosTh, sqrtCosTh;
  double Ne;
  double static r_corm16, r_corm6, r_corm2d5, r_corm17, r_corm7, r_corm3d5;
  double saito_rcor;
  // Linear systems A*a = rhsa and A*b = rhsb
  double static A[4][4], rhsa[4], a[4];
  double rhsb[4], b[4]; 
  
  //
  // Initialization
  //
  if (FirstEntry) {

    FirstEntry = 0;
    radeg = 180.0/pi;
    dtor = pi/180.0;
    r_corm16 = pow(prm->r_corona,-16);
    r_corm6 = pow(prm->r_corona,-6);
    r_corm2d5 = pow(prm->r_corona,-2.5);
    r_corm17 = pow(prm->r_corona,-16);
    r_corm7 = pow(prm->r_corona,-6);
    r_corm3d5 = pow(prm->r_corona,-2.5);


    A[0][0] = 1.; A[0][1] = prm->r_chromo; 
    A[0][2] = pow(prm->r_chromo,2); A[0][3] = pow(prm->r_chromo,3);
    A[1][0] = 1.; A[1][1] = prm->r_corona; 
    A[1][2] = pow(prm->r_corona,2); A[1][3] = pow(prm->r_corona,3);
    A[2][0] = 0.; A[2][1] = 1.; 
    A[2][2] = 2.*prm->r_chromo; A[2][3] = 3.*pow(prm->r_chromo,2);
    A[3][0] = 0.; A[3][1] = 1.;
    A[3][2] = 2.*prm->r_corona; A[3][3] = 3.*pow(prm->r_corona,2);


    minv(4, A);  // A = inv(A), replace A with its inverse
    rhsa[0] = 5.7e11*exp(-7.7e-4*(R0*(prm->r_chromo - 1) - 500.0));
    rhsa[1] = 0.; rhsa[2] = 0.; rhsa[3] = 0.; 
    mvmul(4, A, rhsa, a);  // a = Ainv*rhsa: solution to system A*a = rhsa


  } /* if (FirstEntry) //Initialization */

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

      int flag = 0;

      if(isnan(Ne))
	flag = 1;
      // "Stitching" function 
      //    p(r,th) = a0 + a1*r + a2*r^2 + (r-r1)/(r2-r1)*saito(r_corona,th)
      // The saito function at the point r_corona
      saito_rcor = g1*r_corm16*t1 + g2*r_corm6*t2 + g3*r_corm2d5*t3;

      rhsb[0] = 0.; rhsb[1] = 1.; 
      rhsb[2] = -438.9e6*R0
	*exp(-7.7e-4*(R0*(prm->r_chromo - 1.0) - 500.0))/saito_rcor; 
      rhsb[3] = (- 16.*g1*r_corm17*t1
		 -  6.*g2*r_corm7*t2
		 - 2.5*g3*r_corm3d5*t3)/saito_rcor; 
      mvmul(4, A, rhsb, b);  // b = Ainv*rhsb: solution to system A*b = rhsb

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


  Streamers_calculateDensity(oStreamers, x, y, z, r, &Ne);


  return Ne;
}
