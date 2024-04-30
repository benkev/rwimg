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
// The Newkirk's coronal streamer Ne enhancement is just added to the
// calculated Ne

#include <stdio.h>
#include <math.h>
#include "raytrace.h"

void plasma_density_Saito_streamer(int nRay, 
				   double PhysConst_I[],
				   double Pos_DI[3][nRay], 
				   double Rho_I[nRay], 
				   double GradRho_DI[3][nRay],
				   double DS_I[nRay], 
				   short Flags_I[nRay]) {
  //
  // Coronal density model by Kuniji Saito
  //
  // 
  //
  double DeltaS =      PhysConst_I[0];
  /* double ProtonChargeCGSe = PhysConst_I[1]; */
  double ProtonMass_g =     PhysConst_I[2];
  /* double ElectronMass_g =   PhysConst_I[3]; */
  /* double Rsun_cm =          PhysConst_I[4]; */
  /* double Rsun_km =          PhysConst_I[5];  */
  double RhoSurf =          PhysConst_I[6]; /* = 2.5*3.3452E-16 g/cm^3 */
  /* double h_chromo_km =      PhysConst_I[7]; */
  /* double Te_corona_K =      PhysConst_I[8]; */
  /* double Te_chromo_K =      PhysConst_I[9]; */
  /* double AU_m =             PhysConst_I[10]; */
  /* double AU_Rsun =          PhysConst_I[11]; */
            


  double static const dx = 1e-6, dy = 1e-6, dz = 1e-6;
  double static const ddx = 2e-6, ddy = 2e-6, ddz = 2e-6;

  double static den1, den2;
  double static r2;
  int iRay;
  //double const static a1 = 3.09E8, a2 = 1.58E8, a3 = 0.0251E6;
  double static x, y, z, r;
  double static Ne, Nex, Ney, Nez;

  //double const static r_chr_cor = 1.014388489208633;  // *R0 = R0+1e4, km

  double const static r_chromo = 1.01294964029;  // *R0 = R0 + 9000, km
  double const static r_corona = 1.01582733813;  // *R0 = R0 + 11000, km

  // Bits for Flags_I: PENETR: Eps < 0; Critical surface search

  for (iRay = 0; iRay < nRay; iRay++) {
    x = Pos_DI[0][iRay];
    y = Pos_DI[1][iRay];
    z = Pos_DI[2][iRay];

    r2 = x*x + y*y + z*z;
    r = sqrt(r2);
    
    Ne = density_only(x, y, z);
    
    den1 = density_only(x-dx, y, z);
    den2 = density_only(x+dx, y, z);
    Nex = (den2 - den1)/ddx;

    den1 = density_only(x, y-dy, z);
    den2 = density_only(x, y+dy, z);
    Ney = (den2 - den1)/ddy;

    den1 = density_only(x, y, z-dz);
    den2 = density_only(x, y, z+dz);
    Nez = (den2 - den1)/ddz;

    //
    // Play with the Ne: 
    //

    if (0) { // Ten times Ne in gramms cm^-3
      Rho_I[iRay] = 10.0*Ne*ProtonMass_g;
      GradRho_DI[0][iRay] = 10.0*Nex*ProtonMass_g;
      GradRho_DI[1][iRay] = 10.0*Ney*ProtonMass_g;
      GradRho_DI[2][iRay] = 10.0*Nez*ProtonMass_g;
    }
    else if (1) { // Just Ne in gramms cm^-3
      Rho_I[iRay] = Ne*ProtonMass_g;
      GradRho_DI[0][iRay] = Nex*ProtonMass_g;
      GradRho_DI[1][iRay] = Ney*ProtonMass_g;
      GradRho_DI[2][iRay] = Nez*ProtonMass_g;
    }
    else if (0) {  // Ne in cm^-3
      Rho_I[iRay] = Ne;
      GradRho_DI[0][iRay] = Nex;
      GradRho_DI[1][iRay] = Ney;
      GradRho_DI[2][iRay] = Nez;
    }
   

    // Revert the step to original if it is not in
    // the process of finding the critical surface
    if (!(Flags_I[nRay]&PENETR)) {
      if (r > 2.0) 
        DS_I[iRay] = DeltaS;
      else if (r > (r_corona+DeltaS)) 
        DS_I[iRay] = 0.1*DeltaS;
      else
        DS_I[iRay] = 0.01*DeltaS;
    }
  }
}  // void plasma_density()





double inline density_only(double x, double y, double z) {


  short static FirstEntry = 1;

  double const static pi =  3.1415926535897931;
  double static radeg, dtor;
  double const static R0 = 6.950e5;   // km, Solar radius
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
  double const static r_chromo = 1.01294964029;  // *R0 = R0 + 9000, km
  double const static r_corona = 1.01582733813;  // *R0 = R0 + 11000, km

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
  double static sig, dsig2; // sigma, 2*sigma^2
  double static Cn0, Cn1, Cn2, Cn3;
  double static bet;  // Arc distance b/w streamer axis and the point (x,y,z)

  double static r, r2;
  double static t1, t2, t3;
  double static rm2d5, rm6, rm16;
  double static cosTh, sqrtCosTh;
  double static Ne;
  double static r_corm16, r_corm6, r_corm2d5, r_corm17, r_corm7, r_corm3d5;
  double static saito_rcor;
  // Linear systems A*a = rhsa and A*b = rhsb
  double static A[4][4], rhsa[4], rhsb[4], a[4], b[4]; 
  
  //
  // Initialization
  //
  if (FirstEntry) {
    FirstEntry = 0;
    radeg = 180.0/pi;
    dtor = pi/180.0;
    r_corm16 = pow(r_corona,-16);
    r_corm6 = pow(r_corona,-6);
    r_corm2d5 = pow(r_corona,-2.5);
    r_corm17 = pow(r_corona,-16);
    r_corm7 = pow(r_corona,-6);
    r_corm3d5 = pow(r_corona,-2.5);
    A[0][0] = 1.; A[0][1] = r_chromo; 
    A[0][2] = r_chromo*r_chromo; A[0][3] = pow(r_chromo,3);
    A[1][0] = 1.; A[1][1] = r_corona; 
    A[1][2] = r_corona*r_corona; A[1][3] = pow(r_corona,3);
    A[2][0] = 0.; A[2][1] = 1.; 
    A[2][2] = 2.*r_chromo; A[2][3] = 3.*r_chromo*r_chromo;
    A[3][0] = 0.; A[3][1] = 1.;
    A[3][2] = 2.*r_corona; A[3][3] = 3.*r_corona*r_corona;
    minv(4, A);  // A = inv(A), replace A with its inverse
    rhsa[0] = 5.7e11*exp(-7.7e-4*(R0*(r_chromo-1) - 500.0));
    rhsa[1] = 0.; rhsa[2] = 0.; rhsa[3] = 0.; 
    mvmul(4, A, rhsa, a);  // a = Ainv*rhsa: solution to system A*a = rhsa
    rhsb[0] = 0.; rhsb[1] = 1.; 

    //
    // Streamers definitions
    //
    // Streamer spherical coordinates (theta,phi)
    //th0 = 135.0*dtor;
    //ph0 = 0.0*dtor;
    // 0-Streamer Cartesian coordinates on unit sphere
    //s0x = sin(th0)*cos(ph0);
    //s0y = sin(th0)*sin(ph0);
    //s0z = cos(th0);
    //Cn0 = 4.0; // Strength of the streamer
    //Cn0 = -0.9; // Strength of the streamer

    // 1-Streamer spherical coordinates (theta,phi)
    //th1 = 45.0*dtor;
    //ph1 = 0.0*dtor;
    // Streamer Cartesian coordinates on unit sphere
    //s1x = sin(th1)*cos(ph1);
    //s1y = sin(th1)*sin(ph1);
    //s1z = cos(th1);
    //Cn1 = 2.0; // Strength of the streamer
    //Cn1 = -0.75; // Strength of the streamer

    // 2-Streamer spherical coordinates (theta,phi)
    //th2 = 45.0*dtor;
    //ph2 = -90.0*dtor;
    // Streamer Cartesian coordinates on unit sphere
    //s2x = sin(th2)*cos(ph2);
    //s2y = sin(th2)*sin(ph2);
    //s2z = cos(th2);
    //Cn2 = 1.0; // Strength of the streamer
    //Cn2 = -0.5; // Strength of the streamer

    // 3-Streamer spherical coordinates (theta,phi)
    //th3 = 135.0*dtor;
    //ph3 = -90.0*dtor;
    // Streamer Cartesian coordinates on unit sphere
    //s3x = sin(th3)*cos(ph3);
    //s3y = sin(th3)*sin(ph3);
    //s3z = cos(th3);
    //Cn3 = 0.5; // Strength of the streamer
    //Cn3 = -0.25; // Strength of the streamer
    
    //
    // Hole definitions
    //
    // Hole spherical coordinates (theta,phi) Upper
    //th0 = 180.0*dtor;
    //ph0 = 0.0*dtor;
    // 0-Hole Cartesian coordinates on unit sphere
    //s0x = sin(th0)*cos(ph0);
    //s0y = sin(th0)*sin(ph0);
    //s0z = cos(th0);
    //Cn0 = 4.0; // Strength of the streamer
    //Cn0 = -0.75; // Strength of the streamer

    // 1-Hole spherical coordinates (theta,phi) Lower
    //th1 = 0.0*dtor;
    //ph1 = 0.0*dtor;
    // Streamer Cartesian coordinates on unit sphere
    //s1x = sin(th1)*cos(ph1);
    //s1y = sin(th1)*sin(ph1);
    //s1z = cos(th1);
    //Cn1 = 2.0; // Strength of the streamer
    //Cn1 = -0.50; // Strength of the streamer

    // 2-Hole spherical coordinates (theta,phi) Towards Earth
    //th2 = 100.0*dtor;
    //ph2 = 0.0*dtor;
    // Streamer Cartesian coordinates on unit sphere
    //s2x = sin(th2)*cos(ph2);
    //s2y = sin(th2)*sin(ph2);
    //s2z = cos(th2);
    //Cn2 = 1.0; // Strength of the streamer
    //Cn2 = -0.75; // Strength of the streamer



  }

  //x = Pos_D[0];
  //= Pos_D[1];
  //= Pos_D[2];
  //r2 = dot_product(Pos_D, Pos_D);
  r2 = x*x + y*y + z*z;
  r = sqrt(r2);
  //rm1 = 1.0/r;                 // r^(-1)

  //if ((cnt0--) > 0) printf("x, y, z = %g %g %g, r2 = %10.8f\n", x, y, z,r2);
    
    
  if (r < r_chromo) { 
      
    // Chromospheric exponential density distribution
    Ne = 5.7e11*exp(-7.7e-4*(R0*(r-1) - 500));  // Ne, electron # density

  }

  else { // r >= r_chromo

    rm2d5 = pow(r,-2.5);         // r^(-2.5)
    rm6 = pow(r,-6);             // r^(-6)
    rm16 = pow(r,-16);           // r^(-16)
    cosTh = fabs(z/r);         // cos(theta) = z/r: only abs valyue used
    sqrtCosTh = sqrt(cosTh);
    
    t1 = 1.0 - 0.5*cosTh;          // (1 - 0.5cos(theta))
    t2 = 1.0 - 0.95*cosTh;         // (1 - 0.95cos(theta))
    t3 = 1.0 - sqrtCosTh;          // (1 - sqrt(cos(theta)))

    if (r < r_corona) { // r >= r_chromo, but still r < r_corona:
      
      // "Stitching" function 
      //    p(r,th) = a0 + a1*r + a2*r^2 + (r-r1)/(r2-r1)*saito(r_corona,th)
      // The saito function at the point r_corona
      saito_rcor = g1*r_corm16*t1 + g2*r_corm6*t2 + g3*r_corm2d5*t3;

      rhsb[2] = -438.9e6*R0
	*exp(-7.7e-4*(R0*(r_chromo-1.0) - 500.0))/saito_rcor; 
      rhsb[3] = (- 16.*g1*r_corm17*t1
		 -  6.*g2*r_corm7*t2
		 - 2.5*g3*r_corm3d5*t3)/saito_rcor; 
      mvmul(4, A, rhsb, b);  // b = Ainv*rhsb: solution to system A*b = rhsb
      Ne = a[0] + r*(a[1] + r*(a[2] + r*a[3]))
		     + (b[0] + r*(b[1] + r*(b[2] + r*b[3])))*saito_rcor;
    }

    else { // r > r_corona:

      // The Saito's density distribution
      Ne = g1*rm16*t1 + g2*rm6*t2 + g3*rm2d5*t3;
    }
  }
  
  return Ne;

  // Streamer or hole
  //sig = 0.15*(1.0 + exp(-pow(r - 1.5, 2)/0.22));
  //dsig2 = 2.0*pow(sig,2);
  //bet = acos((x*s0x + y*s0y + z*s0z)/r);
  //Ne = Ne + Ne*Cn0*exp(-pow(bet,2)/dsig2);

  //sig = 0.15*(1.0 + exp(-pow(r - 1.5, 2)/0.22));
  //dsig2 = 2.0*pow(sig,2);
  //bet = acos((x*s1x + y*s1y + z*s1z)/r);
  //Ne = Ne + Ne*Cn1*exp(-pow(bet,2)/dsig2);

  //sig = 0.15*(1.0 + exp(-pow(r - 1.5, 2)/0.22));
  //dsig2 = 2.0*pow(sig,2);
  //bet = acos((x*s2x + y*s2y + z*s2z)/r);
  //Ne = Ne + Ne*Cn2*exp(-pow(bet,2)/dsig2);

  //sig = 0.15*(1.0 + exp(-pow(r - 1.5, 2)/0.22));
  //dsig2 = 2.0*pow(sig,2);
  //bet = acos((x*s3x + y*s3y + z*s3z)/r);
  //Ne = Ne + Ne*Cn3*exp(-pow(bet,2)/dsig2);

 return Ne;
}
