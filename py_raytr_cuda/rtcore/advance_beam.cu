#include "simulation.h"
#include "advance_beam.h"
#include "plasma_parameters.h"
#include "calc.h"
#include <cmath>
#include <stdio.h>
#include <cstdio>
#include <cstdlib>
#include <cuda.h>
#include <curand_kernel.h>
#include <dlfcn.h>
#include <err.h>




using namespace std;

/*
* This function initializes the CUDA random number generator.
*/
__global__ void setup_randgen(curandState *state, int nRay){

    int iRay = blockDim.x*blockIdx.x+threadIdx.x;

    if (iRay >= nRay) return;

    curand_init(1234, iRay, 0, &state[iRay]);

}

/*
* This is a GPU function to calculate the dot product.
*/
__device__ double dot_product(double a[3], double b[3]) {
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
} 

/*
* This is a GPU function to calculate the cross product.
*/
__device__ void cross_product(double a[3], double b[3],double c[3]) {
    c[0] =  a[1]*b[2] - a[2]*b[1];
    c[1] = -a[0]*b[2] + a[2]*b[0];
    c[2] =  a[0]*b[1] - a[1]*b[0];
}

/*
* This is a GPU function to calculate the magnituded of a 3D vector
*/
__device__ double v3magn(double a[]) {
    return sqrt(pow(a[0],2) + pow(a[1],2) +pow(a[2],2));
}


/*
 * Given a mean and sigma, returns a random number from a normal distribution.
 */
//__device__ double normal_random(double mean, double sigma){
//  double random1;
//  double random2;  
//  double V1;
//  double V2;
//  double ret = 2;
//
//  while(ret >= 1){
//    random1 = curand_uniform_double();
//    random2 = curand_uniform_double();  
//    V1 = random1*2.0-1.0;
//    V2 = random2*2.0-1.0;
//    ret = (pow(V1,2)+pow(V2,2));
//  }
//
//  return V1*sqrt((-2.0*pow(sigma,2)*log(ret))/ret);
//}


/*
 * Given the density, gradient of density, step size, plasma parameters,
 * direction, and the ray number, this function alters the direction based on
 * a normal distribution 
 */

__global__ void scatter(Simulation *sim,curandState *state) {
  double b;
  double mu;
  double plasma_freq;
  double sigma;
  double Dir_length;
  double Rho;

  int iRay = blockDim.x*blockIdx.x+threadIdx.x;
 
  if(iRay >= sim->nRay)      return; 

  if (BIT_IS_ON(INACTIVE, iRay)) return;     //------------------------>>

  Rho = sim->Rho_I[iRay];
  Rho = Rho/sim->prm->ProtonMass_g; 
  plasma_freq = 8980*sqrt(Rho);

  mu = sqrt(1.0-pow(plasma_freq,2)/pow(sim->prm->Freq,2));//

  if(isnan(mu))
      return;

  b = .00003927*pow(plasma_freq,4)/(pow(sim->prm->Freq,4)*pow(mu,4));

  sigma = mu*sqrt(b*sim->DS_I[iRay]*sim->prm->Rsun_km);

  /*
   * normal_random(0,sigma);
   */
  sim->Dir_ID[iRay*3+0] += curand_normal_double(&state[iRay])*sigma;
  sim->Dir_ID[iRay*3+1] += curand_normal_double(&state[iRay])*sigma;
  sim->Dir_ID[iRay*3+2] += curand_normal_double(&state[iRay])*sigma;

  Dir_length = v3magn(&sim->Dir_ID[iRay*3+0]);

  sim->Dir_ID[iRay*3+0] = sim->Dir_ID[iRay*3+0]/Dir_length;
  sim->Dir_ID[iRay*3+1] = sim->Dir_ID[iRay*3+1]/Dir_length;
  sim->Dir_ID[iRay*3+2] = sim->Dir_ID[iRay*3+2]/Dir_length;
}


/*
* Calculates the brightness temperature.
*/
__device__ double calc_tbr(
        struct param *prm, 
		double Rho,
		double Te,
		double ds,
		double Tbr,
		double *dtau) {
  /* 
   * Calculate brightness temperature
   */

  double nu_eff, Ne, Te32_inv, omega_p2, znsq;
  double alpha; 

  /*
   * Dimensionless plasma parameters
   */
  Ne =  Rho*(prm->ProtonMassInv);  /* Electron number density */
  Te32_inv = pow(Te,-1.5);            /* Te^(-3/2) to calc nu_eff */
  nu_eff = prm->Cnu*Ne*Te32_inv*log(220.0*Te*pow(Ne,-1./3.)); /*Ginzburg*/
  omega_p2 = (prm->e2_4pi_ovr_m)*Ne;       /* Plasma frequency squared */
  znsq = sqrt(prm->Omega2 - omega_p2);     /* sqrt(w^2 - wp^2) */

  /* 
   * Absorption coefficient 
   */
  alpha = omega_p2*nu_eff/(znsq*(prm->Omega)*(prm->c_light_cms)); 
 
  /*
   * Optical depth inctement
   */
  *dtau = alpha*ds*(prm->Rsun_cm);

  Tbr = Te + exp(-(*dtau))*(Tbr - Te);
  /* printf("alpha = %e, dtau = %e \n", alpha, *dtau); */
 
  return Tbr;
}



/*
* Calculates the Stokes parameters.
*/ 
__device__ void calc_tbriquv(
        struct param *prm,
		double Dir_D[3], 
		double Rho,
		double Bfield_D[3],
		double Te,
		double ds,
		double TbrIQUV_P[4],
		double *dtau,
		int thread_num) {
/* Calculate Stokes I, Q, U, and V components of 
 * brightness temperature.
 *  
 * Compute Mueller matrix, Mue, for ds ray arclength, assuming the
 * plasma parameters are constant over ds.
 * The Mueller matrix is the matrix exponential of the differential
 * Mueller matrix:
 *     M = e^(Mdif*ds).
 * The differential Mueller matrix depends on 7 parameters:
 *            |al  be  ga  de|
 *            |be  al  mu  nu|
 *     Mdif = |ga -mu  al  et|
 *            |de -nu -et  al|
 *The simplified matrix, devoid the ray torque effects, is
 *            |al  be   0  de|,
 *            |be  al  mu   0|
 *     Mdif = |0  -mu  al  et|
 *            |de   0 -et  al|
 */

  double Bfield_abs2, Bfield_abs, Bdir_D[3];
  double costh, sinth2;
  double nu_eff, Ne, Te32_inv;
  double u_plasma, v_plasma, w_plasma;
  double P_plasma, Q_plasma, R_plasma, V_plasma;
  double Mcoef, Mcoefw, McoefV;
  double al, be, de, mu, et, be2, de2, mu2, et2;
  double x2, y2, xy, xy2, xmy, xmy2, xymd, xi, om, xi2, om2, h;
  double bede, bemu, deet, muxi, muom, etmu, etxi, etom;
  double ombe, omde, xibe, xide, xids, omds; 
  double ch, sh, cs, sn, chcs, atn; //, ds_cm;
  double M23_2, ImTe, one_minus_u, one_minus_u2_inv, sqrt_u_plasma;
  double M11, M21, M31, M41;
  double M12, M22, M32, M42;
  double M13, M23, M33, M43;
  double M14, M24, M34, M44;
  double I, Q, U, V;



  //ds_cm = ds*(prm->Rsun_cm); /* Ray arclength increment in centimeters */

  /*
   * Dimensionless plasma parameters
   */
  Bfield_abs2 = dot_product(Bfield_D, Bfield_D);   /* = |Bfield|^2 */
  Bfield_abs = sqrt(Bfield_abs2);              /* Bfield_abs = |Bfield| */
  Bdir_D[0] = Bfield_D[0]/Bfield_abs;
  Bdir_D[1] = Bfield_D[1]/Bfield_abs;
  Bdir_D[2] = Bfield_D[2]/Bfield_abs;
  costh = dot_product(Dir_D, Bdir_D); /* Angle b/w ray and B */
  Ne =  Rho*(prm->ProtonMassInv);     /* Electron number density */
  Te32_inv = pow(Te,-1.5);            /* Te^(-3/2) to calc nu_eff */
  nu_eff = prm->Cnu*Ne*Te32_inv*log(220.0*Te*pow(Ne,-1./3.)); /*Ginzburg*/
  u_plasma = prm->e_ovr_mcw2*Bfield_abs2;        /* u = (e*B/mcw)^2 */  
  v_plasma = prm->e2_4pi_ovr_mw2*Ne;             /* v = 4*pi*e^2*Ne/mw^2 */
  w_plasma = nu_eff/prm->Omega;                  /* w = nu_eff/w */
  one_minus_u = 1.0 - u_plasma;
  one_minus_u2_inv = pow(one_minus_u,-2);  /* 1/(1-u)^2 */
  sqrt_u_plasma = sqrt(u_plasma);
  one_minus_u2_inv = pow(one_minus_u,-2);  /* 1/(1-u)^2 */
  P_plasma = v_plasma*(1.0 + u_plasma)*one_minus_u2_inv;
  Q_plasma = u_plasma*v_plasma*(3.0 - u_plasma)*one_minus_u2_inv;
  R_plasma = 2.0*v_plasma*sqrt_u_plasma*one_minus_u2_inv;
  V_plasma = v_plasma/one_minus_u;
  /* Mcoef = 0.5*prm->k0/sqrt(1.0 - V_plasma); */
  Mcoef = 0.5*prm->k0_Rsun/sqrt(1.0 - V_plasma);
  if(isnan(Mcoef)){
    printf("nan\n");
  }
  Mcoefw = Mcoef*w_plasma;
  McoefV = Mcoef*V_plasma;
    

  /*
   * Differential Mueller matrix M
   */
  sinth2 = 1. - pow(costh,2);

  /* Absorption terms */
  al = -Mcoefw*(2.*P_plasma - Q_plasma*sinth2); /* Diagonal elements */
    
  /* Dichroic terms */
  be = Mcoefw*Q_plasma*sinth2;
  /* ga = Mcoefw*Q*sinth2*sin2ph; ph = 0 */
  de = -2.*Mcoefw*R_plasma*costh; 

  /* Birefringent terms */
  mu = -2.*McoefV*sqrt_u_plasma*costh;
  /* nu =  -McoefV*u_plasma*sinth2*sin2ph; ph = 0 */
  et =  McoefV*u_plasma*sinth2;
    
  /*
   * Mueller matrix - analytical solution by Mosino, Barbosa-Garcia,
   * Starodumov et al., Optics Communications, 2000.
   */
  //    Mue = matrix(empty((4,4),dtype=double))

  be2 = pow(be,2);
  /* ga2 = pow(ga,2); */
  de2 = pow(de,2);
  mu2 = pow(mu,2);
  /* nu2 = pow(nu,2); */
  et2 = pow(et,2);
    /* x = array((be, ga, de)); */
    /* y = array((et, -nu, mu)); */
  xy = be*et + de*mu;  /* -ga*nu */
  xy2 = pow(xy,2);
  /* sign of xy */
  if (xy == 0.) h = 0.;
  else if (xy < 0.) h = -1.;
  else h = 1.;
  x2 = be*be + de*de; /* + ga*ga */
  y2 = et*et + mu*mu; /* + nu*nu */
  xmy = 0.5*(x2 - y2);
  xmy2 = pow(xmy,2);
  xymd = sqrt(xmy2 + xy2); 
  xi2 = xymd + xmy;
  om2 = xymd - xmy;
  xi = sqrt(xi2);
  om = sqrt(om2);
  bede = be*de;
  /* bega = be*ga */
  bemu = be*mu;
  /* benu = be*nu; */
  deet = de*et;
  /* denu = de*nu; */
  /* gade = ga*de
   * gaet = ga*et
   * gamu = ga*mu */
  muxi = mu*xi;
  muom = mu*om;
  /* nuet = nu*et
   * numu = nu*mu
   * nuxi = nu*xi
   * nuom = nu*om */
  etmu = et*mu;
  etxi = et*xi;
  etom = et*om;
  ombe = om*be;
  omde = om*de;
  /* omga = om*ga; */
  xibe = xi*be;
  xide = xi*de;
  /* xiga = xi*ga; */
  /* xids = xi*ds_cm; */
  /* omds = om*ds_cm; */
  xids = xi*ds;
  omds = om*ds;
  ch = cosh(xids);
  sh = sinh(xids);
  sn = sin(omds); 
  cs = cos(omds);
  chcs = ch - cs;

  /* atn = exp(al*ds_cm)/(xi2 + om2); /\* Attenuation factor *\/ */
  atn = exp(al*ds)/(xi2 + om2); /* Attenuation factor */
  
  /*
   * Optical depth inctement
   */
  /* *dtau = -al*ds_cm; */
  *dtau = -al*ds;
  
    
  /* Simplified Mueller matrix without the ray torsion effects */
  M23_2 = (h*omde - muxi)*sh;
  M11 = atn*((om2 + x2)*ch + (xi2 - x2)*cs);                /* M11 */
  M12 = atn*((ombe - h*etxi)*sn + (xibe + h*etom)*sh);      /* M12 = M21 */
  M13 = atn*(deet - bemu)*chcs;                             /* M13 = -M31*/
  M14 = atn*((omde - h*muxi)*sn + (xide + h*muom)*sh);      /* M14 = M41 */
  M21 = M12;                                                /* M21 = M12 */
  M22 = atn*((om2 + be2 - mu2)*ch + (xi2 - be2 + mu2)*cs);  /* M22 */
  M23 = atn*(-h*(xide + muom)*sn + M23_2);                  /* M23 */
  M24 = atn*(bede + etmu)*chcs;                             /* M24 = M42 */
  M31 = -M13;                                               /* M31 = -M13 */
  M32 =  atn*((h*xide + muom)*sn - M23_2);                  /* M32 */
  M33 = atn*((om2 - mu2 - et2)*ch + (xi2 + mu2 + et2)*cs);  /* M33 */
  M34 = atn*(-(h*xibe + etom)*sn + (h*ombe - etxi)*sh);     /* M34 = -M43 */
  M41 = M14;                                                /* M41 = M14*/
  M42 = M24;                                                /* M42 = M24 */
  M43 = -M34;                                               /* M43 = -M34 */
  M44 = atn*((om2 + de2 - et2)*ch + (xi2 - de2 + et2)*cs);  /* M44 */

  /* Full Mueller matrix with the ray torsion effects */
  /* M14_1 = (omde - h*muxi)*sn + (xide + h*muom)*sh; */
  /* M24_2 = (bede + etmu)*chcs; */
  /* M34_1 = (h*xibe + etom)*sn; */
  /* M34_2 = (h*ombe - etxi)*sh; */
  /* M23_2 = (h*omde - muxi)*sh; */
  /* M12_1 = (ombe - h*etxi)*sn + (xibe + h*etom)*sh; */
  /* M13_2 = (deet - bemu)*chcs; */
  /* M14_2 = (benu + gaet)*chcs; */
  /* M24_1 = (h*xiga - nuom)*sn - (h*omga + nuxi)*sh */
  /* M24_2 = (bede + etmu)*chcs; */
  /* M34_3 = (gade - numu)*chcs; */
  /* M23_3 = (bega - nuet)*chcs; */
  /* M12_2 = (gamu + denu)*chcs; */
  /* M13_1 = (omga + h*nuxi)*sn + (xiga - h*nuom)*sh; */
  /* Mue[0,0] = (om2 + x2)*ch + (xi2 - x2)*cs; M11 */
  /* Mue[1,0] = M12_1 + M12_2;  M12 */
  /* Mue[2,0] = M13_1 + M13_2;  M13 */
  /* Mue[3,0] = M14_1 - M14_2   M14 */
  /* Mue[0,1] = M12_1 - M12_2;  M21 */
  /* Mue[1,1] = (om2 + be2 - mu2 - nu2)*ch + (xi2 - be2 + mu2 + nu2)*cs; */
  /* Mue[2,1] = -h*(xide + muom)*sn + M23_2 + M23_3; M23 */
  /* Mue[3,1] = M24_1 + M24_2;  M24 */
  /* Mue[0,2] = M13_1 - M13_2;  M31 */
  /* Mue[1,2] =  (h*xide + muom)*sn - M23_2 + M23_3;  M32 */
  /* Mue[2,2] = (om2 + ga2 - mu2 - et2)*ch + (xi2 - ga2 + mu2 + et2)*cs M33 */
  /* Mue[3,2] = -M34_1 + M34_2 + M34_3;  M34 */
  /* Mue[0,3] =  M14_1 + M14_2;  M41 */
  /* Mue[1,3] = -M24_1 + M24_2;  M42 */
  /* Mue[2,3] =  M34_1 - M34_2 + M34_3;  M43 */
  /* Mue[3,3] = (om2 + de2 - nu2 - et2)*ch + (xi2 - de2 + nu2 + et2)*cs; M44 */

  /* 
   * Calculate new Stokes vector as S(i+1) = Te(i) + M*(S(i) - Te(i)),
   * where S = |i, Q, U, V|, and Te = |Te, 0, 0, 0|.
   */
  I = TbrIQUV_P[0];
  Q = TbrIQUV_P[1];
  U = TbrIQUV_P[2];
  V = TbrIQUV_P[3];

  ImTe = I - Te;
  TbrIQUV_P[0] = M11*ImTe + M21*Q + M31*U + M41*V + Te;
  TbrIQUV_P[1] = M12*ImTe + M22*Q + M32*U + M42*V;
  TbrIQUV_P[2] = M13*ImTe + M23*Q + M33*U + M43*V;
  TbrIQUV_P[3] = M14*ImTe + M24*Q + M34*U + M44*V;


  /* if (thread_num == 0) { */
  /*   printf("Rho = %e, prm->ProtonMassInv = %e, Ne = %g\n", Rho,  */
  /* 	   prm->ProtonMassInv, Ne); */
  /* printf("Bx = %g, By = %g, Bz = %g, B = %g\n", */
  /* 	 Bfield_D[0], Bfield_D[1], Bfield_D[2], */
  /* 	 sqrt(pow(Bfield_D[0],2) + pow(Bfield_D[1],2) + pow(Bfield_D[2],2))); */
  /* printf("Bdirx= %g, Bdiry = %g, Bdirz = %g, |Bdir|^2 = %g\n", */
  /* 	 Bdir_D[0], Bdir_D[1], Bdir_D[2], dot_product(Bdir_D,Bdir_D)); */
  /* printf("Bfield_abs = %g, costh = %g, th = %d, sinth2 = %e, Ne = %g\n", */
  /* 	 Bfield_abs, costh, (int)(acos(costh)*180./pi), sinth2, Ne); */
  /* printf("Te = %e, nu_eff = %g, Ne = %g\n", */
  /* 	 Te, nu_eff, Ne); */
  /* printf("u_plasma = %e, v_plasma = %e, w_plasma = %e\n", */
  /* 	 u_plasma, v_plasma, w_plasma); */
  /* printf("P_plasma = %e, Q_plasma = %e, R_plasma = %e, V_plasma = %e\n", */
  /* 	 P_plasma, Q_plasma, R_plasma, V_plasma); */

  /* printf("dtau = %e, ds = %g, Rsun_cm=%g\n", *dtau, ds, prm->Rsun_cm); */
  /* printf("al=%e, be=%e, de=%e, mu=%e, et=%e \n", al, be, de, mu, et); */
  /* printf("xi=%e, om=%e, xids=%e, omds=%e, chcs=%e\n", xi, om, xids, omds); */
  /* printf("=%e, sh=%e, cs=%e, sn=%e, chcs=%e\n", ch, sh, cs, sn, chcs); */
  /* printf("M11, M21, M31, M41 = %e %e %e %e\n", M11, M21, M31, M41); */
  /* printf("M12, M22, M32, M42 = %e %e %e %e\n", M12, M22, M32, M42); */
  /* printf("M13, M23, M33, M43 = %e %e %e %e\n", M13, M23, M33, M43); */
  /* printf("M14, M24, M34, M44 = %e %e %e %e\n", M14, M24, M34, M44); */
  /* printf("I,Q,U,V   = %e %e %e %e\n", I, Q, U, V); */
  /* printf("ImTe = %e, atn = %e\n", ImTe, atn); */
  /* printf("TbrIQUV_P = %e %e %e %e\n", TbrIQUV_P[0], TbrIQUV_P[1], */
  /*  	 TbrIQUV_P[2], TbrIQUV_P[3]); */
/* } //if (thread_num == 0) */
}


/*
* This is the first part of the simulation, makes the first half step, as well
* as calling scattering.
*/
__global__ void first(Simulation *sim){

    double HalfDS;
    int iRay = blockDim.x*blockIdx.x+threadIdx.x;

    if(iRay >= sim->nRay)     return; 

    if(v3magn(&sim->Pos_ID[iRay*3+0])>sim->rsph+1)
        SET_BIT(INACTIVE, iRay);

    // Do not process the rays that are done or bad
    if (BIT_IS_ON(INACTIVE, iRay)) return;     //------------------------>>

    // Advance r by 1/2 DS 
    HalfDS = 0.5*sim->DS_I[iRay];

    printf("iRay=%d, sim->nRay=%d, HalfDS=%g\n", iRay, sim->nRay, HalfDS);

    // sim->Pos_ID is moved by 1/2 DS in the Dir_DI direction 
    sim->Pos_ID[iRay*3+0] += sim->Dir_ID[iRay*3+0]*HalfDS; 
    sim->Pos_ID[iRay*3+1] += sim->Dir_ID[iRay*3+1]*HalfDS;
    sim->Pos_ID[iRay*3+2] += sim->Dir_ID[iRay*3+2]*HalfDS;

    return;
}


/*
* The second part of the simulation, the main portion with most of the
* simulation. This is everything after the call to plasma_parameters
*/
__global__ void second(Simulation *sim, int iIter) {

  int iTracedPts, itrj, i0, i1, im;
    double const cThird = 1.0/3.0;

    double  Dir1_D[3], Omega_D[3];	
    double  DirVert_D[3];	
    double  StepY_D[3], RelGradRefrInx_D[3];
    double  GradEps_D[3], PosHalfBk_D[3];
    double  Dir_D[3], Xprod_D[3];
    double Te;
    double dtau;

    double HalfDS;            // DS halved
    double Eps, EpsHalfBk, Rho2RhoCr; //, Rho2RhoCr1, Rho2RhoCr2;
    double Coef, Curv, Curv1;
    double LCosAl; // L is inverse grad of \epsilon, Alpha is incidence angle
    double GradEps2, GradEpsDotDir;
    double SolDist;
    //---------------------------
    //double dzeds;            /* dze*DS in cm */
    //double cosh_dzeds, sinh_dzeds;   /* sh(dze*DS), ch(dze*DS) */
    //double exp_alds;                 /* exp(alpha*DS) */


    int iRay = blockDim.x*blockIdx.x+threadIdx.x;


    if(iRay>=sim->nRay)        return; 

//    printf("%x\n", &sim);

//    printf("%f\n", sim.Pos_ID[0]);

//    printf("Here1\n");
    //fprintf(stderr,"scat: %d\n",sim.scattering);
    /*
     * Advance all the rays by 1/2 of the DS step
     */
    /*
       if(sim.scattering)
       scatter(sim.Rho_I[iRay],sim.GradRho_ID[iRay],sim.DS_I[iRay],
               sim.prm,sim.Dir_ID[iRay],iRay);
       */


   //printf("Here ray=%d\n",iRay); 


    

    /*if(sim->rtmode == 3){*/
    /* if (iRay == 1225) { */
    /* 	printf("Call # %g ============================================\n", */
    /* 	       sim->prm->callcount); */
    /* 	printf("          Flag = %d\n", sim->Flags_I[iRay]); */
    /* } */
    /*   /\* if (iRay == 1225) *\/ */
    /*   /\* 	printf("sim->TbrIQUV_IP[0] = %g, TbrIV[1] = %g Entry\n",  *\/ */
    /*   /\* 	       sim->TbrIQUV_IP[iRay*3+0], sim->TbrIQUV_IP[iRay*3+1]);	 *\/ */
    /*}*/

    // Do not process the rays that are done or bad
    if (BIT_IS_ON(INACTIVE, iRay)) return;     //------------------------>>


    HalfDS = 0.5*sim->DS_I[iRay];


    Rho2RhoCr = sim->Rho_I[iRay]*sim->prm->RhoCrInv;


    Eps = 1.0 - Rho2RhoCr;  // Dielectric permittivity

    //GradEps_D = -sim->GradRho_ID(:,iRay)*RhoCrInv


    GradEps_D[0] = -sim->GradRho_ID[iRay*3+0]*sim->prm->RhoCrInv;
    GradEps_D[1] = -sim->GradRho_ID[iRay*3+1]*sim->prm->RhoCrInv;
    GradEps_D[2] = -sim->GradRho_ID[iRay*3+2]*sim->prm->RhoCrInv;




    GradEps2 = dot_product(GradEps_D, GradEps_D);


    //GradEpsDotDir = sum(GradEps_D*Dir_DI(:,iRay))

    GradEpsDotDir = dot_product(GradEps_D, &sim->Dir_ID[iRay*3+0]);


    if (BIT_IS_ON(PENETR, iRay)) { // Condition "Critical surface penetration"

        if (fabs(Eps) > sim->prm->TolEps) { // Another convergence step needed

            //
            // Restore the r shifted along v by DS/2 before the call to 
            //     plasma_params():
            // r = r + v*DS/2
            //
            sim->Pos_ID[iRay*3+0] -= sim->Dir_ID[iRay*3+0]*HalfDS; 
            sim->Pos_ID[iRay*3+1] -= sim->Dir_ID[iRay*3+1]*HalfDS;
            sim->Pos_ID[iRay*3+2] -= sim->Dir_ID[iRay*3+2]*HalfDS;

            //
            // Use the bisection method to reach the critical surface
            //

            // Stop working on the ray if count exceeded maximum
            sim->DirPr_ID[iRay*3+2] -= 1.0;      // Cnt--;
            if (sim->DirPr_ID[iRay*3+2] < 0.0) { // if max count exceeded,
                SET_BIT(BADRAY|INACTIVE, iRay);  // Do not process the ray anymore

                return; //----------------------------------------------------->>
            }

            if (Eps > 0.0) {
                // DS_1 = DS; Eps_1 = Eps
                sim->PosPr_ID[iRay*3+0] = sim->DS_I[iRay];      // DS_1 = DS
                sim->DirPr_ID[iRay*3+0] = Eps;             // Eps_1 = Eps

            }
            else { // Eps <= 0.0:
                // Newton's method is used to improve RHS DS point
                // DS = DS - Eps/GradEpsDotDir
                sim->DS_I[iRay] -= Eps/GradEpsDotDir;
                // DS_2 = DS; Eps_2 = Eps
                sim->PosPr_ID[iRay*3+1] = sim->DS_I[iRay]; // DS_2 = DS
                sim->DirPr_ID[iRay*3+1] = Eps;        // Eps_2 = Eps

            }

            // DS = (DS_1 + DS_2)/2 
            sim->DS_I[iRay] = (sim->PosPr_ID[iRay*3+0] + sim->PosPr_ID[iRay*3+1])*0.5;


            return; //------------------------------------------------------->>

        } // if (fabs(Eps) > sim->prm->TolEps) 

        else { // fabs(Eps) <= sim->prm->TolEps: REFLECTION
            // The critical surface found. Reflecting the ray.

            CLEAR_BIT(PENETR, iRay);   // Clear the penetration condition
            SET_BIT(WASREFL, iRay);  // Mark the ray as "was Snell reflected"

            //sim->DS_I[iRay] = sim->PosPr_ID[iRay*3+2];      // Restore old DS !Not needed!
            HalfDS = 0.5*sim->DS_I[iRay];

            GradEps2 = dot_product(GradEps_D, GradEps_D);
            LCosAl = -GradEpsDotDir/GradEps2;


            DirVert_D[0] = -LCosAl*GradEps_D[0]; // Here v_||
            DirVert_D[1] = -LCosAl*GradEps_D[1]; // Here v_||
            DirVert_D[2] = -LCosAl*GradEps_D[2]; // Here v_||

            //
            // Reflection by the Snell law:
            // v = v - 2*v_||
            //
            sim->Dir_ID[iRay*3+0] = sim->Dir_ID[iRay*3+0] - 2.0*DirVert_D[0];
            sim->Dir_ID[iRay*3+1] = sim->Dir_ID[iRay*3+1] - 2.0*DirVert_D[1];
            sim->Dir_ID[iRay*3+2] = sim->Dir_ID[iRay*3+2] - 2.0*DirVert_D[2];


            return; //------------------------------------------------------>>

        } // fabs(Eps) < TolEps
    } // if (sim->Flags_I[iRay] & PENETR) // Condition "Critical surface penetration"




    if (Eps < 0.0) { // Eps < 0.0: Penetration!
        //if (Eps < -TolEps)  // Eps < 0.0: Penetration!

        SET_BIT(PENETR, iRay); /* Mark the ray as penetrated past crit. surf. */

        //
        // Revert r and v to previous integer point (i.e. before penetration)
        // 
        sim->Dir_ID[iRay*3+0] = sim->DirPr_ID[iRay*3+0];
        sim->Dir_ID[iRay*3+1] = sim->DirPr_ID[iRay*3+1];
        sim->Dir_ID[iRay*3+2] = sim->DirPr_ID[iRay*3+2];

        sim->Pos_ID[iRay*3+0] = sim->PosPr_ID[iRay*3+0];
        sim->Pos_ID[iRay*3+1] = sim->PosPr_ID[iRay*3+1];
        sim->Pos_ID[iRay*3+2] = sim->PosPr_ID[iRay*3+2];
        //
        // Use sim->PosPr_ID and sim->DirPr_ID as memory for Cnt and 
        // for left and right values of DS and Eps
        //
        sim->PosPr_ID[iRay*3+0] = 0.0;             // Store here DS_1
        sim->PosPr_ID[iRay*3+1] = 2.0*sim->DS_I[iRay];  // Store here DS_2
        sim->PosPr_ID[iRay*3+2] = sim->DS_I[iRay];      // Store here old DS
        sim->DirPr_ID[iRay*3+0] = 1.0;             // Store here Eps_1
        sim->DirPr_ID[iRay*3+1] = Eps;             // Store here Eps_2
        sim->DirPr_ID[iRay*3+2] = sim->prm->CntMax;     // Store here iteration counter


        sim->DS_I[iRay] -= Eps/GradEpsDotDir; // Newton's root 1st approximation


        return; //--------------------------------------------------------->>
    } // (Eps < 0.0)  // Penetration!


    //
    // Calculate main ray parameters
    //

    // Restore the original Position (at an integer point): 
    PosHalfBk_D[0] = sim->Pos_ID[iRay*3+0] - sim->Dir_ID[iRay*3+0]*HalfDS; 
    PosHalfBk_D[1] = sim->Pos_ID[iRay*3+1] - sim->Dir_ID[iRay*3+1]*HalfDS;
    PosHalfBk_D[2] = sim->Pos_ID[iRay*3+2] - sim->Dir_ID[iRay*3+2]*HalfDS;


    //GradEps_D = -GradRho_DI(:,iRay)*RhoCrInv
    GradEps_D[0] = -sim->GradRho_ID[iRay*3+0]*sim->prm->RhoCrInv;
    GradEps_D[1] = -sim->GradRho_ID[iRay*3+1]*sim->prm->RhoCrInv;
    GradEps_D[2] = -sim->GradRho_ID[iRay*3+2]*sim->prm->RhoCrInv;






    GradEps2 = dot_product(GradEps_D, GradEps_D);
    GradEpsDotDir = dot_product(GradEps_D, &sim->Dir_ID[iRay*3+0]);

    LCosAl = -GradEpsDotDir/GradEps2;




    DirVert_D[0] = -LCosAl*GradEps_D[0]; // Here v_||
    DirVert_D[1] = -LCosAl*GradEps_D[1]; // Here v_||
    DirVert_D[2] = -LCosAl*GradEps_D[2]; // Here v_||



    EpsHalfBk = Eps - GradEpsDotDir*HalfDS;

    //Curv = (cHalf*HalfDS/EpsHalfBk)**2* &
    //       (GradEps2 - GradEpsDotDir**2)
    /* Explanation of the simplification:
     *  (a x b)^2 = a^2*b^2 - (a dot b)^2
     * if |v|=1, (a x v)^2 = a^2 - (a dot v)^2 */
    Curv = pow(0.5*HalfDS/EpsHalfBk,2) * (GradEps2 - pow(GradEpsDotDir,2));


    if (BIT_IS_OFF(SHARP, iRay)) { /* If the ray is NOT sharp */

        //
        // Check if the trajectory curvature is too sharp to meet the Tolerance
        // If so, reduce the DS step  for the iRay-th ray and leave it until
        // the next call.
        //
        //
        if (BIT_IS_ON(WASREFL, iRay)) {
            CLEAR_BIT(WASREFL, iRay);  // Clear the reflection condition
            // and do not test for curvature!
        }
        else {

            if (Curv >=  sim->prm->Tol2) { // If curvature is too high, decrease step 
                sim->DS_I[iRay] = sim->DS_I[iRay]/(2.0*sqrt(Curv/(sim->prm->Tol2)));

                // r = r_hb
                sim->Pos_ID[iRay*3+0] = PosHalfBk_D[0];
                sim->Pos_ID[iRay*3+1] = PosHalfBk_D[1];
                sim->Pos_ID[iRay*3+2] = PosHalfBk_D[2];
                return;    //---------------------------------------------------->>
            } // if (Curv >=  sim->prm->Tol2) 

        }
        //
        // Test if some of the next points can get into the prohibited part of 
        // space with "negative" dielectric permittivity
        //

        if (GradEpsDotDir*HalfDS <= -cThird*EpsHalfBk) {
            //
            // Mark the ray as steep; memorize the distance to the critical 
            // surface; reduce step
            //
            SET_BIT(SHARP, iRay); /* Mark simulation ray "Sharp" */
            sim->DistToCrSurf_I[iRay] = EpsHalfBk/sqrt(GradEps2);
            sim->DS_I[iRay] =  0.5*(sim->prm->Tol)*sim->DistToCrSurf_I[iRay];

            //Pos_DI(:,iRay) = PosHalfBk_D
            sim->Pos_ID[iRay*3+0] = PosHalfBk_D[0];
            sim->Pos_ID[iRay*3+1] = PosHalfBk_D[1];
            sim->Pos_ID[iRay*3+2] = PosHalfBk_D[2];



            return;    //----------------------------------------------------->>
        } /* if (GradEpsDotDir*HalfDS <= -cThird*EpsHalfBk) */
    } /* if (SHARP & (~sim->Flags_I[iRay]))  // If the ray is NOT sharp */

    //
    // In normal mode (i.e. no penetration)
    // Save previous values of Pos and Dir before they will be changed
    //
    if (BIT_IS_OFF(PENETR, iRay)) { /* i.e. if not "penetration" */
        /* Save Pos */
        sim->PosPr_ID[iRay*3+0] = sim->Pos_ID[iRay*3+0];
        sim->PosPr_ID[iRay*3+1] = sim->Pos_ID[iRay*3+1];
        sim->PosPr_ID[iRay*3+2] = sim->Pos_ID[iRay*3+2];
        /* Save Dir */
        sim->DirPr_ID[iRay*3+0] = sim->Dir_ID[iRay*3+0];
        sim->DirPr_ID[iRay*3+1] = sim->Dir_ID[iRay*3+1];
        sim->DirPr_ID[iRay*3+2] = sim->Dir_ID[iRay*3+2];
    }


    SolDist = v3magn(&sim->Pos_ID[iRay*3+0]); /* Solar Distance */

    //
    // Either:
    // - switch to opposite branch of parabola;
    // - or make a Boris step
    //
    if ((GradEpsDotDir*HalfDS <= -cThird*EpsHalfBk) 
            || (sim->DS_I[iRay] < (sim->prm->AbsMinStep))) {

        // Switch to the opposite branch of parabolic trajectory
        //
        // When a next step can drive the ray into the area with the
        // plasma density greater then its critical value, then a special 
        // technique of "parabolic ray reflection" is employed. 
        // It can be shown that a ray trajectory in the medium with constant
        // gradient of dielectric permittivity is a parabola. If the step sim->DS_I
        // is small enough we can assume the grad \epsilon constant and hence
        // assume that the ray approaches the critical surface in a parabolic 
        // path.
        // We calculate the parameters of the parabola, namely:
        // StepX_D -- vector along the -grad \epsilon, with the length equal 
        //     to the distance from the PosHalfBk_D to the parabola extremum;  
        // StepY_D -- vector perpendicular to StepX_D, lying inside of the
        //     parabola plane, pointing at the opposite parabola branch, and with
        //     the length equal the distance from PosHalfBk_D to the 
        //     opposite branch of parabola.
        // The parabolic reflection is just replacement of the sim->Pos_ID with
        // the symmetric point at the opposite branch of parabola, and changing
        // the "incident" direction sim->Dir_ID for the "departing" ray direction
        // according to Snell law.
        //

        // Here c; |c| = sin \alpha
        StepY_D[0] = sim->Dir_ID[iRay*3+0] - DirVert_D[0];
        StepY_D[1] = sim->Dir_ID[iRay*3+1] - DirVert_D[1];
        StepY_D[2] = sim->Dir_ID[iRay*3+2] - DirVert_D[2];



        sim->Dir_ID[iRay*3+0] = sim->Dir_ID[iRay*3+0] - 2.0*DirVert_D[0];
        sim->Dir_ID[iRay*3+1] = sim->Dir_ID[iRay*3+1] - 2.0*DirVert_D[1];
        sim->Dir_ID[iRay*3+2] = sim->Dir_ID[iRay*3+2] - 2.0*DirVert_D[2];


        //
        //   We need to have |Step_Y| = 4 * sin(\alpha)*cos(\alpha)*EpsHalfBk*L
        //   Now the direction of Step_Y is right, 
        //   the length of it is equal to sin(\alpha)
        //   Multiply it by L*Cos(\alpha)*EpsHalfBk
        //



        StepY_D[0] = 4.0*StepY_D[0]*LCosAl*EpsHalfBk; 
        StepY_D[1] = 4.0*StepY_D[1]*LCosAl*EpsHalfBk; 
        StepY_D[2] = 4.0*StepY_D[2]*LCosAl*EpsHalfBk; 



        sim->Pos_ID[iRay*3+0] = PosHalfBk_D[0] + StepY_D[0];
        sim->Pos_ID[iRay*3+1] = PosHalfBk_D[1] + StepY_D[1];
        sim->Pos_ID[iRay*3+2] = PosHalfBk_D[2] + StepY_D[2];

        //
        //   Step_X is in the direction of -grad \epsilon, whose vector is of 
        //   the length of 1/L
        //   The length of Step_X is cos^2(\alpha)*L*EpsHalfBk 
        //   Thus,
        //



        /*
         * Brightness temperature computation over parabola
         */

        if (sim->rtmode == 2 || sim->rtmode == 3){
            /* Electron temperature Te */ 
            if (SolDist > (sim->prm->r_chro_cor))
                Te = sim->prm->Te_corona_K;
            else
                Te = sim->prm->Te_chromo_K;
        }
        if (sim->rtmode == 2){ /* Calculate the brightness temperature, Tbr */
            sim->Tbr_I[iRay] = calc_tbr(sim->prm, sim->Rho_I[iRay], 
                    Te, sim->DS_I[iRay], sim->Tbr_I[iRay], &dtau);
            sim->OpDepth_I[iRay] += dtau;
        } 
        else if(sim->rtmode==3){ /* Calculate Stokes I, Q, U, and V components of 
                                        * brightness temperature                        */
            calc_tbriquv(sim->prm, &sim->Dir_ID[iRay*3+0], sim->Rho_I[iRay], 
                    &sim->Bfield_ID[iRay*3+0], Te, sim->DS_I[iRay], 
                    &sim->TbrIQUV_IP[iRay*4+0], &dtau, 0);
            sim->OpDepth_I[iRay] += dtau;
        }

    } // if ((GradEpsDotDir*HalfDS <= -cThird*EpsHalfBk) 
    //  || (sim->DS_I[iRay] < (sim->prm->AbsMinStep)))


    else { // Boris' step
        //
        // Make a step using Boris' algorithm
        //

        // grad(n) = grad(eps(i+1/2))/(2*eps(i+1/2))
        Coef = 0.5*HalfDS/(1.0 - Rho2RhoCr);

        RelGradRefrInx_D[0] = Coef*GradEps_D[0];
        RelGradRefrInx_D[1] = Coef*GradEps_D[1]; 
        RelGradRefrInx_D[2] = Coef*GradEps_D[2]; 

        Dir_D[0] = sim->Dir_ID[iRay*3+0];
        Dir_D[1] = sim->Dir_ID[iRay*3+1];
        Dir_D[2] = sim->Dir_ID[iRay*3+2];

        //Omega_D = cross_product(RelGradRefrInx_D, Dir_DI(:,iRay))
        cross_product(RelGradRefrInx_D, Dir_D, Omega_D);

        //Dir1_D = Dir_DI(:,iRay) + cross_product(Dir_DI(:,iRay), Omega_D)
        cross_product(Dir_D, Omega_D, Xprod_D);
        Dir1_D[0] = Dir_D[0] + Xprod_D[0];
        Dir1_D[1] = Dir_D[1] + Xprod_D[1];
        Dir1_D[2] = Dir_D[2] + Xprod_D[2];

        //Omega_D = RelGradRefrInx_D x Dir1_D
        cross_product(RelGradRefrInx_D, Dir1_D, Omega_D);

        // v1 = v + v x Omega
        cross_product(Dir_D, Omega_D, Xprod_D);
        Dir1_D[0] = Dir_D[0] + Xprod_D[0];
        Dir1_D[1] = Dir_D[1] + Xprod_D[1];
        Dir1_D[2] = Dir_D[2] + Xprod_D[2];

        //Curv1 = pow(Omega_D[0],2) + pow(Omega_D[1],2) + pow(Omega_D[2],2);
        Curv1 = dot_product(Omega_D, Omega_D);
        Coef = 2.0/(1.0 + Curv1);
        //Dir_DI(:,iRay) = Dir_DI(:,iRay) + Coef*cross_product(Dir1_D, Omega_D)
        cross_product(Dir1_D, Omega_D, Xprod_D);



        sim->Dir_ID[iRay*3+0] = sim->Dir_ID[iRay*3+0] + Coef*Xprod_D[0];
        sim->Dir_ID[iRay*3+1] = sim->Dir_ID[iRay*3+1] + Coef*Xprod_D[1];
        sim->Dir_ID[iRay*3+2] = sim->Dir_ID[iRay*3+2] + Coef*Xprod_D[2];


        //Pos_DI(:,iRay) = Pos_DI(:,iRay) + Dir_DI(:,iRay)*HalfDS
        sim->Pos_ID[iRay*3+0] += sim->Dir_ID[iRay*3+0]*HalfDS;
        sim->Pos_ID[iRay*3+1] += sim->Dir_ID[iRay*3+1]*HalfDS;
        sim->Pos_ID[iRay*3+2] += sim->Dir_ID[iRay*3+2]*HalfDS;

        /*
         * Brightness temperature computation over Boris step
         */

        if (sim->rtmode == 2 || sim->rtmode == 3){
            /* Electron temperature Te */ 
            if (SolDist > (sim->prm->r_chro_cor))
                Te = sim->prm->Te_corona_K;
            else
                Te = sim->prm->Te_chromo_K;
        }
        if (sim->rtmode == 2){ /* Calculate the brightness temperature, Tbr */        
            sim->Tbr_I[iRay] = calc_tbr(sim->prm, sim->Rho_I[iRay], 
                    Te, sim->DS_I[iRay], sim->Tbr_I[iRay], &dtau);
            sim->OpDepth_I[iRay] += dtau;
        } 
        else if(sim->rtmode==3){ /* Calculate Stokes I, Q, U, and V components of 
                                        * brightness temperature                        */

            /* printf("********* THREAD NUMBER = %d\n", sim->thread_num); */
            /* if (sim->thread_num == 0) { */
            /*   printf("\n"); */
            /*   printf("In adv_beam: \nBx = %g, By = %g, Bz = %g, B = %g\n", */
            /* 	 sim->Bfield_ID[iRay*3+0], sim->Bfield_ID[iRay*3+1], sim->Bfield_ID[iRay*3+2], */
            /* 	 sqrt(pow(sim->Bfield_ID[iRay*3+0],2) + pow(sim->Bfield_ID[iRay*3+1],2) + */
            /* 	      pow(sim->Bfield_ID[iRay*3+2],2))); */
            /* } */

            calc_tbriquv(sim->prm, &sim->Dir_ID[iRay*3+0], sim->Rho_I[iRay], 
                    &sim->Bfield_ID[iRay*3+0], Te, sim->DS_I[iRay], 
                    &sim->TbrIQUV_IP[iRay*4+0], &dtau, 0);
            /* printf("Returned:sim->TbrIQUV_IP = %e %e %e %e\n",  */
            /*       sim->TbrIQUV_IP[iRay*3+0], sim->TbrIQUV_IP[iRay*3+1], */
            /*       sim->TbrIQUV_IP[iRay*3+2], sim->TbrIQUV_IP[iRay*3+3]); */
            /* sim->OpDepth_I[iRay] += dtau; */

        }

    } // else  // Boris' step



    //
    //   The code below makes gradual increases of the DS up to the value
    // specified in DSNew. The smooth step increase is required so as not to
    // get into the space behind the critical surface, stepping with DS that
    // instantly changes from a very little size to the normal DSNew length.
    // DS is usually reduced in a close vicinity of the critical surface,
    // where the ray is travelling along a very sharp curve with high curvature.
    // For many rays it means fractioning of the DS down several orders of 
    // magnitude, therefore the new step trial should sim->start from a bigger step
    // of the same order of magnitude.
    //   This problem is solved using a non-linear difference equation:
    //           Y(i+1) = [2 - Y(i)/X(i)]*Y(i),
    // where X(i) is the desired final DS length from DSNew, and
    // Y(i) is the actual DS length. A simple analysis of the equation shows
    // that, when Y is small compared to X, the next value of Y will be almost 
    // 2*X, so the DS would grow in a geometrical progression. However, as
    // Y approaches X, its growth rate becomes slower. However, Y always reaches
    // X in several steps. One can check that for Y = X the next value of Y is 
    // always that of X.
    // 
    if (BIT_IS_OFF(SHARP, iRay)) { /* If the ray is NOT sharp */
        //
        // For shallow rays the DS is increased unconditionally
        //
        if (sim->DS_New_I[iRay] > sim->DS_I[iRay])
            sim->DS_I[iRay] = (2.0 - sim->DS_I[iRay]/sim->DS_New_I[iRay])*
                               sim->DS_I[iRay];
        else
            sim->DS_I[iRay] = sim->DS_New_I[iRay];


    } /* if (BIT_IS_OFF(SHARP, iRay))  // If the ray is NOT sharp */
    else { // Steep ray
        //
        // If the iRay-th ray is marked as sharp (i.e. "sloping" or 
        // "oblique") then the code below increases its DS only if the 
        // current distance to the critical surface, calculated as 
        //     \epsilon / grad \epsilon, 
        // is greater than simulation distance value saved along with marking the ray 
        // as steep in the DitToCrSurf_I. 
        //   This can happen in two cases: 
        // (1) either the ray was "parabola reflected" and after several steps it
        // went away from the surface by the distance where the parabola 
        // switching occurred; 
        // (2) or the ray is not steep any more because the current DS is so 
        // small, that the next step does not penetrate the critical surface.
        //   The ray is reverted to "gentle" or "shallow"
        //
        if (EpsHalfBk > sim->DistToCrSurf_I[iRay]*sqrt(GradEps2)) {
            CLEAR_BIT(SHARP, iRay);  /* Mark simulation ray NOT sharp */
            if (sim->DS_New_I[iRay] > sim->DS_I[iRay])
                sim->DS_I[iRay] = (2.0 - sim->DS_I[iRay]/sim->DS_New_I[iRay])*
                                   sim->DS_I[iRay];
            else
                sim->DS_I[iRay] = sim->DS_New_I[iRay];

        } // if (EpsHalfBk > sim->DistToCrSurf_I[iRay]*sqrt(GradEps2)) 
    } // else  // Sharp rays

    // if (iRay==0)
    //   printf("blockIdx.x=%d, threadIdx.x=%d, iRay=%d, iIter=%d, itrj=%d\n", \
    // 	     blockIdx.x, threadIdx.x, iRay, iIter, itrj);


    if (BIT_IS_ON(TRACED, iRay)) {
      /*
       * Save positions of particular rays
       *
       * Find the index of the stored ray using binary search 
       * in the ordered table TracedPts_I[nTracedPts]
       */
      i0 = 0;
      i1 = sim->nTracedPts - 1;

      while(i0 < i0) {
	im = (i0 + i1)/2;     /* Point at the middle of table */
	if (sim->TracedPts_I[im] < iRay)
	  i0 = im + 1;
	else 
	  i1 = im;
      }
      iTracedPts = i0;
      /* Since we are sure iRay is in the table, we do not spend valuable
       * time on the "deferred test for equality": */
      /* if ((i0 == i1) and (TracedPts_I[i0] == iRay)) iTracedPts = i0; */

      /*
       * Trajectories_I[iTracedPts,iIter,:] = Pos_ID
       */
      itrj = (sim->nIter*iTracedPts + iIter)*3;
      sim->Trajectories_I[itrj+0] = sim->Pos_ID[iRay*3+0];
      sim->Trajectories_I[itrj+1] = sim->Pos_ID[iRay*3+1];
      sim->Trajectories_I[itrj+2] = sim->Pos_ID[iRay*3+2];
    }

    return;
}





/*
* Multiplies a matrix and a vector
*/
void mvmul1(int N, double *a, double x[], double b[]){
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


/*
 * At compiletime the following preprocessor symbols 
 * in the command line as -D<symbol> can be used:
 * TBR:  compute brightness temperature, no polarization due to magnetic field
 * TBRIV: compute brightness temperature for I and V Stokes parameters,
 *        assuming polarization due to magnetic field
 * TRAJ:  Store trajectory points for specified rays
 * DEB:   Generate debugging output
 */




void advance_beam(Simulation *sim_h, Simulation *sim_d) {

  //
  //   The subroutine beam_path() makes ray tracing and emissivity integration 
  // along ray paths.
  // 
  // Author:
  //   Leonid Benkevitch
  // Revisions:
  // 2007-Jul-17 Initially written in Fortran-90 and tested for frequencies 
  //             around 42.3 MHz
  // 2009-Jan-05 Rewritten in GNU C 
  // 2009-Jan-09 Added ray simulation.stopping at the solar surface for 
  //             frequencies above
  //             ~127 MHz
  // 2009-Jan-13 The array simulation.stopping correction was unsuccessful, 
  //             so then program
  //             reverted to the state before 2009-Jan-09. This requires more 
  //             work.
  // 2009-Jul-01 For calculating the brightness temperature in the Intgrl_I[],
  //             the parameters Freq and simulation.OpDepth_I[] are added 
  // 2011-Mar-15 Added to Python extension module rtcore (part of the 
  //             raytrace pachage). The name beam_path() changed to 
  //             advance_beam().
  // 2011-Apr-01 Introduced dynamic linking for plasma_params() in 
  //             rtcoremodule.c. The pointer to plasma_params() is passed 
  //             as a parameter.
  // 2015-Jun-05 Added saving trajectories of particular rays.
  //             
  
  /*
   * Revision Author:
   *   Mark Benjamin
   * 2011-Jun-24 This file was modified to include the necessary parameters
   * for threading. Two arguments were added, simulation.start and 
   * simulation.stop so that the number
   * of rays to be processed could be split up. All of the for loops that
   * used to go from 0 to simulation.nRay now go from simulation.start to 
   * simulation.stop.
   */
    int iIter=0;
    int nRay;
    int n_blocks;
    int n_threads = 64;

    //int *ray2trj[10] = malloc(sizeof(int[20][10]));

    const int ldbar_len=50;
    int i;

    //struct plasma_constants *p_consts_h, *p_consts_d;

    // Linear systems A*a = rhsa and A*b = rhsb
    double const R0 = 6.955e5;
    curandState *devStates;


    void *handle = dlopen("./plasma_parameters.so",RTLD_LAZY);

    void (*intro)(int,int,Simulation *,Simulation *);
    void (*cleanup)(void);

    if(!handle)
        errx(1,dlerror());

    intro = (void (*)(int,int,Simulation *,Simulation *))dlsym(handle, "intro");
    cleanup = (void (*)(void))dlsym(handle, "cleanup");

    if(!intro)
        errx(1,dlerror());



    nRay = sim_h->nRay;
    
    n_blocks = nRay/n_threads + 1;

    printf("nRay=%d, n_threads=%d, \n", nRay, n_threads, n_blocks);

    // Sets up all the parameters needed for the plasma parameters array
    CUDA_CALL(cudaMalloc((void**)&devStates,sizeof(curandState)*nRay));


    setup_randgen<<<n_blocks,n_threads>>>(devStates, nRay);

    check_call("Kernel Setup Call");


    printf("\n");
    //printf("A[3,3] = %g\n",p_consts_h->A[15]);
    for(iIter=0; iIter<sim_h->nIter; iIter++){

        if(sim_h->scattering){
            scatter<<<n_blocks,n_threads>>>(sim_d, devStates);
            check_call("Scatter Call");
        }

	printf("After call 'scatter'... sim_h->scattering=%d \n", 
	       sim_h->scattering);

        first<<<n_blocks,n_threads>>>(sim_d);
        check_call("Kernel First Call");

	printf("After call 'first'...\n");

        intro(n_blocks, n_threads, sim_d, sim_h);
        check_call("Kernel Plasma Call");

	printf("After call 'intro'...\n");

	second<<<n_blocks,n_threads>>>(sim_d, iIter);
        check_call("Kernel Second Call");

	printf("After call 'second'...\n");

	//printf("iIter=%d\n", iIter);

        // printf("\rCompletion: [");
        // for(i=0;i<ldbar_len*iIter/sim_h->nIter;i++){
        //     printf("=");
        // }
        // printf(">");
        // for(i=0; i<ldbar_len-ldbar_len*iIter/sim_h->nIter-1; i++){
        //     printf(" ");
        // }
        // printf("]  %2.0f%% Complete",100*iIter/(float)sim_h->nIter);

    }

    printf("\n\n");

    cleanup();

    dlclose(handle);

    cudaFree(devStates);
    // cudaFree(p_consts_d);
    // free(p_consts_h);

}

