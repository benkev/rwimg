#include <stdio.h>
//#include <string.h>
//#include <stdlib.h>
#include <math.h>
#include "raytrace.h"
 
inline void calc_tbriquv(
                double Param_I[],
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
 *            |al  be   0  de|
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


  /* Access parameters by structure field names */
  struct param *prm = (void *) Param_I; 

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
    fprintf(stderr,"nan\n");
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


  if (thread_num == 0) {
  /* printf("Bx = %g, By = %g, Bz = %g, B = %g\n", */
  /* 	 Bfield_D[0], Bfield_D[1], Bfield_D[2], */
  /* 	 sqrt(pow(Bfield_D[0],2) + pow(Bfield_D[1],2) + pow(Bfield_D[2],2))); */
  /* printf("Bdirx= %g, Bdiry = %g, Bdirz = %g, |Bdir|^2 = %g\n", */
  /* 	 Bdir_D[0], Bdir_D[1], Bdir_D[2], dot_product(Bdir_D,Bdir_D)); */
  /* printf("Bfield_abs = %g, costh = %g, th = %d, sinth2 = %e, Ne = %g\n", */
  /* 	 Bfield_abs, costh, (int)(acos(costh)*180./pi), sinth2, Ne); */
  /* printf("Te = %e, nu_eff = %g, Ne = %g\n", */
  /* 	 Te, nu_eff, Ne); */
  /* printf("u_plasma = %e, v_plasma = %e, w_plasma = %e\n",  */
  /* 	 u_plasma, v_plasma, plasma); */
  /* printf("P_plasma = %e, Q_plasma = %e, R_plasma = %e, V_plasma = %e\n",  */
  /* 	 P_plasma, Q_plasma, R_plasma, V_plasma); */

  /* printf("dtau = %e, ds = %g, Rsun_cm=%g\n", *dtau, ds, prm->Rsun_cm); */
  printf("al=%e, be=%e, de=%e, mu=%e, et=%e \n", al, be, de, mu, et);
  /* printf("xi=%e, om=%e, xids=%e, omds=%e, chcs=%e\n", xi, om, xids, omds); */
  /* printf("=%e, sh=%e, cs=%e, sn=%e, chcs=%e\n", ch, sh, cs, sn, chcs); */
  printf("M11, M21, M31, M41 = %e %e %e %e\n", M11, M21, M31, M41);
  printf("M12, M22, M32, M42 = %e %e %e %e\n", M12, M22, M32, M42);
  printf("M13, M23, M33, M43 = %e %e %e %e\n", M13, M23, M33, M43);
  printf("M14, M24, M34, M44 = %e %e %e %e\n", M14, M24, M34, M44);
  /* printf("I,Q,U,V   = %e %e %e %e\n", I, Q, U, V); */
  /* printf("ImTe = %e, atn = %e\n", ImTe, atn); */
  /* printf("TbrIQUV_P = %e %e %e %e\n", TbrIQUV_P[0], TbrIQUV_P[1], */
  /*  	 TbrIQUV_P[2], TbrIQUV_P[3]); */
  } /* if (thread_num == 0) */
}
