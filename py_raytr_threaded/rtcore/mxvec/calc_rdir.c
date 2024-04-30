/*
 * The SGI coordinate system: the center is at the sun center; X axis
 * points at the earth; Z is normal to the ecliptic plane; Y makes right-
 * handed system (XYZ). (Roughly speaking, X points to us, Y - to the left,
 * and Z upwards).
 */
/* 
 * For the given observer SGI position (obs) and the target (tau,xi) 
 * coordinates in the image plane, calculates the unity-length direction 
 * vector dir (dir is in SGI). It points from obs to targ.
 */
#include <stdio.h>
#include "raytrace.h"

void calc_rdir(double obs[3], double targ[2], double dir[3]) {
  double nobs[3];       /* obs normalized; nobs is normal to image plane */
  double rpix[3];       /* SGI coordinates of targ */
  double slope[3];      /* vector from obs to rpix */
  double zaxis[3] = {0., 0., 1.}; /* SGI Z axis */
  double tau[3], xi[3]; /* Local basis vectors in the image plane */ 
  double mobs;          /* obs magnitude */
  double mtau;          /* magnitude of the vector collinear to tau */
  double mslope;        /*  magnitude of the vector collinear to dir */
  
  mobs = v3magn(obs);  /* Distance from the SGI origin to the observer */
  nobs[0] = obs[0]/mobs;
  nobs[1] = obs[1]/mobs;
  nobs[2] = obs[2]/mobs;

  cross_product(zaxis, nobs, tau);    /* tau = Z x nobs */

  mtau = v3magn(tau);
  tau[0] = tau[0]/mtau;
  tau[1] = tau[1]/mtau;
  tau[2] = tau[2]/mtau;

  cross_product(nobs, tau, xi);       /* xi = nobs x tau */

  /* 
   * rpix is the target coordinates in SGI
   * Conversion targ -> rpix:
   * rpix = targ[tau]*tau + targ[xi]*xi 
   */
  rpix[0] = targ[0]*tau[0] + targ[1]*xi[0];
  rpix[1] = targ[0]*tau[1] + targ[1]*xi[1];
  rpix[2] = targ[0]*tau[2] + targ[1]*xi[2];
  
  slope[0] = rpix[0] - obs[0];
  slope[1] = rpix[1] - obs[1];
  slope[2] = rpix[2] - obs[2];

  mslope = v3magn(slope);
  dir[0] = slope[0]/mslope;
  dir[1] = slope[1]/mslope;
  dir[2] = slope[2]/mslope;
}
 
