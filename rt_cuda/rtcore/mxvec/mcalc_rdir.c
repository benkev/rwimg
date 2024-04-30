/*
 * The SGI coordinate system: the center is at the sun center; X axis
 * points at the earth; Z is normal to the ecliptic plane; Y makes right-
 * handed system (XYZ). (Roughly speaking, X points to us, Y - to the left,
 * and Z upwards).
 */
/* 
 * For the given observer SGI position (obs) and the target 
 * coordinates in the image plane specified in the input arrays xruler[nx]
 * and yruler[ny], where ij-th pixel coordinates are 
 *    (yruler[i],xruler[j]), and nx*ny=nRay, 
 * the function calculates nRay unity-length direction 
 * vectors dir[nRay][3] (dir is in SGI). They point from obs to the pixels.
 */
#include <stdio.h>
#include "raytrace_cuda.h"

void mcalc_rdir(double obs[3], double xruler[], double yruler[], 
		int nx, int ny,	double dir[][3]) {
  double nobs[3];       /* obs normalized; nobs is normal to image plane */
  double rpix[3];       /* SGI coordinates of targ */
  double slope[3];      /* vector from obs to rpix */
  double zaxis[3] = {0., 0., 1.}; /* SGI Z axis */
  double tau[3], xi[3]; /* Local basis vectors in the image plane */ 
  double mobs;          /* obs magnitude */
  double mtau;          /* magnitude of the vector collinear to tau */
  double mslope;        /* magnitude of the vector collinear to dir */
  double x, y;          /* (x,y) run over all the pixel coordinates */
  int i, j, k;
  
  mobs = v3magn(obs);  /* Distance from the SGI origin to the observer */
  
  /* nobs[3] is a unity-length vector pointing from  
   * the center of image plane (and the origin of the SGI)
   * to the observer */
  nobs[0] = obs[0]/mobs;
  nobs[1] = obs[1]/mobs;
  nobs[2] = obs[2]/mobs;

  /* Calculate in the SGI the unit vectors xi and tau that make up an 
   * ortogonal basis inside of the image plane, so that every its pixel
   * can be expressed in SGI as a linear combination
   *     x*tau + y*xi,
   * where (x,y) are plane pixel coordinates in the image plane. */

  cross_product(zaxis, nobs, tau);    /* tau' = Z x nobs */

  mtau = v3magn(tau);    /* mtau is the magitude of tau' */
  tau[0] = tau[0]/mtau;
  tau[1] = tau[1]/mtau;
  tau[2] = tau[2]/mtau;

  cross_product(nobs, tau, xi);       /* xi = nobs x tau */

  /* k is the index into the dir[nRay][3] array 
   * ordered as 3D [ny][nx][3] C-language array */
  k = 0; 
  for (i = 0; i < ny; i++) {               /* !!!!! Was i < nx  - ERROR!!! */
    y = yruler[i];
    for (j = 0; j < nx; j++) {             /* !!!!! Was j < ny  - ERROR!!! */
      x = xruler[j];
      /* rpix[3] is ij-th pixel coordinates in SGI */
      rpix[0] = x*tau[0] + y*xi[0];
      rpix[1] = x*tau[1] + y*xi[1];
      rpix[2] = x*tau[2] + y*xi[2];
      
      slope[0] = rpix[0] - obs[0];
      slope[1] = rpix[1] - obs[1];
      slope[2] = rpix[2] - obs[2];
      
      mslope = v3magn(slope);
      dir[k][0] = slope[0]/mslope;
      dir[k][1] = slope[1]/mslope;
      dir[k][2] = slope[2]/mslope;
      k++;
    } /* for (j = 0; j < ny; j++) */ 
  } /* for (i = 0; i < nx; i++) */

}
 
