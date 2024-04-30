/*
 * The SGI coordinate system: the center is at the sun center; X axis
 * points at the earth; Z is normal to the ecliptic plane; Y makes right-
 * handed system (XYZ). (Roughly speaking, X points to us, Y - to the left,
 * and Z upwards).
 */

/* 
 * For a given observer SGI position (obs) and nRay direction vectors,
 * dir[nRay][3], of the straight rays, calculates the positions,
 * pos[nRay][3] of the ray intersections with the sphere of the radius rsph
 * drawn around the SGI center (0,0,0).
 * Of the two line-sphere intersection points, the one closer to the
 * observer position (obs) is stored in pos[][3]. 
 * If rsph = 1., it is the solar surface.
 *
 * Upon return, i-th element of the array isec[nRay] contains the number of 
 * intersections of the i-th ray with the integration sphere.
 * Also, a zero triplet in obs[i][:] indicates that the i-th ray 
 * does not penetrate the sphere. The tangential point coordinats are saved 
 * in pos[][:], but for the same ray index isec[] contains 1.
 * Raytracing makes sense for only those rays whose isec == 2.  
 * The function returns the total number of intersection points as the
 * sum of those over all the rays. 
 * If the ray does not intersect the sphere, 0 is added to the total. 
 * In case of a tangential ray, 1 is added. For two intersection points, 
 * 2 is added. Only if the returned number is 2*nRay, all the rays 
 * penetrate the integration sphere.
 *
 * The problem is reduced to solving the quadratic equation:
 *
 *  dir^2*dist^2 + 2(dir.obs)*dist + (obs^2 - rsph^2) = 0
 *  
 */

#include <math.h>
#include <stdio.h>
#include "raytrace_cuda.h"

int mcalc_rsph(double obs[3], double dir[][3], double rsph, int nRay, 
	       short isec[], double pos[][3]) {
  double rsph2 = rsph*rsph;
  double vxo[3], vdo;
  int iRay;    /* Ray index */
  int isec_total = 0;   /* Returned value */
  double discr;         /* Discriminant of quadratic equation */
  double sqrt_discr;    /* sqrt(discr) */
  double dist;          /* Distance from observer to sphere along the ray */
  

  for (iRay = 0; iRay < nRay; iRay++) {
    cross_product(&dir[iRay][0], obs, vxo);          /* vxo = dir x obs */
    vdo = dot_product(&dir[iRay][0], obs);           /* vdo = dir . obs */
    discr = rsph2 - dot_product(vxo, vxo); /* discriminant */
    if (discr < 0) { /* No intersections */
      pos[iRay][0] = pos[iRay][1] = pos[iRay][2] = 0.0;
      isec[iRay] = 0;
    }
    else { /* Either 1 (tangent) or 2 (penetr) intersection points */ 
      sqrt_discr = sqrt(discr);
      /* Distance from observer to sphere along the ray */
      dist = -vdo - sqrt_discr; 
      /* The closest to observer point on the surface of the sphere */
      pos[iRay][0] = obs[0] + dir[iRay][0]*dist;
      pos[iRay][1] = obs[1] + dir[iRay][1]*dist;
      pos[iRay][2] = obs[2] + dir[iRay][2]*dist;
      if (discr == 0.0) { /* Unlikely */
	isec[iRay] = 1;   /* Tangential ray */
	isec_total += 1;
      }
      else {              /* i.e. discr > 0.0 */
	isec[iRay] = 2;   /* Penetrating ray */
	isec_total += 2;
      }
    }

  } /* for (iRay = 0; iRay < nRay; i++) */

  return isec_total; /* 1: tangential; or 2: two-point intersection */
}
