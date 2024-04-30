/*
 * The SGI coordinate system: the center is at the sun center; X axis
 * points at the earth; Z is normal to the ecliptic plane; Y makes right-
 * handed system (XYZ). (Roughly speaking, X points to us, Y - to the left,
 * and Z upwards).
 */

/* 
 * For a given observer SGI position (obs) and the direction (dir) of the
 * straight ray, calculates the position (pos) of the ray intersection with
 * the sphere of the radius (rsph) drawn around the SGI center (0,0,0).
 * Of the two line-sphere intersection points, the one closer to the
 * observer position (obs) is returned. 
 * If rsph = 1., it is the solar surface.
 *
 * The function returns the number of intersection points. 
 * If the ray does not intersect the sphere, 0 is returned. In case of 
 * tangential ray, 1 is returned. For two intersection points, 2 is returned.
 *
 * The problem is reduced to solving the quadratic equation:
 *
 *  dir^2*dist^2 + 2(dir.obs)*dist + (obs^2 - rsph^2) = 0
 *  
 */

#include <math.h>
#include <stdio.h>
#include "raytrace_cuda.h"

int calc_rpos(double obs[3], double dir[3], double rsph, double pos[3]) {
  double rsph2 = rsph*rsph;
  double vxo[3], vdo;
  int retval = 0;  /* Returned value */
  double discr;    /* Discriminant of quadratic equation */
  double dist;     /* Distance from observer to sphere along the ray */

  cross_product(dir, obs, vxo);          /* vxo = dir x obs */
  vdo = dot_product(dir, obs);           /* vdo = dir . obs */
  discr = rsph2 - dot_product(vxo, vxo); /* discriminant */
  if (discr < 0) {
    pos[0] = pos[1] = pos[2] = 0.0;
    return 0;
  }
  if (discr == 0) retval = 1; /* A single intersection point - tangential */
  else            retval = 2; /* Two intersection points - intersection */
  discr = sqrt(discr);
  dist = -vdo - discr; /* Distance from observer to sphere along the ray */
  
  /* The closest to observer point on the surface of the sphere */
  pos[0] = obs[0] + dir[0]*dist;
  pos[1] = obs[1] + dir[1]*dist;
  pos[2] = obs[2] + dir[2]*dist;
  
  return retval; /* 1: tangential; or 2: two-point intersection */
}
