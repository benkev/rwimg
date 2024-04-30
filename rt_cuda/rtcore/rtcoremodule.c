/************************************************************
 * rtcoremodule.c                                           *
 *                                                          *
 * The Python extension module implementing the low level   *
 * ray tracing algorithm and service routines.              *
 * Makes use of the NumPy arrays.                           *
 *                                                          *
 * Created 23 March 2011 by Leonid Benkevitch               *
 * 2011-Apr-01 Introduced dynamic linking for               *
 *             plasma_parameters() sibroutine. The pointer  *
 *             is passed as a parameter to advance_beam.    * 
 *                                                          *
 ************************************************************/
#include <stdio.h>
#include <dlfcn.h>
#include <math.h>
#include <string.h>
#include <Python.h>
#include "numpy/arrayobject.h"
#include "simulation.h"
#include <time.h>
//#include <pthread.h>



/* 
 * This defintion is the number of threads the program will use 
 * in its calculation.
 * If set to 1 this will be completely serial.
 */
#define NUM_THREADS 1 

#define OBJ_TO_ARRAY(obj, npy_type, arr)			  	    \
  if (obj == NULL) return NULL;						    \
  arr = (PyArrayObject *) PyArray_FROM_OTF(obj, npy_type, NPY_INOUT_ARRAY); \
  if (arr == NULL) return NULL;

//
//#define OBJ_TO_DOUBLE_ARRAY(obj,arr)					\
//  if (obj == NULL) return NULL;					\
//  arr = (PyArrayObject *) PyArray_FROM_OTF(obj, NPY_DOUBLE, NPY_INOUT_ARRAY);\
//  if (arr == NULL) return NULL;
//
//#define OBJ_TO_INT_ARRAY(obj,arr)					\
//  if (obj == NULL) return NULL;					\
//  arr = (PyArrayObject *) PyArray_FROM_OTF(obj, NPY_INT, NPY_INOUT_ARRAY); \
//  if (arr == NULL) return NULL;
//

#define CHKSIZE_AND_GETPTR(arr, nrq, ty, ptr)				\
  {									\
    int nel = PyArray_SIZE(arr);					\
    if (nel != nrq) {							\
      PyErr_Format(PyExc_ValueError, "Wrong size of array " #arr ": %d; " \
		   "must be %d.", nel, nrq);				\
      return NULL;							\
    }									\
    ptr = (ty *) PyArray_DATA(arr);					\
  }






/*
 * Call from Python:
 * trace_beam(...)
 */

/*
 * This function performs the raytrace from start to finish, 
 * calling advance_beam at every step.
 * In this threaded version this acts only as a shell function, 
 * calling the threadFunc for every thread.
 */
PyObject *trace_beam(PyObject *self, PyObject *args) {


  /* Python objects; arrays to be converted to ndarrays 
   *
   * 25 pyobjects, 25 pyarrays, 25 arrays
   */


  PyObject *pyobjParam_I = NULL; 
  PyObject *pyobjPos_ID = NULL; 
  PyObject *pyobjDir_ID = NULL; 
  PyObject *pyobjDS_I = NULL; 
  PyObject *pyobjTbr_I = NULL; 
  PyObject *pyobjTbrIQUV_IP = NULL; 
  PyObject *pyobjOpDepth_I = NULL; 
  PyObject *pyobjFlags_I = NULL; 
  PyObject *pyobjRho_I = NULL; 
  PyObject *pyobjGradRho_ID = NULL; 
  PyObject *pyobjBfield_ID = NULL; 
  PyObject *pyobjPosPr_ID = NULL; 
  PyObject *pyobjDirPr_ID = NULL; 
  PyObject *pyobjDS_New_I = NULL; 
  PyObject *pyobjDistToCrSurf_I = NULL; 
  PyObject *pyobjTracedPts_I = NULL;
  PyObject *pyobjTrajectories_I = NULL;
  PyObject *pyobjLastStep_I = NULL;

  PyObject *pyobjTheta_I = NULL; 
  PyObject *pyobjPhi_I = NULL; 
  PyObject *pyobjOrientation_I = NULL; 
  PyObject *pyobjDensity_I = NULL; 
  PyObject *pyobjBaseStrength_I = NULL;
  PyObject *pyobjStalkStrength_I = NULL;
  PyObject *pyobjScale_I = NULL;

  /* 
   * 25 pyarrays
   */
  PyArrayObject *pyarParam_I = NULL;
  PyArrayObject *pyarPos_ID = NULL;
  PyArrayObject *pyarDir_ID = NULL;
  PyArrayObject *pyarDS_I = NULL;
  PyArrayObject *pyarTbr_I = NULL;
  PyArrayObject *pyarTbrIQUV_IP = NULL;
  PyArrayObject *pyarOpDepth_I = NULL;
  PyArrayObject *pyarFlags_I = NULL;
  PyArrayObject *pyarRho_I = NULL;
  PyArrayObject *pyarGradRho_ID = NULL;
  PyArrayObject *pyarBfield_ID = NULL;
  PyArrayObject *pyarPosPr_ID = NULL;
  PyArrayObject *pyarDirPr_ID = NULL;
  PyArrayObject *pyarDS_New_I = NULL;
  PyArrayObject *pyarDistToCrSurf_I = NULL;
  PyArrayObject *pyarTracedPts_I = NULL;
  PyArrayObject *pyarTrajectories_I = NULL;
  PyArrayObject *pyarLastStep_I = NULL;

  PyArrayObject *pyarTheta_I = NULL; 
  PyArrayObject *pyarPhi_I = NULL; 
  PyArrayObject *pyarOrientation_I = NULL; 
  PyArrayObject *pyarDensity_I = NULL; 
  PyArrayObject *pyarBaseStrength_I = NULL;
  PyArrayObject *pyarStalkStrength_I = NULL;
  PyArrayObject *pyarScale_I = NULL;

  /* Constants */
  /* double const pi = 3.1415926535897931; // -- Now pi as a macros */

  /* Elemental input variables */
  double  Tol, RIntgSph;

  /* 
   * Array parameters: 25 pyarrays
   */
  double *Pos_ID = NULL, *Dir_ID = NULL, *DS_I = NULL, *Tbr_I = NULL;
  double *OpDepth_I = NULL, *Param_I = NULL;
  short  *Flags_I = NULL;
  double *Rho_I = NULL, *GradRho_ID = NULL;
  double *PosPr_ID = NULL, *DirPr_ID = NULL, *DS_New_I = NULL;
  double *TbrIQUV_IP = NULL, *Bfield_ID = NULL;
  double *DistToCrSurf_I = NULL;

  int    *TracedPts_I = NULL, *LastStep_I = NULL;
  double *Trajectories_I = NULL;

  double *Theta_I = NULL, *Phi_I = NULL; 
  double *Orientation_I = NULL, *Density_I = NULL;
  double *BaseStrength_I = NULL, *StalkStrength_I = NULL, *Scale_I = NULL;


  int nRay;
  int nRay3;   /* = 3*nRay, size of the coordinate-triplet arrays */
  int nStreamers;
  int i, iRay, iIter, scattering;
  int nactive;    /* Number of active rays (max. nRay) */

  int nTracedPts; /*Number of points and max_trajectory size */

  double percent; /* Percent of active rays = (nactive/nRay)*100% */
  double dblnRay;
  int icnt = 0;
  int icall;
  double rsph2; 
  double minsd; /* minimal solar distance */
  /* Ray tracing mode: 1:just trajectories, 2:Tbr, 3:TbrIQUV */
  int rtmode = 0;
  int nIter;
  struct param *prm; /* Access parameters by structure field names */


  /* Parse the (void *) array into the standard arguments 
   * required for the raytracing alogrithm. */


  if (!PyArg_ParseTuple(args, "OOOOdiiiOOOOOOOOOOOOOOOOOOOOO:trace_beam",
			&pyobjParam_I,
			&pyobjPos_ID,
			&pyobjDir_ID,
			&pyobjDS_I,
			&RIntgSph,
			&nIter,
			&rtmode,
			&scattering,
			&pyobjFlags_I,
			&pyobjTbr_I,
			&pyobjTbrIQUV_IP,
			&pyobjOpDepth_I,
			&pyobjRho_I,
			&pyobjGradRho_ID,
			&pyobjBfield_ID,
			&pyobjPosPr_ID,
			&pyobjDirPr_ID,
			&pyobjDS_New_I,
			&pyobjDistToCrSurf_I,
			&pyobjTracedPts_I,
			&pyobjTrajectories_I,
			&pyobjLastStep_I,
			&pyobjTheta_I,
			&pyobjPhi_I,
			&pyobjOrientation_I,
			&pyobjDensity_I,
			&pyobjBaseStrength_I,
			&pyobjStalkStrength_I,
			&pyobjScale_I))
    return NULL;


  /* Convert Python objects to Numpy array objects */

  OBJ_TO_ARRAY(pyobjParam_I, NPY_DOUBLE,        pyarParam_I);
  OBJ_TO_ARRAY(pyobjPos_ID,  NPY_DOUBLE,        pyarPos_ID);
  OBJ_TO_ARRAY(pyobjDir_ID,  NPY_DOUBLE,        pyarDir_ID);
  OBJ_TO_ARRAY(pyobjDS_I,    NPY_DOUBLE,        pyarDS_I);

  if (rtmode == 2)    /* Brightness temperature */
    OBJ_TO_ARRAY(pyobjTbr_I,      NPY_DOUBLE,   pyarTbr_I);
  if (rtmode == 3) {  
    /* Brightness temperature in Stokes param. I,Q,U,V */
    OBJ_TO_ARRAY(pyobjTbrIQUV_IP, NPY_DOUBLE,   pyarTbrIQUV_IP);
    /* Magnetic field for Tb in Stokes param. I,Q,U,V */
    OBJ_TO_ARRAY(pyobjBfield_ID,  NPY_DOUBLE,   pyarBfield_ID);  }
  if (rtmode == 2 || rtmode == 3)
    OBJ_TO_ARRAY(pyobjOpDepth_I,  NPY_DOUBLE,   pyarOpDepth_I);

  OBJ_TO_ARRAY(pyobjRho_I,          NPY_DOUBLE,     pyarRho_I);
  OBJ_TO_ARRAY(pyobjGradRho_ID,     NPY_DOUBLE,     pyarGradRho_ID);
  OBJ_TO_ARRAY(pyobjPosPr_ID,       NPY_DOUBLE,     pyarPosPr_ID);
  OBJ_TO_ARRAY(pyobjDirPr_ID,       NPY_DOUBLE,     pyarDirPr_ID);
  OBJ_TO_ARRAY(pyobjDS_New_I,       NPY_DOUBLE,     pyarDS_New_I);
  OBJ_TO_ARRAY(pyobjDistToCrSurf_I, NPY_DOUBLE,     pyarDistToCrSurf_I);

  OBJ_TO_ARRAY(pyobjTrajectories_I, NPY_DOUBLE,     pyarTrajectories_I);
  OBJ_TO_ARRAY(pyobjLastStep_I,     NPY_INT,        pyarLastStep_I);
  OBJ_TO_ARRAY(pyobjTracedPts_I,    NPY_INT,        pyarTracedPts_I);

  nTracedPts = PyArray_SIZE(pyarTracedPts_I);

  OBJ_TO_ARRAY(pyobjTheta_I,         NPY_DOUBLE, pyarTheta_I);
  OBJ_TO_ARRAY(pyobjPhi_I,           NPY_DOUBLE, pyarPhi_I);
  OBJ_TO_ARRAY(pyobjOrientation_I,   NPY_DOUBLE, pyarOrientation_I);
  OBJ_TO_ARRAY(pyobjDensity_I,       NPY_DOUBLE, pyarDensity_I);
  OBJ_TO_ARRAY(pyobjBaseStrength_I,  NPY_DOUBLE, pyarBaseStrength_I);
  OBJ_TO_ARRAY(pyobjStalkStrength_I, NPY_DOUBLE, pyarStalkStrength_I);
  OBJ_TO_ARRAY(pyobjScale_I,         NPY_DOUBLE, pyarScale_I);

  /* The only non-double and non-int32 array: individual treatment */
  pyarFlags_I  = (PyArrayObject *) 
    PyArray_FROM_OTF(pyobjFlags_I, NPY_SHORT, NPY_INOUT_ARRAY);
  if (pyarFlags_I == NULL) return NULL;

  /* extract sizes, check if correct, and get pointers to arrays */
  nRay3 = PyArray_SIZE(pyarPos_ID); /* Size of Pos_ID is used as reference */ 
  nRay = nRay3/3; /* All arrays with _D are composed of coord. triplets */ 


  /* Check array sizes and get pointers to array data from PyArrayObject-s */ 

  CHKSIZE_AND_GETPTR(pyarParam_I,        NPARAM,  double,  Param_I);
  CHKSIZE_AND_GETPTR(pyarPos_ID,         nRay3,   double,  Pos_ID);
  CHKSIZE_AND_GETPTR(pyarDir_ID,         nRay3,   double,  Dir_ID);
  CHKSIZE_AND_GETPTR(pyarDS_I,           nRay,    double,  DS_I);
  if (rtmode == 2)   { /* Brightness temperature */
    CHKSIZE_AND_GETPTR(pyarTbr_I,        nRay,    double,  Tbr_I);
  }
  if (rtmode == 3) { 
    /* Brightness temperature in Stokes param. I,Q,U,V */
    CHKSIZE_AND_GETPTR(pyarTbrIQUV_IP,   4*nRay,  double,  TbrIQUV_IP);
    /* For TbrIQUV, polarization due to magnetic field Bfield */	
    CHKSIZE_AND_GETPTR(pyarBfield_ID,    nRay3,   double,  Bfield_ID);  
  }
  if (rtmode == 2 || rtmode == 3) 
    CHKSIZE_AND_GETPTR(pyarOpDepth_I,    nRay,    double,  OpDepth_I);

  CHKSIZE_AND_GETPTR(pyarFlags_I,        nRay,    short,   Flags_I);
  CHKSIZE_AND_GETPTR(pyarRho_I,          nRay,    double,  Rho_I);
  CHKSIZE_AND_GETPTR(pyarGradRho_ID,     nRay3,   double,  GradRho_ID);
  CHKSIZE_AND_GETPTR(pyarPosPr_ID,       nRay3,   double,  PosPr_ID);
  CHKSIZE_AND_GETPTR(pyarDirPr_ID,       nRay3,   double,  DirPr_ID);
  CHKSIZE_AND_GETPTR(pyarDS_New_I,       nRay,    double,  DS_New_I);
  CHKSIZE_AND_GETPTR(pyarDistToCrSurf_I, nRay,    double,  DistToCrSurf_I);
  CHKSIZE_AND_GETPTR(pyarTrajectories_I, nTracedPts*nIter*3, double,   
		                                        Trajectories_I);
  CHKSIZE_AND_GETPTR(pyarTracedPts_I, nTracedPts, int,  TracedPts_I);
  CHKSIZE_AND_GETPTR(pyarLastStep_I,  nTracedPts, int,  LastStep_I);

  nStreamers = PyArray_SIZE(pyarTheta_I);

  CHKSIZE_AND_GETPTR(pyarTheta_I,         nStreamers, double, Theta_I);
  CHKSIZE_AND_GETPTR(pyarPhi_I,           nStreamers, double, Phi_I);
  CHKSIZE_AND_GETPTR(pyarOrientation_I,   nStreamers, double, Orientation_I);
  CHKSIZE_AND_GETPTR(pyarDensity_I,       nStreamers, double, Density_I);
  CHKSIZE_AND_GETPTR(pyarBaseStrength_I,  nStreamers, double, BaseStrength_I);
  CHKSIZE_AND_GETPTR(pyarStalkStrength_I, nStreamers, double, StalkStrength_I);
  CHKSIZE_AND_GETPTR(pyarScale_I,         nStreamers, double, Scale_I);

  prm = (struct param *) Param_I; /* Access parameters by struct field names */


  /*
   * The main loop, where the advance_beam() function is called repeatedly.
   * Every call advances the rays by some increments, generally different
   * for the rays in different plasma density.
   * This loop runs until either all the rays are moved outside the
   * integration sphere, or the maximum number of steps, nIter, is reached.
   * 
   * 
   */

  /* Slightly more than int. sphere radius^2 */
  rsph2 = pow(RIntgSph,2) + 0.01;
  dblnRay = (double) nRay;
  /* nactive = nRay; */
  icnt = 0; /* Step counter */



  make_simulation(nIter,
		  nRay,
		  nTracedPts,
		  prm,
		  Pos_ID, 
		  Dir_ID, 
		  DS_I, 
		  Flags_I,
		  Tbr_I,
		  TbrIQUV_IP,
		  OpDepth_I,
		  Rho_I,
		  GradRho_ID,
		  Bfield_ID,
		  PosPr_ID,
		  DirPr_ID,
		  DS_New_I,
		  DistToCrSurf_I,
		  TracedPts_I,
		  Trajectories_I,
		  LastStep_I,
		  rtmode,
		  scattering,
		  RIntgSph,
		  Theta_I,
		  Phi_I,
		  Orientation_I,
		  Density_I,
		  BaseStrength_I,
		  StalkStrength_I,
		  Scale_I);


  Py_INCREF(Py_None);
  return Py_None;

} /* PyObject *trace_beam(PyObject *self, PyObject *args) */



/************************************************************************* 
 * For a given observer SGI position (obs) and the target (tau,xi)       *
 * coordinates in the image plane, calculates the unity-length direction *
 * vectors dir (dir is in SGI), pointing from obs to targ.               *
 *************************************************************************/
/*
 * Call from Python: 
 * dir = raydir(obs, targ)
 */

PyObject *raydir(PyObject *self, PyObject *args) {

  /* Python objects; arrays to be converted to ndarrays */
  PyObject *pyobjObs_D = NULL; 
  PyObject *pyobjXruler_I = NULL; 
  PyObject *pyobjYruler_I = NULL; 
  PyObject *pyobjDir_ID = NULL; 

  PyArrayObject *pyarObs_D = NULL; 
  PyArrayObject *pyarXruler_I = NULL; 
  PyArrayObject *pyarYruler_I = NULL; 
  PyArrayObject *pyarDir_ID = NULL; 

  int nx, ny, nRay, nRay2, nRay3;

  /* Array data pointers */
  double *Obs_D = NULL, *Xruler_I = NULL, *Yruler_I = NULL, *Dir_ID = NULL;

  if (!PyArg_ParseTuple(args, "OOOO:raydir",
			&pyobjObs_D,
			&pyobjXruler_I,
			&pyobjYruler_I,
			&pyobjDir_ID))
    return NULL;

  /* Convert Python objects to Numpy array objects */
  OBJ_TO_ARRAY(pyobjObs_D, NPY_DOUBLE,    pyarObs_D);
  OBJ_TO_ARRAY(pyobjXruler_I, NPY_DOUBLE, pyarXruler_I);
  OBJ_TO_ARRAY(pyobjYruler_I, NPY_DOUBLE, pyarYruler_I);
  OBJ_TO_ARRAY(pyobjDir_ID, NPY_DOUBLE,   pyarDir_ID);

  /* Sizes of Xruler_I and Yruler_I, nx and ny, are used as reference */
  nx = PyArray_DIM(pyarXruler_I,0); 
  ny = PyArray_DIM(pyarYruler_I,0);
  nRay = nx*ny;


  /* Check array sizes and get pointers to array data from PyArrayObject-s */ 
  CHKSIZE_AND_GETPTR(pyarObs_D, 3, double, Obs_D);
  CHKSIZE_AND_GETPTR(pyarDir_ID, 3*nRay, double, Dir_ID);
  Xruler_I = (double *) PyArray_DATA(pyarXruler_I);
  Yruler_I = (double *) PyArray_DATA(pyarYruler_I);

  /* Calculate the direction vectors and save them in Dir_ID */  
  mcalc_rdir(Obs_D, Xruler_I, Yruler_I, nx, ny, Dir_ID);

  Py_INCREF(Py_None);
  return Py_None;

} /* PyObject *raydir() */

/* 
 * Get ray positions on the surface of integration sphere
 * Call from Python:
 *  nisec = raysph(Obs_D, Dir_ID, rsph, Isec_I, Pos_ID)
 * Input:
 *   Obs_D[3] - observer position vector (the earth position in SGI) 
 *   Dir_ID[nRay][3] - direction vectors of the rays
 *   rsph - radius of the integration sphere
 * Output:
 *   Isec_I[nRay] - for each ray, number intersections with the sphere: 0,1,2
 *   Pos_ID[nRay][3] - ray positions on the surface of integration sphere
 *
 * The problem is reduced to solving the quadratic equation:
 *
 *  dir^2*dist^2 + 2(dir.obs)*dist + (obs^2 - rsph^2) = 0
 *  
 */
/*
 * Call from Python: 
 * raysph(obs, targ)
 */
PyObject *raysph(PyObject *self, PyObject *args) {

  /* Python objects; arrays to be converted to ndarrays */
  PyObject *pyobjObs_D = NULL; 
  PyObject *pyobjDir_ID = NULL; 
  PyObject *pyobjIsec_I = NULL; 
  PyObject *pyobjPos_ID = NULL; 

  PyArrayObject *pyarObs_D = NULL; 
  PyArrayObject *pyarDir_ID = NULL; 
  PyArrayObject *pyarIsec_I = NULL; 
  PyArrayObject *pyarPos_ID = NULL; 

  /* Elemental input variables */
  double rsph; /* Radius (in solar radii) of the integration sphere */
  int nRay; /* Total number of rays */
  int nIsec; /* Total number of ray intersections with the sphere */
  int ndim, nx, ny;

  /* Array data pointers */
  double *Obs_D = NULL, *Dir_ID = NULL, *Pos_ID = NULL;
  short  *Isec_I  = NULL;


  if (!PyArg_ParseTuple(args, "OOdOO:raysph",
			&pyobjObs_D,
			&pyobjDir_ID,
			&rsph,
			&pyobjIsec_I,
			&pyobjPos_ID))
    return NULL;
  /* Convert Python objects to Numpy array objects */
  OBJ_TO_ARRAY(pyobjObs_D, NPY_DOUBLE,  pyarObs_D);
  OBJ_TO_ARRAY(pyobjDir_ID, NPY_DOUBLE, pyarDir_ID);
  OBJ_TO_ARRAY(pyobjIsec_I, NPY_DOUBLE, pyarIsec_I);
  OBJ_TO_ARRAY(pyobjPos_ID, NPY_DOUBLE, pyarPos_ID);

  /* extract sizes, check if correct, and get pointers to arrays */
  /* Size(s) of pyarDir_ID are used as reference */
  ndim =  PyArray_NDIM(pyarDir_ID);
  if      (ndim == 1)   nRay = PyArray_DIM(pyarDir_ID,0)/3;
  else if (ndim == 2)   nRay = PyArray_DIM(pyarDir_ID,0); 
  else if (ndim == 3) {
    ny = PyArray_DIM(pyarDir_ID,0); 
    nx = PyArray_DIM(pyarDir_ID,1);
    nRay = nx*ny;
  } 

  /* Check array sizes and get pointers to array data from PyArrayObject-s */ 
  CHKSIZE_AND_GETPTR(pyarObs_D,       3, double, Obs_D);
  /* CHKSIZE_AND_GETPTR(pyarDir_ID, 3*nRay, double, Dir_ID); // no need! */
  CHKSIZE_AND_GETPTR(pyarIsec_I,   nRay, short,  Isec_I);
  CHKSIZE_AND_GETPTR(pyarPos_ID, 3*nRay, double, Pos_ID);

  /* Calculate the position vectors and save them in Pos_ID */  
  nIsec = mcalc_rsph(Obs_D, Dir_ID, rsph, nRay, Isec_I, Pos_ID);

  return Py_BuildValue("i", nIsec);

} /* PyObject *raysph() */


/* 
 * Get ray directions toward the image plane pixels and their positions 
 * on the surface of integration sphere
 * Call from Python:
 *  nisec = raypd(Obs_D, Xruler, Yruler, rsph, Isec_I, Dir_ID, Pos_ID)
 * Input:
 *   Obs_D[3] - observer position vector (the earth position in SGI)
 *   Xruler_I[nx] - array of pixel positions in X direction
 *   Yruler_I[ny] - array of pixel positions in Y direction
 *   rsph - radius of the integration sphere
 * Output:
 *   Isec_I[nRay] - for each ray, number intersections with the sphere: 0,1,2
 *   Dir_ID[nRay][3] - direction vectors of the rays
 *   Pos_ID[nRay][3] - ray positions on the surface of integration sphere
 *
 * After the directions are found, the problem is reduced to solving the 
 * quadratic equation:
 *
 *  dir^2*dist^2 + 2(dir.obs)*dist + (obs^2 - rsph^2) = 0
 *  
 */
/*
 * Call from Python: 
 * raypd(obs, targ)
 */
PyObject *raypd(PyObject *self, PyObject *args) {

  /* Python objects; arrays to be converted to ndarrays */
  PyObject *pyobjObs_D = NULL; 
  PyObject *pyobjXruler_I = NULL; 
  PyObject *pyobjYruler_I = NULL; 
  PyObject *pyobjIsec_I = NULL; 
  PyObject *pyobjDir_ID = NULL; 
  PyObject *pyobjPos_ID = NULL; 

  PyArrayObject *pyarObs_D = NULL; 
  PyArrayObject *pyarXruler_I = NULL; 
  PyArrayObject *pyarYruler_I = NULL; 
  PyArrayObject *pyarIsec_I = NULL; 
  PyArrayObject *pyarDir_ID = NULL; 
  PyArrayObject *pyarPos_ID = NULL; 

  /* Elemental input variables */
  double rsph; /* Radius (in solar radii) of the integration sphere */
  int nRay; /* Total number of rays */
  int nIsec = 0; /* Total number of ray intersections with the sphere */
  int nx, ny;

  /* Array data pointers */
  double *Obs_D = NULL, *Dir_ID = NULL, *Pos_ID = NULL;
  double *Xruler_I = NULL, *Yruler_I = NULL;
  short  *Isec_I  = NULL;


  if (!PyArg_ParseTuple(args, "OOOdOOO:raypd",
			&pyobjObs_D,
			&pyobjXruler_I,
			&pyobjYruler_I,
			&rsph,
			&pyobjIsec_I,
			&pyobjDir_ID,
			&pyobjPos_ID))
    return NULL;

  /* Convert Python objects to Numpy array objects */
  OBJ_TO_ARRAY(pyobjObs_D, NPY_DOUBLE,    pyarObs_D);
  OBJ_TO_ARRAY(pyobjXruler_I, NPY_DOUBLE, pyarXruler_I);
  OBJ_TO_ARRAY(pyobjYruler_I, NPY_DOUBLE, pyarYruler_I);
  OBJ_TO_ARRAY(pyobjIsec_I, NPY_DOUBLE,   pyarIsec_I);
  OBJ_TO_ARRAY(pyobjDir_ID, NPY_DOUBLE,   pyarDir_ID);
  OBJ_TO_ARRAY(pyobjPos_ID, NPY_DOUBLE,   pyarPos_ID);

  /* extract sizes, check if correct, and get pointers to arrays */
  /* Sizes of Xruler_I and Yruler_I, nx and ny, are used as reference */
  nx = PyArray_DIM(pyarXruler_I,0); 
  ny = PyArray_DIM(pyarYruler_I,0);
  nRay = nx*ny;

  /* Check array sizes and get pointers to array data from PyArrayObject-s */ 
  Xruler_I = (double *) PyArray_DATA(pyarXruler_I);
  Yruler_I = (double *) PyArray_DATA(pyarYruler_I);
  CHKSIZE_AND_GETPTR(pyarObs_D,       3, double, Obs_D);
  CHKSIZE_AND_GETPTR(pyarIsec_I,   nRay, short,  Isec_I);
  CHKSIZE_AND_GETPTR(pyarDir_ID, 3*nRay, double, Dir_ID);
  CHKSIZE_AND_GETPTR(pyarPos_ID, 3*nRay, double, Pos_ID);

  /* Calculate the direction vectors and save them in Dir_ID */  
  mcalc_rdir(Obs_D, Xruler_I, Yruler_I, nx, ny, Dir_ID);

  /* Calculate the position vectors and save them in Pos_ID */  
  nIsec = mcalc_rsph(Obs_D, Dir_ID, rsph, nRay, Isec_I, Pos_ID);

  return Py_BuildValue("i", nIsec);

} /* PyObject *raypd() */




/* ==== methods table ====================== */
static PyMethodDef rtcore_methods[] = {
  {"trace_beam", trace_beam, METH_VARARGS, "Make several steps of ray tracing"},
  {"raydir", raydir,  METH_VARARGS, "Calculate initial ray directions"},
  {"raysph", raysph,  METH_VARARGS, 
   "Calculate initial ray positions on the surface of a sphere"},
  {"raypd", raypd,  METH_VARARGS, 
   "Calculate initial ray positions and directions on the surface of a sphere"},
  {NULL, NULL, 0, NULL}
};


/* ==== Initialize ====================== */
PyMODINIT_FUNC initrtcore()  {
  Py_InitModule("rtcore", rtcore_methods);
  import_array();  // for NumPy
}

