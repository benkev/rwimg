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
//#include								\
//    "/usr/lib/python2.6/site-packages/numpy/core/include/numpy/arrayobject.h"
#include "raytrace.h"
#include <time.h>
#include <pthread.h>

/* 
 * This defintion is the number of threads the program will use 
 * in its calculation.
 * If set to 1 this will be completely serial.
 */
#define NUM_THREADS 10

#define OBJ_TO_INOUT_ARRAY(obj,arr)     \
    if (obj == NULL) return NULL;	        \
arr = (PyArrayObject *) PyArray_FROM_OTF(obj, NPY_DOUBLE, NPY_INOUT_ARRAY); \
if (arr == NULL) return NULL;

#define CHKSIZE_AND_GETPTR(arr,nrq,ty,ptr)	\
{					\
    int i, ndim, nel = 1;		\
    ndim =  PyArray_NDIM(arr);		\
    for (i = 0; i < ndim; i++)		\
    nel = nel*PyArray_DIM(arr,i);     \
    if (nel != nrq) {			\
        PyErr_Format(PyExc_ValueError, "Wrong size of array " #arr ": %d; " \
                "must be %d.", nel, nrq); \
        return NULL;			     \
    }					     \
    ptr =    (ty *) PyArray_DATA(arr);   \
}

/*
 * Declaration of the threadFunction
 */
PyObject *threadFunc(void *args[]);

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
PyObject *trace_beam(PyObject *self, PyObject *args){
{

    /* Python objects; arrays to be converted to ndarrays */
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


    /* Constants */
    /* double const pi = 3.1415926535897931; // -- Now pi as a macros */

    /* Elemental input variables */
    double  Tol, RIntgSph;
    /* Array parameters */
    double *Pos_ID = NULL, *Dir_ID = NULL, *DS_I = NULL, *Tbr_I = NULL;
    double *OpDepth_I = NULL, *Param_I = NULL;
    short  *Flags_I = NULL;
    double *Rho_I = NULL, *GradRho_ID = NULL;
    double *PosPr_ID = NULL, *DirPr_ID = NULL, *DS_New_I = NULL;
    double *TbrIQUV_IP = NULL, *Bfield_ID = NULL;
    double *DistToCrSurf_I = NULL;
    double *TracedPts_I = NULL, *Trajectories_I = NULL, *LastStep_I = NULL;
    int nRay;
    int nRay3;
    int i, iRay, iIter, scattering;
    int nactive;    /* Number of active rays (max. nRay) */

    int nTracedPts; /*Number of points and max_trajectory size */

    double percent; /* Percent of active rays = (nactive/nRay)*100% */
    double dblnRay;
    int icnt = 0;
    int icall;
    double rsph2; 
    double minsd; /* minimal solar distance */
    int rtmode = 0; /* Ray tracing mode: 1:just trajectories, 2:Tbr, 3:TbrIQUV */

    struct param *prm; /* Access parameters by structure field names */


    /* Parse the (void *) array into the standard arguments 
     * required for the raytracing alogrithm. */


    if (!PyArg_ParseTuple(args, "OOOOsdiiiOOOOOOOOOOOOOO:trace_beam",
                &pyobjParam_I,
                &pyobjPos_ID,
                &pyobjDir_ID,
                &pyobjDS_I,
                &plfname,
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
                &pyobjLastStep_I))
                    return NULL;





    /* Convert Python objects to Numpy array objects */

    OBJ_TO_INOUT_ARRAY(pyobjParam_I,        pyarParam_I);
    OBJ_TO_INOUT_ARRAY(pyobjPos_ID,         pyarPos_ID);
    OBJ_TO_INOUT_ARRAY(pyobjDir_ID,         pyarDir_ID);
    OBJ_TO_INOUT_ARRAY(pyobjDS_I,           pyarDS_I);
    if (rtmode == 2)    /* Brightness temperature */
        OBJ_TO_INOUT_ARRAY(pyobjTbr_I,          pyarTbr_I);
    if (rtmode == 3) {  /* Brightness temperature in Stokes param. I,Q,U,V */
        OBJ_TO_INOUT_ARRAY(pyobjTbrIQUV_IP,          pyarTbrIQUV_IP);
    }
    if (rtmode == 2 || rtmode == 3)
        OBJ_TO_INOUT_ARRAY(pyobjOpDepth_I,      pyarOpDepth_I);
    OBJ_TO_INOUT_ARRAY(pyobjRho_I,          pyarRho_I);
    OBJ_TO_INOUT_ARRAY(pyobjGradRho_ID,     pyarGradRho_ID);
    if (rtmode == 3)    /* Magnetic field for Tb in Stokes param. I,Q,U,V */
        OBJ_TO_INOUT_ARRAY(pyobjBfield_ID,      pyarBfield_ID);
    OBJ_TO_INOUT_ARRAY(pyobjPosPr_ID,       pyarPosPr_ID);
    OBJ_TO_INOUT_ARRAY(pyobjDirPr_ID,       pyarDirPr_ID);
    OBJ_TO_INOUT_ARRAY(pyobjDS_New_I,       pyarDS_New_I);
    OBJ_TO_INOUT_ARRAY(pyobjDistToCrSurf_I, pyarDistToCrSurf_I);
    OBJ_TO_INOUT_ARRAY(pyobjTracedPts_I, pyarTracedPts_I);
    OBJ_TO_INOUT_ARRAY(pyobjTrajectories_I, pyarTrajectories_I);
    OBJ_TO_INOUT_ARRAY(pyobjLastStep_I, pyarLastStep_I);

    /* The only non-double array: individual treatment */
    pyarFlags_I  = (PyArrayObject *) 
        PyArray_FROM_OTF(pyobjFlags_I, NPY_SHORT, NPY_INOUT_ARRAY);
    if (pyarFlags_I == NULL) return NULL;

    /* extract sizes, check if correct, and get pointers to arrays */
    nRay3 = PyArray_SIZE(pyarPos_ID); /* Size of Pos_ID is used as reference */ 
    nRay = nRay3/3; /* All arrays with _D are composed of coord. triplets */ 

    /* Check array sizes and get pointers to array data from PyArrayObject-s */ 

    if(thread_num == 0)
        printf("NPARAM = %d\n", NPARAM);

    CHKSIZE_AND_GETPTR(pyarParam_I,    NPARAM, double, Param_I);
    CHKSIZE_AND_GETPTR(pyarPos_ID,     nRay3, double,  Pos_ID);
    CHKSIZE_AND_GETPTR(pyarDir_ID,     nRay3, double,  Dir_ID);
    CHKSIZE_AND_GETPTR(pyarDS_I,       nRay, double,   DS_I);
    if (rtmode == 2)    /* Brightness temperature */
        CHKSIZE_AND_GETPTR(pyarTbr_I,      nRay, double,   Tbr_I);
    if (rtmode == 3) {  /* Brightness temperature in Stokes param. I,Q,U,V */
        CHKSIZE_AND_GETPTR(pyarTbrIQUV_IP,  4*nRay, double,   TbrIQUV_IP);
    }
    if (rtmode == 2 || rtmode == 3)
        CHKSIZE_AND_GETPTR(pyarOpDepth_I,  nRay, double,   OpDepth_I);
    CHKSIZE_AND_GETPTR(pyarFlags_I,    nRay, short,    Flags_I);
    CHKSIZE_AND_GETPTR(pyarRho_I,      nRay, double,   Rho_I);
    CHKSIZE_AND_GETPTR(pyarGradRho_ID, nRay3, double,  GradRho_ID);
    if (rtmode == 3)  /* For TbrIQUV, polarization due to magnetic field Bfield */
        CHKSIZE_AND_GETPTR(pyarBfield_ID, nRay3, double,  Bfield_ID);
    CHKSIZE_AND_GETPTR(pyarPosPr_ID,   nRay3, double,  PosPr_ID);
    CHKSIZE_AND_GETPTR(pyarDirPr_ID,   nRay3, double,  DirPr_ID);
    CHKSIZE_AND_GETPTR(pyarDS_New_I,   nRay, double,   DS_New_I);
    CHKSIZE_AND_GETPTR(pyarDistToCrSurf_I, nRay, double,   DistToCrSurf_I);

    nTracedPts = PyArray_SIZE(pyarTracedPts_I);

    CHKSIZE_AND_GETPTR(pyarTracedPts_I, nTracedPts, double,  TracedPts_I);
    CHKSIZE_AND_GETPTR(pyarTrajectories_I, nTracedPts*3*nIter, double,   Trajectories_I);
    CHKSIZE_AND_GETPTR(pyarLastStep_I, nTracedPts, double,  LastStep_I);


    prm = (void *) Param_I; /* Access parameters by structure field names */

    /* printf("Param_I[0..%d] = \n", NPARAM); */
    /* for (i = 0; i < NPARAM; i++) */
    /*   printf("Param_I[%d] = %20.12e\n", i, Param_I[i]); */

    /* printf("ProtonMassInv = %20.12e\n", prm->ProtonMassInv); */
    /* printf("Boltzmannk_ergK = %20.12e\n", prm->Boltzmannk_ergK); */
    /* printf("Freq = %20.12e\n", prm->Freq); */
    /* printf("Freq2 = %20.12e\n", prm->Freq2); */
    /* printf("RhoCrInv = %20.12e\n", prm->RhoCrInv); */
    /* printf("TolEps = %20.12e\n", prm->TolEps); */
    /* printf("CntMax = %20.12e\n", prm->CntMax); */
    /* printf("Tol2 = %20.12e\n", prm->Tol2); */
    /* printf("AbsMinStep = %20.12e\n", prm->AbsMinStep); */
    /* printf("r_chro_cor = %20.12e\n", prm->r_chro_cor); */
    /* printf("Te_chromo_K = %20.12e\n", prm->Te_chromo_K); */
    /* printf("Rsun_cm = %20.12e\n", prm->Rsun_cm); */
    /* printf("lnLambda_13p7 = %20.12e\n", prm->lnLambda_13p7); */
    /* printf(" = %20.12e\n", prm->); */
    /* printf(" = %20.12e\n", prm->); */
    /* printf(" = %20.12e\n", prm->); */
    /* printf(" = %20.12e\n", prm->); */

    /*
     * The main loop, where the advance_beam() function is called repeatedly.
     * Every call advances the rays by some increments, generally different
     * for the rays in different plasma density.
     * This loop runs until either all the rays are moved outside the
     * integration sphere, or the maximum number of steps, nIter, is reached.
     * 
     * 
     */

    rsph2 = pow(RIntgSph,2) + 0.01; /* Slightly more than int. sphere radius^2 */
    dblnRay = (double) nRay;
    /* nactive = nRay; */
    icnt = 0; /* Step counter */



    make_simulation(nIter,
                nRay,
                Param_I,
                Pos_ID, 
                Dir_ID, 
                DS_I, 
                plasma_parameters,
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
                rtmode,
                scattering);


    Py_INCREF(Py_None);
    return Py_None;

}



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
    OBJ_TO_INOUT_ARRAY(pyobjObs_D,    pyarObs_D);
    OBJ_TO_INOUT_ARRAY(pyobjXruler_I, pyarXruler_I);
    OBJ_TO_INOUT_ARRAY(pyobjYruler_I, pyarYruler_I);
    OBJ_TO_INOUT_ARRAY(pyobjDir_ID,   pyarDir_ID);

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
    OBJ_TO_INOUT_ARRAY(pyobjObs_D,  pyarObs_D);
    OBJ_TO_INOUT_ARRAY(pyobjDir_ID, pyarDir_ID);
    OBJ_TO_INOUT_ARRAY(pyobjIsec_I, pyarIsec_I);
    OBJ_TO_INOUT_ARRAY(pyobjPos_ID, pyarPos_ID);

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
    OBJ_TO_INOUT_ARRAY(pyobjObs_D,    pyarObs_D);
    OBJ_TO_INOUT_ARRAY(pyobjXruler_I, pyarXruler_I);
    OBJ_TO_INOUT_ARRAY(pyobjYruler_I, pyarYruler_I);
    OBJ_TO_INOUT_ARRAY(pyobjIsec_I,   pyarIsec_I);
    OBJ_TO_INOUT_ARRAY(pyobjDir_ID,   pyarDir_ID);
    OBJ_TO_INOUT_ARRAY(pyobjPos_ID,   pyarPos_ID);

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


/*
 * The python frontend passes an array of position values each entry has
 * 3 coordinate. It also passes in empty arrays to hold the return values for
 * rho, gradrho, and bfield. Those values are calculated at the various points
 * and returned.
 */
PyObject *plprofile(PyObject *self, PyObject *args) {

    int rtmode = 0; /* Ray tracing mode: 1:just trajectories, 2:Tbr, 3:TbrIQUV */
    int nPoints, nPoints3;

    const char *plfname; /* Name of plasma parameters function */

    /* Dynamic linking data */
    void *dlh;            /* DL library handle */
    char plfsoname[300]; /* [Path]name of plasma parameters DL library */
    void (* plasma_parameters)();

    double *Pos_ID = NULL, *DS_I = NULL;
    double *Param_I = NULL;
    short  *Flags_I = NULL;
    double *Rho_I = NULL, *GradRho_ID = NULL;
    double *Bfield_ID = NULL;


    /* Python objects; arrays to be converted to ndarrays */
    PyObject *pyobjParam_I = NULL; 
    PyObject *pyobjPos_ID = NULL; 
    PyObject *pyobjRho_I = NULL; 
    PyObject *pyobjGradRho_ID = NULL; 
    PyObject *pyobjBfield_ID = NULL; 
    PyObject *pyobjDS_I = NULL; 
    PyObject *pyobjFlags_I = NULL; 

    PyArrayObject *pyarParam_I = NULL;
    PyArrayObject *pyarPos_ID = NULL;
    PyArrayObject *pyarRho_I = NULL;
    PyArrayObject *pyarGradRho_ID = NULL;
    PyArrayObject *pyarBfield_ID = NULL;
    PyArrayObject *pyarDS_I = NULL;
    PyArrayObject *pyarFlags_I = NULL;



    if (!PyArg_ParseTuple(args, "OOOOOisOO:plprofile",
                &pyobjParam_I,
                &pyobjPos_ID,
                &pyobjRho_I,
                &pyobjGradRho_ID,
                &pyobjBfield_ID,
                &rtmode,
                &plfname,
                &pyobjDS_I,
                &pyobjFlags_I))
        return NULL;

    /* Convert Python objects to Numpy array objects */
    OBJ_TO_INOUT_ARRAY(pyobjParam_I,        pyarParam_I);
    OBJ_TO_INOUT_ARRAY(pyobjPos_ID,         pyarPos_ID);
    OBJ_TO_INOUT_ARRAY(pyobjRho_I,          pyarRho_I);
    OBJ_TO_INOUT_ARRAY(pyobjGradRho_ID,     pyarGradRho_ID);
    if (rtmode == 3)    /* Magnetic field for Tb in Stokes param. I and V */
        OBJ_TO_INOUT_ARRAY(pyobjBfield_ID,      pyarBfield_ID);
    OBJ_TO_INOUT_ARRAY(pyobjDS_I,           pyarDS_I);



    /* The only non-double array: individual treatment */
    pyarFlags_I  = (PyArrayObject *) 
        PyArray_FROM_OTF(pyobjFlags_I, NPY_SHORT, NPY_INOUT_ARRAY);
    if (pyarFlags_I == NULL) return NULL;

    /* extract sizes, check if correct, and get pointers to arrays */
    nPoints = PyArray_SIZE(pyarPos_ID); /* Size of Pos_ID is used as reference */ 
    nPoints3 = nPoints/3;
    /* Check array sizes and get pointers to array data from PyArrayObject-s */ 

    printf("NPARAM = %d\n", NPARAM);

    CHKSIZE_AND_GETPTR(pyarParam_I,    NPARAM, double, Param_I);
    CHKSIZE_AND_GETPTR(pyarPos_ID,     nPoints, double,  Pos_ID);
    CHKSIZE_AND_GETPTR(pyarRho_I,      nPoints3, double,   Rho_I);
    CHKSIZE_AND_GETPTR(pyarGradRho_ID, nPoints, double,  GradRho_ID);
    if (rtmode == 3)  /* For TbrIQUV, polarization due to magnetic field Bfield */
        CHKSIZE_AND_GETPTR(pyarBfield_ID, nPoints, double,  Bfield_ID);
    CHKSIZE_AND_GETPTR(pyarDS_I,       nPoints3, double,   DS_I);
    CHKSIZE_AND_GETPTR(pyarFlags_I,    nPoints3, short,    Flags_I);



    /*
     * Dynamically link the plasma density & magnetic field function
     */
    /* Make "library.so" name as plfsoname = plfname + ".so" */
    //plfsoname = (char *) malloc((strlen(plfname) + 4)*sizeof(char));
    getcwd(plfsoname, 255); /* The DL name must include all the path; here cwd */
    strcat(plfsoname, "/");
    strcat(plfsoname, plfname);
    strcat(plfsoname, ".so");

    //printf("Library name: %s\n\n", plfsoname);
    /* Get the library handle */
    dlh = dlopen(plfsoname, RTLD_LAZY);
    if (!dlh) {
        fputs (dlerror(), stderr);
        return NULL;
    }

    plasma_parameters = dlsym(dlh, plfname);

    plasma_parameters(0,
            nPoints3, 
            nPoints3, 
            Param_I, 
            Pos_ID, 
            Rho_I, 
            GradRho_ID,
            Bfield_ID,
            DS_I, 
            Flags_I,
            rtmode);

    Py_INCREF(Py_None);
    return Py_None;

}

/*
 * The python frontend passes an array of parameters to make a streamer
 */
PyObject *make_streamer(PyObject *self, PyObject *args) {

    double dRadius, dTheta, dPhi, dDensity, dBaseStrength, dStalkStrength;
    double dOrientation,dScaleX,dScaleY,dScaleZ;

    const char *plfname; /* Name of plasma parameters function */

    /* Dynamic linking data */
    void *dlh;            /* DL library handle */
    char plfsoname[300]; /* [Path]name of plasma parameters DL library */
    void (* add_streamer)();


    if (!PyArg_ParseTuple(args, "dddddddddds:makes_streamer",
                &dRadius,
                &dTheta,
                &dPhi,
                &dOrientation,
                &dDensity,
                &dBaseStrength,
                &dStalkStrength,
                &dScaleX,
                &dScaleY,
                &dScaleZ,
                &plfname))
        return NULL;

    /*
     * Dynamically link the plasma density & magnetic field function
     */
    /* Make "library.so" name as plfsoname = plfname + ".so" */
    //plfsoname = (char *) malloc((strlen(plfname) + 4)*sizeof(char));
    getcwd(plfsoname, 255); /* The DL name must include all the path; here cwd */
    strcat(plfsoname, "/");
    strcat(plfsoname, plfname);
    strcat(plfsoname, ".so");

    //printf("Library name: %s\n\n", plfsoname);
    /* Get the library handle */
    dlh = dlopen(plfsoname, RTLD_LAZY);
    if (!dlh) {
        fputs (dlerror(), stderr);
        return NULL;
    }


    add_streamer = dlsym(dlh, "add_streamer");

    add_streamer(dRadius,
            dTheta,
            dPhi,
            dOrientation,
            dDensity,
            dBaseStrength,
            dStalkStrength,
            dScaleX,
            dScaleY,
            dScaleZ);


    Py_INCREF(Py_None);
    return Py_None;

}

/*
 * The python frontend asks to destroy all the streamers
 */
PyObject *remove_streamers(PyObject *self, PyObject *args) {

    const char *plfname; /* Name of plasma parameters function */

    /* Dynamic linking data */
    void *dlh;            /* DL library handle */
    char plfsoname[300]; /* [Path]name of plasma parameters DL library */
    void (* remove_streamers)();


    if (!PyArg_ParseTuple(args, "s:remove_streamers",
                &plfname))
        return NULL;

    /*
     * Dynamically link the plasma density & magnetic field function
     */
    /* Make "library.so" name as plfsoname = plfname + ".so" */
    //plfsoname = (char *) malloc((strlen(plfname) + 4)*sizeof(char));
    getcwd(plfsoname, 255); /* The DL name must include all the path; here cwd */
    strcat(plfsoname, "/");
    strcat(plfsoname, plfname);
    strcat(plfsoname, ".so");

    //printf("Library name: %s\n\n", plfsoname);
    /* Get the library handle */
    dlh = dlopen(plfsoname, RTLD_LAZY);
    if (!dlh) {
        fputs (dlerror(), stderr);
        return NULL;
    }

    remove_streamers = dlsym(dlh, "remove_streamers");

    remove_streamers();


    Py_INCREF(Py_None);
    return Py_None;

}


/* ==== methods table ====================== */
static PyMethodDef rtcore_methods[] = {
    {"trace_beam", trace_beam, METH_VARARGS, "Make several steps of ray tracing"},
    {"raydir", raydir,  METH_VARARGS, "Calculate initial ray directions"},
    {"raysph", raysph,  METH_VARARGS, 
        "Calculate initial ray positions on the surface of a sphere"},
    {"raypd", raypd,  METH_VARARGS, 
        "Calculate initial ray positions and directions on the surface of a sphere"},
    {"plprofile", plprofile,  METH_VARARGS, 
        "Get a profile of the plasma"},
    {"make_streamer", make_streamer,  METH_VARARGS, 
        "Make a streamer"},
    {"remove_streamers", remove_streamers,  METH_VARARGS, 
        "Make a streamer"},
    {NULL, NULL, 0, NULL}
};


/* ==== Initialize ====================== */
PyMODINIT_FUNC initrtcore()  {
    Py_InitModule("rtcore", rtcore_methods);
    import_array();  // for NumPy
}

