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
#include <stdint.h>
#include <math.h>
#include <float.h>
#include <string.h>
#include <Python.h>

//#include	\
//  "/usr/lib/python2.6/site-packages/numpy/core/include/numpy/arrayobject.h"

#include "numpy/arrayobject.h"
#include "raytrace.h"
#include <time.h>
#include <pthread.h>

/* 
 * This defintion is the number of threads the program will use 
 * in its calculation.
 * If set to 1 this will be completely serial.
*/
#define NUM_THREADS 1

// NPY_DOUBLE
#define OBJ_TO_INOUT_ARRAY(obj,npy_ty,arr)	\
  if (obj == NULL) return NULL;	        \
  arr = (PyArrayObject *) PyArray_FROM_OTF(obj, npy_ty, NPY_INOUT_ARRAY); \
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
      
/* Global variables */
char str_pos[] = "pos"; 
char str_dir[] = "dir"; 

/*
 * Declaration of the threadFunction
 */
PyObject *threadFunc(void *args[]);

/*
 * This function performs the raytrace from start to finish, 
 * calling advance_beam at every step.
 * In this threaded version this acts only as a shell function, 
 * calling the threadFunc for every thread.
 */

/*
 * Call from Python:
 * trace_beam(...)
 */
PyObject *trace_beam(PyObject *self, PyObject *args) {

  /* An array to hold the arguments to be passed to the threads */
  void *thread_args[NUM_THREADS][3];
  int i;
  /* An array to hold the thread numbers */
  pthread_t threads[NUM_THREADS];


/* create all threads */
  for(i = 0;i < NUM_THREADS; i++){
    /* Put the necessary arguments in the array and then pass that */
    thread_args[i][0] = self;
    thread_args[i][1] = args;
    thread_args[i][2] = i;
    pthread_create(&threads[i], NULL, threadFunc, (void *) thread_args[i]);
  }

/* wait for all threads to complete */
  for (i=0; i<NUM_THREADS; ++i) {
    pthread_join(threads[i], NULL);
  }

  Py_INCREF(Py_None);
  return Py_None;

} /* PyObject *trace_beam(PyObject *self, PyObject *args) */

/*
 * This function is called for every thread to perform its own beam tracing.
 *
 * Only the root thread (thread 0) prints out messages.
 *
 * parameters: Param_I,Pos_ID,Dir_ID,Arclen_I,Dist_I,DS_I,plfname,RIntgSph,
 *   nIter,rtmode,scattering,Flags_I,
 *   Tbr_I,TbrIQUV_IP,OpDepth_I,Rho_I,GradRho_ID,
 *   Bfield_ID,PosPr_ID,DirPr_ID,DS_New_I,
 *   DistToCrSurf_I,Trace
 */
PyObject *threadFunc(void *thread_args[]){

  PyObject *self = NULL;
  PyObject *args = NULL;

  /* Python objects; arrays to be converted to ndarrays */
  PyObject *pyobjParam_I = NULL; 
  PyObject *pyobjPos_ID = NULL; 
  PyObject *pyobjDir_ID = NULL; 
  PyObject *pyobjArclen_I = NULL; 
  PyObject *pyobjDist_I = NULL; 
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
  PyObject *pyobjTrace = NULL;
  PyObject *pyobjTrcPos_ID = NULL;
  PyObject *pyobjTrcDir_ID = NULL;
  PyObject *pyobjTrcArclen_I = NULL;
  PyObject *pyobjTrcDist_I = NULL;
  PyObject *pyobjTrcRho_I = NULL;
  PyObject *pyobjTrcGradRho_ID = NULL;
  PyObject *pyobjTrcBfield_ID = NULL;
  PyObject *pyobjTrcTbr_I = NULL;
  PyObject *pyobjTrcTbrIQUV_IP = NULL;
  PyObject *pyobjTrcLast_I = NULL;
  PyObject *pyobjTrcIrays_I = NULL;
  PyObject *pyobjTrcNtraj = NULL;
  PyObject *pyobjTrcNpmax = NULL;

  PyArrayObject *pyarParam_I = NULL;
  PyArrayObject *pyarPos_ID = NULL;
  PyArrayObject *pyarDir_ID = NULL;
  PyArrayObject *pyarArclen_I = NULL;
  PyArrayObject *pyarDist_I = NULL;
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
  PyArrayObject *pyarTrcPos_ID = NULL;
  PyArrayObject *pyarTrcDir_ID = NULL;
  PyArrayObject *pyarTrcArclen_I = NULL;
  PyArrayObject *pyarTrcDist_I = NULL;
  PyArrayObject *pyarTrcRho_I = NULL;
  PyArrayObject *pyarTrcGradRho_ID = NULL;
  PyArrayObject *pyarTrcBfield_ID = NULL;
  PyArrayObject *pyarTrcTbr_I = NULL;
  PyArrayObject *pyarTrcTbrIQUV_IP = NULL;
  PyArrayObject *pyarTrcLast_I = NULL;
  PyArrayObject *pyarTrcIrays_I = NULL;


  /* Constants */
  /* double const pi = 3.1415926535897931; // -- Now pi as a macros */

  /* Elemental input variables */
  double  Tol, RIntgSph;
  /* Array parameters */
  double *Pos_ID = NULL, *Dir_ID = NULL, *DS_I = NULL, *Arclen_I = NULL;
  double *Dist_I, *Tbr_I = NULL, *OpDepth_I = NULL, *Param_I = NULL;
  short  *Flags_I = NULL;
  double *Rho_I = NULL, *GradRho_ID = NULL;
  double *PosPr_ID = NULL, *DirPr_ID = NULL, *DS_New_I = NULL;
  double *TbrIQUV_IP = NULL, *Bfield_ID = NULL;
  double *DistToCrSurf_I = NULL;

  int32_t *TrcIrays_I = NULL, *TrcLast_I = NULL;
  double *TrcPos_ID = NULL, *TrcDir_ID = NULL, *TrcArclen_I = NULL;
  double *TrcRho_I = NULL,  *TrcGradRho_ID = NULL, *TrcDist_I = NULL;
  double *TrcBfield_ID = NULL;
  double *TrcTbrIQUV_IP = NULL, *TrcTbr_I = NULL;
  /* double * = NULL, * = NULL, * = NULL; */
  /* double * = NULL, * = NULL, * = NULL; */
  
  int ndimPos, ndimDir, ndimTb, ndimOpDepth, ndimFlags, nRay;
  int ndimDS, ndimRho, ndimGradRho, ndimPosPr, ndimDirPr, ndimDS_New;
  int ndimDistToCrSurf, nIter, nRay3;
  int i, iRay, iter, scattering, itraj;
  int nactive;    /* Number of active rays in a thread */
  /* int nactive_global;    // Number of active rays */

  int ntraj; /* Number of trajectories traced */
  int npmax; /* Max number of points - max trajectory length */

  struct TraceBits { /* Bit flags of the traced parameters */
    unsigned int Pos    : 1;
    unsigned int Dir    : 1;
    unsigned int Arclen : 1;
    unsigned int Dist   : 1;
    unsigned int Rho    : 1;
    unsigned int gRho   : 1;
    unsigned int Bfield : 1;
    unsigned int Tbr    : 1;
    unsigned int Stokes : 1;
  } toTrace = {0,0,0,0,0,0,0,0,0}; 

  double percent; /* Percent of active rays = (nactive/nRay)*100% */
  double dblnRay;
  int icnt = 0;
  int icall;
  double rsph2; 
  double minsd; /* minimal solar distance */
  int rtmode = 0; /* Ray tracing mode: 1:just trajectories, 2:Tbr, 3:TbrIQUV */
  void (* plasma_parameters)();
 const char *plfname; /* Name of plasma parameters function */
  char plfsoname[300]; /* [Path]name of plasma parameters DL library */

  /* Dynamic linking data */
  void *dlh;            /* DL library handle */
  char *error;

  struct param *prm; /* Access parameters by structure field names */

  /* The thread number, and what rays to begin and end calculating at */
  int thread_num, start, stop;

  char s_pos[16]; 
  char s_dir[16]; 

  strcpy((char *)s_pos, str_pos);
  strcpy((char *)s_dir, str_dir);

  /* Parse the (void *) array into the standard arguments 
   * required for the raytracing alogrithm. */
  self = (PyObject *)thread_args[0];
  args = (PyObject *)thread_args[1];
  thread_num = (int)thread_args[2];



  if (!PyArg_ParseTuple(args, "OOOOOOsdiiiOOOOOOOOOOOO:trace_beam",
  			&pyobjParam_I,
  			&pyobjPos_ID,
  			&pyobjDir_ID,
  			&pyobjArclen_I,
  			&pyobjDist_I,
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
			&pyobjTrace))
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
  if(thread_num == 0)
    printf("Library name: %s\n\n", plfsoname);

  /* Get the plasma_parameters pointer to the library function */
  plasma_parameters = dlsym(dlh, plfname);
  if ((error = dlerror()) != NULL)  {
    fputs(error, stderr);
    return NULL;
  }

  /* Convert Python objects to Numpy array objects */

  OBJ_TO_INOUT_ARRAY(pyobjParam_I,  NPY_DOUBLE,  pyarParam_I);
  OBJ_TO_INOUT_ARRAY(pyobjPos_ID,   NPY_DOUBLE,  pyarPos_ID);
  OBJ_TO_INOUT_ARRAY(pyobjDir_ID,   NPY_DOUBLE,  pyarDir_ID);
  OBJ_TO_INOUT_ARRAY(pyobjArclen_I, NPY_DOUBLE,  pyarArclen_I);
  OBJ_TO_INOUT_ARRAY(pyobjDist_I,   NPY_DOUBLE,  pyarDist_I);
  OBJ_TO_INOUT_ARRAY(pyobjDS_I,     NPY_DOUBLE,  pyarDS_I);
  if (rtmode == 2)    /* Brightness temperature */
    OBJ_TO_INOUT_ARRAY(pyobjTbr_I, NPY_DOUBLE, pyarTbr_I);
  if (rtmode == 3) {  /* Brightness temperature in Stokes param. I,Q,U,V */
    OBJ_TO_INOUT_ARRAY(pyobjTbrIQUV_IP, NPY_DOUBLE, pyarTbrIQUV_IP);
  }
  if (rtmode == 2 || rtmode == 3)
    OBJ_TO_INOUT_ARRAY(pyobjOpDepth_I, NPY_DOUBLE, pyarOpDepth_I);
  OBJ_TO_INOUT_ARRAY(pyobjRho_I,       NPY_DOUBLE, pyarRho_I);
  OBJ_TO_INOUT_ARRAY(pyobjGradRho_ID,  NPY_DOUBLE, pyarGradRho_ID);
  if (rtmode == 3)    /* Magnetic field for Tb in Stokes param. I,Q,U,V */
    OBJ_TO_INOUT_ARRAY(pyobjBfield_ID, NPY_DOUBLE, pyarBfield_ID);
  OBJ_TO_INOUT_ARRAY(pyobjPosPr_ID, NPY_DOUBLE, pyarPosPr_ID);
  OBJ_TO_INOUT_ARRAY(pyobjDirPr_ID, NPY_DOUBLE, pyarDirPr_ID);
  OBJ_TO_INOUT_ARRAY(pyobjDS_New_I, NPY_DOUBLE, pyarDS_New_I);
  OBJ_TO_INOUT_ARRAY(pyobjDistToCrSurf_I, NPY_DOUBLE, pyarDistToCrSurf_I);
  OBJ_TO_INOUT_ARRAY(pyobjFlags_I, NPY_SHORT, pyarFlags_I);

 /* If the Trace parameter is not None: */

  printf("pyobjTrace =  %p at %p (address), before PyObject_RichCompareBool()\n"
	 , pyobjTrace,  &pyobjTrace);

  if (PyObject_RichCompareBool(pyobjTrace, Py_None, Py_NE)) {
    /* Mandatory fields of the traj class:
     * Check if the 'pos', 'iray', and 'last' attributes are present
     * in the Trace (a class instance) argument. */
    
    int hasattr = 1, ntnp, ntnp3, ntnp4;
    hasattr &= PyObject_HasAttrString(pyobjTrace, s_pos);
    hasattr &= PyObject_HasAttrString(pyobjTrace, "irays");
    hasattr &= PyObject_HasAttrString(pyobjTrace, "last");
    hasattr &= PyObject_HasAttrString(pyobjTrace, "ntraj");
    hasattr &= PyObject_HasAttrString(pyobjTrace, "npmax");
    if (!hasattr) {
      if(thread_num == 0) printf("Error: The Trace parameter must have " \
	     "'pos', 'iray', 'last', 'ntraj', and 'npmax' attributes.\n");
      return NULL;
    }
    toTrace.Pos = 1; /* Always trace positions, if Trace != None */
    
    if(thread_num == 0) printf("hasattr done\n");

    pyobjTrcPos_ID = PyObject_GetAttrString(pyobjTrace, s_pos);
    if(thread_num == 0) printf("pyobjTrcPos_ID is at %p\n", pyobjTrcPos_ID);
    pyobjTrcIrays_I = PyObject_GetAttrString(pyobjTrace, "irays");
    if(thread_num == 0) printf("pyobjTrcIrays_I is at %p\n", pyobjTrcIrays_I);
    pyobjTrcLast_I = PyObject_GetAttrString(pyobjTrace, "last");
    if(thread_num == 0) printf("pyobjLast_I is at %p\n", pyobjTrcLast_I);
    pyobjTrcNtraj = PyObject_GetAttrString(pyobjTrace, "ntraj");
    if(thread_num == 0) printf("pyobjTrcNtraj is at %p\n", pyobjTrcNtraj);
    pyobjTrcNpmax = PyObject_GetAttrString(pyobjTrace, "npmax");
    if(thread_num == 0) printf("pyobjTrcNpmax is at %p\n", pyobjTrcNpmax);
    OBJ_TO_INOUT_ARRAY(pyobjTrcPos_ID, NPY_DOUBLE, pyarTrcPos_ID);
    if(thread_num == 0) printf("pyarTrcPos_ID is at %p\n", pyarTrcPos_ID);
    OBJ_TO_INOUT_ARRAY(pyobjTrcIrays_I, NPY_INT, pyarTrcIrays_I);

    if(thread_num == 0) printf("pyarTrcIrays_I is at %p\n", pyarTrcIrays_I);
    OBJ_TO_INOUT_ARRAY(pyobjTrcLast_I, NPY_INT, pyarTrcLast_I);
    if(thread_num == 0) printf("pyarTrcLast_I is at %p\n", pyarTrcLast_I);

    ntraj =  (int)PyInt_AS_LONG(pyobjTrcNtraj);  /* # of trajectories traced */
    npmax = (int)PyInt_AS_LONG(pyobjTrcNpmax); /* max# of pts in each traj.*/
    ntnp = ntraj*npmax;
    ntnp3 = ntnp*3;
    ntnp4 = ntnp3 + ntnp;
    CHKSIZE_AND_GETPTR(pyarTrcPos_ID,  ntnp3, double, TrcPos_ID);
    CHKSIZE_AND_GETPTR(pyarTrcIrays_I, ntraj, int,    TrcIrays_I);
    CHKSIZE_AND_GETPTR(pyarTrcLast_I,  ntraj, int,    TrcLast_I);
    if(thread_num == 0) {
      printf("ntraj = %d, npmax = %d, ntnp = %d\n", ntraj, npmax, ntnp);
      for (i = 0; i < ntraj; i++) 
	printf("Irays[%d] = %d\n", i, TrcIrays_I[i]); 
    }

    /*
     * Get pointers to the data of arrays requested for tracking by the user
     */
    if (PyObject_HasAttrString(pyobjTrace, s_dir)) {
      pyobjTrcDir_ID  = PyObject_GetAttrString(pyobjTrace, s_dir);
      OBJ_TO_INOUT_ARRAY(pyobjTrcDir_ID, NPY_DOUBLE, pyarTrcDir_ID);
      CHKSIZE_AND_GETPTR(pyarTrcDir_ID,  ntnp3, double, TrcDir_ID);
      if(thread_num == 0) printf("TrcDir_ID is at %p\n", TrcDir_ID);
      toTrace.Dir = 1; }
    if (PyObject_HasAttrString(pyobjTrace, "arclen")) {
      pyobjTrcArclen_I  = PyObject_GetAttrString(pyobjTrace, "arclen");
      OBJ_TO_INOUT_ARRAY(pyobjTrcArclen_I, NPY_DOUBLE, pyarTrcArclen_I);
      CHKSIZE_AND_GETPTR(pyarTrcArclen_I,  ntnp, double, TrcArclen_I);
      if(thread_num == 0) printf("TrcArclen_I is at %p\n", TrcArclen_I);
      toTrace.Arclen = 1; }
    if (PyObject_HasAttrString(pyobjTrace, "dist")) {
      pyobjTrcDist_I  = PyObject_GetAttrString(pyobjTrace, "dist");
      OBJ_TO_INOUT_ARRAY(pyobjTrcDist_I, NPY_DOUBLE, pyarTrcDist_I);
      CHKSIZE_AND_GETPTR(pyarTrcDist_I,  ntnp, double, TrcDist_I);
      if(thread_num == 0) printf("TrcDist_ID is at %p\n", TrcDist_I);
      toTrace.Dist = 1; }
    if (PyObject_HasAttrString(pyobjTrace, "rho")) {
      pyobjTrcRho_I = PyObject_GetAttrString(pyobjTrace, "rho");
      OBJ_TO_INOUT_ARRAY(pyobjTrcRho_I, NPY_DOUBLE, pyarTrcRho_I);
      CHKSIZE_AND_GETPTR(pyarTrcRho_I,  ntnp, double, TrcRho_I);
      if(thread_num == 0) printf("TrcRho_I is at %p\n", TrcRho_I);
      toTrace.Rho = 1; }
    if (PyObject_HasAttrString(pyobjTrace, "grho")) { 
      pyobjTrcGradRho_ID = PyObject_GetAttrString(pyobjTrace, "grho");
      OBJ_TO_INOUT_ARRAY(pyobjTrcGradRho_ID, NPY_DOUBLE, pyarTrcGradRho_ID);
      CHKSIZE_AND_GETPTR(pyarTrcGradRho_ID,  ntnp3, double, TrcGradRho_ID);
      if(thread_num == 0) printf("TrcGradRho_ID is at %p\n", TrcGradRho_ID);
      toTrace.gRho = 1; }
    if (PyObject_HasAttrString(pyobjTrace, "bfield")) {
      pyobjTrcBfield_ID = PyObject_GetAttrString(pyobjTrace, "bfield");
      OBJ_TO_INOUT_ARRAY(pyobjTrcBfield_ID, NPY_DOUBLE, pyarTrcBfield_ID);
      CHKSIZE_AND_GETPTR(pyarTrcBfield_ID,  ntnp3, double, TrcBfield_ID);
      if(thread_num == 0) printf("TrcBfield_ID is at %p\n", TrcBfield_ID);
      toTrace.Bfield = 1; }
    if (PyObject_HasAttrString(pyobjTrace, "tbr")) {
      pyobjTrcTbr_I = PyObject_GetAttrString(pyobjTrace, "tbr");
      OBJ_TO_INOUT_ARRAY(pyobjTrcTbr_I, NPY_DOUBLE, pyarTrcTbr_I);
      CHKSIZE_AND_GETPTR(pyarTrcTbr_I,  ntnp, double, TrcTbr_I);
      if(thread_num == 0) printf("TrcTbr_I is at %p\n", TrcTbr_I);
      toTrace.Tbr = 1; }
    if (PyObject_HasAttrString(pyobjTrace, "tbriquv")) { 
      pyobjTrcTbrIQUV_IP = PyObject_GetAttrString(pyobjTrace, "tbriquv");
      OBJ_TO_INOUT_ARRAY(pyobjTrcTbrIQUV_IP, NPY_DOUBLE, pyarTrcTbrIQUV_IP);
      CHKSIZE_AND_GETPTR(pyarTrcTbrIQUV_IP,  ntnp4, double, TrcTbrIQUV_IP);
      if(thread_num == 0) printf("TrcTbrIQUV_IP is at %p\n", TrcTbrIQUV_IP);
      toTrace.Stokes = 1; }
  } /* if (PyObject_RichCompareBool(pyobjTrace, Py_None, Py_NE)) */
  else 
    printf("pyobjTrace = None\n");

  printf("1 pyobjTrace:     %d\n", PyArray_REFCOUNT(pyobjTrace));

  /* extract sizes, check if correct, and get pointers to arrays */
  nRay3 = PyArray_SIZE(pyarPos_ID); /* Size of Pos_ID is used as reference */ 
  nRay = nRay3/3; /* All arrays with _D are composed of coord. triplets */ 

    /* Check array sizes and get pointers to array data from PyArrayObject-s */ 

  if(thread_num == 0)
    printf("NPARAM = %d\n", NPARAM);

  CHKSIZE_AND_GETPTR(pyarParam_I,    NPARAM, double,  Param_I);
  CHKSIZE_AND_GETPTR(pyarPos_ID,     nRay3,  double,  Pos_ID);
  CHKSIZE_AND_GETPTR(pyarDir_ID,     nRay3,  double,  Dir_ID);
  CHKSIZE_AND_GETPTR(pyarArclen_I,   nRay,   double,  Arclen_I);
  CHKSIZE_AND_GETPTR(pyarDist_I,     nRay,   double,  Dist_I);
  CHKSIZE_AND_GETPTR(pyarDS_I,       nRay,   double,  DS_I);
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
    CHKSIZE_AND_GETPTR(pyarBfield_ID, nRay3, double,    Bfield_ID);
  CHKSIZE_AND_GETPTR(pyarPosPr_ID,   nRay3, double,     PosPr_ID);
  CHKSIZE_AND_GETPTR(pyarDirPr_ID,   nRay3, double,     DirPr_ID);
  CHKSIZE_AND_GETPTR(pyarDS_New_I,   nRay,  double,     DS_New_I);
  CHKSIZE_AND_GETPTR(pyarDistToCrSurf_I, nRay, double,  DistToCrSurf_I);

  prm = (void *) Param_I; /* Access parameters by structure field names */

  if(thread_num == 0)
    printf("Freq = %g, Tol = %g, RIntgSph = %g, nIter = %d\n",
  	 prm->Freq, prm->Tol, RIntgSph, nIter);

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
  icnt = 0; /* Step counter */

  /* Determine which rays this thread needs to trace */
  start = nRay/NUM_THREADS*thread_num;
  stop = nRay/NUM_THREADS*(thread_num + 1);

  /* If this is the last thread make sure to finish the caluclation
     i.e. don't leave rays out because of a rounding error */
  if (thread_num == NUM_THREADS - 1)
    stop = nRay; 

  for (iter = 0; iter < nIter; iter++) {

    /* Calculate this step of the beam */
    advance_beam(start,
		 stop,
		 thread_num,
		 nRay, 
		 Param_I,
		 Pos_ID, 
		 Dir_ID, 
		 Arclen_I,
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

    /*
     * Find solar distances for the current thread,
     * mark ray INACTIVE if it is outside the integration sphere,
     * and count the inactive rays
     */
    nactive = stop-start;
    for (iRay = start; iRay < stop; iRay++) {
      double r2; /* Solar distance squared */
      size_t i3;
      if (BIT_IS_ON(INACTIVE,iRay)) {
	--nactive;
	continue; //================>>>
      }
      i3 = 3*iRay;
      r2 = dot_product(&Pos_ID[i3], &Pos_ID[i3]); 
      Dist_I[iRay] = sqrt(r2);             /* Solar distance */
      if (r2 > rsph2) {
	SET_BIT(INACTIVE,iRay); /* The i-th ray is done */
	--nactive;
      }
    }

    if (toTrace.Pos) 
      /* 
       * For the rays we wish to track, go through and keep track of them 
       */
      for(itraj = 0; itraj < ntraj; itraj++) {
	iRay = TrcIrays_I[itraj];
	if (BIT_IS_ON(INACTIVE, iRay)) continue; //=======================>>>
	if (TrcLast_I[itraj] >= npmax) continue; //=======================>>>
	if (iRay >= start && iRay < stop) {
	  int npt = TrcLast_I[itraj];   /* Last saved point position */
	  size_t i = itraj*npmax + npt; /* Index into Trc* arrays */
	  size_t i3 = 3*i;              /* Points at xyz components in Trc* */
	  size_t j3 = 3*iRay;           /* Index into source's xyz components */
	  TrcPos_ID[i3] =   Pos_ID[j3];
	  TrcPos_ID[i3+1] = Pos_ID[j3+1];
	  TrcPos_ID[i3+2] = Pos_ID[j3+2];

	  if (toTrace.Dir) {
	    TrcDir_ID[i3] =   Dir_ID[j3];
	    TrcDir_ID[i3+1] = Dir_ID[j3+1];
	    TrcDir_ID[i3+2] = Dir_ID[j3+2];
	  }
	  if (toTrace.Arclen) {
	    TrcArclen_I[i] = Arclen_I[iRay];
	  }
	  if (toTrace.Dist) {
	    TrcDist_I[i] = Dist_I[iRay];
	  }
	  if (toTrace.Rho) {
	    TrcRho_I[i] = Rho_I[iRay];
	  }
	  if (toTrace.gRho) {
	    TrcGradRho_ID[i3] =   GradRho_ID[j3];
	    TrcGradRho_ID[i3+1] = GradRho_ID[j3+1];
	    TrcGradRho_ID[i3+2] = GradRho_ID[j3+2];
	  }
	  if (toTrace.Bfield) {
	    TrcBfield_ID[i3] =   Bfield_ID[j3];
	    TrcBfield_ID[i3+1] = Bfield_ID[j3+1];
	    TrcBfield_ID[i3+2] = Bfield_ID[j3+2];
	  }
	  if (toTrace.Tbr) {
	    TrcTbr_I[i] =   Tbr_I[iRay];
	  }
	  if (toTrace.Stokes) {
	    size_t i4 = 4*i;    /* Points at iquv Stokes components in Trc* */
	    size_t j4 = 4*iRay; /* Points at iquv Stokes components in srces */
	    TrcTbrIQUV_IP[i4] =   TbrIQUV_IP[j4];
	    TrcTbrIQUV_IP[i4+1] = TbrIQUV_IP[j4+1];
	    TrcTbrIQUV_IP[i4+2] = TbrIQUV_IP[j4+2];
	    TrcTbrIQUV_IP[i4+3] = TbrIQUV_IP[j4+3];
	  }

	  TrcLast_I[itraj]++;
	} /* if(iRay >= start && iRay < stop) */
      } /* for(itraj = 0; itraj < ntraj; itraj++) */
    

    /* If this thread is the root thread, then print out the diagnostics */
    if (thread_num == 0) {
      /* 
       * Find the ray which is closest to the sun 
       * and determine if there are active rays
       * The previous nactive counter specific for 
       * current thread will be recalculated to represent
       * total number of inactive rays.
       * Therefore, the control will not exit thread 0
       * until all the threads are done.
       */
      nactive = nRay;
      minsd = DBL_MAX;  /* Just a very big value */
      for (iRay = 0; iRay < nRay; iRay++) {
	double r; /* Solar distance */
	r = Dist_I[iRay];
	if ((minsd > r) && (BIT_IS_OFF(PENETR,iRay))) minsd = r;
	if (BIT_IS_ON(INACTIVE,iRay)) nactive--;
      }

      ////         icall = (int)prm->callcount;

      /* Here nactive counts for all the rays */
      if ((icnt % (int)(10.0/(prm->DeltaS))) == 0) {
	percent = 100.0*(double)nactive/dblnRay;
	printf("%5i: =====(%5i) %8.4f%% ====== Min Dist from Sun center: %g\n",
	       icnt, nactive, percent, minsd);
      }
      icnt++;   
    } /* if (thread_num == 0) */

    /* Terminate this thread if all the rays given have been calculated 
     * (unless you are the root thread in which case do not terminate 
     * until all the rays from all the threads have finished */

    if (nactive == 0) break; //--------All rays in thread are done --------->>


  } /* for (iter = 0; iter < nIter; iter++) */
 
  /* Determine the time it took to complete and print it out. */
  if(thread_num==0){
    printf("LOOP COUNTER = %d\n", icnt);
    printf("RAYS REMAINING = %d\n",nactive);
  }

  /*
   * Close the dynamically linked (DL) library
   */
  dlclose(dlh); 


  /* printf("Reference counts at exit:\n"); */
  /* printf("pyobjPos_ID:    %d\n", PyArray_REFCOUNT(pyobjPos_ID)); */
  /* printf("pyarPos_ID:     %d\n", PyArray_REFCOUNT(pyarPos_ID)); */
  /* printf("pyobjTrace:     %d\n", PyArray_REFCOUNT(pyobjTrace)); */
  /* printf("pyobjTrcPos_ID: %d\n", PyArray_REFCOUNT(pyobjTrcPos_ID)); */
  /* printf("pyarTrcPos_ID:  %d\n", PyArray_REFCOUNT(pyarTrcPos_ID)); */
  /* printf("pyobjRho_I:     %d\n", PyArray_REFCOUNT(pyobjRho_I)); */
  /* printf("pyarRho_I:      %d\n", PyArray_REFCOUNT(pyarRho_I)); */
  /* printf("%d\n", PyArray_REFCOUNT()); */

  Py_DECREF(pyarParam_I);
  Py_DECREF(pyarPos_ID);
  Py_DECREF(pyarDir_ID);
  Py_DECREF(pyarArclen_I);
  Py_DECREF(pyarDist_I);
  Py_DECREF(pyarDS_I);
  if (rtmode == 2) Py_DECREF(pyarTbr_I);
  if (rtmode == 3) { Py_DECREF(pyarTbrIQUV_IP); Py_DECREF(pyarBfield_ID); }
  Py_DECREF(pyarOpDepth_I);
  Py_DECREF(pyarFlags_I);
  Py_DECREF(pyarRho_I);
  Py_DECREF(pyarGradRho_ID);
  Py_DECREF(pyarPosPr_ID);
  Py_DECREF(pyarDirPr_ID);
  Py_DECREF(pyarDS_New_I);
  Py_DECREF(pyarDistToCrSurf_I);
  /* Py_DECREF(); */

  if (toTrace.Pos) {
    Py_DECREF(pyobjTrcNtraj); 
    Py_DECREF(pyobjTrcNpmax); 
    Py_DECREF(pyarTrcLast_I); Py_DECREF(pyobjTrcLast_I); 
    Py_DECREF(pyarTrcIrays_I); Py_DECREF(pyobjTrcIrays_I); 
    Py_DECREF(pyarTrcPos_ID); Py_DECREF(pyobjTrcPos_ID); 
  }
  if (toTrace.Dir) { Py_DECREF(pyarTrcDir_ID); Py_DECREF(pyobjTrcDir_ID); }
  if (toTrace.Arclen) 
    { Py_DECREF(pyarTrcArclen_I); Py_DECREF(pyobjTrcArclen_I); }
  if (toTrace.Dist) { Py_DECREF(pyarTrcDist_I); Py_DECREF(pyobjTrcDist_I); }
  if (toTrace.Rho) { Py_DECREF(pyarTrcRho_I); Py_DECREF(pyobjTrcRho_I); }
  if (toTrace.gRho) 
    { Py_DECREF(pyarTrcGradRho_ID); Py_DECREF(pyobjTrcGradRho_ID); }
  if (toTrace.Bfield) 
    { Py_DECREF(pyarTrcBfield_ID); Py_DECREF(pyobjTrcBfield_ID); }
  if (toTrace.Tbr) { Py_DECREF(pyarTrcTbr_I); Py_DECREF(pyobjTrcTbr_I); }
  if (toTrace.Stokes) 
    { Py_DECREF(pyarTrcTbrIQUV_IP); Py_DECREF(pyobjTrcTbrIQUV_IP); }
  /* if (toTrace.) { Py_DECREF(pyarTrc); Py_DECREF(pyobjTrc); } */

  Py_INCREF(Py_None);
  return Py_None;

} /* PyObject *threadFunc(void *thread_args[]) */




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
  OBJ_TO_INOUT_ARRAY(pyobjObs_D, NPY_DOUBLE,    pyarObs_D);
  OBJ_TO_INOUT_ARRAY(pyobjXruler_I, NPY_DOUBLE, pyarXruler_I);
  OBJ_TO_INOUT_ARRAY(pyobjYruler_I, NPY_DOUBLE, pyarYruler_I);
  OBJ_TO_INOUT_ARRAY(pyobjDir_ID, NPY_DOUBLE,   pyarDir_ID);

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
  OBJ_TO_INOUT_ARRAY(pyobjObs_D, NPY_DOUBLE,  pyarObs_D);
  OBJ_TO_INOUT_ARRAY(pyobjDir_ID, NPY_DOUBLE, pyarDir_ID);
  OBJ_TO_INOUT_ARRAY(pyobjIsec_I, NPY_DOUBLE, pyarIsec_I);
  OBJ_TO_INOUT_ARRAY(pyobjPos_ID, NPY_DOUBLE, pyarPos_ID);

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
  OBJ_TO_INOUT_ARRAY(pyobjObs_D, NPY_DOUBLE,    pyarObs_D);
  OBJ_TO_INOUT_ARRAY(pyobjXruler_I, NPY_DOUBLE, pyarXruler_I);
  OBJ_TO_INOUT_ARRAY(pyobjYruler_I, NPY_DOUBLE, pyarYruler_I);
  OBJ_TO_INOUT_ARRAY(pyobjIsec_I, NPY_DOUBLE,   pyarIsec_I);
  OBJ_TO_INOUT_ARRAY(pyobjDir_ID, NPY_DOUBLE,   pyarDir_ID);
  OBJ_TO_INOUT_ARRAY(pyobjPos_ID, NPY_DOUBLE,   pyarPos_ID);

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

  /* Adjust reference counters, incremented by PyArray_FROM_OTF() */
  Py_DECREF(pyarObs_D);
  Py_DECREF(pyarXruler_I);
  Py_DECREF(pyarYruler_I);
  Py_DECREF(pyarIsec_I);
  Py_DECREF(pyarPos_ID);
  Py_DECREF(pyarDir_ID);

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
  OBJ_TO_INOUT_ARRAY(pyobjParam_I, NPY_DOUBLE,        pyarParam_I);
  OBJ_TO_INOUT_ARRAY(pyobjPos_ID, NPY_DOUBLE,         pyarPos_ID);
  OBJ_TO_INOUT_ARRAY(pyobjRho_I, NPY_DOUBLE,          pyarRho_I);
  OBJ_TO_INOUT_ARRAY(pyobjGradRho_ID, NPY_DOUBLE,     pyarGradRho_ID);
  if (rtmode == 3)    /* Magnetic field for Tb in Stokes param. I and V */
    OBJ_TO_INOUT_ARRAY(pyobjBfield_ID, NPY_DOUBLE,    pyarBfield_ID);
  OBJ_TO_INOUT_ARRAY(pyobjDS_I, NPY_DOUBLE,           pyarDS_I);

  

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

  /* Adjust reference counters, incremented by PyArray_FROM_OTF() */

  Py_DECREF(pyarParam_I);
  Py_DECREF(pyarPos_ID);
  Py_DECREF(pyarRho_I);
  Py_DECREF(pyarGradRho_ID);
  if (rtmode == 3) Py_DECREF(pyarBfield_ID); 
  Py_DECREF(pyarDS_I);
  Py_DECREF(pyarFlags_I);
 

  Py_INCREF(Py_None);
  return Py_None;

}

/*
 * The python frontend passes an array of parameters to make a streamer
 * Written by Mark D. Benjamin, July 2011
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

