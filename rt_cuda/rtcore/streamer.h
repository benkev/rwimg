/****************************/
/* Streamer.h               */
/* Author: Mark Benjamin    */
/****************************/


#ifndef STREAMER_INCLUDED
#define STREAMER_INCLUDED

typedef struct Streamers *Streamers_T;
typedef struct Streamer *Streamer_T;

int Streamers_makeStreamer(Streamers_T oStreamers, 
			   double dRadius, 
			   double dTheta, 
			   double dPhi, 
			   double dOrientation,
			   double dDensity, 
			   double dBaseStrength, 
			   double dStalkStrength,
			   double dScaleX,
			   double dScaleY,
			   double dScaleZ);

void Streamers_free(Streamers_T oStreamers);
Streamers_T Streamers_new();
int Streamers_length(Streamers_T oStreamers);
double *Streamers_getSpherical(Streamers_T oStreamers,int iIndex);
double *Streamers_getCartesean(Streamers_T oStreamers,int iIndex);
__device__ void Streamers_calculateBField(Streamers_T oStreamers, 
			       double x,
			       double y,
			       double z, 
			       double B[3]);

__device__ void Streamers_calculateDensity(Streamers_T oStreamers,
				  double x,
				  double y,
				  double z,
				  double r,
				  double *Ne);

#endif
