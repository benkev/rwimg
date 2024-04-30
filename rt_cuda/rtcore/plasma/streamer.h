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
			   double dStalkStrength);
Streamers_T Streamers_new();
int Streamers_length(Streamers_T oStreamers);
double *Streamers_getSpherical(Streamers_T oStreamers,int iIndex);
double *Streamers_getCartesean(Streamers_T oStreamers,int iIndex);
void Streamers_calculateBField(Streamers_T oStreamers, 
			       double x,
			       double y,
			       double z, 
			       double B[3]);

double Streamers_CalculateDensity(Streamers_T oStreamers,
				  double x,
				  double y,
				  double z,
				  double r,
				  double Ne);

#endif
