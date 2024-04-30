/****************************/
/* Streamer.h               */
/* Author: Mark Benjamin    */
/****************************/


#ifndef STREAMER_INCLUDED
#define STREAMER_INCLUDED

typedef struct Streamers *Streamers_T;
typedef struct Streamer *Streamer_T;

Streamer_T Streamers_addStreamer(double dRadius, double dTheta, double dPhi);
Streamers_T Streamers_new(char *filename);
int Streamers_length(Streamers_T oStreamers);
double *Streamers_getSpherical(Streamers_T oStreamers,int iIndex);
double *Streamers_getCartesean(Streamers_T oStreamers,int iIndex);

#endif
