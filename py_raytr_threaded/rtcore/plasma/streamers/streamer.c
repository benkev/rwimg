/****************************/
/* Streamer.c               */
/* Author: Mark Benjamin    */
/****************************/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "streamer.h"


struct Streamer{

  double pSpherical[3];//dRadius,dTheta,dPhi;
  double pCartesean[3];//dX,dY,dZ;

};

struct Streamers{

  int iLenth;
  Streamer_T *pArray;

};

Streamer_T Streamers_addStreamer(double dRadius, double dTheta, double dPhi){
  
  Streamer_T  oStreamer = (struct Streamer*)malloc(sizeof(struct Streamer));
  if(oStreamer == NULL)
    return NULL;
  
  oStreamer->pSpherical[0] = dRadius;
  oStreamer->pSpherical[1] = dTheta;
  oStreamer->pSpherical[2] = dPhi;

  oStreamer->pCartesean[0] = sin(dTheta)*cos(dPhi);
  oStreamer->pCartesean[1] = sin(dTheta)*sin(dPhi);
  oStreamer->pCartesean[2] = cos(dTheta);

  /*oStreamer->dTheta = dTheta;
  oStreamer->dPhi = dPhi;

  oStreamer->dX = sin(dTheta)*cos(dPhi);
  oStreamer->dY = sin(dTheta)*sin(dPhi);
  oStreamer->dZ = cos(dTheta);

  oStreamer->dRadius = dRadius;
  */
  return oStreamer;

}

Streamers_T Streamers_new(char *filename){

  Streamers_T oStreamers;
  int iNumStreamers;
  double dRadius,dTheta,dPhi;
  Streamer_T *pArray;
  int iIndex;
  FILE *fp;


  oStreamers = (struct Streamers*)malloc(sizeof(struct Streamers));
  if(oStreamers == NULL)
    return NULL;

  fp = fopen(filename,"r");
  fscanf(fp, "%d", &iNumStreamers);

  oStreamers->iLenth = iNumStreamers;
  
  pArray = (Streamer_T *)malloc(sizeof(double)*iNumStreamers);
  if(oStreamers == NULL)
    return NULL;

  iIndex = 0;

  while(fscanf(fp, "%f %f %f", &dRadius, &dTheta, &dPhi) != EOF){
    pArray[iIndex] = Streamers_addStreamer(dRadius, dTheta, dPhi);
    iIndex++;
  }

  oStreamers->pArray = pArray;


}

int Streamers_length(Streamers_T oStreamers){

  return oStreamers->iLenth;

}

double *Streamers_getSpherical(Streamers_T oStreamers,int iIndex){

  return oStreamers->pArray[iIndex]->pSpherical;

}

double *Streamers_getCartesian(Streamers_T oStreamers, int iIndex){

  return oStreamers->pArray[iIndex]->pCartesean;

}





