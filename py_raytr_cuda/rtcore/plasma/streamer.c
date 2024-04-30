/****************************/
/* Streamer.c               */
/* Author: Mark Benjamin    */
/****************************/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "streamer.h"


struct Streamer{
  
  Streamer_T oNextStreamer;
  double pSpherical[3];//dRadius,dTheta,dPhi;
  double pCartesean[3];//dX,dY,dZ;
  double dDensity,dBaseStrength,dStalkStrength;
  double dOrientation;

};

struct Streamers{

  int iLength;
  Streamer_T oFirstStreamer;

};


int Streamers_makeStreamer(Streamers_T oStreamers, 
			   double dRadius, 
			   double dTheta, 
			   double dPhi,
			   double dOrientation,
			   double dDensity, 
			   double dBaseStrength, 
			   double dStalkStrength){
  
  Streamer_T  oStreamer = (struct Streamer*)malloc(sizeof(struct Streamer));
  if(oStreamer == NULL)
    return 1;
  

  
  oStreamer->pSpherical[0] = dRadius;
  oStreamer->pSpherical[1] = dTheta;
  oStreamer->pSpherical[2] = dPhi;

  oStreamer->pCartesean[0] = sin(dTheta)*cos(dPhi)*dRadius;
  oStreamer->pCartesean[1] = sin(dTheta)*sin(dPhi)*dRadius;
  oStreamer->pCartesean[2] = cos(dTheta)*dRadius;

  oStreamer->dDensity = dDensity;
  oStreamer->dBaseStrength = dBaseStrength;
  oStreamer->dStalkStrength = dStalkStrength;
  
  oStreamer->dOrientation = dOrientation;


  /* Set up the linked list */
  oStreamer->oNextStreamer = oStreamers->oFirstStreamer;
  oStreamers->oFirstStreamer = oStreamer;

  oStreamers->iLength++;
  return 0;
  
}

Streamers_T Streamers_new(){

  Streamers_T oStreamers;

  oStreamers = (struct Streamers*)malloc(sizeof(struct Streamers));
  if(oStreamers == NULL)
    return NULL;

  oStreamers->iLength = 0;
  oStreamers->oFirstStreamer = NULL;

  return oStreamers;

}

int Streamers_length(Streamers_T oStreamers){

  return oStreamers->iLength;

}



void Streamers_calculateDensity(Streamers_T oStreamers, double x,double y,double z, double r, double *Ne){

  double Cn0;
  double sig,dsig2;
  double bet;

  Streamer_T oCurrentStreamer;

  

  if(oStreamers==NULL)
    return;

  oCurrentStreamer = oStreamers->oFirstStreamer;

  while(oCurrentStreamer!=NULL){
    Cn0 = oCurrentStreamer->dDensity;
    sig = 0.15*(1.0 + exp(-pow(r - 1.5, 2)/0.22));
    dsig2 = 2.0*pow(sig,2);
    bet = acos((x*oCurrentStreamer->pCartesean[0] + 
		y*oCurrentStreamer->pCartesean[1] + 
		z*oCurrentStreamer->pCartesean[2])/r);
    *Ne = *Ne + *Ne*Cn0*exp(-pow(bet,2)/dsig2);
    oCurrentStreamer = oCurrentStreamer->oNextStreamer;
  }

  return;

}


/*
 * Revision Author:
 *   Mark Benjamin
 * 2011-Jun-24: This function to calculate the Bfield due to coronal streamers 
 * has been added.
 */
void Streamers_calculateBField(Streamers_T oStreamers, double x,double y,double z, double B[3]){

  double static B0m;          // ~ field at the parabola vertex
  double static B0s;          // ~ field in the stalk flux tubes
  double static C2 = 2;            // Parabola vertex height    
  double static yshm = 0.5;
  double static yshs = 0.125;     // Distance of the stalk's B extremum from 0
  double static ccy = 5;          // Y-compression coefficient

  double static ysi;
    
  double static C1;


    
  /* rho, theta, phi */
  int static num_streamers = 1;
  double static spherical_coord_location[3];
  double rotation_matrix_1[2][2];
  double rotation_matrix_2[2][2];
  double rotation_matrix_3[2][2];

  double unprimed[3];
  double hold, hold1[2], hold2[2];
  double double_primed[3];

  /* rho, theta, phi */
  double spherical[3];

  double delta_phi,delta_theta;
  double new_x,new_y,new_z;

  double yr,yl,z2,dif12;

  double r1; 
  double r2;
  double rqub1;
  double rqub2;

  double f1;
  double f2;
  double f;

  double B_length = 0;

  int i = 0;

  Streamer_T oCurrentStreamer;


  ysi = 1.0/yshs;     // Inverse of distance of the stalk's B extremum from 0
    
  C1 = C2/(yshm*yshm);


  B[0] = 0;
  B[1] = 0;
  B[2] = 0;

  if(oStreamers==NULL)
    return;

  oCurrentStreamer = oStreamers->oFirstStreamer;
  
  while(oCurrentStreamer!=NULL){

    B0m = oCurrentStreamer->dBaseStrength;
    B0s = oCurrentStreamer->dStalkStrength;

    i++;


    spherical_coord_location[0] = oCurrentStreamer->pSpherical[0]; 
    spherical_coord_location[1] = oCurrentStreamer->pSpherical[1];
    spherical_coord_location[2] = oCurrentStreamer->pSpherical[2];


    rotation_matrix_1[0][0] = cos(spherical_coord_location[2]);
    rotation_matrix_1[0][1] = sin(spherical_coord_location[2]);
    rotation_matrix_1[1][0] = -sin(spherical_coord_location[2]);
    rotation_matrix_1[1][1] = cos(spherical_coord_location[2]);

    rotation_matrix_2[0][0] = cos(spherical_coord_location[1]);
    rotation_matrix_2[0][1] = sin(spherical_coord_location[1]);
    rotation_matrix_2[1][0] = -sin(spherical_coord_location[1]);
    rotation_matrix_2[1][1] = cos(spherical_coord_location[1]);

    rotation_matrix_3[0][0] = cos(oCurrentStreamer->dOrientation);
    rotation_matrix_3[0][1] = -sin(oCurrentStreamer->dOrientation);
    rotation_matrix_3[1][0] = sin(oCurrentStreamer->dOrientation);
    rotation_matrix_3[1][1] = cos(oCurrentStreamer->dOrientation);

    hold1[0] = x;
    hold1[1] = y;

    mvmul(2, rotation_matrix_1, hold1, hold2); 
  
    new_y = hold2[1];

    hold1[1] = hold2[0];
    hold1[0] = z;

    mvmul(2, rotation_matrix_2, hold1, hold2);
    new_z = hold2[0];
    new_x = hold2[2];


    if(new_z < 1.2){
      return;
    }

    /* Rotate the streamer based on the orientation */
    hold1[0] = new_x;
    hold1[1] = new_y;

    mvmul(2, rotation_matrix_3, hold1, hold2);

    new_x = hold2[0];
    new_y = hold2[1];

    f1 =  exp(-(pow((ysi*new_x),2) + ccy*pow(((ysi*new_y) - 1.),2)));
    f2 = -exp(-(pow((ysi*new_x),2) + ccy*pow(((ysi*new_y) + 1.),2)));
    f = (f1 + f2);

    if (new_z < C2){
      f = f*exp(new_z-C2); // Rapidly fall off below the flux tube 
    }else{
      f = f*exp(-0.05*(new_z-C2)); // Gradually fall off with height 
    }

    z2 = new_z + C1*pow(new_y,2)- C2; // The space along Y is bent parabolically

    
    /*
     * Field in the streamer base
     */
    yr = new_y - yshm;
    yl = new_y + yshm;
    r1 = 0.5 + sqrt(pow(new_x,2) + pow(yr,2) + pow(z2,2));
    r2 = 0.5 + sqrt(pow(new_x,2) + pow(yl,2) + pow(z2,2));
    rqub1 = pow(r1,3);
    rqub2 = pow(r2,3);

    dif12 = (1./rqub1 - 1./rqub2);
    B[1] = B0m*(yr/rqub1 - yl/rqub2);
    B[0] = B0m*new_x*dif12;
    B[2] = B0m*new_z*dif12;

    B[2] = B[2] + f*B0s;
    oCurrentStreamer = oCurrentStreamer->oNextStreamer;

  }
}

