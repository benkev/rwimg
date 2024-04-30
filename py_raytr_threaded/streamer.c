/****************************/
/* Streamer.c               */
/* Author: Mark Benjamin    */
/****************************/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "streamer.h"

/*
 * A definition of a streamer.
 */
struct Streamer{
  // A pointer to the next streamer
  Streamer_T oNextStreamer;
  // An array of spherical coords of the streamer on the sun
  double pSpherical[3];//dRadius,dTheta,dPhi
  // An array of cartesian coords of the streamer
  double pCartesean[3];//dX,dY,dZ
  // An array of x,y, and z scales
  double pScale[3];//ScaleX,ScaleY,ScaleZ
  // Strengths of the streamer
  double dDensity,dBaseStrength,dStalkStrength;
  // Orientation of the streamer between 0 and 2pi
  double dOrientation;

};

/*
 * A definition of a collection of streamers.
 */
struct Streamers{

  // The number of streamers
  int iLength;
  // A pointer to the first streamer in the list
  Streamer_T oFirstStreamer;

};

/* 
 * Creates a streamer and adds it to the list headed by oStreamers
 */
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
			   double dScaleZ){
  
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

  oStreamer->pScale[0] = dScaleX;
  oStreamer->pScale[1] = dScaleY;
  oStreamer->pScale[2] = dScaleZ;

  /* Set up the linked list */
  oStreamer->oNextStreamer = oStreamers->oFirstStreamer;
  oStreamers->oFirstStreamer = oStreamer;

  oStreamers->iLength++;
  return 0;
  
}

/*
 * This makes a new streamer list and returns it.
 */
Streamers_T Streamers_new(){

  Streamers_T oStreamers;


  oStreamers = (struct Streamers*)malloc(sizeof(struct Streamers));
  if(oStreamers == NULL)
    return NULL;

  oStreamers->iLength = 0;
  oStreamers->oFirstStreamer = NULL;

  return oStreamers;

}

/*
 * Frees all the memory assosiated with the streamer list oStreamers
 */
void Streamers_free(Streamers_T oStreamers){

  Streamer_T oCurrentStreamer,oNextStreamer;

  if(oStreamers == NULL)
    return;

  oCurrentStreamer = oStreamers->oFirstStreamer;

  while(oCurrentStreamer!=NULL){
    oNextStreamer = oCurrentStreamer->oNextStreamer;
    free(oCurrentStreamer);
    oCurrentStreamer = oNextStreamer;
  }

  free(oStreamers);
  

}

/*
 * Returns the number of streamers in the list oStreamers
 */
int Streamers_length(Streamers_T oStreamers){

  return oStreamers->iLength;

}


/*
 * Given a list of streamers, oStreamers, as well as a coordinate in space
 * and the plasma density, this function adds the density of the streamers
 * to the given density.
 */
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
 * Given a list of streamers, oStreamers, and a coordinate in space as well
 * as a Bfield vector, this function calculates the Bfield due to the streamer, 
 * and adds it to the given B field.
 */
void Streamers_calculateBField(Streamers_T oStreamers, double x, double y, 
			       double z, double B[3]){

  double static B0m;          // ~ field at the parabola vertex
  double static B0s;          // ~ field in the stalk flux tubes
  double static C2 = 1;            // Parabola vertex height    
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

  double Bx,By,Bz;

  int i = 0;

  Streamer_T oCurrentStreamer;

  B[0] = 0;
  B[1] = 0;
  B[2] = 0;

  if(oStreamers==NULL)
    return;


  ysi = 1.0/yshs;     // Inverse of distance of the stalk's B extremum from 0
    
  C1 = C2/(yshm*yshm);


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
    new_x = hold2[1];


    if(new_z < 1){
      oCurrentStreamer = oCurrentStreamer->oNextStreamer;
      continue;
     }

    /* Rotate the streamer based on the orientation */
    hold1[0] = new_x;
    hold1[1] = new_y;

    mvmul(2, rotation_matrix_3, hold1, hold2);

    new_x = hold2[0]/oCurrentStreamer->pScale[0];
    new_y = hold2[1]/oCurrentStreamer->pScale[1];

    /* We want to the streamer to be at the surface, not the center */
    new_z = (new_z - spherical_coord_location[0])/oCurrentStreamer->pScale[2]; 

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
    By = B0m*(yr/rqub1 - yl/rqub2);
    Bx = B0m*new_x*dif12;
    Bz = B0m*new_z*dif12;

    Bz = Bz + f*B0s;

    /* Rotate the BField Back */
    rotation_matrix_1[0][1] = -rotation_matrix_1[0][1];
    rotation_matrix_1[1][0] = -rotation_matrix_1[1][0];

    rotation_matrix_2[0][1] = -rotation_matrix_2[0][1];
    rotation_matrix_2[1][0] = -rotation_matrix_2[1][0];

    rotation_matrix_3[0][1] = -rotation_matrix_3[0][1];
    rotation_matrix_3[1][0] = -rotation_matrix_3[1][0];

    
    hold1[0] = Bx;
    hold1[1] = By;

    mvmul(2, rotation_matrix_3, hold1, hold2);

    Bx = hold2[0];
    By = hold2[1];
    


    hold1[0] = Bz;
    hold1[1] = Bx;

    mvmul(2, rotation_matrix_2, hold1, hold2); 
  
    Bz = hold2[0];

    hold1[0] = hold2[1];
    hold1[1] = By;

    mvmul(2, rotation_matrix_1, hold1, hold2);
    Bx = hold2[0];
    By = hold2[1];



    B[1] = B[1] + By;
    B[0] = B[0] + Bx;
    B[2] = B[2] + Bz;



    oCurrentStreamer = oCurrentStreamer->oNextStreamer;



  }
}

