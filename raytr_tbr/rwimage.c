#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "raytrace.h"

//
// beam_intens_traj()
//
//   Calculation of the sun image in the radio frequencies 10 - 127 MHz
//
//   This subprogram accepts the observer's position in the solar system in SGI
// coordinates in solar radii, radio wave frequency in Hertz, coordinates of 
// the image plane's lower left and upper right corners like (x0,y0), (x1,y1), 
// radius of integration sphere, dimensions of the image in pixels and returns 
// the radiotelescope image of the sun (or any object located inside of the 
// integration sphere) in the units of intensity.    
//   The beam_intens_traj() is an interface to the beam_path() subprogram. 
// The beam_path() makes wave raytracing and emissivity integration for each 
// ray intensity, but it works with just linear arrays of ray parameters. It is
// the beam_intens_traj() that gathers the integration data in the form of 
// 2-dimensional image.  
//   
// Selected ray trajectores are stored in Traj_DII[3][nTraj][lenTrajMax].
//
//   Inputs:
//
// nTrajMax: number of array dimensions for trajectories to store
// nTraj: actual number of trajectories to store
// lenTrajMax: maximum length of the stored trajectories
// SelTraj_II[2][nTraj]: (y,x)-indices for each trajectory to be stored. These 
//   numbers are [i][j]-indices of the pixels in the image 
//   Intens_II[nYPixel][nXPixel].
//
//   Output:
//
// Traj_DII[3][nTraj][lenTrajMax]: 3D coordinates for each of nTraj 
//   trajectories, lenTrajMax points maximum for each trajectory. 
// LenTraj_I[nTraj]: Length (number of points) for each trajectory. Since some 
//   rays are longer than the other ones, the lengths are required to know 
//   which values in the returned Traj_DII are valid. Say, for some trajectory 
//   iTraj in [0:nTraj-1], valid are Traj_DII[0:2][iTraj][0:lenTrajMax-1].
//   In the course of the ray calculation LenTraj_I[] contains current 
//   lengths-1, so that for iTraj ray LenTraj_I[iTraj] is the index k into 
//   Traj_DII[:][iTraj][k].
//
// Author:
//   Leonid Benkevitch
// Revisions:
// 2007-Jul-17 Initially written in Fortran-90 and tested for frequencies 
//             around 42.3 MHz
// 2009-Jan-05 Rewritten in GNU C
// 2009-Jan-15 Added saving of selected ray trajectories for subsequent 
//             plotting; renamed from beam_intensity() into beam_intens_traj()
// 2009-Jul-01 For calculating the brightness temperature in the Intens_I[],
//             the parameters RadioFrequency and OpDepth_I[nRay] are added 
//

void beam_intens_traj(double rEarth_D[3], double Freq, 
		      double ImageRange_I[4], double rISph, 
		      int nXPixel, int nYPixel, 
		      double Intens_II[nYPixel][nXPixel], 
		      int const nTrajMax, int const lenTrajMax, int nTraj, 
		      int SelTraj_II[2][nTrajMax], 
		      int LenTraj_I[nTrajMax], 
		      double Traj_DII[3][nTrajMax][lenTrajMax],
		      double Tb_II[nTrajMax][lenTrajMax])

{
  double extern const cPi;              // = 3.1415926535897931;
  double extern const ProtonChargeSGSe; // = 4.8e-10; // StatCoulombs, SGSe
  double extern const cProtonMass;      // = 1.672621636E-24; //g
  double extern const cElectronMass;    // = 9.10938215E-28; //g
  double extern const cSolarRadius;     // = 6.955E10 cm
  double extern const DeltaS;           // R_sun
  short static const FALSE = 0, TRUE = 1;

  int nRay = nXPixel*nYPixel;
  short Excl_I[nRay]; // A ray is excluded from processing if it is true
  // XPixel, YPixel - Pixel coordinates INSIDE of the image plane
  double XPixel_II[nYPixel][nXPixel], YPixel_II[nYPixel][nXPixel];   
  double Normal_D[3];       // Unity vector normal to image plane
  double static ZAxisOrt_D[3] = {0, 0, 1};
  double Tau_D[3], Xi_D[3]; // Image plane inner Cartesian orts;
  double rEarthLen, SlopeUnscaledLen;
  double XyzPixel_DII[3][nYPixel][nXPixel]; // SGI pixel coordinates
  double XyzPixel_D[3];     // SGI pixel coordinates
  // SGI unity slope vectors for all the line-of-sights pointing at the pixels:
  double Dir_DII[3][nYPixel][nXPixel];    
  // Distance from the radiotelescope to the integration sphere
  double dEarth2ISph_II[nYPixel][nXPixel]; 
  // Squared distance from the radiotelescope to the integration sphere
  double EarthToIntSphere2_II[nYPixel][nXPixel];  
  double Pos_DII[3][nYPixel][nXPixel]; 
  double XPosition_II[nYPixel][nXPixel], YPosition_II[nYPixel][nXPixel];
  double ZPosition_II[nYPixel][nXPixel], SolarDistSqr_II[nYPixel][nXPixel];
  double SolarDistSqr_I[nYPixel*nXPixel];
  double XPixelDel, YPixelDel;    
  double SlopeUnscaled_D[3];
  double Dir_D[3], DirXrEarth_D[3], Pos_D[3];
  double DiscrmSqr, SlopeDotRearth;
  double XLower, XUpper, YLower, YUpper;
  double Pos_DI[3][nRay], Dir_DI[3][nRay];
  double Intens_I[nYPixel*nXPixel], DS_I[nYPixel*nXPixel]; 
  double RayPath_I[nYPixel*nXPixel];
  short Flags_I[nYPixel*nXPixel]; 
  short NewEntry; // Must be set to .true. before a call with new value of nRay
  double OneAU = 215.0, Tol = 0.01;
  double MaxRayPath = 60.;
  double XPixel, YPixel, SolarDistMin, MinRayPath, PercentRayLeft;
  double rISphSqr;
  double DensityCr;
  int nIter, i, j, k, iRay, nRayInISph, iTraj;
  short RayInISph_I[nRay];
  double minSolarDistance;
  short deb = FALSE;  //TRUE;
  short ActiveRaysExist;
  double OpDepth_I[nYPixel*nXPixel];
  FILE *fh; // Text file to save some parameters during beam_path() calls
  // Bits for Flags_I
  short const static Penetr = 0x0001;  // Eps < 0; Critical surface search
  short const static WasRefl = 0x0002; // First step after reflection

  //
  // Calculate the critical density from the frequency
  //
  DensityCr = cPi*cProtonMass*cElectronMass*pow(Freq/ProtonChargeSGSe,2);

  if (deb) {
    printf("\nDensityCr = %g\n", DensityCr);
  }

  //
  // Determine the image plane inner coordinates of pixel centers
  //
  XLower = ImageRange_I[0];
  YLower = ImageRange_I[1];
  XUpper = ImageRange_I[2];
  YUpper = ImageRange_I[3];
  XPixelDel = (XUpper - XLower)/nXPixel;
  YPixelDel = (YUpper - YLower)/nYPixel;

  for (j = 0; j < nXPixel; j++) {
    XPixel_II[0][j] = XLower + ((double)j+0.5)*XPixelDel;
  }
  for (i = 1; i < nYPixel; i++)
    for (j = 0; j < nXPixel; j++) {
      XPixel_II[i][j] = XPixel_II[0][j];
    }

  for (i = 0; i < nYPixel; i++) {
    YPixel_II[i][0] = YLower + ((double)i+0.5)*YPixelDel;
  }
  for (j = 1; j < nXPixel; j++) 
    for (i = 0; i < nYPixel; i++) {
      YPixel_II[i][j] = YPixel_II[i][0];
    }

  if (deb) {
    printf("\nXPixelDel = %g, YPixelDel = %g\n", XPixelDel, YPixelDel);
    printf("\nXPixel:\n");
    print2d(nYPixel, nXPixel, XPixel_II);
    printf("\nYPixel:\n");
    print2d(nYPixel, nXPixel, YPixel_II);
  }

  //
  // Determune the orts, Tau and Xi, of the inner coordinate
  //  system of the image plane
  //
  rEarthLen = sqrt(dot_product(rEarth_D,rEarth_D));
  Normal_D[0] = rEarth_D[0]/rEarthLen;
  Normal_D[1] = rEarth_D[1]/rEarthLen;
  Normal_D[2] = rEarth_D[2]/rEarthLen;
  cross_product(ZAxisOrt_D, Normal_D, Tau_D); // Tau = ZAxisOrt x Normal
  cross_product(Normal_D, Tau_D, Xi_D);       // Xi = Normal x Tau

  //if (deb) {
    printf("\nNormal_D = ");
    print1d(3, Normal_D);
    printf("\nTau_D = ");
    print1d(3, Tau_D);
    printf("\nXi_D = ");
    print1d(3, Xi_D);
    //}

  //
  // Calculate coordinates of all the pixels in the SGI
  //
  for (i = 0; i < nYPixel; i++)
    for (j = 0; j < nXPixel; j++) {
      XyzPixel_D[0] = XPixel_II[i][j]*Tau_D[0] + YPixel_II[i][j]*Xi_D[0];
      XyzPixel_D[1] = XPixel_II[i][j]*Tau_D[1] + YPixel_II[i][j]*Xi_D[1];
      XyzPixel_D[2] = XPixel_II[i][j]*Tau_D[2] + YPixel_II[i][j]*Xi_D[2];
      SlopeUnscaled_D[0] = XyzPixel_D[0] - rEarth_D[0];
      SlopeUnscaled_D[1] = XyzPixel_D[1] - rEarth_D[1];
      SlopeUnscaled_D[2] = XyzPixel_D[2] - rEarth_D[2];
      SlopeUnscaledLen = sqrt(dot_product(SlopeUnscaled_D, SlopeUnscaled_D));
      Dir_DII[0][i][j] = SlopeUnscaled_D[0]/SlopeUnscaledLen;  // v
      Dir_DII[1][i][j] = SlopeUnscaled_D[1]/SlopeUnscaledLen;  // v
      Dir_DII[2][i][j] = SlopeUnscaled_D[2]/SlopeUnscaledLen;  // v
      printf("Dir_D = %g %g %g\n", Dir_DII[0][i][j],  Dir_DII[1][i][j],  
      	     Dir_DII[2][i][j]);
    }

  //
  // Find the points on the integration sphere where it intersects 
  // with the straight "rays" 
  //

  for (i = 0; i < nYPixel; i++)
    for (j = 0; j < nXPixel; j++) {
      Dir_D[0] = Dir_DII[0][i][j];
      Dir_D[1] = Dir_DII[1][i][j];
      Dir_D[2] = Dir_DII[2][i][j];
      cross_product(Dir_D, rEarth_D, DirXrEarth_D);
      SlopeDotRearth = dot_product(Dir_D, rEarth_D);
      DiscrmSqr = sqrt(dot_product(Dir_D,Dir_D)*rISph*rISph 
			   - dot_product(DirXrEarth_D,DirXrEarth_D));
      dEarth2ISph_II[i][j] = -SlopeDotRearth - DiscrmSqr;
      EarthToIntSphere2_II[i][j] =    -SlopeDotRearth + DiscrmSqr;
      Pos_DII[0][i][j] = rEarth_D[0] + Dir_DII[0][i][j]*dEarth2ISph_II[i][j];
      Pos_DII[1][i][j] = rEarth_D[1] + Dir_DII[1][i][j]*dEarth2ISph_II[i][j];
      Pos_DII[2][i][j] = rEarth_D[2] + Dir_DII[2][i][j]*dEarth2ISph_II[i][j];
      printf("Pos_D = %g %g %g\n", Pos_DII[0][i][j],  Pos_DII[1][i][j],  
      	     Pos_DII[2][i][j]);
    }

  if (1) {
    printf("dEarth2ISph_II = \n");
    print2d(nYPixel, nXPixel, dEarth2ISph_II);
  }

  //
  // Do emissivity integration inside of the integration sphere 
  //
  for (i = 0; i < nYPixel; i++)
    for (j = 0; j < nXPixel; j++) {
      k = i*nXPixel + j;
      Pos_DI[0][k] = Pos_DII[0][i][j];
      Pos_DI[1][k] = Pos_DII[1][i][j];
      Pos_DI[2][k] = Pos_DII[2][i][j];
      if (deb) {
	printf("Pos = %g %g %g\n", Pos_DI[0][k], Pos_DI[1][k], Pos_DI[2][k]);
      }
      Dir_DI[0][k] = Dir_DII[0][i][j];
      Dir_DI[1][k] = Dir_DII[1][i][j];
      Dir_DI[2][k] = Dir_DII[2][i][j];
    }

  if (1) {
    printf("Dir_DI = \n");
    print2d(3, nRay, Dir_DI);
  }

  //
  // Store the initial trajectories' coordinates in Traj_DII
  //
  // int nTraj, int lenTrajMax, int SelTraj_II[2][nTraj], 
  // int LenTraj_I[nTraj], double Traj_DII[3][nTraj][lenTrajMax]
  //
  printf("\nRay indices and starting positions\n");
  for (iTraj = 0; iTraj < nTraj; iTraj++) {
    i = SelTraj_II[0][iTraj];
    j = SelTraj_II[1][iTraj];
    Traj_DII[0][iTraj][0] = Pos_DII[0][i][j];
    Traj_DII[1][iTraj][0] = Pos_DII[1][i][j];
    Traj_DII[2][iTraj][0] = Pos_DII[2][i][j];
    Tb_II[iTraj][0] = 0.0;   // Brightness temperature at the initial point
    LenTraj_I[iTraj] = 1;
    k = nXPixel*i + j;  // ray index into the linear array Pos_DI[:][k] DEBUG!
    printf("i = %3i, j = %3i, k = %6i, Traj: %g %g %g\n", 
	   i, j, k, Traj_DII[0][iTraj][0], Traj_DII[1][iTraj][0], 
	   Traj_DII[2][iTraj][0]);

  }
  
  for (i = 0; i < nRay; i++) {
    Intens_I[i] = 0.0;
    Flags_I[i] = 0x0000;
    Excl_I[i] = FALSE;
    RayInISph_I[i] = 1;
    DS_I[i] = DeltaS;
    RayPath_I[i] = 0.0;
    OpDepth_I[i] = 0.0;
  }
  Tol = 0.005;
  nIter = 0;
  rISphSqr = rISph*rISph + 0.01;
  NewEntry = 1;
  MinRayPath = 0.0;
  nRayInISph = nRay;
  ActiveRaysExist = TRUE;
  if (deb) {
    printf("rISphSqr = %g\n",rISphSqr );
    printf("NewEntry = %i\n", NewEntry);
  }

  printf("\n\nNe_CR = %g\n\n", DensityCr/cProtonMass); 


  fh = fopen("beam.txt", "w");


  while (minSolarDistance < rISph) {
    //int ii, i1;
    //for (ii = 0; ii < 2000; ii++) {
    nIter = nIter + 1;
    if (deb) { printf("nIter = %i ============================\n", nIter); }
    minSolarDistance = 1e16;  // Just a very big value
    for (i = 0; i < nRay; i++) {
      double r2;
      Pos_D[0] = Pos_DI[0][i];
      Pos_D[1] = Pos_DI[1][i];
      Pos_D[2] = Pos_DI[2][i];
      r2 = SolarDistSqr_I[i] = dot_product(Pos_D, Pos_D);  //sum(Pos_DI**2,1);
      if (SolarDistSqr_I[i] > rISphSqr) { 
	RayInISph_I[i] = 0;
	Excl_I[i] = TRUE;
      } // if (SolarDistSqr_I[i] > rISphSqr)
      if ((minSolarDistance > r2) && (!(Flags_I[i]&Penetr)))
	minSolarDistance = r2;
    }
    minSolarDistance = sqrt(minSolarDistance);

    if (deb) {
      printf("SolarDistSqr_I = "); print1d(nRay, SolarDistSqr_I); 
      printf("SolarDistSqr_I[0] = %g\n", SolarDistSqr_I[0]);
      printf("nRayInISph = %i\n", nRayInISph);
    }

    if (nIter % (int)(10.0/DeltaS) == 0) {
    
      PercentRayLeft = 100.0*(double)nRayInISph/(double)nRay;
      //minSolarDistance = minval(SolarDistSqr_I, nRay);
      printf("%5i: =====(%5i) %8.4f\% ====== Min Dist from Sun center: %g\n", 
	     nIter, nRayInISph, PercentRayLeft, minSolarDistance);
    }
    



    beam_path(plasma_density, nRay, Excl_I, Pos_DI, Dir_DI, DS_I,
              Tol, DensityCr, Intens_I, Flags_I, &NewEntry,
	      Freq, OpDepth_I, fh);



    //
    // Store calculated ray positions  and brightness temperatures 
    // for the selected rays
    //
    // double Tb_II[nTrajMax][lenTrajMax]
    // double Traj_DII[3][nTraj][lenTrajMax]) 
    //
    for (iTraj = 0; iTraj < nTraj; iTraj++) {
      int p;
      i = SelTraj_II[0][iTraj];
      j = SelTraj_II[1][iTraj];
      k = nXPixel*i + j; // ray index into the linear array Pos_DI[:][k]
      if ((!Excl_I[k]) && (LenTraj_I[iTraj] < lenTrajMax)) {
	// points where to put the calculated coordinates into Traj_DII
	p = LenTraj_I[iTraj];    
	Traj_DII[0][iTraj][p] = Pos_DI[0][k];
	Traj_DII[1][iTraj][p] = Pos_DI[1][k];
	Traj_DII[2][iTraj][p] = Pos_DI[2][k];
	Tb_II[iTraj][p] = Intens_I[k];
	LenTraj_I[iTraj] += 1;
      }
    }
    
    MinRayPath = RayPath_I[0];
    for (i = 0; i < nRay; i++) {
      RayPath_I[i] = RayPath_I[i] + DS_I[i];
      if (MinRayPath > RayPath_I[i]) MinRayPath = RayPath_I[i];
    }
    //
    // Get sum of all RayInISph_I[] elements and 
    // find if there are active rays 
    //
    nRayInISph = 0; 
    ActiveRaysExist = TRUE;
    for(i = 0; i < nRay; i++) {
      nRayInISph += RayInISph_I[i];
      //ActiveRaysExist = ActiveRaysExist && (Excl_I[i] || Flags_I[i]);
      ActiveRaysExist = ActiveRaysExist && (Excl_I[i]);
      //printf("%i ", (int)Excl_I[i]);
    } 
    ActiveRaysExist = !ActiveRaysExist;
    //printf("ActiveRaysExist = %i;\n", ActiveRaysExist);
  } // while (ActiveRaysExist > 0)

  for (iRay = 0; iRay < nRay; iRay++) {
    if (Flags_I[iRay]&Penetr) printf("+++ Snell ray %i", iRay);
  }


  fclose(fh);


  //
  // Free internal memory
  //
  NewEntry = -1;
  //beam_path(plasma_density, nRay, Excl_I, Pos_DI, Dir_DI, DS_I,
  //	    Tol, DensityCr, Intens_I, Flags_I, &NewEntry);
  beam_path(plasma_density, nRay, Excl_I, Pos_DI, Dir_DI, DS_I,
	    Tol, DensityCr, Intens_I, Flags_I, &NewEntry,
	    Freq, OpDepth_I, fh);


  //
  // Restore 2-Dim image from linear array
  //
  for (i = 0; i < nYPixel; i++)
    for (j = 0; j < nXPixel; j++) {
      k = nXPixel*i + j;
      Intens_II[i][j] = Intens_I[k];     // *exp(-OpDepth_I[k]);
    }

} // void beam_intens_traj()


