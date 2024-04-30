#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "raytrace.h"

//
// beam_intens_traj()
//
//   Calculation of the sun image in the radio frequencies 10 - 127 MHz
//
//   This subprogram accepts the observer's position in the solar system in SGI coordinates
// in solar radii, radio wave frequency in Hertz, coordinates of the image plane's lower left
// and upper right corners like (x0,y0), (x1,y1), radius of integration sphere, dimensions
// of the image in pixels and returns the radiotelescope image of the sun (or any object
// located inside of the integration sphere) in the units of intensity.    
//   The beam_intens_traj() is an interface to the beam_path() subprogram. The beam_path() makes
// wave raytracing and emissivity integration for each ray intensity, but it works with just
// linear arrays of ray parameters. It is the beam_intens_traj() that gathers the integration
// data in the form of 2-dimensional image.  
//   
// Selected ray trajectores are stored in Traj_DII[3][nTraj][lenTrajMax].
//   Inputs:
// nTrajMax: number of array dimensions for trajectories to store
// nTraj: actual number of trajectories to store
// lenTrajMax: maximum length of the stored trajectories
// SelTraj_II[2][nTraj]: (y,x)-indices for each trajectory to be stored. These numbers are 
//   [i][j]-indices of the pixels in the image Intensity_II[nYPixel][nXPixel]. 
//   Output:
// Traj_DII[3][nTraj][lenTrajMax]: 3D coordinates for each of nTraj trajectories, lenTrajMax 
//   points maximum for each trajectory. 
// LenTraj_I[nTraj]: Length (number of points) for each trajectory. Since some rays 
//   are longer than the other ones, the lengths are required to know which values in 
//   the returned Traj_DII are valid. Say, for some trajectory iTraj in [0:nTraj-1],
//   valid are Traj_DII[0:2][iTraj][0:lenTrajMax-1].
//   In the course of the ray calculation LenTraj_I[] contains current lengths-1,
//   so that for iTraj ray LenTraj_I[iTraj] is the index k into Traj_DII[:][iTraj][k].
//
// Author:
//   Leonid Benkevitch
// Revisions:
// 2007-Jul-17 Initially written in Fortran-90 and tested for frequencies around 42.3 MHz
// 2009-Jan-05 Rewritten in GNU C
// 2009-Jan-15 Added saving of selected ray trajectories for subsequent plotting;
//             renamed from beam_intensity() into beam_intens_traj()
//

void beam_intens_traj(double XyzEarth_D[3], double RadioFrequency, double ImageRange_I[4], 
		      double rIntegration, int nXPixel, int nYPixel, double Intensity_II[nYPixel][nXPixel], 
		      int const nTrajMax, int const lenTrajMax, int nTraj, int SelTraj_II[2][nTrajMax], 
		      int LenTraj_I[nTrajMax], double Traj_DII[3][nTrajMax][lenTrajMax])

{
  double extern const cPi;                         // = 3.1415926535897931;
  double extern const ProtonChargeSGSe;            // = 4.8e-10; // StatCoulombs, SGSe
  double extern const cProtonMass;                 // = 1.672621636E-24; //g
  double extern const cElectronMass;               // = 9.10938215E-28; //g
  short static const FALSE = 0, TRUE = 1;

  int nRay = nXPixel*nYPixel;
  //short ExcludeRay_I[nYPixel*nXPixel];   // A ray is excluded from processing if it is .true.
  short ExcludeRay_I[nRay];   // A ray is excluded from processing if it is .true.
  double XPixel_II[nYPixel][nXPixel], YPixel_II[nYPixel][nXPixel];   // Pixel coordinates INSIDE of the image plane
  double Normal_D[3];                                          // Unity vector normal to image plane
  double static ZAxisOrt_D[3] = {0, 0, 1};
  double Tau_D[3], Xi_D[3];                                   // Image plane inner Cartesian orts;
  double XyzEarthLen, SlopeUnscaledLen;
  double XyzPixel_DII[3][nYPixel][nXPixel];         // SGI pixel coordinates
  double XyzPixel_D[3];         // SGI pixel coordinates
  double Slope_DII[3][nYPixel][nXPixel];    // SGI unity slope vectors for all the line-of-sights pointing at the pixels 
  double EarthToIntSphereDist_II[nYPixel][nXPixel];  // Distance from the radiotelescope to the integration sphere
  double EarthToIntSphere2_II[nYPixel][nXPixel];  // Squared distance from the radiotelescope to the integration sphere
  double Position_DII[3][nYPixel][nXPixel]; 
  double XPosition_II[nYPixel][nXPixel], YPosition_II[nYPixel][nXPixel];
  double ZPosition_II[nYPixel][nXPixel], SolarDistSqr_II[nYPixel][nXPixel];
  double SolarDistSqr_I[nYPixel*nXPixel];
  double XPixelDel, YPixelDel;    
  double SlopeUnscaled_D[3];
  double Slope_D[3], SlopeCrossREarth_D[3], Pos_D[3];
  double DiscrmSquared, SlopeDotRearth;
  double XLower, XUpper, YLower, YUpper;
  double Position_DI[3][nRay], Slope_DI[3][nRay];
  double Intensity_I[nYPixel*nXPixel], DeltaS_I[nYPixel*nXPixel], RayPath_I[nYPixel*nXPixel];
  short RayFlag_I[nYPixel*nXPixel];         // 1 if a ray is OK, 0 otherwise
  short NewEntry;           // Must be set to .true. before a call with new value of nRay
  double OneAU = 215.0, Tolerance = 0.01, DeltaS = 1.0;
  double MaxRayPath = 60.;
  double XPixel, YPixel, SolarDistMin, MinRayPath, PercentRayLeft, rIntegrationSqr;
  double DensityCr;
  int nIter, i, j, k, iRay, nRayInsideIntSphere, iTraj;
  short RayInsideIntSphere_I[nRay];
  double minSolarDistance;
  //double Pos_D[3];
  short deb = FALSE;  //TRUE;
  short ActiveRaysExist;
  //
  // Calculate the critical density from the frequency
  //
  DensityCr = cPi*cProtonMass*cElectronMass*pow(RadioFrequency/ProtonChargeSGSe,2);

  if (deb) {
    printf("\nDensityCr = %g\n", DensityCr);
  }
  //write(*,*) 'DensityCr = ', DensityCr
  //
  // Determine the image plane inner coordinates of pixel centers
  //
  XLower = ImageRange_I[0];
  YLower = ImageRange_I[1];
  XUpper = ImageRange_I[2];
  YUpper = ImageRange_I[3];
  XPixelDel = (XUpper - XLower)/nXPixel;
  YPixelDel = (YUpper - YLower)/nYPixel;

  // XPixel_II(1,:) = (/ (XLower + (real(j)-0.5)*XPixelDel, j = 1, nXPixel) /)
  //do i = 2, nXPixel
  //   XPixel_II(i,:) = XPixel_II(1,:)
  //end do
  for (j = 0; j < nXPixel; j++) {
    XPixel_II[0][j] = XLower + ((double)j+0.5)*XPixelDel;
  }
  for (i = 1; i < nYPixel; i++)
    for (j = 0; j < nXPixel; j++) {
      XPixel_II[i][j] = XPixel_II[0][j];
    }

  //YPixel_II(:,1) = (/ (YLower + (real(j)-0.5)*YPixelDel, j = 1, nYPixel) /)
  //do i = 2, nYPixel
  //   YPixel_II(:,i) = YPixel_II(:,1)
  //end do
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
  // Determune the orts, Tau and Xi, of the inner coordinate system of the image plane
  //
  XyzEarthLen = sqrt(dot_product(XyzEarth_D,XyzEarth_D));
  Normal_D[0] = XyzEarth_D[0]/XyzEarthLen;
  Normal_D[1] = XyzEarth_D[1]/XyzEarthLen;
  Normal_D[2] = XyzEarth_D[2]/XyzEarthLen;
  cross_product(ZAxisOrt_D, Normal_D, Tau_D);    // Tau = ZAxisOrt x Normal
  cross_product(Normal_D, Tau_D, Xi_D);          // Xi = Normal x Tau

  //if (deb) {
    printf("\nNormal_D = ");
    print1d(3, Normal_D);
    printf("\nTau_D = ");
    print1d(3, Tau_D);
    printf("\nXi_D = ");
    print1d(3, Xi_D);
    //}
  //write(*,*)
  //write(*,'(a,3f8.5)') 'Normal_D = ', Normal_D
  //write(*,'(a,3f8.5)') 'Tau_D = ', Tau_D
  //write(*,'(a,3f8.5)') 'Xi_D = ', Xi_D
  //write(*,*)

  //
  // Calculate coordinates of all the pixels in the SGI
  //
  //do i = 1, nYPixel
  //   do j = 1, nXPixel
  //      XyzPixel_DII(:,i,j) = XPixel_II(i,j)*Tau_D + YPixel_II(i,j)*Xi_D
  //      SlopeUnscaled_D = XyzPixel_DII(:,i,j) - XyzEarth_D
  //      Slope_DII(:,i,j) = SlopeUnscaled_D/sqrt(sum(SlopeUnscaled_D**2))             // v
  //   end do
  //end do
  for (i = 0; i < nYPixel; i++)
    for (j = 0; j < nXPixel; j++) {
      //XyzPixel_DII[0][i][j] = XPixel_II[i][j]*Tau_D[0] + YPixel_II[i][j]*Xi_D[0];
      //XyzPixel_DII[1][i][j] = XPixel_II[i][j]*Tau_D[1] + YPixel_II[i][j]*Xi_D[1];
      //XyzPixel_DII[2][i][j] = XPixel_II[i][j]*Tau_D[2] + YPixel_II[i][j]*Xi_D[2];
      XyzPixel_D[0] = XPixel_II[i][j]*Tau_D[0] + YPixel_II[i][j]*Xi_D[0];
      XyzPixel_D[1] = XPixel_II[i][j]*Tau_D[1] + YPixel_II[i][j]*Xi_D[1];
      XyzPixel_D[2] = XPixel_II[i][j]*Tau_D[2] + YPixel_II[i][j]*Xi_D[2];
      SlopeUnscaled_D[0] = XyzPixel_D[0] - XyzEarth_D[0];
      SlopeUnscaled_D[1] = XyzPixel_D[1] - XyzEarth_D[1];
      SlopeUnscaled_D[2] = XyzPixel_D[2] - XyzEarth_D[2];
      SlopeUnscaledLen = sqrt(dot_product(SlopeUnscaled_D, SlopeUnscaled_D));
      Slope_DII[0][i][j] = SlopeUnscaled_D[0]/SlopeUnscaledLen;             // v
      Slope_DII[1][i][j] = SlopeUnscaled_D[1]/SlopeUnscaledLen;             // v
      Slope_DII[2][i][j] = SlopeUnscaled_D[2]/SlopeUnscaledLen;             // v
      //if (deb) {
      //printf("XyzPixel_D = "); print1d(3,XyzPixel_D);
      //printf("Slope_DII = %g %g %g\n", Slope_DII[0][i][j], Slope_DII[1][i][j], Slope_DII[2][i][j]);
	//}
    }

  //
  // Find the points on the integration sphere where it intersects with the straight "rays" 
  //
  //do i = 1, nYPixel
  //   do j = 1, nXPixel
  //      EarthToIntSphereDist_II(i,j) = -sum(XyzEarth_D*Slope_DII(:,i,j)) - &
  //           sqrt(sum((Slope_DII(:,i,j)*rIntegration)**2) - sum(cross_product(Slope_DII(:,i,j), XyzEarth_D)**2))
  //      EarthToIntSphere2_II(i,j) = -sum(XyzEarth_D*Slope_DII(:,i,j)) + &
  //           sqrt(sum((Slope_DII(:,i,j)*rIntegration)**2) - sum(cross_product(Slope_DII(:,i,j), XyzEarth_D)**2))
  //      Position_DII(:,i,j) = XyzEarth_D + Slope_DII(:,i,j)*EarthToIntSphereDist_II(i,j)
  //   end do
  //end do

  for (i = 0; i < nYPixel; i++)
    for (j = 0; j < nXPixel; j++) {
      Slope_D[0] = Slope_DII[0][i][j];
      Slope_D[1] = Slope_DII[1][i][j];
      Slope_D[2] = Slope_DII[2][i][j];
      cross_product(Slope_D, XyzEarth_D, SlopeCrossREarth_D);
      SlopeDotRearth = dot_product(Slope_D, XyzEarth_D);
      DiscrmSquared = sqrt(dot_product(Slope_D,Slope_D)*rIntegration*rIntegration 
			   - dot_product(SlopeCrossREarth_D,SlopeCrossREarth_D));
      EarthToIntSphereDist_II[i][j] = -SlopeDotRearth - DiscrmSquared;
      EarthToIntSphere2_II[i][j] =    -SlopeDotRearth + DiscrmSquared;
      Position_DII[0][i][j] = XyzEarth_D[0] + Slope_DII[0][i][j]*EarthToIntSphereDist_II[i][j];
      Position_DII[1][i][j] = XyzEarth_D[1] + Slope_DII[1][i][j]*EarthToIntSphereDist_II[i][j];
      Position_DII[2][i][j] = XyzEarth_D[2] + Slope_DII[2][i][j]*EarthToIntSphereDist_II[i][j];
      if (deb) {
	printf("Pos_D = %g %g %g\n", Position_DII[0][i][j],  Position_DII[1][i][j],  Position_DII[2][i][j]);
      }
    }

  if (deb) {
    printf("EarthToIntSphereDist_II = \n");
    print2d(nYPixel, nXPixel, EarthToIntSphereDist_II);
    //write(*,*)
    //write(*,*) 'EarthToIntSphere2 = '
    //write(*,'(10f7.2)') (EarthToIntSphere2_II(i,:), i = 1, nYPixel)
    //write(*,*)
    //================
    //write(*,*) 'Slope_X = ' 
    //write(*,'(10f10.5)') (Slope_DII(1,i,:), i = 1, nYPixel)
    //write(*,*)
    //write(*,*) 'Slope_Y = ' 
    //write(*,'(10f10.5)') (Slope_DII(2,i,:), i = 1, nYPixel)
    //write(*,*)
    //write(*,*) 'Slope_Z = ' 
    //write(*,'(10f10.5)') (Slope_DII(3,i,:), i = 1, nYPixel)
    //write(*,*)
    //end if

  }

  //
  // Do emissivity integration inside of the integration sphere 
  //
  //nRay = nXPixel*nYPixel;

  //Position_DI = reshape(Position_DII, (/3, nRay/))
  //Slope_DI = reshape(Slope_DII, (/3, nRay/))
  for (i = 0; i < nYPixel; i++)
    for (j = 0; j < nXPixel; j++) {
      k = i*nXPixel + j;
      Position_DI[0][k] = Position_DII[0][i][j];
      Position_DI[1][k] = Position_DII[1][i][j];
      Position_DI[2][k] = Position_DII[2][i][j];
      if (deb) {
	printf("Pos = %g %g %g\n", Position_DI[0][k], Position_DI[1][k], Position_DI[2][k]);
      }
      Slope_DI[0][k] = Slope_DII[0][i][j];
      Slope_DI[1][k] = Slope_DII[1][i][j];
      Slope_DI[2][k] = Slope_DII[2][i][j];
    }

  if (deb) {
    printf("Slope_DI = \n");
    print2d(3, nRay, Slope_DI);
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
    Traj_DII[0][iTraj][0] = Position_DII[0][i][j];
    Traj_DII[1][iTraj][0] = Position_DII[1][i][j];
    Traj_DII[2][iTraj][0] = Position_DII[2][i][j];
    LenTraj_I[iTraj] = 1;
    k = nXPixel*i + j;       // ray index into the linear array Position_DI[:][k] DEBUG!
    printf("i = %3i, j = %3i, k = %6i, Traj: %g %g %g\n", i, j, k, 
	   Traj_DII[0][iTraj][0], Traj_DII[1][iTraj][0], Traj_DII[2][iTraj][0]);

  }
  


  for (i = 0; i < nRay; i++) {
    Intensity_I[i] = 0.0;
    RayFlag_I[i] = FALSE;
    ExcludeRay_I[i] = FALSE;
    RayInsideIntSphere_I[i] = 1;
    DeltaS_I[i] = DeltaS;
    RayPath_I[i] = 0.0;
  }
  Tolerance = 0.01;
  nIter = 0;
  rIntegrationSqr = rIntegration*rIntegration + 0.01;
  NewEntry = 1;
  MinRayPath = 0.0;
  nRayInsideIntSphere = nRay;
  ActiveRaysExist = TRUE;
  if (deb) {
    printf("rIntegrationSqr = %g\n",rIntegrationSqr );
    printf("NewEntry = %i\n", NewEntry);
  }

  
  //write(*,*) 'nRayInsideIntSphere = ', nRayInsideIntSphere
  //printf("iter = %i; DeltaS_I[70] = %g; Position_DI[*][70] = %g %g %g; Slope_DI[*][70] = %g %g %g;\n",
  //	 0 , DeltaS_I[70], 
  //	 Position_DI[0][70], Position_DI[1][70], Position_DI[2][70], Slope_DI[0][70], Slope_DI[1][70], Slope_DI[2][70]);

  //Slope_DI[2][5100] -= 1e-14;

  //while (nRayInsideIntSphere) {
  while (ActiveRaysExist) {
  int ii, i1;
  //for (ii = 0; ii < 2; ii++) {
    nIter = nIter + 1;
    if (deb) { printf("nIter = %i ==================================\n", nIter); }
    for (i = 0; i < nRay; i++) {
      double r2;
      Pos_D[0] = Position_DI[0][i];
      Pos_D[1] = Position_DI[1][i];
      Pos_D[2] = Position_DI[2][i];
      r2 = SolarDistSqr_I[i] = dot_product(Pos_D, Pos_D);  //sum(Position_DI**2,1);
      //if (i == 70) printf("SolarDist_I[70] = %g; rIntegration = %g\n", sqrt(SolarDistSqr_I[i]), rIntegration);
      //if (SolarDistSqr_I[i] < rIntegrationSqr) printf("iRay inside: %i; SolarDist_I[i] = %g;\n", i, sqrt(SolarDistSqr_I[i]));
      if (SolarDistSqr_I[i] > rIntegrationSqr) { 
	RayInsideIntSphere_I[i] = 0;
	ExcludeRay_I[i] = TRUE;
      } // if (SolarDistSqr_I[i] > rIntegrationSqr)
    }
    //printf("ExclRay[] = "); for(i1 = 0; i1 < nRay; i1++) {printf("%i ", ExcludeRay_I[i1]);}; printf("\n");

    if (deb) {
      printf("SolarDistSqr_I = "); print1d(nRay, SolarDistSqr_I); 
      printf("SolarDistSqr_I[0] = %g\n", SolarDistSqr_I[0]);
      printf("nRayInsideIntSphere = %i\n", nRayInsideIntSphere);
    }

    if (nIter % (int)(10.0/DeltaS) == 0) {
    
      PercentRayLeft = 100.0*(double)nRayInsideIntSphere/(double)nRay;
      //write(*,'(i5,a,f8.4,a,f11.6)') nIter,': ======= ', PercentRayLeft, &
      //     ' % ====== Minimum Distance from Sun: ', minval(sqrt(SolarDistSqr_I)) 
      minSolarDistance = minval(SolarDistSqr_I, nRay);
      //printf("minSolarDistance^2 = %g", minSolarDistance);
      minSolarDistance = sqrt(minSolarDistance);
      printf("%5i: ======= %8.4f% ====== Minimum Distance from the center of Sun: %g\n", nIter, PercentRayLeft, minSolarDistance);
      //for (k = 0; k < nRay; k++) {
      //  if (sqrt(SolarDistSqr_I[k]) <= 1.0) printf("i = %i, j = %i, R = %g;\n", k/nXPixel, k%nXPixel, sqrt(SolarDistSqr_I[k]));
      //}
    }
    
    //beam_path(void (* Get_Plasma_Density)(), int nRay, short ExcludeRay_I[nRay], double Position_DI[3][nRay],
    //          double Slope_DI[3][nRay], double DeltaS_I[nRay], double ToleranceInit, double DensityCr, 
    //          double Intensity_I[nRay], short RayFlag_I[nRay], short *NewEntry)
    beam_path(plasma_density, nRay, ExcludeRay_I, Position_DI, Slope_DI, DeltaS_I,
	      Tolerance, DensityCr, Intensity_I, RayFlag_I, &NewEntry);

    //Pos_D[0] = Position_DI[0][70];
    //Pos_D[1] = Position_DI[1][70];
    //Pos_D[2] = Position_DI[2][70];

    //
    // Store calculated ray positions for the selected rays
    //
    // double Traj_DII[3][nTraj][lenTrajMax]) 
    //
    for (iTraj = 0; iTraj < nTraj; iTraj++) {
      int p;
      i = SelTraj_II[0][iTraj];
      j = SelTraj_II[1][iTraj];
      k = nXPixel*i + j;       // ray index into the linear array Position_DI[:][k]
      if ((!ExcludeRay_I[k]) && (LenTraj_I[iTraj] < lenTrajMax)) {
	p = LenTraj_I[iTraj];    // points where to put the calculated coordinates into Traj_DII
	Traj_DII[0][iTraj][p] = Position_DI[0][k];
	Traj_DII[1][iTraj][p] = Position_DI[1][k];
	Traj_DII[2][iTraj][p] = Position_DI[2][k];
	LenTraj_I[iTraj] += 1;
      }
    }
    

    //printf("iter = %i; Excl[70] = %i; DeltaS_I[70] = %g; Position_DI[*][70] = %g %g %g; Slope_DI[*][70] = %g %g %g;\n",
    //	   nIter , ExcludeRay_I[70], DeltaS_I[70], 
    //     Pos_D[0], Pos_D[1], Pos_D[2], Slope_DI[0][70], Slope_DI[1][70], Slope_DI[2][70]);

    MinRayPath = RayPath_I[0];
    for (i = 0; i < nRay; i++) {
      RayPath_I[i] = RayPath_I[i] + DeltaS_I[i];
      if (MinRayPath > RayPath_I[i]) MinRayPath = RayPath_I[i];
    }
    //
    // Get sum of all RayInsideIntSphere_I[] elements and 
    // find if there are active rays 
    //
    nRayInsideIntSphere = 0; 
    ActiveRaysExist = TRUE;
    for(i = 0; i < nRay; i++) {
      nRayInsideIntSphere += RayInsideIntSphere_I[i];
      ActiveRaysExist = ActiveRaysExist && ExcludeRay_I[i];
      //printf("%i ", (int)ExcludeRay_I[i]);
    } 
    ActiveRaysExist = !ActiveRaysExist;
    //printf("ActiveRaysExist = %i;\n", ActiveRaysExist);
  } // while (ActiveRaysExist > 0)

  for (iRay = 0; iRay < nRay; iRay++) {
    if (RayFlag_I[iRay]) printf("+++ Bad ray %i", iRay);
  }

  //
  // Free internal memory
  //
  NewEntry = -1;
  beam_path(plasma_density, nRay, ExcludeRay_I, Position_DI, Slope_DI, DeltaS_I,
	    Tolerance, DensityCr, Intensity_I, RayFlag_I, &NewEntry);


  //
  // Restore 2-Dim image from linear array
  //
  for (i = 0; i < nYPixel; i++)
    for (j = 0; j < nXPixel; j++) {
      k = nXPixel*i + j;
      Intensity_II[i][j] = Intensity_I[k];
    }

} // void beam_intens_traj()


//================================================

void plasma_density_inv_squares(int nRay, double Position_DI[3][nRay], double Density_I[nRay], double GradDensity_DI[3][nRay],
		    double DeltaS_I[nRay], short RayFlag_I[nRay]) {
  //
  // at a specified locations, Position_DI
  //
  // DensityAtSolarSurface = Ne(@SolarSurface)*ProtonMass = 2x10^8(cm^-3)*1.6726x10^-24(g) = 3.3452e-16(g/cm^3)
  //
  double extern const cPi;                         // = 3.1415926535897931;
  double extern const ProtonChargeSGSe;            // = 4.8e-10; // StatCoulombs, SGSe
  double extern const cProtonMass;                 // = 1.672621636E-24; //g
  double extern const cElectronMass;               // = 9.10938215E-28; //g

  double const static OneAU = 215.0; // solar radii
  double const static DensityAtSolarSurface = 3.3452E-16;    // g/cm^3
  double Pos_D[3], SolarDistSqr, SolarDistQuad;
  int iRay, i;
  short deb = 1;

  //printf("DensityAtSolarSurface = %g\n", DensityAtSolarSurface);

  //SolarDistSqr_I = sum(Position_DI**2,1)
  //Density_I = DensityAtSolarSurface/SolarDistSqr_I
  //SolarDistQuad_I = SolarDistSqr_I**2
  //do i = 1, 3
  //  GradDensity_DI(i,:) = -2.*DensityAtSolarSurface*Position_DI(i,:)/SolarDistQuad_I
  //end do

  for (iRay = 0; iRay < nRay; iRay++) {
    Pos_D[0] = Position_DI[0][iRay];
    Pos_D[1] = Position_DI[1][iRay];
    Pos_D[2] = Position_DI[2][iRay];
    SolarDistSqr = dot_product(Pos_D, Pos_D);
    Density_I[iRay] = DensityAtSolarSurface/SolarDistSqr;
    SolarDistQuad = SolarDistSqr*SolarDistSqr;
    GradDensity_DI[0][iRay] = -2.*DensityAtSolarSurface*Pos_D[0]/SolarDistQuad;
    GradDensity_DI[1][iRay] = -2.*DensityAtSolarSurface*Pos_D[1]/SolarDistQuad;
    GradDensity_DI[2][iRay] = -2.*DensityAtSolarSurface*Pos_D[2]/SolarDistQuad;
    DeltaS_I[iRay] = 1.0;       // Revert the step to original
  }
  if (deb) {
    printf("GradDensity_DI = \n"); print2d(3, nRay, GradDensity_DI);
  }

}  // void plasma_density()



//================================================

void plasma_density(int nRay, double Position_DI[3][nRay], double Density_I[nRay], double GradDensity_DI[3][nRay],
		    double DeltaS_I[nRay], short RayFlag_I[nRay]) {
  //
  // Coronal density model by Kuniji Saito
  //
  // 
  //
  double extern const cPi;                         // = 3.1415926535897931;
  double extern const ProtonChargeSGSe;            // = 4.8e-10; // StatCoulombs, SGSe
  double extern const cProtonMass;                 // = 1.672621636E-24; //g
  double extern const cElectronMass;               // = 9.10938215E-28; //g

  double const static DensityAtSolarSurface = 3.3452E-16;    // g/cm^3
  double static Pos_D[3], SolarDistSqr;
  int iRay, i;
  double const static a1 = 3.09E8, a2 = 1.58E8, a3 = 0.0251E6;
  double static t1, t2, t3;
  double static x, y, z, r, rhoCyl2, rhoCyl, rhoCylm1;
  double static rm1, rm2d5, rm3d5, rm6, rm7, rm16, rm17;
  double static cosTh, sinTh, cosPh, sinPh, absCosTh;
  double static Ne, Ner, Neth, Nerho, Nex, Ney, Nez;
  short deb = 0;

  for (iRay = 0; iRay < nRay; iRay++) {
    x = Pos_D[0] = Position_DI[0][iRay];
    y = Pos_D[1] = Position_DI[1][iRay];
    z = Pos_D[2] = Position_DI[2][iRay];
    SolarDistSqr = dot_product(Pos_D, Pos_D);
    r = sqrt(SolarDistSqr);
    rhoCyl2 = pow(x,2) + pow(y,2);
    rhoCyl = sqrt(rhoCyl2);      // Cylindrical coordinate rho = sqrt(x^2 + Y^2)
    rhoCylm1 = 1.0/rhoCyl;
    rm1 = 1/r;                   // r^(-1)
    rm2d5 = pow(r,-2.5);         // r^(-2.5)
    rm3d5 = rm2d5*rm1;           // r^(-3.5)
    rm6 = pow(r,-6);             // r^(-6)
    rm7 = rm6*rm1;               // r^(-7)
    rm16 = pow(r,-16);           // r^(-16)
    rm17 = rm16*rm1;             // r^(-17)
    cosTh = z*rm1;               // cos(theta) = z/r
    sinTh = rhoCyl*rm1;          // sin(theta) = rho/r
    cosPh = x*rhoCylm1;           // cos(phi) = x/rho
    sinPh = y*rhoCylm1;           // sin(phi) = y/rho
    absCosTh = fabs(cosTh);        // abs(cos(theta)) = abs(z/r)

    t1 = 1.0 - 0.5*absCosTh;          // (1 - 0.5cos(theta))
    t2 = 1.0 - 0.95*absCosTh;         // (1 - 0.95cos(theta))
    t3 = 1.0 - sqrt(absCosTh); 
    Ne = a1*rm16*t1 + a2*rm6*t2 + a3*rm2d5*t3;                   // Ne, electron number density
    Ner = -16.0*a1*rm17 - 6.0*a2*rm7*t2 - 2.5*a3*rm3d5*t3;       // Ne derivative by r
    Neth = 0.5*a1*rm16*sinTh + 0.95*a2*rm6*sinTh + 0.5*a3*rm2d5*sinTh/sqrt(absCosTh); //  Ne derivative by theta<90
    if (z < 0.0) Neth = -Neth; //  Ne derivative by theta>90
    Nerho = rm1*Neth*cosTh + Ner*sinTh;  // Ne derivative by rho = sqrt(x^2+y^2)
    Nex = Nerho*cosPh;                   // Ne derivative by x
    Ney = Nerho*sinPh;                   // Ne derivative by y
    Nez = -rm1*Neth*sinTh + Ner*cosTh;   // Ne derivative by z
    Density_I[iRay] = Ne*cProtonMass;
    GradDensity_DI[0][iRay] = Nex*cProtonMass;
    GradDensity_DI[1][iRay] = Ney*cProtonMass;
    GradDensity_DI[2][iRay] = Nez*cProtonMass;
    DeltaS_I[iRay] = 1.0;       // Revert the step to original
    //printf("Density = %g, GradDensity_DI = %g, %g, %g\n", Density_I[iRay], 
    //   GradDensity_DI[0][iRay], GradDensity_DI[1][iRay], GradDensity_DI[2][iRay]);
  }
  if (deb) {
    int i;
    i = 11*5 + 5;
    printf("Density = %g, GradDensity_DI = %g, %g, %g\n", Density_I[i], 
	   GradDensity_DI[0][i], GradDensity_DI[1][i], GradDensity_DI[2][i]);
    printf("cos(th) = %g,\n", cosTh); 
  }

}  // void plasma_density()

