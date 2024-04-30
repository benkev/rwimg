#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "raytrace.h"


//
// Calculates plasma density, Rho_I, and its gradient, 
// GradRho_DI(3,nRay), at specified locations Pos_DI(3,nRay)
// Also, it provides appropriate step, DS_New_I, conceivably dependent
// on the numeric grid size
//
void Get_Plasma_Density(int nRay, double Pos_DI[3][nRay], 
			double Rho_I[nRay], double GradRho_DI[3][nRay],
			double DS_New_I[nRay], short Flags_I[nRay]);

void beam_path(void (* Get_Plasma_Density)(), int nRay, 
	       short ExclRay_I[nRay], double Pos_DI[3][nRay],
	       double Dir_DI[3][nRay], double DS_I[nRay], 
	       double TolInit, double RhoCr, 
	       double Intens_I[nRay], short Flags_I[nRay], 
	       short *NewEntry, 
	       double Freq, double OpDepth_I[nRay], FILE *fh) 
{
  //
  //   The subroutine beam_path() makes raytracing and emissivity integration 
  // along ray paths.
  // It works on a group of rays (a beam). Each ray is specified by its 
  // Cartesian Pos_DI and its direction cosines Dir_DI. The subroutine 
  // calculates new Pos_DI, which is DS_I away from the initial 
  // position. It calculates the intensity along the step as emissivity by 
  // DS and adds this to the Intens_I. Thus, after series of calls
  // to beam_path(), the Intens_I contains the result of integration along 
  // the paths of every ray. 
  //   The electromagnetic rays are refracted in plasma due to its non-zero 
  // refractive index, which is the square root of the dielectric permittivity
  // \epsilon. The \epsilon, in turn, is a function of the plasma density. The 
  // external subprogram, Get_Plasma_Density(), provides the plasma Rho_I 
  // along with its gradient, GradRho_DI at the Pos_DI.
  // The \epsilon can be calculated as 1 - Rho_I/RhoCr, where RhoCr
  // is the "critical" plasma density at which the dielectric permittivity
  // \epsilon falls down to zero. The value of RhoCr is proportional to the
  // square of the wave frequency. For example, for the 42.3 MHz radiowaves the
  // critical density is ~3.71x10^(-17) g/cm^3. 
  // A surface where the plasma density achieves the critical value acts like a
  // mirror. No radiation penetrates the critical surface. The radiowaves can 
  // only travel in the regions with lower density.
  //   The emissivity w of plasma at a selected plasma frequency is calculated 
  // as a polynomial function 
  //               w = (Rho2RhoCr)^2*(0.5 - Rho2RhoCr)^2, 
  // where Rho2RhoCr is the quotient Rho_I/RhoCr. So, here the plasma
  // density is used for calculation of both emissivity and dielectric
  // permittyvity.
  //   The parameters of the beam_path() are briefly described below.
  //
  // Get_Plasma_Density():  external subroutine that returns the plasma density
  //   Rho_I and its gradient GradRho_DI. It also provides the recommended step 
  //   size DS_New_I and asserts the Flags_I (assigns TRUE) for a ray should it
  //   (illegally) penetrate into the region with "negative" dielectric 
  //   permittivity.
  // nRay:           number of rays being processed.
  // ExclRay_I:  the caller program can use this logical array to stop 
  //   processing of any individual ray when it already finished travelling 
  //   inside of a specified part of space. Set the corresponding element of 
  //   ExclRay_I to TRUE to leave it unprocessed during the subsequent calls to
  //   beam_path(). Before the first call to beam_path() all the elements of 
  //   ExclRay_I should be set to FALSE; then all the rays will be processed.
  // Pos_DI:    Cartesian position vectors of the rays.
  // Dir_DI:    Direction cosines of the rays.
  // DS_I:      Current step values for the rays. Set the elements of Dir_DI to
  //   some reasonable value, say, 1.0. The DS_I elements are individually 
  //   modified by beam_path() to satisfy the precision requirements set by 
  //   the TolInit.
  // TolInit:   determines the precision of ray paths calculation. TolInit is 
  //   the inverse of the minimum number of ray trajectory points per one 
  //   radian of the ray curvature. If this requirement is not satisfied, the 
  //   corresponding element of DS_I is decreased. Do not set the TolInit to 
  //   any value greater than 0.1 (which means 10 points per curvature radian): 
  //   it will be internally set to 0.1 anyway.   
  // RhoCr: the plasma density at which its dielectric permittivity becomes 
  //   zero for chosen wave frequency.
  // Intens_I:    the intensities of each ray calculated as the integral of 
  //   emissivity along the ray path. During each call to ray_path(), each 
  //   element of the Intens_I is incremented by the integral along the step 
  //   DS_I. Set all the elements of Intens_I to 0.0 before the first call to 
  //   the beam_path().
  // Flags_I:      the TRUE elements of this logical array indicate that the 
  //   corresponding ray penetrated the "prohibited" region of space with the 
  //   plasma density above its critical value. Normally, it should never 
  //   happen. However, in case the algorithm made such an error, the flagged 
  //   rays should be considered as "bad" and thrown away from the resultant 
  //   Intens_I. Set all the elements of Flags_I to FALSE before calling 
  //   ray_path() for the first time.
  // NewEntry:   A pointer to a short integer variable that controls the 
  //   internal array allocation and deallocation. Set this variable to 1 
  //   before the first call to beam_path(). This value forces the beam_path() 
  //   to allocate internal dynamic arrays and take several initial actions. 
  //   During subsequent calls to beam_path() the *NewEntry will keep the 
  //   value 0, which leaves the allocated arrays intact. Setting the *NewEntry
  //   to -1 during the same run will free the array memory, after which 
  //   another allocation of the internal allocatables, with possibly different
  //   dimensions, can be made via setting the *NewEntry to 1 again. This 
  //   variable was introduced out of the speed considerations: the automatic
  //   arrays work slower than the allocatable ones, because the allocation 
  //   occurs only once at the first call.

  // 
  // Author:
  //   Leonid Benkevitch
  // Revisions:
  // 2007-Jul-17 Initially written in Fortran-90 and tested for frequencies 
  //             around 42.3 MHz
  // 2009-Jan-05 Rewritten in GNU C 
  // 2009-Jan-09 Added ray stopping at the solar surface for frequencies above
  //             ~127 MHz
  // 2009-Jan-13 The array stopping correction was unsuccessful, so the program
  //             reverted to the state before 2009-Jan-09. This requires more 
  //             work.
  // 2009-Jul-01 For calculating the brightness temperature in the Intens_I[],
  //             the parameters Freq and OpDepth_I[] are added 
  //


  double extern const cPi;              // = 3.1415926535897931;
  double extern const ProtonChargeSGSe; // = 4.8e-10; // StatCoulombs, SGSe
  double extern const cProtonMass;      // = 1.672621636E-24; //g
  double extern const cElectronMass;    // = 9.10938215E-28; //g
  double extern const cSolarRadius;     // = 6.955E10cm

  short static const FALSE = 0, TRUE = 1;
  const static double cZero = 0.0, cOne = 1.0, cTwo = 2.0, cHalf = 0.5;
  const static double cThird = 0.33333333333333333333, cFour = 4.0;;

  int i, j, iRay, n2Ry = 2*nRay;
  double static Dir1_D[3], Omega_D[3];	
  double static DirVert_D[3];	
  double static StepX_D[3], StepY_D[3], RelGradRefrInx_D[3];
  double static GradEps_D[3], PosHalfBk_D[3];
  double static Dir_D[3], Xprod_D[3];
  double static Tol, TolSqr, RhoCrInv, AbsMinStep;
  double static StepXSqr, StepYSqr, StepXdbl_D[3];
  double Te, dTb;
  double Ne, CAbsorp, dOpDepth;
  double xi;   // 1.4 in Chromosphere, 2.0 in Corona,  cm^6 K^1.5 Hz^2
  double const static Te_corona = 1.0e6; // K
  double const static Te_chromo = 3.0e4; // K
  double const static R0 = 6.955e5;   // km, Solar radius
  double const static h_chromo = 10000.0; // km, Chromosphere height
  double static r_chromo; // km, Chromosphere radius

  double HalfDS;            // DS halved
  double Pos_D[3];          // Position vector for current iRay
  double PosBefBoris_D[3];  // Position vector before Boris step 
  double Eps, EpsHalfBk, Rho2RhoCr, Rho2RhoCr1, Rho2RhoCr2;
  double Coef, Curv, Curv1;
  double LCosAl; // L is inverse grad of \epsilon, Alpha is incidence angle
  double LCosAl2;    // = LCosAl^2
  double GradEpsSqr, GradEpsDotDir;
  double LPrb, GradEps, dDS;
  short deb = 0;
  double SolDist;
  double DSEpsPrec;
  double const static EpsPrec = 0.0001;
  double static DistCr = 0.0;
  double const static SolarRadius = 1.0;
  int const static TrcRay = 21500; //19500; //21100; //22700; //25100;
  double const static TolEps = 1e-6;
  short const static CntMax = 50;
  // Bits for Flags_I
  short const static Penetr = 0x0001;  // Eps < 0; Critical surface search
  short const static WasRefl = 0x0002; // First step after reflection

  // Allocatables:
  double static *GradRho_DI = 0;	// [3,nRay]
  double static *PosPr_DI = 0;	        // [3,nRay], previous Pos_DI 
  double static *DirPr_DI = 0;	        // [3,nRay], previous Dir_DI 
  double static *Rho_I = 0, *DS_New_I = 0, *DistToCrSurf_I = 0;	// [nRay]
  short static *GentleRay_I = 0;        // [nRay] TRUE for shallow rays
  //
  // Initialization at first entry
  // or at 
  //
  if (*NewEntry) {
    *NewEntry = FALSE;
    if (GradRho_DI) {free(GradRho_DI); GradRho_DI = 0;}
    if (PosPr_DI) {free(PosPr_DI); PosPr_DI = 0;}
    if (DirPr_DI) {free(DirPr_DI); DirPr_DI = 0;}
    if (Rho_I) {free(GradRho_DI); GradRho_DI = 0;}
    if (DS_New_I) {free(DS_New_I); DS_New_I = 0;}
    if (DistToCrSurf_I) {free(DistToCrSurf_I); DistToCrSurf_I = 0;}
    if (GentleRay_I) {free(GentleRay_I); GentleRay_I = 0;}
    if (*NewEntry == -1) {
      return;    // Free memory and exit
    }
    printf("\n*NewEntry = %i\n\n", (int)(*NewEntry));

    GradRho_DI = (double *)calloc(3*nRay, sizeof(double)); 
    PosPr_DI = (double *)calloc(3*nRay, sizeof(double)); 
    DirPr_DI = (double *)calloc(3*nRay, sizeof(double)); 
    Rho_I = (double *)calloc(nRay, sizeof(double));
    DS_New_I = (double *)calloc(nRay, sizeof(double));
    DistToCrSurf_I = (double *)calloc(nRay, sizeof(double));
    GentleRay_I = (short *)calloc(nRay, sizeof(short));
    for (i = 0; i < nRay; i++) {
      DistToCrSurf_I[i] = cZero;
      GentleRay_I[i] = TRUE;
    }
    RhoCrInv = cOne/RhoCr;
    // Allow minimum ten points between a vacuum and a critical surface and
    // minimum 10 points over 1 rad of the curvature:
    Tol = min(TolInit,0.1);  
    TolSqr = pow(Tol,2);
    // One ten-thousandth of average step
    AbsMinStep = 1e-4*sum(DS_I, nRay)/nRay; 
    r_chromo = 1.0 + h_chromo/R0;  /* Chromospheric radius in solar radii */
    printf("\nAbsMinStep = %g\n\n", AbsMinStep);
    if (1) { 
      printf("TolInit = %g, Tol = %g\n", TolInit, Tol);
      printf("RhoCr = %g\n", RhoCr);
      printf("RhoCrInv = %g\n", RhoCrInv);
      printf("r_chromo = %g\n", r_chromo);
      //printf("Dir_DI = \n");
      //print2d(3, nRay, Dir_DI);
    }
  } // if (NewEntry) 


  //
  // Advance all the rays by 1/2 of the DS step
  //
  for (iRay=0; iRay < nRay; iRay++) {
    
    // Do not process the rays that are done or bad
    if (ExclRay_I[iRay])
      continue;     //---------------------------------------------------->>
    
    // Advance r by 1/2 DS 
    HalfDS = cHalf*DS_I[iRay];
    // Pos_DI is moved by 1/2 DS in the Dir_DI direction 
    Pos_DI[0][iRay] += Dir_DI[0][iRay]*HalfDS; 
    Pos_DI[1][iRay] += Dir_DI[1][iRay]*HalfDS;
    Pos_DI[2][iRay] += Dir_DI[2][iRay]*HalfDS;
  }
  


  Get_Plasma_Density(nRay, Pos_DI, Rho_I, GradRho_DI, DS_New_I, Flags_I);



  /* printf("Rho_I = \n"); */
  /* print2d(4, 4, Rho_I); */
  /* printf("\n"); */




  for (iRay = 0; iRay < nRay; iRay++) {
    
    // Do not process the rays that are done or bad
    if (ExclRay_I[iRay]) continue; //------------------------------------>>

    HalfDS = cHalf*DS_I[iRay];

    Rho2RhoCr = Rho_I[iRay]*RhoCrInv;

    Eps = cOne - Rho2RhoCr;  // Dielectric permittivity

    //GradEps_D = -GradRho_DI(:,iRay)*RhoCrInv
    GradEps_D[0] = -GradRho_DI[iRay]*RhoCrInv;
    GradEps_D[1] = -GradRho_DI[nRay+iRay]*RhoCrInv;
    GradEps_D[2] = -GradRho_DI[n2Ry+iRay]*RhoCrInv;
    
    GradEpsSqr = dot_product(GradEps_D, GradEps_D);
    
    //GradEpsDotDir = sum(GradEps_D*Dir_DI(:,iRay))
    GradEpsDotDir = GradEps_D[0]*Dir_DI[0][iRay] + 
      GradEps_D[1]*Dir_DI[1][iRay] +
      GradEps_D[2]*Dir_DI[2][iRay];
    
    
    
    
    if (iRay == TrcRay) {
      printf("===============================================\n");
      SolDist = 0.0;
      printf("AFTER Get_Pl_Den: Pos = "); 
      for (i = 0; i < 3; i++) {
	printf("%g ", Pos_DI[i][iRay]);
	SolDist += pow(Pos_DI[i][iRay],2);
      }
      printf("\n"); 
      SolDist = sqrt(SolDist);
      printf("AFTER Get_Pl_Den: Ray [%i,%i]\n", iRay%200, iRay/200);

      printf("AFTER Get_Pl_Den: SolDist = %g, DSnew = %g\n", SolDist, 
	     DS_New_I[iRay]);
      printf("AFTER Get_Pl_Den: Rho_I = %g, RhoCr = %g, Flags_I = %i\n", 
	     Rho_I[iRay], RhoCr, Flags_I[iRay]);
      printf("                : Eps = %g\n", Eps); 

      fprintf(fh,"|r|=%12.4e r=[%12.4e %12.4e %12.4e] " \
	      "v=[%12.4e %12.4e %12.4e] "			\
	      "rho=%12.4e grRho=[%12.4e %12.4e %12.4e]\n", 
	      SolDist, Pos_DI[0][iRay], Pos_DI[1][iRay], Pos_DI[2][iRay],
	      Dir_DI[0][iRay], Dir_DI[1][iRay], Dir_DI[2][iRay],
	      Rho_I[iRay], 
	      GradRho_DI[iRay], GradRho_DI[nRay+iRay],
	      GradRho_DI[n2Ry+iRay]);
    }




    
    if (Flags_I[iRay]&Penetr) { // Condition "Critical surface penetration"

      if (iRay ==  TrcRay) { 
	  printf("FLAGS: [%i,%i]\n", iRay/200, iRay%200);
      }


      if (fabs(Eps) > TolEps) { // Another convergence step needed
	
	//
	// Restore the r shifted along v by DS/2 before the call to 
	//     get_plasma_density():
	// r = r + v*DS/2
	//
	Pos_DI[0][iRay] -= Dir_DI[0][iRay]*HalfDS; 
	Pos_DI[1][iRay] -= Dir_DI[1][iRay]*HalfDS;
	Pos_DI[2][iRay] -= Dir_DI[2][iRay]*HalfDS;

	if (iRay ==  TrcRay) { 
	  printf("CONVERGENCE STEP: Cnt = %g\n", 
		 DirPr_DI[n2Ry+iRay]);
	}

	//
	// Use the bisection method to reach the critical surface
	//
       
	// Stop working on the ray if count exceeded maximum
	DirPr_DI[n2Ry+iRay] -= 1.0;  // Cnt--;
	if (DirPr_DI[n2Ry+iRay] < 0.0) { // if max count exceeded
	  ExclRay_I[iRay] = TRUE; // Bad Ray: do not process it anymore
	  continue; //----------------------------------------------------->>
	}
       
	if (Eps > 0.0) {
	  // DS_1 = DS; Eps_1 = Eps
	  PosPr_DI[iRay] = DS_I[iRay];      // DS_1 = DS
	  DirPr_DI[iRay] = Eps;             // Eps_1 = Eps
	}
	else { // Eps <= 0.0:
	  if (iRay == TrcRay) { 
	    printf("CONVERGENCE STEP CALC: NEWTON'S IMPROVEMENT!\n)");
	  }
	  // Newton's method is used to improve RHS DS point
	  // DS = DS - Eps/GradEpsDotDir
	  DS_I[iRay] -= Eps/GradEpsDotDir;
	  // DS_2 = DS; Eps_2 = Eps
	  PosPr_DI[nRay+iRay] = DS_I[iRay]; // DS_2 = DS
	  DirPr_DI[nRay+iRay] = Eps;        // Eps_2 = Eps
	}
       
	// DS = (DS_1 + DS_2)/2 
	DS_I[iRay] = (PosPr_DI[iRay] + PosPr_DI[nRay+iRay])*cHalf;
       
	if (iRay ==  TrcRay) { 
	  printf("CONVERGENCE STEP CALC: DS_1 = %g, DS_2 = %g\n", 
		 PosPr_DI[iRay], PosPr_DI[nRay+iRay]);
	  printf("CONVERGENCE STEP CALC: Eps_1 = %g, Eps_2 = %g\n", 
		 DirPr_DI[iRay], DirPr_DI[nRay+iRay]);
	  printf("CONVERGENCE STEP CALC: DS_2 - DS_1 = %g\n", 
		 PosPr_DI[nRay+iRay] - PosPr_DI[iRay]);
	}
       
	continue; //------------------------------------------------------->>
       
      } // if (fabs(Eps) > TolEps) 
     
      else { // fabs(Eps) <= TolEps: REFLECTION
	
	// The critical surface found. Reflecting the ray.
       
	//Flags_I[iRay] = FALSE;
	Flags_I[iRay] &= ~Penetr;  // Reset the penetration condition
	Flags_I[iRay] |= WasRefl;  // Mark the ray as "was Snell reflected"

	//DS_I[iRay] = PosPr_DI[n2Ry+iRay];      // Restore old DS !Not needed!
	HalfDS = cHalf*DS_I[iRay];
	
	//GradEps_D = -GradRho_DI(:,iRay)*RhoCrInv
	//GradEps_D[0] = -GradRho_DI[iRay]*RhoCrInv;
	//GradEps_D[1] = -GradRho_DI[nRay+iRay]*RhoCrInv;
	//GradEps_D[2] = -GradRho_DI[n2Ry+iRay]*RhoCrInv;
	
	GradEpsSqr = dot_product(GradEps_D, GradEps_D);
	
	//GradEpsDotDir = sum(GradEps_D*Dir_DI(:,iRay))
	//GradEpsDotDir = GradEps_D[0]*Dir_DI[0][iRay] + 
	//  GradEps_D[1]*Dir_DI[1][iRay] +
	//  GradEps_D[2]*Dir_DI[2][iRay];
	
	LCosAl = -GradEpsDotDir/GradEpsSqr;
	
	DirVert_D[0] = -LCosAl*GradEps_D[0]; // Here v_||
	DirVert_D[1] = -LCosAl*GradEps_D[1]; // Here v_||
	DirVert_D[2] = -LCosAl*GradEps_D[2]; // Here v_||
	
	if (iRay ==  TrcRay) { 
	  printf("REFLECTION: v =    (BEFORE CHANGE):");
	  printf("%g %g %g\n", Dir_DI[0][iRay], Dir_DI[1][iRay], 
		 Dir_DI[2][iRay]);
	  printf("REFLECTION: v_|| = ");
	  printf("%g %g %g\n", DirVert_D[0], DirVert_D[1], 
		 DirVert_D[2]);
	}

	//
	// Reflection by the Snell law:
	// v = v - 2*v_||
	//
	Dir_DI[0][iRay] = Dir_DI[0][iRay] - cTwo*DirVert_D[0];
	Dir_DI[1][iRay] = Dir_DI[1][iRay] - cTwo*DirVert_D[1];
	Dir_DI[2][iRay] = Dir_DI[2][iRay] - cTwo*DirVert_D[2];

	continue; //------------------------------------------------------>>

      } // fabs(Eps) < TolEps
    } // if (Flags_I[iRay]&Penetr) // Condition "Critical surface penetration"
    




    if (Eps < 0.0) { // Eps < 0.0: Penetration!
      //if (Eps < -TolEps) { // Eps < 0.0: Penetration!

      if (iRay ==  TrcRay) { 
	printf("PENETRATION!: Eps = %g\n", Eps);
      }

      //Flags_I[iRay] = TRUE;   // Mark "bad ray"
      Flags_I[iRay] |= Penetr;  // Mark "penetration"

      //
      // Revert r and v to previous integer point
      // 
      Dir_DI[0][iRay] = DirPr_DI[iRay];
      Dir_DI[1][iRay] = DirPr_DI[nRay+iRay];
      Dir_DI[2][iRay] = DirPr_DI[n2Ry+iRay];

      Pos_DI[0][iRay] = PosPr_DI[iRay];
      Pos_DI[1][iRay] = PosPr_DI[nRay+iRay];
      Pos_DI[2][iRay] = PosPr_DI[n2Ry+iRay];

      //
      // Use PosPr_DI and DirPr_DI as memory for Cnt and 
      // for left and right values of DS and Eps
      //
      PosPr_DI[iRay] = 0.0;                  // Store here DS_1
      PosPr_DI[nRay+iRay] = 2.0*DS_I[iRay];  // Store here DS_2
      PosPr_DI[n2Ry+iRay] = DS_I[iRay];      // Store here old DS
      DirPr_DI[iRay] = 1.0;                  // Store here Eps_1
      DirPr_DI[nRay+iRay] = Eps;             // Store here Eps_2
      DirPr_DI[n2Ry+iRay] = (double)CntMax;  // Store here iteration counter

      DS_I[iRay] -= Eps/GradEpsDotDir; // Newton's root approximation

      continue; //--------------------------------------------------------->>
    } // (Eps < 0.0) { // Penetration!





    //
    // Calculate main ray parameters
    //

    // Original Position (at an integer point): 
    PosHalfBk_D[0] = Pos_DI[0][iRay] - Dir_DI[0][iRay]*HalfDS; 
    PosHalfBk_D[1] = Pos_DI[1][iRay] - Dir_DI[1][iRay]*HalfDS;
    PosHalfBk_D[2] = Pos_DI[2][iRay] - Dir_DI[2][iRay]*HalfDS;
 
    //Rho2RhoCr = Rho_I[iRay]*RhoCrInv;

    //Eps = cOne - Rho2RhoCr;

    //GradEps_D = -GradRho_DI(:,iRay)*RhoCrInv
    GradEps_D[0] = -GradRho_DI[iRay]*RhoCrInv;
    GradEps_D[1] = -GradRho_DI[nRay+iRay]*RhoCrInv;
    GradEps_D[2] = -GradRho_DI[n2Ry+iRay]*RhoCrInv;

    //GradEpsSqr = sum(GradEps_D**2);
    //GradEpsSqr = pow(GradEps_D[0],2) + pow(GradEps_D[1],2) + 
    //             pow(GradEps_D[2],2);
    GradEpsSqr = dot_product(GradEps_D, GradEps_D);

    //GradEpsDotDir = sum(GradEps_D*Dir_DI(:,iRay))
    GradEpsDotDir = GradEps_D[0]*Dir_DI[0][iRay] + GradEps_D[1]*Dir_DI[1][iRay]
      + GradEps_D[2]*Dir_DI[2][iRay];

    LCosAl = -GradEpsDotDir/GradEpsSqr;

    DirVert_D[0] = -LCosAl*GradEps_D[0]; // Here v_||
    DirVert_D[1] = -LCosAl*GradEps_D[1]; // Here v_||
    DirVert_D[2] = -LCosAl*GradEps_D[2]; // Here v_||

    EpsHalfBk = Eps - GradEpsDotDir*HalfDS;

    //Curv = (cHalf*HalfDS/EpsHalfBk)**2* &
    //       (GradEpsSqr - GradEpsDotDir**2)
    // (a x b)^2 = a^2*b^2 - (a dot b)^2
    // if |v|=1, (a x v)^2 = a^2 - (a dot v)^2
    Curv = pow(cHalf*HalfDS/EpsHalfBk,2) * (GradEpsSqr - pow(GradEpsDotDir,2));

    if (iRay == TrcRay) {
      DistCr = EpsHalfBk/sqrt(GradEpsSqr);
      printf("DISTANCE TO CRIT SURF: Before = %g, Now = %g, RelDiff = %g\n", 
	     DistToCrSurf_I[iRay], DistCr, 
	     (DistToCrSurf_I[iRay]-DistCr)/DistCr);
      DistToCrSurf_I[iRay] = DistCr;
      printf("GradEps = %g %g %g\n", 
	     GradEps_D[0], GradEps_D[1], GradEps_D[2]);
      printf("DirVert = %g %g %g\n", 
	     DirVert_D[0], DirVert_D[1], DirVert_D[2]);
      printf("GradEpsDotDir = %g, Curv = %g, DS = %g\n", GradEpsDotDir, 
	     Curv, DS_I[iRay]);
      printf("Eps = %g, EpsHalfBk = %g\n", Eps, EpsHalfBk);
      printf("EpsHalfBk = %g, LCosAl = %g\n", 
	     EpsHalfBk, LCosAl);
    }


    
    if (GentleRay_I[iRay]) {
      
      //
      // Check if the trajectory curvature is too sharp to meet the Tolerance
      // If so, reduce the DS step  for the iRay-th ray and leave it until
      // the next call.
      //
      //
      if (Flags_I[iRay] & WasRefl) {
	Flags_I[iRay] &= ~WasRefl;  // Reset the reflection condition
	// and do not test for curvature!
      }
      else {

	if (Curv >=  TolSqr) { // Testing curvature 
	  DS_I[iRay] = DS_I[iRay]/(cTwo*sqrt(Curv/TolSqr));
	  // r = r_hb
	  Pos_DI[0][iRay] = PosHalfBk_D[0];
	  Pos_DI[1][iRay] = PosHalfBk_D[1];
	  Pos_DI[2][iRay] = PosHalfBk_D[2];
	  if (iRay == TrcRay) {
	    printf("Curv = %g, TolSqr = %g\n", Curv, TolSqr);
	  }
	  continue;    //---------------------------------------------------->>
	} // if (Curv >=  TolSqr) 

      }
      //
      // Test if some of the next points can get into the prohibited part of 
      // space with "negative" dielectric permittivity
      //
      
      if (GradEpsDotDir*HalfDS <= -cThird*EpsHalfBk) {
	//
	// Mark the ray as steep; memorize the distance to the critical 
	// surface; reduce step
	//
	GentleRay_I[iRay] = FALSE;
	DistToCrSurf_I[iRay] = EpsHalfBk/sqrt(GradEpsSqr);
	DS_I[iRay] =  cHalf*Tol*DistToCrSurf_I[iRay];
	
	if (iRay == TrcRay) {
	  printf("DistanceToCritSurf = %g, DS = %g\n", 
		 DistToCrSurf_I[iRay], DS_I[iRay]);
	}
	
	
	//Pos_DI(:,iRay) = PosHalfBk_D
	Pos_DI[0][iRay] = PosHalfBk_D[0];
	Pos_DI[1][iRay] = PosHalfBk_D[1];
	Pos_DI[2][iRay] = PosHalfBk_D[2];
	continue;    //----------------------------------------------------->>
      } // if (GradEpsDotDir*HalfDS <= -cThird*EpsHalfBk)
      
    } // if (GentleRay_I[iRay])
    
    






    //
    // In normal mode (i.e. no penetration)
    // Save previous values of Pos and Dir before they will be changed
    //
    //if (Flags_I[iRay] == FALSE) {
    if (!(Flags_I[iRay] & Penetr)) { // i.e. if not "penetration"

      PosPr_DI[iRay] =      Pos_DI[0][iRay];
      PosPr_DI[nRay+iRay] = Pos_DI[1][iRay];
      PosPr_DI[n2Ry+iRay] = Pos_DI[2][iRay];
      
      DirPr_DI[iRay] =      Dir_DI[0][iRay];
      DirPr_DI[nRay+iRay] = Dir_DI[1][iRay];
      DirPr_DI[n2Ry+iRay] = Dir_DI[2][iRay];

    }
 
    
    //Solar Distance
    SolDist = 0.0;
    for (i = 0; i < 3; i++) {
      SolDist += pow(Pos_DI[i][iRay],2);
    }
    SolDist = sqrt(SolDist);
    

    //
    // Either:
    // - switch to opposite branch of parabola;
    // - or make a Boris step
    //
    if ((GradEpsDotDir*HalfDS <= -cThird*EpsHalfBk) 
	|| (DS_I[iRay] < AbsMinStep)) {

      // Switch to the opposite branch of parabolic trajectory
      //
      // When a next step can drive the ray into the area with the
      // plasma density greater then its critical value, then a special 
      // technique of "parabolic ray reflection" is employed. 
      // It can be shown that a ray trajectory in the medium with constant
      // gradient of dielectric permittivity is a parabola. If the step DS_I
      // is small enough we can assume the grad \epsilon constant and hence
      // assume that the ray approaches the critical surface in a parabolic 
      // path.
      // We calculate the parameters of the parabola, namely:
      // StepX_D -- vector along the -grad \epsilon, with the length equal 
      //     to the distance from the PosHalfBk_D to the parabola extremum;  
      // StepY_D -- vector perpendicular to StepX_D, lying inside of the
      //     parabola plane, pointing at the opposite parabola branch, and with
      //     the length equal the distance from PosHalfBk_D to the 
      //     opposite branch of parabola.
      // The parabolic reflection is just replacement of the Pos_DI with
      // the symmetric point at the opposite branch of parabola, and changing
      // the "incident" direction Dir_DI for the "departing" ray direction
      // according to Snell law.
      //

      //LCosAl = -GradEpsDotDir/GradEpsSqr;

      if (iRay == TrcRay) {
	fprintf(fh, "PARABOLA\n");
	printf("PARABOLA SWITCHING, [%i,%i]=iRay = %i, DS = %g\n", 
	       iRay/200, iRay%200, iRay, DS_I[iRay]);
 	printf("PARABOLA SWITCHING, x,y,z    = %g %g %g\n", 
	       Pos_DI[0][iRay], Pos_DI[1][iRay], Pos_DI[2][iRay]);
	printf("PARABOLA SWITCHING, Dir_DI = %g %g %g\n", 
	       Dir_DI[0][iRay], Dir_DI[1][iRay], Dir_DI[2][iRay]);
	printf("PARABOLA SWITCHING, LCosAl=%g, gradEps*Dir=%g, GradEps^2=%g\n", 
	       LCosAl, GradEpsDotDir, GradEpsSqr);
      }

      //DirVert_D[0] = -LCosAl*GradEps_D[0]; // Here v_proj
      //DirVert_D[1] = -LCosAl*GradEps_D[1]; // Here v_proj
      //DirVert_D[2] = -LCosAl*GradEps_D[2]; // Here v_proj

      //StepY_D = Dir_DI(:,iRay) - DirVert_D       // Here c; |c| = sin \alpha
      StepY_D[0] = Dir_DI[0][iRay] - DirVert_D[0]; // Here c; |c| = sin \alpha
      StepY_D[1] = Dir_DI[1][iRay] - DirVert_D[1]; // Here c; |c| = sin \alpha
      StepY_D[2] = Dir_DI[2][iRay] - DirVert_D[2]; // Here c; |c| = sin \alpha

      //Dir_DI(:,iRay) = Dir_DI(:,iRay) - cTwo*DirVert_D       
      Dir_DI[0][iRay] = Dir_DI[0][iRay] - cTwo*DirVert_D[0];
      Dir_DI[1][iRay] = Dir_DI[1][iRay] - cTwo*DirVert_D[1];
      Dir_DI[2][iRay] = Dir_DI[2][iRay] - cTwo*DirVert_D[2];

      //
      //   We need to have |Step_Y| = 4 * sin(\alpha)*cos(\alpha)*EpsHalfBk*L
      //   Now the direction of Step_Y is right, 
      //   the length of it is equal to sin(\alpha)
      //   Multiply it by L*Cos(\alpha)*EpsHalfBk
      //

      StepY_D[0] = cFour*StepY_D[0]*LCosAl*EpsHalfBk; 
      StepY_D[1] = cFour*StepY_D[1]*LCosAl*EpsHalfBk; 
      StepY_D[2] = cFour*StepY_D[2]*LCosAl*EpsHalfBk; 

      Pos_DI[0][iRay] = PosHalfBk_D[0] + StepY_D[0];
      Pos_DI[1][iRay] = PosHalfBk_D[1] + StepY_D[1];
      Pos_DI[2][iRay] = PosHalfBk_D[2] + StepY_D[2];
	
      //
      //   Step_X is in the direction of -grad \epsilon, whose vector is of 
      //   the length of 1/L
      //   The length of Step_X is cos^2(\alpha)*L*EpsHalfBk 
      //   Thus,
      //

      LCosAl2 = LCosAl*LCosAl;
      StepX_D[0] = (EpsHalfBk*LCosAl2)*GradEps_D[0];
      StepX_D[1] = (EpsHalfBk*LCosAl2)*GradEps_D[1];
      StepX_D[2] = (EpsHalfBk*LCosAl2)*GradEps_D[2];

      if (iRay == TrcRay) {
	printf("PARABOLA SWITCHING, EpsHalfBk = %g, LCosAl = %g\n", 
	       EpsHalfBk, LCosAl);
	printf("PARABOLA SWITCHING, StepX_D = "); print1d(3, StepX_D);
	printf("PARABOLA SWITCHING, StepY_D = "); print1d(3, StepY_D);
	printf("\n"); 
      }

      ////LPrb = cTwo*sqrt(sum(StepX_D**2))
      //LPrb = sqrt(sum((cTwo*StepX_D)**2) + sum(StepY_D**2))
      //StepXSqr = pow(cTwo*StepX_D[0],2) + pow(cTwo*StepX_D[1],2) + 
      //           pow(cTwo*StepX_D[2],2);
      //StepYSqr = pow(StepY_D[0],2) + pow(StepY_D[1],2) + pow(StepY_D[2],2);
      StepXdbl_D[0] = cTwo*StepX_D[0];
      StepXdbl_D[1] = cTwo*StepX_D[1];
      StepXdbl_D[2] = cTwo*StepX_D[2];
      StepXSqr = dot_product(StepXdbl_D, StepXdbl_D); 
      StepYSqr = dot_product(StepY_D, StepY_D); 
      LPrb = sqrt(StepXSqr + StepYSqr);  // the parabola segment length
      if (deb) { printf("LPrb = %g\n", LPrb);}
      //Intens_I(iRay) = Intens_I(iRay) 
      //                  + LPrb*(Rho2RhoCr**2)*(cHalf - Rho2RhoCr)**2
      ////-->Intens_I[iRay] += LPrb*pow(Rho2RhoCr,2)*pow(cHalf - Rho2RhoCr,2);


      // In the Intensity[] the brightness temperature is saved!
      // Calculate brightness temperature as integral along ray path
      // of electron temperature times exp(-tau), tau is optical depth: ??????
      // The tau is integral along ray path of the absorption coefficient kappa
      //
      //Te = (SolDist > r_chromo) ? Te_corona : Te_chromo;
      //z = (a > b) ? a : b;
      if (SolDist > r_chromo) {
	Te = Te_corona;
	//xi = 2.0; // cm^6 K^1.5 Hz^2 -- from Sheridan
	xi = 0.0765766935; // cm^6 K^1.5 Hz^2 -- from Smerd
      }
      else {
	Te = Te_chromo;
	//xi = 1.4; // cm^6 K^1.5 Hz^2 -- from Sheridan
	xi = 0.1148650403; // cm^6 K^1.5 Hz^2 -- from Smerd
      }
      Ne =  Rho_I[iRay]/cProtonMass;  // Electron density
      CAbsorp = xi*pow(Ne,2)/(pow(Te,1.5)*
			      pow(Freq,2)*sqrt(Eps));
      dOpDepth = CAbsorp*LPrb*cSolarRadius;
      OpDepth_I[iRay] += dOpDepth;
      // Actually, Intens_I here holds the  Tb (brightness temperature)
      //Intens_I[iRay] += Te*exp(-OpDepth_I[iRay])*dOpDepth;  
      dTb =  Te*exp(-OpDepth_I[iRay])*dOpDepth;
      Intens_I[iRay] += dTb; 
 
      if (iRay ==  TrcRay) { 
	/*
	fprintf(fh,"%12.4e %12.4e %12.4e %12.4e %12.4e %12.4e %12.4e " \
		"%12.4e %12.4e\n", 
		SolDist, LPrb, Ne, CAbsorp, dOpDepth, OpDepth_I[iRay], 
		dTb, dTb/(LPrb*R0), Intens_I[iRay]);
	*/
	printf("PARABOLA: r = ");
	printf("%g %g %g\n", Pos_DI[0][iRay], Pos_DI[1][iRay], 
	       Pos_DI[2][iRay]);
	printf("PARABOLA: v = ");
	printf("%g %g %g\n", Dir_DI[0][iRay], Dir_DI[1][iRay], 
	       Dir_DI[2][iRay]);
	printf("PARA INTEGR: Te = %g, Ne = %g\n", Te, Ne);
	printf("PARA INTEGR: CAbsorp = %g, dOpDepth = %g, OpDepth_I = %g\n", 
	       CAbsorp, dOpDepth, OpDepth_I[iRay]);
	printf("PARA INTEGR: exp(-tau) = %g, dTb = %g, Tb = %g\n", 
	       exp(-OpDepth_I[iRay]), dTb, Intens_I[iRay]);
      }



    } // if ((GradEpsDotDir*HalfDS <= -cThird*EpsHalfBk) 
      //  || (DS_I[iRay] < AbsMinStep))


    else { // Boris' step
      //
      // Make a step using Boris' algorithm
      //
	
      if (deb) {printf("// Make a step using Boris' algorithm...;\n");}

      PosBefBoris_D[0] = Pos_DI[0][iRay];
      PosBefBoris_D[1] = Pos_DI[1][iRay];
      PosBefBoris_D[2] = Pos_DI[2][iRay];

      // grad(n) = grad(eps(i+1/2))/(2*eps(i+1/2))
      Coef = cHalf*HalfDS/(cOne - Rho2RhoCr);
      RelGradRefrInx_D[0] = Coef*GradEps_D[0];
      RelGradRefrInx_D[1] = Coef*GradEps_D[1]; 
      RelGradRefrInx_D[2] = Coef*GradEps_D[2]; 

      Dir_D[0] = Dir_DI[0][iRay];
      Dir_D[1] = Dir_DI[1][iRay];
      Dir_D[2] = Dir_DI[2][iRay];

      //Omega_D = cross_product(RelGradRefrInx_D, Dir_DI(:,iRay))
      cross_product(RelGradRefrInx_D, Dir_D, Omega_D);

      //Dir1_D = Dir_DI(:,iRay) + cross_product(Dir_DI(:,iRay), Omega_D)
      cross_product(Dir_D, Omega_D, Xprod_D);
      Dir1_D[0] = Dir_D[0] + Xprod_D[0];
      Dir1_D[1] = Dir_D[1] + Xprod_D[1];
      Dir1_D[2] = Dir_D[2] + Xprod_D[2];
	
      //Omega_D = RelGradRefrInx_D x Dir1_D
      cross_product(RelGradRefrInx_D, Dir1_D, Omega_D);

      // v1 = v + v x Omega
      cross_product(Dir_D, Omega_D, Xprod_D);
      Dir1_D[0] = Dir_D[0] + Xprod_D[0];
      Dir1_D[1] = Dir_D[1] + Xprod_D[1];
      Dir1_D[2] = Dir_D[2] + Xprod_D[2];

      //Curv1 = pow(Omega_D[0],2) + pow(Omega_D[1],2) + pow(Omega_D[2],2);
      Curv1 = dot_product(Omega_D, Omega_D);
      Coef = cTwo/(cOne + Curv1);
      //Dir_DI(:,iRay) = Dir_DI(:,iRay) + Coef*cross_product(Dir1_D, Omega_D)
      cross_product(Dir1_D, Omega_D, Xprod_D);
      Dir_DI[0][iRay] = Dir_DI[0][iRay] + Coef*Xprod_D[0];
      Dir_DI[1][iRay] = Dir_DI[1][iRay] + Coef*Xprod_D[1];
      Dir_DI[2][iRay] = Dir_DI[2][iRay] + Coef*Xprod_D[2];

      //Pos_DI(:,iRay) = Pos_DI(:,iRay) + Dir_DI(:,iRay)*HalfDS
      Pos_DI[0][iRay] += Dir_DI[0][iRay]*HalfDS;
      Pos_DI[1][iRay] += Dir_DI[1][iRay]*HalfDS;
      Pos_DI[2][iRay] += Dir_DI[2][iRay]*HalfDS;


      if (iRay == TrcRay) {
	printf("BORIS' STEP\n");
	printf("Dir = "); 
	for (i = 0; i < 3; i++) {
	  printf("%g ", Dir_DI[i][iRay]);
	}
	printf("\n"); 
	printf("Pos = "); 
	for (i = 0; i < 3; i++) {
	  printf("%g ", Pos_DI[i][iRay]);
	}
	printf("\n"); 
	printf("SolDist = %g\n", SolDist);
      }

      // In the Intensity[] the brightness temperature is saved!
      // Calculate brightness temperature as integral along ray path
      // of electron temperature times exp(-tau), tau is optical depth:
      // The tau is integral along ray path of the absorption coefficient kappa
      //
      //Te = (SolDist > r_chromo) ? Te_corona : Te_chromo;
      //z = (a > b) ? a : b;
     if (SolDist > r_chromo) {
	Te = Te_corona;
	//xi = 2.0; // cm^6 K^1.5 Hz^2 -- from Sheridan
	xi = 0.0765766935; // cm^6 K^1.5 Hz^2 -- from Smerd
      }
      else {
	Te = Te_chromo;
	//xi = 1.4; // cm^6 K^1.5 Hz^2 -- from Sheridan
	xi = 0.1148650403; // cm^6 K^1.5 Hz^2 -- from Smerd
      }
      Ne =  Rho_I[iRay]/cProtonMass;  // Electron density
      CAbsorp = xi*pow(Ne,2)/(pow(Te,1.5)*
			      pow(Freq,2)*sqrt(Eps));
      dOpDepth = CAbsorp*DS_I[iRay]*cSolarRadius;
      OpDepth_I[iRay] += dOpDepth;
      // Actually, Intens_I here holds the  Tb (brightness temperature)
      //Intens_I[iRay] += Te*exp(-OpDepth_I[iRay])*dOpDepth;  
      dTb =  Te*exp(-OpDepth_I[iRay])*dOpDepth;
      Intens_I[iRay] += dTb;  


      if (iRay ==  TrcRay) { 
	/*
	fprintf(fh,"%12.4e %12.4e %12.4e %12.4e %12.4e %12.4e %12.4e " \
		"%12.4e %12.4e\n", 
		SolDist, DS_I[iRay], Ne, CAbsorp, dOpDepth, OpDepth_I[iRay], 
		dTb, dTb/(DS_I[iRay]*R0), Intens_I[iRay]);
	*/
	printf("BORIS: r = ");
	printf("%g %g %g\n", Pos_DI[0][iRay], Pos_DI[1][iRay], 
	       Pos_DI[2][iRay]);
	printf("BORIS: v = ");
	printf("%g %g %g\n", Dir_DI[0][iRay], Dir_DI[1][iRay], 
	       Dir_DI[2][iRay]);
	printf("BORIS INTEGR: Te = %g, Ne = %g\n", Te, Ne);
	printf("BORIS INTEGR: CAbsorp = %g, dOpDepth = %g, OpDepth_I = %g\n", 
	       CAbsorp, dOpDepth, OpDepth_I[iRay]);
	printf("BORIS INTEGR: exp(-tau) = %g, dTb = %g, Tb = %g\n", 
	       exp(-OpDepth_I[iRay]), dTb, Intens_I[iRay]);

      }


      // Calculate Itensity as integral along ray path of Volumetric Emissivity 
      // Emissivity ~ (rho/rho_cr)^2*(1/2 - rho/rho_cr)^2  (Sokolov, 2006):
      //--> Intens_I[iRay] = Intens_I[iRay] + 
      //            DS_I[iRay]*pow(Rho2RhoCr,2)*pow(cHalf - Rho2RhoCr, 2);
      // Alternative Emissivity calcylation (Smerd, 1949):
      // Emissivity ~ (rho/rho_cr)*(1 - rho/rho_cr)
      //--> Intens_I[iRay] += DS_I[iRay]*Rho2RhoCr*sqrt(cOne - Rho2RhoCr);

    } // else { // Boris' step




    //
    //   The code below makes gradual increases of the DS up to the value
    // specified in DSNew. The smooth step increase is required so as not to
    // get into the space behind the critical surface, stepping with DS that
    // instantly changes from a very little size to the normal DSNew length.
    // DS is usually reduced in a close vicinity of the critical surface,
    // where the ray is travelling along a very sharp curve with high curvature.
    // For many rays it means fractioning of the DS down several orders of 
    // magnitude, therefore the new step trial should start from a bigger step
    // of the same order of magnitude.
    //   This problem is solved using a non-linear difference equation:
    //           Y(i+1) = [2 - Y(i)/X(i)]*Y(i),
    // where X(i) is the desired final DS length from DSNew, and
    // Y(i) is the actual DS length. A simple analysis of the equation shows
    // that, when Y is small compared to X, the next value of Y will be almost 
    // 2*X, so the DS would grow in a geometrical progression. However, as
    // Y approaches X, its growth rate becomes slower. However, Y always reaches
    // X in several steps. One can check that for Y = X the next value of Y is 
    // always that of X.
    // 
    if (GentleRay_I[iRay]) {
      //
      // For shallow rays the DS is increased unconditionally
      //
      if (DS_New_I[iRay] > DS_I[iRay])
	DS_I[iRay] = (cTwo - DS_I[iRay]/DS_New_I[iRay])*DS_I[iRay];
      else
	DS_I[iRay] = DS_New_I[iRay];
    } // if (GentleRay_I[iRay])
    else { // Steep ray
      //
      // If the iRay-th ray is marked as steep (i.e. "not gentle" or 
      // "not shallow") then the code below increases its DS only if the 
      // current distance to the critical surface, calculated as 
      //     \epsilon / grad \epsilon, 
      // is greater than this distance value saved along with marking the ray 
      // as steep in the DistToCrSurf_I. 
      //   This can happen in two cases: 
      // (1) either the ray was "parabola reflected" and after several steps it
      // went away from the surface by the distance where the parabola 
      // switching occurred; 
      // (2) or the ray is not steep any more because the current DS is so 
      // small, that the next step does not penetrate the critical surface.
      //   The ray is reverted to "gentle" or "shallow"
      //
      if (EpsHalfBk > DistToCrSurf_I[iRay]*sqrt(GradEpsSqr)) {
	GentleRay_I[iRay] = TRUE;
	if (DS_New_I[iRay] > DS_I[iRay])
	  DS_I[iRay] = (cTwo - DS_I[iRay]/DS_New_I[iRay])*DS_I[iRay];
	else
	  DS_I[iRay] = DS_New_I[iRay];
      } // if (EpsHalfBk > DistToCrSurf_I[iRay]*sqrt(GradEpsSqr)) 
    } // else { // Steep ray
  } // for (iRay = 1; iRay < nRay; iRay++) 
} // beam_path()
