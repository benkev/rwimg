#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "raytrace.h"



void advance_beam_stokesIV(int nRay, 
			   double Param_I[],
			   double Pos_ID[nRay][3], 
			   double Dir_ID[nRay][3], 
			   double DS_I[nRay], 
			   void (*plasma_params)(),
			   short Flags_I[nRay], 
			   double Stokes_I[nRay,2], 
			   double OpDepth_I[nRay],
			   double Rho_I[nRay],
			   double GradRho_ID[nRay][3],
			   double Bfield_ID[nRay][3],
			   double PosPr_ID[nRay][3],
			   double DirPr_ID[nRay][3],
			   double DS_New_I[nRay],
			   double DistToCrSurf_I[nRay]) {
  //
  //   The subroutine beam_path() makes raytracing and emissivity integration 
  // along ray paths.
  // It works on a group of rays (a beam). Each ray is specified by its 
  // Cartesian Pos_ID and its direction cosines Dir_ID. The subroutine 
  // calculates new Pos_ID, which is DS_I away from the initial 
  // position. It calculates the intensity along the step as emissivity by 
  // DS and adds this to the Intgrl_I. Thus, after series of calls
  // to beam_path(), the Intgrl_I contains the result of integration along 
  // the paths of every ray. 
  //   The electromagnetic rays are refracted in plasma due to its non-zero 
  // refractive index, which is the square root of the dielectric permittivity
  // \epsilon. The \epsilon, in turn, is a function of the plasma density. The 
  // external subprogram, plasma_density(), provides the plasma Rho_I 
  // along with its gradient, GradRho_ID at the Pos_ID.
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
  // plasma_density():  external subroutine that returns the plasma density
  //   Rho_I and its gradient GradRho_ID. It also provides the recommended step 
  //   size DS_New_I and asserts the Flags_I (assigns TRUE) for a ray should it
  //   (illegally) penetrate into the region with "negative" dielectric 
  //   permittivity.
  // nRay:           number of rays being processed.
  // INACTIVE:  the caller program can use this bit in Flags_I to stop 
  //   processing of any individual ray when it already finished travelling 
  //   inside of a specified part of space. Set the INACTIVE bit in  
  //   Flags_I to 1 to leave it unprocessed during the subsequent calls to
  //   beam_path(). Before the first call to beam_path() all the elements of 
  //   Flags_I should be set to 0; then all the rays will be processed.
  // Pos_ID:    Cartesian position vectors of the rays.
  // Dir_ID:    Direction cosines of the rays.
  // DS_I:      Current step values for the rays. Set the elements of Dir_ID to
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
  // Intgrl_I:    the intensities of each ray calculated as the integral of 
  //   emissivity along the ray path. During each call to ray_path(), each 
  //   element of the Intgrl_I is incremented by the integral along the step 
  //   DS_I. Set all the elements of Intgrl_I to 0.0 before the first call to 
  //   the beam_path().
  // Flags_I:      the TRUE elements of this logical array indicate that the 
  //   corresponding ray penetrated the "prohibited" region of space with the 
  //   plasma density above its critical value. Normally, it should never 
  //   happen. However, in case the algorithm made such an error, the flagged 
  //   rays should be considered as "bad" and thrown away from the resultant 
  //   Intgrl_I. Set all the elements of Flags_I to FALSE before calling 
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
  // 2009-Jul-01 For calculating the brightness temperature in the Intgrl_I[],
  //             the parameters Freq and OpDepth_I[] are added 
  // 2011-Mar-15 Added to Python extension module rtcore (part of the 
  //             raytrace pachage). The name beam_path() changed to 
  //             advance_beam().
  // 2011-Apr-01 Introduced dynamic linking for plasma_density() in 
  //             rtcoremodule.c. The pointer to plasma_density() is passed 
  //             as a parameter.
  //             


  /* self.param = numpy.array(( */
  /*     self.DeltaS,     # Initial ray path increment in solar radii  */
  /*     self.freq,            # Hz, radio wave frequency */
  /*     self.rhocr,           # Critical density at given self.freq */
  /*     self.RhoCrInv,        # 1/self.RhoCr */
  /*     self.tol,        # Maximum radians per one step  */
  /*     self.tol2,       # Tolerance self.tol squared (tol2 = tol^2) */
  /*     self.AbsMinStep, # Limit to adaptive step decrease */
  /*     self.TolEps,     # Tolerance for dielectric permittivity */
  /*     self.CntMax      # Upper limit to Newton iterations number */
  /*     self.h_chromo_km,     # km, Chromosphere thickness */
  /*     self.r_chromo,   # in Rsun units: "top" of chromosphere */
  /*     self.r_chro_cor, # in Rsun units: chromosphere-corona  "border" */
  /*     self.r_corona,   # in Rsun units: "bottom" of corona */
  /*     self.ProtonChargeCGSe,   # StatCoulombs, CGSe */
  /*     self.ProtonMass_g,    # g */
  /*     self.ProtonMassInv,   # 1/g */
  /*     self.ElectronMass_g,  # g */
  /*     self.Rsun_cm,         # cm, Solar radius */
  /*     self.Rsun_km,         # km, Solar radius */
  /*     self.RhoSurf,         # g/cm^3, density at solar surface */
  /*     self.Te_corona_K,       # K */
  /*     self.Te_chromo_K,       # K */
  /*     self.AU_m,              # Sun-earth distance (1 AU) in meters */
  /*     self.AU_Rsun,           # Sun-earth distance (1 AU) in solar radii */
  /*     ), dtype=numpy.double) */
  double Freq =             Param_I[1];
  double RhoCr =            Param_I[2];
  double RhoCrInv =         Param_I[3];
  double Tol =              Param_I[4];
  double Tol2 =             Param_I[5];
  double AbsMinStep =       Param_I[6]; 
  double TolEps =           Param_I[7]; 
  double CntMax =           Param_I[8]; 
  double r_chro_cor =       Param_I[11]; /* Chromosphere radius in Rsun */
  double ProtonMassInv =    Param_I[15];
  double Rsun_cm =          Param_I[17];
  double Te_corona_K =      Param_I[20];
  double Te_chromo_K =      Param_I[21];

  double const static cThird = 1.0/3.0;

  int i, iRay;
  double static Dir1_D[3], Omega_D[3];	
  double static DirVert_D[3];	
  double static StepX_D[3], StepY_D[3], RelGradRefrInx_D[3];
  double static GradEps_D[3], PosHalfBk_D[3];
  double static Dir_D[3], Xprod_D[3];
  double static StepXSqr, StepYSqr, StepXdbl_D[3];
  double Te, dTb;
  double Ne, CAbsorp, dOpDepth;
  double xi;   // 1.4 in Chromosphere, 2.0 in Corona,  cm^6 K^1.5 Hz^2

  double HalfDS;            // DS halved
  double PosBefBoris_D[3];  // Position vector before Boris step 
  double Eps, EpsHalfBk, Rho2RhoCr; //, Rho2RhoCr1, Rho2RhoCr2;
  double Coef, Curv, Curv1;
  double LCosAl; // L is inverse grad of \epsilon, Alpha is incidence angle
  double LCosAl2;    // = LCosAl^2
  double GradEpsSqr, GradEpsDotDir;
  double LPrb;
  double SolDist;
  double static DistCr = 0.0;
  int const static TrcRay = 22281; // for 200 MHz

  /*
   * Advance all the rays by 1/2 of the DS step
   */
  for (iRay=0; iRay < nRay; iRay++) {
    
    // Do not process the rays that are done or bad
    if (BIT_IS_ON(INACTIVE, iRay)) continue;     //------------------------>>
    
    // Advance r by 1/2 DS 
    HalfDS = 0.5*DS_I[iRay];
    // Pos_ID is moved by 1/2 DS in the Dir_DI direction 
    Pos_ID[iRay][0] += Dir_ID[iRay][0]*HalfDS; 
    Pos_ID[iRay][1] += Dir_ID[iRay][1]*HalfDS;
    Pos_ID[iRay][2] += Dir_ID[iRay][2]*HalfDS;
  }
  

  plasma_params(nRay, Param_I, Pos_ID, Rho_I, GradRho_ID, 
		Bfield_ID, DS_New_I, Flags_I);


  for (iRay = 0; iRay < nRay; iRay++) {
    
    // Do not process the rays that are done or bad
    if (BIT_IS_ON(INACTIVE, iRay)) continue;     //------------------------>>

    HalfDS = 0.5*DS_I[iRay];

    Rho2RhoCr = Rho_I[iRay]*RhoCrInv;

    Eps = 1.0 - Rho2RhoCr;  // Dielectric permittivity

    GradEps_D[0] = -GradRho_ID[iRay][0]*RhoCrInv;
    GradEps_D[1] = -GradRho_ID[iRay][1]*RhoCrInv;
    GradEps_D[2] = -GradRho_ID[iRay][2]*RhoCrInv;
    
    GradEpsSqr = dot_product(GradEps_D, GradEps_D);
    
    GradEpsDotDir = GradEps_D[0]*Dir_ID[iRay][0] + 
      GradEps_D[1]*Dir_ID[iRay][1] +
      GradEps_D[2]*Dir_ID[iRay][2];
    /* GradEpsDotDir = dot_product(GradEps_D, &Dir_ID[iRay][0]); */
    
    
    
    if (iRay == TrcRay) {
      printf("===============================================\n");
      SolDist = 0.0;
      printf("AFTER Get_Pl_Den: Pos = "); 
      for (i = 0; i < 3; i++) {
	printf("%g ", Pos_ID[iRay][i]);
	SolDist += pow(Pos_ID[iRay][i],2);
      }
      printf("\n"); 
      SolDist = sqrt(SolDist);
      printf("AFTER Get_Pl_Den: Ray[%i,%i] = Ray [%i]\n", 
	     iRay/200, iRay%200, iRay);

      printf("AFTER Get_Pl_Den: SolDist = %g, DSnew = %g\n", SolDist, 
	     DS_New_I[iRay]);
      printf("AFTER Get_Pl_Den: Rho_I = %g, RhoCr = %g, Flags_I = %i\n", 
	     Rho_I[iRay], RhoCr, Flags_I[iRay]);
      printf("                : Eps = %g\n", Eps); 

    }

    
    if (BIT_IS_ON(PENETR, iRay)) { // Condition "Critical surface penetration"
      if (iRay ==  TrcRay) { 
	printf("FLAGS: [%i,%i] = %d\n", iRay/200, iRay%200, Flags_I[iRay]);
      }

      if (fabs(Eps) > TolEps) { // Another convergence step needed
	
	//
	// Restore the r shifted along v by DS/2 before the call to 
	//     plasma_density():
	// r = r + v*DS/2
	//
	Pos_ID[iRay][0] -= Dir_ID[iRay][0]*HalfDS; 
	Pos_ID[iRay][1] -= Dir_ID[iRay][1]*HalfDS;
	Pos_ID[iRay][2] -= Dir_ID[iRay][2]*HalfDS;

	if (iRay ==  TrcRay) { 
	  printf("CONVERGENCE STEP: Cnt = %g\n", 
		 DirPr_ID[iRay][2]);
	}

	//
	// Use the bisection method to reach the critical surface
	//

       	// Stop working on the ray if count exceeded maximum
	DirPr_ID[iRay][2] -= 1.0;  // Cnt--;
	if (DirPr_ID[iRay][2] < 0.0) { // if max count exceeded
	  SET_BIT(BADRAY|INACTIVE, iRay);  // Do not process the ray anymore
	  continue; //----------------------------------------------------->>
	}
       
	if (Eps > 0.0) {
	  // DS_1 = DS; Eps_1 = Eps
	  PosPr_ID[iRay][0] = DS_I[iRay];      // DS_1 = DS
	  DirPr_ID[iRay][0] = Eps;             // Eps_1 = Eps
	}
	else { // Eps <= 0.0:
	  if (iRay == TrcRay) { 
	    printf("CONVERGENCE STEP CALC: NEWTON'S IMPROVEMENT!\n)");
	  }
	  // Newton's method is used to improve RHS DS point
	  // DS = DS - Eps/GradEpsDotDir
	  DS_I[iRay] -= Eps/GradEpsDotDir;
	  // DS_2 = DS; Eps_2 = Eps
	  PosPr_ID[iRay][1] = DS_I[iRay]; // DS_2 = DS
	  DirPr_ID[iRay][1] = Eps;        // Eps_2 = Eps
	}
       
	// DS = (DS_1 + DS_2)/2 
	DS_I[iRay] = (PosPr_ID[iRay][0] + PosPr_ID[iRay][1])*0.5;
       
	if (iRay ==  TrcRay) { 
	  printf("CONVERGENCE STEP CALC: DS_1 = %g, DS_2 = %g\n", 
		 PosPr_ID[iRay][0], PosPr_ID[iRay][1]);
	  printf("CONVERGENCE STEP CALC: Eps_1 = %g, Eps_2 = %g\n", 
		 DirPr_ID[iRay][0], DirPr_ID[iRay][1]);
	  printf("CONVERGENCE STEP CALC: DS_2 - DS_1 = %g\n", 
		 PosPr_ID[iRay][1] - PosPr_ID[iRay][0]);
	}
       
	continue; //------------------------------------------------------->>
       
      } // if (fabs(Eps) > TolEps) 
     
      else { // fabs(Eps) <= TolEps: REFLECTION
	
	// The critical surface found. Reflecting the ray.
       
	CLEAR_BIT(PENETR, iRay);   // Clear the penetration condition
	SET_BIT(WASREFL, iRay);  // Mark the ray as "was Snell reflected"

	//DS_I[iRay] = PosPr_ID[iRay][2];      // Restore old DS !Not needed!
	HalfDS = 0.5*DS_I[iRay];
	
	GradEpsSqr = dot_product(GradEps_D, GradEps_D);
	
	LCosAl = -GradEpsDotDir/GradEpsSqr;
	
	DirVert_D[0] = -LCosAl*GradEps_D[0]; // Here v_||
	DirVert_D[1] = -LCosAl*GradEps_D[1]; // Here v_||
	DirVert_D[2] = -LCosAl*GradEps_D[2]; // Here v_||
	
	if (iRay ==  TrcRay) { 
	  printf("REFLECTION: v =    (BEFORE CHANGE):");
	  printf("%g %g %g\n", Dir_ID[iRay][0], Dir_ID[iRay][1], 
		 Dir_ID[iRay][2]);
	  printf("REFLECTION: v_|| = ");
	  printf("%g %g %g\n", DirVert_D[0], DirVert_D[1], 
		 DirVert_D[2]);
	}

	//
	// Reflection by the Snell law:
	// v = v - 2*v_||
	//
	Dir_ID[iRay][0] = Dir_ID[iRay][0] - 2.0*DirVert_D[0];
	Dir_ID[iRay][1] = Dir_ID[iRay][1] - 2.0*DirVert_D[1];
	Dir_ID[iRay][2] = Dir_ID[iRay][2] - 2.0*DirVert_D[2];

	continue; //------------------------------------------------------>>

      } // fabs(Eps) < TolEps
    } // if (Flags_I[iRay] & PENETR) // Condition "Critical surface penetration"
    

    if (Eps < 0.0) { // Eps < 0.0: Penetration!
      //if (Eps < -TolEps) { // Eps < 0.0: Penetration!

      if (iRay ==  TrcRay) { 
	printf("PENETRATION!: Eps = %g\n", Eps);
      }

      SET_BIT(PENETR, iRay);  /* Mark "penetration" */

      //
      // Revert r and v to previous integer point
      // 
      Dir_ID[iRay][0] = DirPr_ID[iRay][0];
      Dir_ID[iRay][1] = DirPr_ID[iRay][1];
      Dir_ID[iRay][2] = DirPr_ID[iRay][2];

      Pos_ID[iRay][0] = PosPr_ID[iRay][0];
      Pos_ID[iRay][1] = PosPr_ID[iRay][1];
      Pos_ID[iRay][2] = PosPr_ID[iRay][2];

      //
      // Use PosPr_ID and DirPr_ID as memory for Cnt and 
      // for left and right values of DS and Eps
      //
      PosPr_ID[iRay][0] = 0.0;             // Store here DS_1
      PosPr_ID[iRay][1] = 2.0*DS_I[iRay];  // Store here DS_2
      PosPr_ID[iRay][2] = DS_I[iRay];      // Store here old DS
      DirPr_ID[iRay][0] = 1.0;             // Store here Eps_1
      DirPr_ID[iRay][1] = Eps;             // Store here Eps_2
      DirPr_ID[iRay][2] = CntMax;  // Store here iteration counter

      DS_I[iRay] -= Eps/GradEpsDotDir; // Newton's root approximation

      continue; //--------------------------------------------------------->>
    } // (Eps < 0.0) { // Penetration!

    //
    // Calculate main ray parameters
    //
    
    // Original Position (at an integer point): 
    PosHalfBk_D[0] = Pos_ID[iRay][0] - Dir_ID[iRay][0]*HalfDS; 
    PosHalfBk_D[1] = Pos_ID[iRay][1] - Dir_ID[iRay][1]*HalfDS;
    PosHalfBk_D[2] = Pos_ID[iRay][2] - Dir_ID[iRay][2]*HalfDS;

    GradEps_D[0] = -GradRho_ID[iRay][0]*RhoCrInv;
    GradEps_D[1] = -GradRho_ID[iRay][1]*RhoCrInv;
    GradEps_D[2] = -GradRho_ID[iRay][2]*RhoCrInv;

    GradEpsSqr = dot_product(GradEps_D, GradEps_D);

    GradEpsDotDir = GradEps_D[0]*Dir_ID[iRay][0] + GradEps_D[1]*Dir_ID[iRay][1]
      + GradEps_D[2]*Dir_ID[iRay][2];
    /* GradEpsDotDir = dot_product(GradEps_D, &Dir_ID[iRay][0]) ??? */

    LCosAl = -GradEpsDotDir/GradEpsSqr;

    DirVert_D[0] = -LCosAl*GradEps_D[0]; // Here v_||
    DirVert_D[1] = -LCosAl*GradEps_D[1]; // Here v_||
    DirVert_D[2] = -LCosAl*GradEps_D[2]; // Here v_||

    EpsHalfBk = Eps - GradEpsDotDir*HalfDS;

    /* Simplification in Curv calculation: 
     * (a x b)^2 = a^2*b^2 - (a dot b)^2
     * if |v|=1, (a x v)^2 = a^2 - (a dot v)^2 */
    Curv = pow(0.5*HalfDS/EpsHalfBk,2) * (GradEpsSqr - pow(GradEpsDotDir,2));

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


    
    if (BIT_IS_OFF(SHARP, iRay)) { /* If the ray is NOT sharp */
    
      //
      // Check if the trajectory curvature is too sharp to meet the Tolerance
      // If so, reduce the DS step  for the iRay-th ray and leave it until
      // the next call.
      //
      //
      if (BIT_IS_ON(WASREFL, iRay)) {
	CLEAR_BIT(WASREFL, iRay);  // Clear the reflection condition
	// and do not test for curvature!
      }
      else {

	if (Curv >=  Tol2) { // Testing curvature 
	  DS_I[iRay] = DS_I[iRay]/(2.0*sqrt(Curv/Tol2));
	  // r = r_hb
	  Pos_ID[iRay][0] = PosHalfBk_D[0];
	  Pos_ID[iRay][1] = PosHalfBk_D[1];
	  Pos_ID[iRay][2] = PosHalfBk_D[2];

	  if (iRay == TrcRay) {
	    printf("Curv = %g, Tol2 = %g\n", Curv, Tol2);
	  }

	  continue;    //---------------------------------------------------->>
	} // if (Curv >=  Tol2) 

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
	SET_BIT(SHARP, iRay);
	DistToCrSurf_I[iRay] = EpsHalfBk/sqrt(GradEpsSqr);
	DS_I[iRay] =  0.5*Tol*DistToCrSurf_I[iRay];
	
	if (iRay == TrcRay) {
	  printf("DistanceToCritSurf = %g, DS = %g\n", 
		 DistToCrSurf_I[iRay], DS_I[iRay]);
	}
		
	Pos_ID[iRay][0] = PosHalfBk_D[0];
	Pos_ID[iRay][1] = PosHalfBk_D[1];
	Pos_ID[iRay][2] = PosHalfBk_D[2];
	continue;    //----------------------------------------------------->>
      } /* if (GradEpsDotDir*HalfDS <= -cThird*EpsHalfBk) */
      
    } /* if (SHARP & (~Flags_I[iRay])) { // If the ray is NOT sharp */
    

    //
    // In normal mode (i.e. no penetration)
    // Save previous values of Pos and Dir before they will be changed
    //
    if (BIT_IS_OFF(PENETR, iRay)) { // i.e. if not "penetration"
      /* save previous position, r */
      PosPr_ID[iRay][0] = Pos_ID[iRay][0];
      PosPr_ID[iRay][1] = Pos_ID[iRay][1];
      PosPr_ID[iRay][2] = Pos_ID[iRay][2];
      
      /* save previous direction, v */
      DirPr_ID[iRay][0] = Dir_ID[iRay][0];
      DirPr_ID[iRay][1] = Dir_ID[iRay][1];
      DirPr_ID[iRay][2] = Dir_ID[iRay][2];
    }
 
    
    /* Solar Distance */
    SolDist = 0.0;
    for (i = 0; i < 3; i++) {
      SolDist += pow(Pos_ID[iRay][i],2);
    }
    SolDist = sqrt(SolDist);
    /* SolDist = v3magn(&Pos_ID[iRay][0]) ??? */
    /* SolDist = sqrt(pow(Pos_ID[iRay][0],2) + pow(Pos_ID[iRay][1],2)) + */
    /*                pow(Pos_ID[iRay][2],2));                           */
    

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
      // The parabolic reflection is just replacement of the Pos_ID with
      // the symmetric point at the opposite branch of parabola, and changing
      // the "incident" direction Dir_ID for the "departing" ray direction
      // according to Snell law.
      //

      if (iRay == TrcRay) {
	printf("PARABOLA SWITCHING, [%i,%i]=iRay = %i, DS = %g\n", 
	       iRay/200, iRay%200, iRay, DS_I[iRay]);
	printf("PARABOLA SWITCHING, x,y,z    = %g %g %g\n", 
	       Pos_ID[iRay][0], Pos_ID[iRay][1], Pos_ID[iRay][2]);
	printf("PARABOLA SWITCHING, Dir_ID = %g %g %g\n", 
	       Dir_ID[iRay][0], Dir_ID[iRay][1], Dir_ID[iRay][2]);
	printf("PARABOLA SWITCHING, LCosAl=%g, gradEps*Dir=%g, GradEps^2=%g\n", 
	       LCosAl, GradEpsDotDir, GradEpsSqr);
      }


      //StepY_D = Dir_ID[iRay,:] - DirVert_D       // Here c; |c| = sin \alpha
      StepY_D[0] = Dir_ID[iRay][0] - DirVert_D[0]; // Here c; |c| = sin \alpha
      StepY_D[1] = Dir_ID[iRay][1] - DirVert_D[1]; // Here c; |c| = sin \alpha
      StepY_D[2] = Dir_ID[iRay][2] - DirVert_D[2]; // Here c; |c| = sin \alpha

      //Dir_DI[iRay,:] = Dir_DI[iRay,:] - 2*DirVert_D       
      Dir_ID[iRay][0] = Dir_ID[iRay][0] - 2.0*DirVert_D[0];
      Dir_ID[iRay][1] = Dir_ID[iRay][1] - 2.0*DirVert_D[1];
      Dir_ID[iRay][2] = Dir_ID[iRay][2] - 2.0*DirVert_D[2];

      //
      //   We need to have |Step_Y| = 4 * sin(\alpha)*cos(\alpha)*EpsHalfBk*L
      //   Now the direction of Step_Y is right, 
      //   the length of it is equal to sin(\alpha)
      //   Multiply it by L*Cos(\alpha)*EpsHalfBk
      //

      StepY_D[0] = 4.0*StepY_D[0]*LCosAl*EpsHalfBk; 
      StepY_D[1] = 4.0*StepY_D[1]*LCosAl*EpsHalfBk; 
      StepY_D[2] = 4.0*StepY_D[2]*LCosAl*EpsHalfBk; 

      Pos_ID[iRay][0] = PosHalfBk_D[0] + StepY_D[0];
      Pos_ID[iRay][1] = PosHalfBk_D[1] + StepY_D[1];
      Pos_ID[iRay][2] = PosHalfBk_D[2] + StepY_D[2];
	
      //
      //   Step_X is in the direction of -grad \epsilon, whose vector is of 
      //   the length of 1/L
      //   The length of Step_X is cos^2(\alpha)*L*EpsHalfBk 
      //   Thus,
      //

      LCosAl2 = pow(LCosAl,2);
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

      StepXdbl_D[0] = 2.0*StepX_D[0];
      StepXdbl_D[1] = 2.0*StepX_D[1];
      StepXdbl_D[2] = 2.0*StepX_D[2];
      StepXSqr = dot_product(StepXdbl_D, StepXdbl_D); 
      StepYSqr = dot_product(StepY_D, StepY_D); 
      LPrb = sqrt(StepXSqr + StepYSqr);  // the parabola segment length

      // In Stokes_II[iRay,0]] the brightness temperature is saved!
      // Calculate brightness temperature as integral along ray path
      // of electron temperature times exp(-tau), tau is optical depth: ??????
      // The tau is integral along ray path of the absorption coefficient kappa
      //
      if (SolDist > r_chro_cor) {
	Te = Te_corona_K;
	//xi = 2.0; // cm^6 K^1.5 Hz^2 -- from Sheridan
	xi = 0.0765766935; // cm^6 K^1.5 Hz^2 -- from Smerd
      }
      else {
	Te = Te_chromo_K;
	//xi = 1.4; // cm^6 K^1.5 Hz^2 -- from Sheridan
	xi = 0.1148650403; // cm^6 K^1.5 Hz^2 -- from Smerd
      }
      /* Ne =  Rho_I[iRay]/ProtonMass_g;  // Electron density */
      Ne =  Rho_I[iRay]*ProtonMassInv;  // Electron density
      CAbsorp = xi*pow(Ne,2)/(pow(Te,1.5)*pow(Freq,2)*sqrt(Eps));
      dOpDepth = CAbsorp*LPrb*Rsun_cm;
      OpDepth_I[iRay] += dOpDepth;
      // Actually, Stokes_II here holds the  Tb (brightness temperature)
      dTb =  Te*exp(-OpDepth_I[iRay])*dOpDepth;
      Stokes_II[iRay][0] += dTb; 
 
      if (iRay ==  TrcRay) { 
	printf("PARABOLA: r = ");
	printf("%g %g %g\n", Pos_ID[iRay][0], Pos_ID[iRay][1], 
	       Pos_ID[iRay][2]);
	printf("PARABOLA: v = ");
	printf("%g %g %g\n", Dir_ID[iRay][0], Dir_ID[iRay][1], 
	       Dir_ID[iRay][2]);
	printf("PARA INTEGR: Te = %g, Ne = %g\n", Te, Ne);
	printf("PARA INTEGR: CAbsorp = %g, dOpDepth = %g, OpDepth_I = %g\n", 
	       CAbsorp, dOpDepth, OpDepth_I[iRay]);
	printf("PARA INTEGR: exp(-tau) = %g, dTb = %g, Tb = %g\n", 
	       exp(-OpDepth_I[iRay]), dTb, Stokes_II[iRay][0]);
      }



    } // if ((GradEpsDotDir*HalfDS <= -cThird*EpsHalfBk) 
      //  || (DS_I[iRay] < AbsMinStep))


    else { // Boris' step
      //
      // Make a step using Boris' algorithm
      //

      PosBefBoris_D[0] = Pos_ID[iRay][0];
      PosBefBoris_D[1] = Pos_ID[iRay][1];
      PosBefBoris_D[2] = Pos_ID[iRay][2];

      // grad(n) = grad(eps(i+1/2))/(2*eps(i+1/2))
      Coef = 0.5*HalfDS/(1.0 - Rho2RhoCr);
      RelGradRefrInx_D[0] = Coef*GradEps_D[0];
      RelGradRefrInx_D[1] = Coef*GradEps_D[1]; 
      RelGradRefrInx_D[2] = Coef*GradEps_D[2]; 

      Dir_D[0] = Dir_ID[iRay][0];
      Dir_D[1] = Dir_ID[iRay][1];
      Dir_D[2] = Dir_ID[iRay][2];

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
      Coef = 2.0/(1.0 + Curv1);
      //Dir_DI(:,iRay) = Dir_DI(:,iRay) + Coef*cross_product(Dir1_D, Omega_D)
      cross_product(Dir1_D, Omega_D, Xprod_D);
      Dir_ID[iRay][0] = Dir_ID[iRay][0] + Coef*Xprod_D[0];
      Dir_ID[iRay][1] = Dir_ID[iRay][1] + Coef*Xprod_D[1];
      Dir_ID[iRay][2] = Dir_ID[iRay][2] + Coef*Xprod_D[2];

      //Pos_DI(:,iRay) = Pos_DI(:,iRay) + Dir_DI(:,iRay)*HalfDS
      Pos_ID[iRay][0] += Dir_ID[iRay][0]*HalfDS;
      Pos_ID[iRay][1] += Dir_ID[iRay][1]*HalfDS;
      Pos_ID[iRay][2] += Dir_ID[iRay][2]*HalfDS;


      if (iRay == TrcRay) {
	printf("BORIS' STEP\n");
	printf("Dir = "); 
	for (i = 0; i < 3; i++) {
	  printf("%g ", Dir_ID[iRay][i]);
	}
	printf("\n"); 
	printf("Pos = "); 
	for (i = 0; i < 3; i++) {
	  printf("%g ", Pos_ID[iRay][i]);
	}
	printf("\n"); 
	printf("SolDist = %g\n", SolDist);
      }

      // In Stokes_II[iRay][0] the brightness temperature is saved!
      // Calculate brightness temperature as integral along ray path
      // of electron temperature times exp(-tau), tau is optical depth:
      // The tau is integral along ray path of the absorption coefficient kappa
      //
      //Te = (SolDist > r_chro_cor) ? Te_corona_K : Te_chromo_K;
      //z = (a > b) ? a : b;
      if (SolDist > r_chro_cor) {
	Te = Te_corona_K;
	//xi = 2.0; // cm^6 K^1.5 Hz^2 -- from Sheridan
	xi = 0.0765766935; // cm^6 K^1.5 Hz^2 -- from Smerd
      }
      else {
	Te = Te_chromo_K;
	//xi = 1.4; // cm^6 K^1.5 Hz^2 -- from Sheridan
	xi = 0.1148650403; // cm^6 K^1.5 Hz^2 -- from Smerd
      }
      /* Ne =  Rho_I[iRay]/ProtonMass_g;  // Electron density */
      Ne =  Rho_I[iRay]*ProtonMassInv;  // Electron density
      CAbsorp = xi*pow(Ne,2)/(pow(Te,1.5)*pow(Freq,2)*sqrt(Eps));
      dOpDepth = CAbsorp*DS_I[iRay]*Rsun_cm;
      OpDepth_I[iRay] += dOpDepth;
      // Actually, Intgrl_I here holds the  Tb (brightness temperature)
      //Intgrl_I[iRay] += Te*exp(-OpDepth_I[iRay])*dOpDepth;  
      dTb =  Te*exp(-OpDepth_I[iRay])*dOpDepth;
      Stokes_II[iRay][0] += dTb;  


      if (iRay ==  TrcRay) { 
	printf("BORIS: r = ");
	printf("%g %g %g\n", Pos_ID[iRay][0], Pos_ID[iRay][1], 
	       Pos_ID[iRay][2]);
	printf("BORIS: v = ");
	printf("%g %g %g\n", Dir_ID[iRay][0], Dir_ID[iRay][1], 
	       Dir_ID[iRay][2]);
	printf("BORIS INTEGR: Te = %g, Ne = %g\n", Te, Ne);
	printf("BORIS INTEGR: CAbsorp = %g, dOpDepth = %g, OpDepth_I = %g\n", 
	       CAbsorp, dOpDepth, OpDepth_I[iRay]);
	printf("BORIS INTEGR: exp(-tau) = %g, dTb = %g, Tb = %g\n", 
	       exp(-OpDepth_I[iRay]), dTb, Stokes_II[iRay][0]);

      }


      // Itensity as integral along ray path of Volumetric Emissivity 
      //   Emissivity ~ (rho/rho_cr)^2*(1/2 - rho/rho_cr)^2  (Sokolov, 2006):
      // Alternative Emissivity calculation (Smerd, 1949):
      //   Emissivity ~ (rho/rho_cr)*(1 - rho/rho_cr)
      //--> Intgrl_I[iRay] += DS_I[iRay]*Rho2RhoCr*sqrt(cOne - Rho2RhoCr);

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
    if (BIT_IS_OFF(SHARP, iRay)) { /* If the ray is NOT sharp */
      //
      // For shallow rays the DS is increased unconditionally
      //
      if (DS_New_I[iRay] > DS_I[iRay])
	DS_I[iRay] = (2.0 - DS_I[iRay]/DS_New_I[iRay])*DS_I[iRay];
      else
	DS_I[iRay] = DS_New_I[iRay];
    } /* if (BIT_IS_OFF(SHARP, iRay)) { // If the ray is NOT sharp */
    else { // Steep ray
      //
      // If the iRay-th ray is marked as sharp (i.e. "sloping" or 
      // "oblique") then the code below increases its DS only if the 
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
	CLEAR_BIT(SHARP, iRay);  /* Mark this ray NOT sharp */
	if (DS_New_I[iRay] > DS_I[iRay])
	  DS_I[iRay] = (2.0 - DS_I[iRay]/DS_New_I[iRay])*DS_I[iRay];
	else
	  DS_I[iRay] = DS_New_I[iRay];
      } // if (EpsHalfBk > DistToCrSurf_I[iRay]*sqrt(GradEpsSqr)) 
    } // else { // Sharp ray
  } // for (iRay = 1; iRay < nRay; iRay++) 
} // advance_beam_stokesIV(...)
