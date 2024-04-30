#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "raytrace.h"
#include <pthread.h>



void scatter(double Rho,
	     double GradRho_I[3],
	     double DS,
	     struct param *prm,
	     double ScatterAngle_I[3]);

/*
 * At compiletime the following preprocessor symbols 
 * in the command line as -D<symbol> can be used:
 * TBR:  compute brightness temperature, no polarization due to magnetic field
 * TBRIV: compute brightness temperature for I and V Stokes parameters,
 *        assuming polarization due to magnetic field
 * TRAJ:  Store trajectory points for specified rays
 * DEB:   Generate debugging output
 */


void advance_beam(
	int start,
	int stop,
	int nRay,
        double Param_I[],
        double Pos_ID[nRay][3], 
	double Dir_ID[nRay][3], 
	double DS_I[nRay], 
        void (*plasma_params)(),
        short Flags_I[nRay],
        double Tbr_I[nRay], 
        double TbrIV_II[nRay][2], 
        double TprIV_II[nRay][2],  /* Previous step Tb I and V values */
        double OpDepth_I[nRay],
        double Rho_I[nRay],
        double GradRho_ID[nRay][3],
        double Bfield_ID[nRay][3],
        double PosPr_ID[nRay][3],
        double DirPr_ID[nRay][3],
        double DS_New_I[nRay],
        double DistToCrSurf_I[nRay],
	int rtmode) {
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
  // external subprogram, plasma_params(), provides the plasma Rho_I 
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
  // plasma_params():  external subroutine that returns the plasma density
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
  // 2011-Apr-01 Introduced dynamic linking for plasma_params() in 
  //             rtcoremodule.c. The pointer to plasma_params() is passed 
  //             as a parameter.
  //             

  /*
   * Revision Author:
   *   Mark Benjamin
   * 2011-Jun-24 This file was modified to include the necessary parameters
   * for threading. Two arguments were added, start and stop so that the number
   * of rays to be processed could be split up. All of the for loops that
   * used to go from 0 to nRay now go from start to stop.
   */


  /* Access parameters by structure field names */
  struct param *prm = (void *) Param_I; 

  double const static cThird = 1.0/3.0;

  int iRay;
  double  Dir1_D[3], Omega_D[3];	
  double  DirVert_D[3];	
  double  StepX_D[3], StepY_D[3], RelGradRefrInx_D[3];
  double  GradEps_D[3], PosHalfBk_D[3];
  double  Dir_D[3], Xprod_D[3];
  double  StepXSqr, StepYSqr, StepXdbl_D[3];
  double Te, dTb;
  double Ne, CAbsorp, dtau;
  double xi;   // 1.4 in Chromosphere, 2.0 in Corona,  cm^6 K^1.5 Hz^2

  double HalfDS;            // DS halved
  double PosBefBoris_D[3];  // Position vector before Boris step 
  double Eps, EpsHalfBk, Rho2RhoCr; //, Rho2RhoCr1, Rho2RhoCr2;
  double Coef, Curv, Curv1;
  double LCosAl; // L is inverse grad of \epsilon, Alpha is incidence angle
  double LCosAl2;    // = LCosAl^2
  double GradEps2, GradEpsDotDir;
  double LPrb;
  double SolDist;
  /* Start - Variables for polarization calculations */
  double Bfield_abs; /* Magnitude of magnetic field vector */
  double Bfield_abs2; /* Magnitude of magnetic field vector squared*/
  double *pBfield;    /* Pointer at current B vector */
  double Bdir_D[3];  /* Magnetic field direction vector b */
  double cos_th;  /* Cosine of angle between magnetic field and ray */
  double cos2_th, sin2_th;  /* Functions of (v,b) angle */
  double nu_eff; /* Effective electron-neutral collision frequency, rad/s */
  double u_plasma, v_plasma, w_plasma; /* Plasma parameters u, v, w */
  double P_plasma, Q_plasma, R_plasma, V_plasma; /* Plasma parameters */
  double one_minus_u; /* 1 - u */
  double one_minus_u2_inv; /* 1/(1 - u)^2 */
  double alpha;            /* Absorption coefficient */
  double dze;              /* dze^2 = beta^2 + gamma^2 + delta^2 */
  double chi;              /* eta = dze/alpha */
  //-------------------------
  double Tdif, Tsum, Ted, Tes, Term1, Term2;
  double chim, chip, chim2;
  //---------------------------
  //double dzeds;            /* dze*DS in cm */
  //double cosh_dzeds, sinh_dzeds;   /* sh(dze*DS), ch(dze*DS) */
  //double exp_alds;                 /* exp(alpha*DS) */
  double I_pr, V_pr;  /* Previous-step I and V Stokes in terms of Tb */
  double Te32_inv;    /* Te^(-3/2) */
  double DS_cm;    /* Te^(-3/2) */
  double Mcoef;    /* k0*w / 2sqrt(1 - V) */
  /* End - Variables for polarization calculations */


  prm->callcount += 1.0; /* Count calls to advance_beam() */

  /*
   * Advance all the rays by 1/2 of the DS step
   */
  for (iRay=start; iRay < stop; iRay++) {

    double sdist = v3magn(Pos_ID[iRay]);
    //if(sdist<5)
    scatter(Rho_I[iRay],GradRho_ID[iRay],DS_I[iRay],prm,Dir_ID[iRay]);

   // Do not process the rays that are done or bad
    if (BIT_IS_ON(INACTIVE, iRay)) continue;     //------------------------>>

    // Advance r by 1/2 DS 
    HalfDS = 0.5*DS_I[iRay];

  


    // Pos_ID is moved by 1/2 DS in the Dir_DI direction 
    Pos_ID[iRay][0] += Dir_ID[iRay][0]*HalfDS; 
    Pos_ID[iRay][1] += Dir_ID[iRay][1]*HalfDS;
    Pos_ID[iRay][2] += Dir_ID[iRay][2]*HalfDS;

 
  }
  
  

  
  plasma_params(start,
		stop,
		nRay, 
		Param_I, 
		Pos_ID, 
		Rho_I, 
		GradRho_ID, 
		Bfield_ID,
		DS_New_I, 
		Flags_I,
		rtmode);
  



  for (iRay = start; iRay < stop; iRay++) {

    /*if(rtmode == 3){*/
    /* if (iRay == 1225) { */
    /* 	printf("Call # %g ============================================\n", */
    /* 	       prm->callcount); */
    /* 	printf("          Flag = %d\n", Flags_I[iRay]); */
    /* } */
    /*   /\* if (iRay == 1225) *\/ */
    /*   /\* 	printf("TbrIV_II[0] = %g, TbrIV[1] = %g Entry\n",  *\/ */
    /*   /\* 	       TbrIV_II[iRay][0], TbrIV_II[iRay][1]);	 *\/ */
    /*}*/

    // Do not process the rays that are done or bad
    if (BIT_IS_ON(INACTIVE, iRay)) continue;     //------------------------>>

    
    HalfDS = 0.5*DS_I[iRay];
    Rho2RhoCr = Rho_I[iRay]*prm->RhoCrInv;


    Eps = 1.0 - Rho2RhoCr;  // Dielectric permittivity
    

    //GradEps_D = -GradRho_ID(:,iRay)*RhoCrInv


    GradEps_D[0] = -GradRho_ID[iRay][0]*prm->RhoCrInv;
  
    GradEps_D[1] = -GradRho_ID[iRay][1]*prm->RhoCrInv;
    GradEps_D[2] = -GradRho_ID[iRay][2]*prm->RhoCrInv;




    GradEps2 = dot_product(GradEps_D, GradEps_D);
    

    //GradEpsDotDir = sum(GradEps_D*Dir_DI(:,iRay))

    GradEpsDotDir = dot_product(GradEps_D, &Dir_ID[iRay][0]);

  
    if (BIT_IS_ON(PENETR, iRay)) { // Condition "Critical surface penetration"

      if (fabs(Eps) > prm->TolEps) { // Another convergence step needed

	//
	// Restore the r shifted along v by DS/2 before the call to 
	//     plasma_params():
	// r = r + v*DS/2
	//
	Pos_ID[iRay][0] -= Dir_ID[iRay][0]*HalfDS; 
	Pos_ID[iRay][1] -= Dir_ID[iRay][1]*HalfDS;
	Pos_ID[iRay][2] -= Dir_ID[iRay][2]*HalfDS;

	//
	// Use the bisection method to reach the critical surface
	//
       
	// Stop working on the ray if count exceeded maximum
	DirPr_ID[iRay][2] -= 1.0;      // Cnt--;
	if (DirPr_ID[iRay][2] < 0.0) { // if max count exceeded,
	  SET_BIT(BADRAY|INACTIVE, iRay);  // Do not process the ray anymore

	  continue; //----------------------------------------------------->>
	}
       
	if (Eps > 0.0) {
	  // DS_1 = DS; Eps_1 = Eps
	  PosPr_ID[iRay][0] = DS_I[iRay];      // DS_1 = DS
	  DirPr_ID[iRay][0] = Eps;             // Eps_1 = Eps

	}
	else { // Eps <= 0.0:
	  // Newton's method is used to improve RHS DS point
	  // DS = DS - Eps/GradEpsDotDir
	  DS_I[iRay] -= Eps/GradEpsDotDir;
	  // DS_2 = DS; Eps_2 = Eps
	  PosPr_ID[iRay][1] = DS_I[iRay]; // DS_2 = DS
	  DirPr_ID[iRay][1] = Eps;        // Eps_2 = Eps

	}
       
	// DS = (DS_1 + DS_2)/2 
	DS_I[iRay] = (PosPr_ID[iRay][0] + PosPr_ID[iRay][1])*0.5;


	continue; //------------------------------------------------------->>
       
      } // if (fabs(Eps) > prm->TolEps) 
     
      else { // fabs(Eps) <= prm->TolEps: REFLECTION
	// The critical surface found. Reflecting the ray.
       
	CLEAR_BIT(PENETR, iRay);   // Clear the penetration condition
	SET_BIT(WASREFL, iRay);  // Mark the ray as "was Snell reflected"

	//DS_I[iRay] = PosPr_ID[iRay][2];      // Restore old DS !Not needed!
	HalfDS = 0.5*DS_I[iRay];
	
	GradEps2 = dot_product(GradEps_D, GradEps_D);
	LCosAl = -GradEpsDotDir/GradEps2;

  
	DirVert_D[0] = -LCosAl*GradEps_D[0]; // Here v_||
	DirVert_D[1] = -LCosAl*GradEps_D[1]; // Here v_||
	DirVert_D[2] = -LCosAl*GradEps_D[2]; // Here v_||

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

      SET_BIT(PENETR, iRay); /* Mark the ray as penetrated past crit. surf. */

      //
      // Revert r and v to previous integer point (i.e. before penetration)
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
      DirPr_ID[iRay][2] = prm->CntMax;     // Store here iteration counter

    
      DS_I[iRay] -= Eps/GradEpsDotDir; // Newton's root 1st approximation


      continue; //--------------------------------------------------------->>
    } // (Eps < 0.0) { // Penetration!



    //
    // Calculate main ray parameters
    //
    
    // Restore the original Position (at an integer point): 
    PosHalfBk_D[0] = Pos_ID[iRay][0] - Dir_ID[iRay][0]*HalfDS; 
    PosHalfBk_D[1] = Pos_ID[iRay][1] - Dir_ID[iRay][1]*HalfDS;
    PosHalfBk_D[2] = Pos_ID[iRay][2] - Dir_ID[iRay][2]*HalfDS;


    //GradEps_D = -GradRho_DI(:,iRay)*RhoCrInv
    GradEps_D[0] = -GradRho_ID[iRay][0]*prm->RhoCrInv;
    GradEps_D[1] = -GradRho_ID[iRay][1]*prm->RhoCrInv;
    GradEps_D[2] = -GradRho_ID[iRay][2]*prm->RhoCrInv;






    GradEps2 = dot_product(GradEps_D, GradEps_D);
    GradEpsDotDir = dot_product(GradEps_D, &Dir_ID[iRay][0]);

    LCosAl = -GradEpsDotDir/GradEps2;




    DirVert_D[0] = -LCosAl*GradEps_D[0]; // Here v_||
    DirVert_D[1] = -LCosAl*GradEps_D[1]; // Here v_||
    DirVert_D[2] = -LCosAl*GradEps_D[2]; // Here v_||



    EpsHalfBk = Eps - GradEpsDotDir*HalfDS;

    //Curv = (cHalf*HalfDS/EpsHalfBk)**2* &
    //       (GradEps2 - GradEpsDotDir**2)
    /* Explanation of the simplification:
     *  (a x b)^2 = a^2*b^2 - (a dot b)^2
     * if |v|=1, (a x v)^2 = a^2 - (a dot v)^2 */
    Curv = pow(0.5*HalfDS/EpsHalfBk,2) * (GradEps2 - pow(GradEpsDotDir,2));

 

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
	
	if (Curv >=  prm->Tol2) { // If curvature is too high, decrease step 
	  DS_I[iRay] = DS_I[iRay]/(2.0*sqrt(Curv/(prm->Tol2)));
	  // r = r_hb
	  Pos_ID[iRay][0] = PosHalfBk_D[0];
	  Pos_ID[iRay][1] = PosHalfBk_D[1];
	  Pos_ID[iRay][2] = PosHalfBk_D[2];
	  continue;    //---------------------------------------------------->>
	} // if (Curv >=  prm->Tol2) 

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
	SET_BIT(SHARP, iRay); /* Mark this ray "Sharp" */
	DistToCrSurf_I[iRay] = EpsHalfBk/sqrt(GradEps2);
	DS_I[iRay] =  0.5*(prm->Tol)*DistToCrSurf_I[iRay];
	
	//Pos_DI(:,iRay) = PosHalfBk_D
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
    if (BIT_IS_OFF(PENETR, iRay)) { /* i.e. if not "penetration" */
      /* Save Pos */
      PosPr_ID[iRay][0] = Pos_ID[iRay][0];
      PosPr_ID[iRay][1] = Pos_ID[iRay][1];
      PosPr_ID[iRay][2] = Pos_ID[iRay][2];
      /* Save Dir */
      DirPr_ID[iRay][0] = Dir_ID[iRay][0];
      DirPr_ID[iRay][1] = Dir_ID[iRay][1];
      DirPr_ID[iRay][2] = Dir_ID[iRay][2];
    }
 


    SolDist = v3magn(&Pos_ID[iRay][0]); /* Solar Distance */

    //
    // Either:
    // - switch to opposite branch of parabola;
    // - or make a Boris step
    //
    if ((GradEpsDotDir*HalfDS <= -cThird*EpsHalfBk) 
	|| (DS_I[iRay] < (prm->AbsMinStep))) {

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

      StepY_D[0] = Dir_ID[iRay][0] - DirVert_D[0]; // Here c; |c| = sin \alpha
      StepY_D[1] = Dir_ID[iRay][1] - DirVert_D[1]; // Here c; |c| = sin \alpha
      StepY_D[2] = Dir_ID[iRay][2] - DirVert_D[2]; // Here c; |c| = sin \alpha
      


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

      StepXdbl_D[0] = 2.0*StepX_D[0];
      StepXdbl_D[1] = 2.0*StepX_D[1];
      StepXdbl_D[2] = 2.0*StepX_D[2];
      StepXSqr = dot_product(StepXdbl_D, StepXdbl_D); 
      StepYSqr = dot_product(StepY_D, StepY_D); 
      LPrb = sqrt(StepXSqr + StepYSqr);  // the parabola segment length

      //Intgrl_I(iRay) = Intgrl_I(iRay) 
      //                  + LPrb*(Rho2RhoCr**2)*(cHalf - Rho2RhoCr)**2
      ////-->Intgrl_I[iRay] += LPrb*pow(Rho2RhoCr,2)*pow(cHalf - Rho2RhoCr,2);

      if(rtmode==2 || rtmode==3){
	/* Electron temperature Te */ 
	if (SolDist > (prm->r_chro_cor)) {
	  Te = prm->Te_corona_K;
	  //xi = 2.0; // cm^6 K^1.5 Hz^2 -- from Sheridan
	  xi = 0.0765766935; // cm^6 K^1.5 Hz^2 -- from Smerd
	}
	else {
	  Te = prm->Te_chromo_K;
	  //xi = 1.4; // cm^6 K^1.5 Hz^2 -- from Sheridan
	  xi = 0.1148650403; // cm^6 K^1.5 Hz^2 -- from Smerd
	}
	Te32_inv = pow(Te,-1.5);        /* Te^(-3/2) to calc nu_eff */
	/* Electron density Ne */
	Ne =  Rho_I[iRay]*(prm->ProtonMassInv);  // Electron density
	/* Optical depth, calculated as integral of the absorption 
	 * coefficient CAbsorp over parabolic segment LPrb. */
	CAbsorp = xi*pow(Ne,2)/(pow(Te,1.5)*(prm->Freq2)*sqrt(Eps));
	dtau = CAbsorp*LPrb*(prm->Rsun_cm);
	OpDepth_I[iRay] += dtau;
	//dTb =  Te*exp(-OpDepth_I[iRay])*dOpDepth;
      }
      if(rtmode==2){
	/* Calculate the brightness temperature, Tbr, as integral 
	 *    Tbr = \int_{\tau_1}^{\tau_2}{ T_e e^{-\tau} } d \tau,
	 * where \tau is the optical depth */
	//dTb = (Te - Tbr_I[iRay])*dtau;
	//Tbr_I[iRay] += dTb; 
	Tbr_I[iRay] = (Tbr_I[iRay] - Te)*exp(-dtau) + Te;
      }else if(rtmode==3){

	/* Calculate Stokes I and V components of brightness temperature 
	 * TbrIV_II[iRay][0] is I Stokes component (intensity)
	 * TbrIV_II[iRay][1] is V Stokes component (circular polarization)
	 */
	pBfield = &Bfield_ID[iRay][0];
	Bfield_abs2 = dot_product(pBfield, pBfield); /* = |Bfield|^2 */
	Bfield_abs = sqrt(Bfield_abs2);              /* Bfield_abs = |Bfield| */
	Bdir_D[0] = Bfield_ID[iRay][0]/Bfield_abs;
	Bdir_D[1] = Bfield_ID[iRay][1]/Bfield_abs;
	Bdir_D[2] = Bfield_ID[iRay][2]/Bfield_abs;
	cos_th = dot_product(&Dir_ID[iRay][0], Bdir_D); /* Angle b/w ray and B */
	cos2_th = pow(cos_th,2);
	sin2_th = sqrt(1.0 - cos2_th);
	/* if (iRay == 1225) { */
	/* 	printf("Bfield = %g, %g, %g\n",  */
	/* 	       Bfield_ID[iRay][0], Bfield_ID[iRay][1], 
		       Bfield_ID[iRay][2]); */
	/* 	printf("Bdir = %g, %g, %g\n",  */
	/* 	       Bdir_D[0], Bdir_D[1], Bdir_D[2]); */
	/* 	printf("Bfield_abs2 = %g, Bfield_abs = %g\n", 
		Bfield_abs2, Bfield_abs); */
	/* } */
	//nu_eff = prm->lnLambda_13p7*Ne*Te32_inv; /* Eff. collision freq. */
	//nu_eff = 0.1*prm->lnLambda_13p7*Ne*Te32_inv; /* Eff. collision freq. */
	//nu_eff = prm->lnLambda_13p7*Ne*Te32_inv; /* by Melrose */
	nu_eff = prm->Cnu*Ne*Te32_inv*log(220.0*Te*pow(Ne,-1./3.)); /*Ginzburg*/
	//nu_eff = 3.0*Ne*Te32_inv*log(220.0*Te*pow(Ne,-1./3.)); /*by Ginzburg'*/
	u_plasma = prm->e_ovr_mcw2*Bfield_abs2;        /* u = (e*B/mcw)^2 */
	v_plasma = prm->e2_4pi_ovr_mw2*Ne;             /* v = 4pie^2Ne/mw^2 */
	w_plasma = nu_eff/prm->Omega;
	one_minus_u = 1.0 - u_plasma;
	one_minus_u2_inv = pow(one_minus_u,-2);  /* 1/(1-u)^2 */
	P_plasma = v_plasma*(1.0 + u_plasma)*one_minus_u2_inv;
	Q_plasma = u_plasma*v_plasma*(3.0 - u_plasma)*one_minus_u2_inv;
	R_plasma = 2.0*v_plasma*sqrt(u_plasma)*one_minus_u2_inv;
	V_plasma = v_plasma/one_minus_u;
	Mcoef = 0.5*prm->k0*w_plasma/sqrt(1.0 - V_plasma);
	/* alpha = - Mcoef*(2.0*P_plasma - Q_plasma*sin2_th); */
	alpha = Mcoef*(2.0*P_plasma - Q_plasma*sin2_th);
	dtau = alpha*DS_cm;
	dze = Mcoef*sqrt(pow(Q_plasma*sin2_th,2) + 
			 4.0*pow(R_plasma*cos_th,2));
	//if (cos_th > 0.0) dze = -dze;
	if (cos_th < 0.0) dze = -dze;
	chi = dze/alpha;
	chim = 1.0 - chi;
	chip = 1.0 + chi;
	chim2 = chim*chip;  /* 1 - chi^2 */
	Tdif = TbrIV_II[iRay][0] - TbrIV_II[iRay][1];
	Tsum = TbrIV_II[iRay][0] + TbrIV_II[iRay][1];
	Ted = Te - Tdif*(1.0 - chi);
	Tes = Te - Tsum*(1.0 + chi);
	Term1 = 0.5*exp(-chim*dtau)*Ted/chim;
	Term2 = 0.5*exp(-chip*dtau)*Tes/chip;
	TbrIV_II[iRay][0] = -Term1 - Term2 + Te/chim2;
	TbrIV_II[iRay][1] =  Term1 - Term2 - Te*chi/chim2;
	/* dzeds = dze*DS_cm; */
	/* cosh_dzeds = cosh(dze*DS_cm); */
	/* sinh_dzeds = sinh(dze*DS_cm); */
	/* exp_alds = exp(alpha*DS_cm); */

	/* if (iRay == 1225){ */
	/* 	printf("exp_alds = %g, cosh_dzeds = %g, sinh_dzeds = %g\n", */
	/* 	       exp_alds, cosh_dzeds, sinh_dzeds); */
	/* 	printf("TbrIV_II[0] = %g, TbrIV[1] = %g BEFORE\n",  */
	/* 	       TbrIV_II[iRay][0], TbrIV_II[iRay][1]); */
	/* 	printf("I_pr = %g, V_pr = %g\n", I_pr, V_pr); */
	/* } */
	//TbrIV_II[iRay][0] = exp_alds*(I_pr*cosh_dzeds + V_pr*sinh_dzeds);
	//TbrIV_II[iRay][1] = exp_alds*(I_pr*sinh_dzeds + V_pr*cosh_dzeds);
	/* if (iRay == 1225) */
	/* 	printf("TbrIV_II[0] = %g, TbrIV[1] = %g AFTER\n",  */
	/* 	       TbrIV_II[iRay][0], TbrIV_II[iRay][1]);	 */
	//TprIV_II[iRay][0] =  I_pr; /* Save previous I value */
	//TprIV_II[iRay][1] =  V_pr; /* Save previous V value */
	/* if (iRay == 1225){ */
	/* 	printf("SolDist = %g\n", SolDist); */
	/* 	printf("dTb = %g, k0 w/2sqrt(1-V) = %g\n", dTb, Mcoef); */
	/* 	printf("TbrI = %g, TbrV = %g\n", TbrIV_II[iRay][0], 
		TbrIV_II[iRay][1]); */
	/* 	printf("alpha = %g, dze = %g, DS = %g, DS_cm = %g\n",  */
	/* 	       alpha, dze, DS_I[iRay], DS_cm); */
	/* 	printf("al*DS = %g, dze*DS = %g\n", alpha*DS_cm, dze*DS_cm); */
	/* 	printf("|B|=%g, cos_th=%g, nu_eff=%g\n", Bfield_abs, 
		cos_th, nu_eff); */
	/* 	printf("u = %g, v = %g, w = %g\n", u_plasma, v_plasma, 
		w_plasma); */
	/* 	printf("P = %g, Q = %g, R = %g, V = %g\n",  */
	/* 	       P_plasma, Q_plasma, R_plasma, V_plasma); */
	/* 	printf("dTb = %g\n", dTb); */
	/* } */


      }
      
    } // if ((GradEpsDotDir*HalfDS <= -cThird*EpsHalfBk) 
      //  || (DS_I[iRay] < (prm->AbsMinStep)))


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

      if(rtmode==2 || rtmode==3){
	/* Electron temperature Te */
	if (SolDist > (prm->r_chro_cor)) {
	  Te = prm->Te_corona_K;
	  //xi = 2.0; // cm^6 K^1.5 Hz^2 -- from Sheridan
	  xi = 0.0765766935; // cm^6 K^1.5 Hz^2 -- from Smerd
	}
	else {
	  Te = prm->Te_chromo_K;
	  //xi = 1.4; // cm^6 K^1.5 Hz^2 -- from Sheridan
	  xi = 0.1148650403; // cm^6 K^1.5 Hz^2 -- from Smerd
	}
	Te32_inv = pow(Te,-1.5);        /* Te^(-3/2) to calc nu_eff */
	/* Electron density Ne */
	Ne =  Rho_I[iRay]*(prm->ProtonMassInv);  // Electron density
	/* Optical depth, calculated as integral of the absorption 
	 * coefficient CAbsorp over the ray path element DS */
	CAbsorp = xi*pow(Ne,2)/(pow(Te,1.5)*(prm->Freq2)*sqrt(Eps));
	DS_cm = DS_I[iRay]*(prm->Rsun_cm);
	dtau = CAbsorp*DS_cm;
	OpDepth_I[iRay] += dtau;
	//dTb =  Te*exp(-OpDepth_I[iRay])*dOpDepth; /* Unpolarized increment */
	//dTb = (Te - Tbr_I[iRay])*DS_cm;
      }
      if(rtmode==2){
	/* Calculate the brightness temperature, Tbr, as integral 
	 *    Tbr = \int_{\tau_1}^{\tau_2}{ T_e e^{-\tau} } d \tau,
	 * where \tau is the optical depth */
	/* if (iRay == 53) { */
	/* 	printf("dTb = %g, Tbr_I[iRay] = %g\n", dTb, Tbr_I[iRay]); */
	/* } */
	Tbr_I[iRay] = (Tbr_I[iRay] - Te)*exp(-dtau) + Te;
 
      }
      if(rtmode==3){
	/* Calculate Stokes I and V components of brightness temperature 
	 * TbrIV_II[iRay][0] is I Stokes component (intensity)
	 * TbrIV_II[iRay][1] is V Stokes component (circular polarization)
	 */
	pBfield = &Bfield_ID[iRay][0];
	Bfield_abs2 = dot_product(pBfield, pBfield); /* = |Bfield|^2 */
	Bfield_abs = sqrt(Bfield_abs2);              /* Bfield_abs = |Bfield| */
	Bdir_D[0] = Bfield_ID[iRay][0]/Bfield_abs;
	Bdir_D[1] = Bfield_ID[iRay][1]/Bfield_abs;
	Bdir_D[2] = Bfield_ID[iRay][2]/Bfield_abs;
	cos_th = dot_product(&Dir_ID[iRay][0], Bdir_D); /* Angle b/w ray and B */
	cos2_th = pow(cos_th,2);
	sin2_th = sqrt(1.0 - cos2_th);
	/* if (iRay == 1225) { */
	/* 	printf("Bfield = %g, %g, %g\n",  */
	/* 	       Bfield_ID[iRay][0], Bfield_ID[iRay][1], 
		       Bfield_ID[iRay][2]); */
	/* 	printf("Bdir = %g, %g, %g\n",  */
	/* 	       Bdir_D[0], Bdir_D[1], Bdir_D[2]); */
	/* 	printf("Bfield_abs2 = %g, Bfield_abs = %g\n", 
		Bfield_abs2, Bfield_abs); */
	/* } */
	//nu_eff = prm->lnLambda_13p7*Ne*Te32_inv; /* Eff. collision freq. */
	//nu_eff = 0.1*prm->lnLambda_13p7*Ne*Te32_inv; /* Eff. collision freq. */
	//nu_eff = prm->lnLambda_13p7*Ne*Te32_inv; /* by Melrose */
	nu_eff = prm->Cnu*Ne*Te32_inv*log(220.0*Te*pow(Ne,-1./3.)); /*Ginzburg*/
	//nu_eff = 3.0*Ne*Te32_inv*log(220.0*Te*pow(Ne,-1./3.)); /*by Ginzburg'*/
	u_plasma = prm->e_ovr_mcw2*Bfield_abs2;        /* u = (e*B/mcw)^2 */
	v_plasma = prm->e2_4pi_ovr_mw2*Ne;             /* v = 4pie^2Ne/mw^2 */
	w_plasma = nu_eff/prm->Omega;
	one_minus_u = 1.0 - u_plasma;
	one_minus_u2_inv = pow(one_minus_u,-2);  /* 1/(1-u)^2 */
	P_plasma = v_plasma*(1.0 + u_plasma)*one_minus_u2_inv;
	Q_plasma = u_plasma*v_plasma*(3.0 - u_plasma)*one_minus_u2_inv;
	R_plasma = 2.0*v_plasma*sqrt(u_plasma)*one_minus_u2_inv;
	V_plasma = v_plasma/one_minus_u;
	Mcoef = 0.5*prm->k0*w_plasma/sqrt(1.0 - V_plasma);
	/* alpha = - Mcoef*(2.0*P_plasma - Q_plasma*sin2_th); */
	alpha = Mcoef*(2.0*P_plasma - Q_plasma*sin2_th);

	dtau = alpha*DS_cm;
	dze = Mcoef*sqrt(pow(Q_plasma*sin2_th,2) + 
			 4.0*pow(R_plasma*cos_th,2));
	//if(isnan(alpha) || isnan(dze))
	//  fprintf(stderr,"IRay %d Alpha: %g Dze: %g\n",iRay,alpha,dze);
	//if (cos_th > 0.0) dze = -dze;
	if (cos_th < 0.0) dze = -dze;
	chi = dze/alpha;
	chi = dze/alpha;
	chim = 1.0 - chi;
	chip = 1.0 + chi;
	chim2 = chim*chip;  /* 1 - chi^2 */
	Tdif = TbrIV_II[iRay][0] - TbrIV_II[iRay][1];
	Tsum = TbrIV_II[iRay][0] + TbrIV_II[iRay][1];
	Ted = Te - Tdif*(1.0 - chi);
	Tes = Te - Tsum*(1.0 + chi);
	Term1 = 0.5*exp(-chim*dtau)*Ted/chim;
	Term2 = 0.5*exp(-chip*dtau)*Tes/chip;
	TbrIV_II[iRay][0] = -Term1 - Term2 + Te/chim2;
	TbrIV_II[iRay][1] =  Term1 - Term2 - Te*chi/chim2;
	
	/* dzeds = dze*DS_cm; */
	/* cosh_dzeds = cosh(dze*DS_cm); */
	/* sinh_dzeds = sinh(dze*DS_cm); */
	/* exp_alds = exp(alpha*DS_cm); */

	/* if (iRay == 1225){ */
	/* 	printf("exp_alds = %g, cosh_dzeds = %g, sinh_dzeds = %g\n", */
	/* 	       exp_alds, cosh_dzeds, sinh_dzeds); */
	/* 	printf("TbrIV_II[0] = %g, TbrIV[1] = %g BEFORE\n",  */
	/* 	       TbrIV_II[iRay][0], TbrIV_II[iRay][1]); */
	/* 	printf("I_pr = %g, V_pr = %g\n", I_pr, V_pr); */
	/* } */
	//TbrIV_II[iRay][0] = exp_alds*(I_pr*cosh_dzeds + V_pr*sinh_dzeds);
	//TbrIV_II[iRay][1] = exp_alds*(I_pr*sinh_dzeds + V_pr*cosh_dzeds);
	/* if (iRay == 1225) */
	/* 	printf("TbrIV_II[0] = %g, TbrIV[1] = %g AFTER\n",  */
	/* 	       TbrIV_II[iRay][0], TbrIV_II[iRay][1]);	 */
	//TprIV_II[iRay][0] =  I_pr; /* Save previous I value */
	//TprIV_II[iRay][1] =  V_pr; /* Save previous V value */
	/* if (iRay == 1225){ */
	/* 	printf("SolDist = %g\n", SolDist); */
	/* 	printf("dTb = %g, k0 w/2sqrt(1-V) = %g\n", dTb, Mcoef); */
	/* 	printf("TbrI = %g, TbrV = %g\n", TbrIV_II[iRay][0], 
		TbrIV_II[iRay][1]); */
	/* 	printf("alpha = %g, dze = %g, DS = %g, DS_cm = %g\n",  */
	/* 	       alpha, dze, DS_I[iRay], DS_cm); */
	/* 	printf("al*DS = %g, dze*DS = %g\n", alpha*DS_cm, dze*DS_cm); */
	/* 	printf("|B|=%g, cos_th=%g, nu_eff=%g\n", Bfield_abs, 
		cos_th, nu_eff); */
	/* 	printf("u = %g, v = %g, w = %g\n", u_plasma, v_plasma, 
		w_plasma); */
	/* 	printf("P = %g, Q = %g, R = %g, V = %g\n",  */
	/* 	       P_plasma, Q_plasma, R_plasma, V_plasma); */
	/* 	printf("dTb = %g\n", dTb); */
	/* } */
      }
      // Calculate Itensity as integral along ray path of Volumetric Emissivity 
      // Emissivity ~ (rho/rho_cr)^2*(1/2 - rho/rho_cr)^2  (Sokolov, 2006):
      //--> Intgrl_I[iRay] = Intgrl_I[iRay] + 
      //            DS_I[iRay]*pow(Rho2RhoCr,2)*pow(cHalf - Rho2RhoCr, 2);
      // Alternative Emissivity calcylation (Smerd, 1949):
      // Emissivity ~ (rho/rho_cr)*(1 - rho/rho_cr)
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
      // as steep in the DitToCrSurf_I. 
      //   This can happen in two cases: 
      // (1) either the ray was "parabola reflected" and after several steps it
      // went away from the surface by the distance where the parabola 
      // switching occurred; 
      // (2) or the ray is not steep any more because the current DS is so 
      // small, that the next step does not penetrate the critical surface.
      //   The ray is reverted to "gentle" or "shallow"
      //
      if (EpsHalfBk > DistToCrSurf_I[iRay]*sqrt(GradEps2)) {
	CLEAR_BIT(SHARP, iRay);  /* Mark this ray NOT sharp */
	if (DS_New_I[iRay] > DS_I[iRay])
	  DS_I[iRay] = (2.0 - DS_I[iRay]/DS_New_I[iRay])*DS_I[iRay];
	else
	  DS_I[iRay] = DS_New_I[iRay];
      } // if (EpsHalfBk > DistToCrSurf_I[iRay]*sqrt(GradEps2)) 
    } // else { // Sharp ray
    /*if(rtmode==3){*/
      /* if (iRay == 1225) { */			
      /* 	printf("TbrIV_II[0] = %g, TbrIV[1] = %g Exit\n",  */
      /* 	       TbrIV_II[iRay][0], TbrIV_II[iRay][1]);	 */
      /* 	printf("Exit # %g =========\n", */
      /* 	       prm->callcount); */
      /* } */
    /*}*/

  } // for (iRay = 1; iRay < nRay; iRay++) 


}

inline double normal_random(double mean, double sigma){
  double random1 = random()/(float)RAND_MAX;
  double random2 = random()/(float)RAND_MAX;  

  return sqrt(-2*pow(sigma,2)*log(random1))*cos(2*sigma*pi*random2);
}

void scatter(double Rho,
	     double GradRho_I[3],
	     double DS,
	     struct param *prm,
	     double Dir_I[3]){

  double b;
  double ep2;
  double mu;
  double plasma_freq;
  double li,l0;
  double sigma;
  double Dir_length;

  ep2 = pow(.1,2);/*((GradRho_I[0]*GradRho_I[0]+
	   GradRho_I[1]*GradRho_I[1]+
	   GradRho_I[2]*GradRho_I[2])/
	   (Rho*Rho));*/

  mu = 1.0 - Rho*prm->RhoCrInv;
  plasma_freq = sqrt(80600000*Rho);

  li = pow(.08/prm->Rsun_km,.3333333);
  l0 = pow(li*1000000,.6666666);

  b = ((pi*pow(plasma_freq,4)*ep2)/
       (pow(prm->Freq,4)*pow(mu,4)*li*l0));

  sigma = mu*sqrt(b*DS);


  //  fprintf(stderr,"Sigma %g\n", sigma);

  Dir_I[0] += normal_random(0,sigma);
  Dir_I[1] += normal_random(0,sigma);
  Dir_I[2] += normal_random(0,sigma);

  Dir_length = v3magn(Dir_I);

  //fprintf(stdout,"Length %g Before: %g,%g,%g\n",Dir_length,Dir_I[0],Dir_I[1],Dir_I[2]);

  Dir_I[0] = Dir_I[0]/Dir_length;
  Dir_I[1] = Dir_I[1]/Dir_length;
  Dir_I[2] = Dir_I[2]/Dir_length;

  //fprintf(stdout,"After: %g,%g,%g\n",Dir_I[0],Dir_I[1],Dir_I[2]);


}


