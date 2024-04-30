#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "raytrace.h"


//
// Calculates plasma density, Density_I, and its gradient, 
// GradDensity_DI(3,nRay), at specified locations Position_DI(3,nRay)
// Also, it provides appropriate step, DeltaSNew_I, conceivably dependent
// on the numeric grid size
//
void Get_Plasma_Density(int nRay, double Position_DI[3][nRay], double Density_I[nRay], double GradDensity_DI[3][nRay],
			double DeltaSNew_I[nRay], short RayFlag_I[nRay]);

void beam_path(void (* Get_Plasma_Density)(), int nRay, short ExcludeRay_I[nRay], double Position_DI[3][nRay],
	       double Slope_DI[3][nRay], double DeltaS_I[nRay], double ToleranceInit, double DensityCr, 
	       double Intensity_I[nRay], short RayFlag_I[nRay], short *NewEntry) 
{
  //
  //   The subroutine beam_path() makes raytracing and emissivity integration along ray paths.
  // It works on a group of rays (a beam). Each ray is specified by its Cartesian Position_DI
  // and its direction cosines Slope_DI. The subroutine calculates new Position_DI, which is 
  // DeltaS_I away from the initial position. It calculates the intensity along the step
  // as emissivity by DeltaS and adds this to the Intensity_I. Thus, after series of calls
  // to beam_path(), the Intensity_I contains the result of integration along the paths of 
  // every ray. 
  //   The electromagnetic rays are refracted in plasma due to its non-zero refractive index,
  // which is the square root of the diellectric permittivity \epsilon. The \epsilon, in turn,
  // is a function of the plasma density. The external subprogram, Get_Plasma_Density(),
  // provides the plasma Density_I along with its gradient, GradDensity_DI at the Position_DI.
  // The \epsilon can be calculated as 1 - Density_I/DensityCr, where DensityCr is the "critical"
  // plasma density at which the dielectric permittivity \epsilon falls down to zero. 
  // The value of DensityCr is proportional to the square of the wave frequency. For example,
  // for the 42.3 MHz radiowaves the critical density is ~3.71x10^(-17) g/cm^3. 
  // A surface where the plasma density achieves the critical value acts like a mirror. 
  // No radiation penetrates the critical surface. The radiowaves can only travel in the regions
  // with lower density.
  //   The emissivity w of plasma at a selected plasma frequency is calculated as a polynomial
  // function w = (Dens2DensCr)^2*(0.5 - Dens2DensCr)^2, where Dens2DensCr is the quotient
  // Density_I/DensityCr. So, here the plasma density is used for calculation of both emissivity
  // and dielectric permittyvity.
  //   The parameters of the beam_path() are briefly described below.
  //
  // Get_Plasma_Density():  external subroutine that returns the plasma Density_I and its gradient
  //     GradDensity_DI. It also provides the recommended step size DeltaSNew_I and asserts
  //     the RayFlag_I (assigns TRUE) for a ray should it (illegally) penetrate into the region
  //     with "negative" dielectric permittivity.
  // nRay:           number of rays being processed.
  // ExcludeRay_I:  the caller program can use this logical array to stop processing of any individual 
  //     ray when it already finished travelling inside of a specified part of space. Set the 
  //     corresponding element of ExcludeRay_I to TRUE to leave it unprocessed during the subsequent
  //     calls to beam_path(). Before the first call to beam_path() all the elements of ExcludeRay_I
  //     should be set to FALSE; then all the rays will be processed.
  // Position_DI:    Cartesian position vectors of the rays.
  // Slope_DI:       Direction cosines of the rays.
  // DeltaS_I:       Current step values for the rays. Set the elements of Slope_DI to some reasonable 
  //     value, say, 1.0. The DeltaS_I elements are individually modified by 
  //     beam_path() to satisfy the precision requirements set by ToleranceInit.
  // ToleranceInit:   determines the precision of ray paths calculation. ToleranceInit is the inverse
  //     of the minimum number of ray trajectory points per one radian of the ray curvature. If this
  //     requirement is not satisfied, the corresponding element of DeltaS_I is decreased. Do not set
  //     the ToleranceInit to any value greater than 0.1 (which means 10 points per curvature radian): 
  //     it will be internally set to 0.1 anyway.   
  // DensityCr: the plasma density at which its dielectric permittivity becomes zero for chosen wave
  //     frequency.
  // Intensity_I:    the intensities of each ray calculated as the integral of emissivity along the ray 
  //     path. During each call to ray_path(), each element of the Intensity_I is incremented by the
  //     integral along the step DeltaS_I. Set all the elements of Intensity_I to 0.0 before the 
  //     first call to the beam_path().
  // RayFlag_I:      the TRUE elements of this logical array indicate that the corresponding rays
  //     penetrated into the "prohibited" region of space with the plasma density above its critical 
  //     value. Normally, it should never happen. However, in case the algorithm made such an error,
  //     the flagged rays should be considered as "bad" and thrown away from the resultant Intensity_I.
  //     Set all the elements of RayFlag_I to FALSE before calling ray_path() for the first time.
  // NewEntry:   A pointer to a short integer variable that controls the internal array allocation
  //     and deallocation. Set this variable to 1 before the first call to beam_path(). This value forces
  //     the beam_path() to allocate internal dynamic arrays and take several initial actions. 
  //     During subsequent calls to beam_path() the *NewEntry will keep the value 0, which leaves
  //     the allocated arrays intact. Setting the *NewEntry to -1 during the same run will free 
  //     the array memory, after which another allocation of the internal allocatables, with 
  //     possibly different dimensions is possible via setting the *NewEntry to 1 again. This variable
  //     was introduced out of the speed considerations: the automatic arrays work slower than the 
  //     allocatable ones, because the allocation occurs only once at the first call.

  // 
  // Author:
  //   Leonid Benkevitch
  // Revisions:
  // 2007-Jul-17 Initially written in Fortran-90 and tested for frequencies around 42.3 MHz
  // 2009-Jan-05 Rewritten in GNU C 
  // 2009-Jan-09 Added ray stopping at the solar surface for frequencies above ~127 MHz
  // 2009-Jan-13 The array stopping correction was unsuccessful, so the program reverted to 
  //             the state before 2009-Jan-09. This requires more work.
  //


  //    external Get_Plasma_Density
  //integer, intent(in) :: nRay                                  // # of pixels in the raster
  //real, intent(inout), dimension(3,nRay) :: Position_DI, Slope_DI
  //real, intent(inout), dimension(nRay) :: Intensity_I, DeltaS_I
  //real, intent(in) ::  ToleranceInit, DensityCr
  //logical, intent(inout), dimension(nRay) :: RayFlag_I      // TRUE if a ray is OK, FALSE otherwise
  //logical, intent(inout) :: NewEntry                        // Must be set to TRUE before a call with new value of nRay
  //logical, intent(inout), dimension(nRay) :: ExcludeRay_I   // A ray is excluded from processing if it is TRUE
  

  double extern const cPi;                         // = 3.1415926535897931;
  double extern const ProtonChargeSGSe;            // = 4.8e-10; // StatCoulombs, SGSe
  double extern const cProtonMass;                 // = 1.672621636E-24; //g
  double extern const cElectronMass;               // = 9.10938215E-28; //g

  short static const FALSE = 0, TRUE = 1;
  const static double cZero = 0.0, cOne = 1.0, cTwo = 2.0, cHalf = 0.5, cFour = 4.0;
  const static double cThird = 0.33333333333333333333;

  int i, j, iRay, n2Ry = 2*nRay;
  double static Slope1_D[3], Omega_D[3];	
  double static ProjSlopeOnMinusGradEps_D[3];	
  double static StepX_D[3], StepY_D[3], RelGradRefrInx_D[3];
  double static GradDielPerm_D[3], PositionHalfBack_D[3];
  double static Slope_D[3], Xprod_D[3];
  double static Tolerance, ToleranceSqr, DensityCrInv, AbsoluteMinimumStep;
  double static StepXSqr, StepYSqr, StepXdbl_D[3];

  double HalfDeltaS;                        // DeltaS halved
  double Pos_D[3];                          // Position vector for current iRay
  double PosBeforeBorisStep_D[3];           // Position vector before Boris step for current iRay
  double DielPerm, DielPermHalfBack, Dens2DensCr, Dens2DensCr1, Dens2DensCr2;
  double Coef, Curv, Curv1;
  double LCosAl;    //Where L is inverse grad of \epsilon, Alpha is the incidence angle
  double LCosAl2;    // = LCosAl^2
  double GradDielPermSqr, GradEpsDotSlope;
  double ParabLen, GradDielPerm;
  short deb = 0;
  double const static SolarRadius = 1.0;

  // Allocatable:
  double static *GradDensity_DI = 0;	// [3,nRay]
  double static *Density_I = 0, *DeltaSNew_I = 0, *DistanceToCritSurf_I = 0;	// [nRay]
  short static *GentleRay_I = 0;       // [nRay] TRUE for shallow rays; is set to FALSE for steep rays.
  //double GradDensity_DI[3][nRay];
  //double Density_I[nRay], DeltaSNew_I[nRay], DistanceToCritSurf_I[nRay];
  //short GentleRay_I[nRay];       //  TRUE for shallow rays; is set to FALSE for steep rays.


  //real, save, dimension(3)       :: Slope1_D, Omega_D
  //real, save, dimension(3)       :: ProjSlopeOnMinusGradEps_D
  //real, save, dimension(3)       :: StepX_D, StepY_D, RelGradRefrInx_D 
  //real, save, dimension(3)       :: GradDielPerm_D, PositionHalfBack_D

  //real, save, pointer, dimension(:,:) :: GradDensity_DI
  //real, save, pointer, dimension(:)   :: Density_I
  //real, save, pointer, dimension(:)   :: DeltaSNew_I
  //real, save, pointer, dimension(:)   :: DistanceToCritSurf_I
  //logical, save, pointer, dimension(:)   :: GentleRay_I       // TRUE for shallow rays; is set to FALSE for steep rays.

  //real :: HalfDeltaS                        // DeltaS halved
  //real :: DielPerm, DielPermHalfBack, Dens2DensCr, Dens2DensCr1, Dens2DensCr2
  //real :: Coef, Curv, Curv1
  //real :: LCosAl    //Where L is inverse grad of \epsilon, Alpha is the incidence angle
  //real :: GradDielPermSqr, GradEpsDotSlope
  //real :: ParabLen, GradDielPerm
  //real, save :: Tolerance, ToleranceSqr, DensityCrInv, AbsoluteMinimumStep
  //integer :: i, j, iRay 

  //
  // Initialization at first entry
  // or at 
  //
  if (*NewEntry) {
    *NewEntry = FALSE;
    if (GradDensity_DI) {free(GradDensity_DI); GradDensity_DI = 0;}
    if (Density_I) {free(GradDensity_DI); GradDensity_DI = 0;}
    if (DeltaSNew_I) {free(DeltaSNew_I); DeltaSNew_I = 0;}
    if (DistanceToCritSurf_I) {free(DistanceToCritSurf_I); DistanceToCritSurf_I = 0;}
    if (GentleRay_I) {free(GentleRay_I); GentleRay_I = 0;}
    if (*NewEntry == -1) return;    // Free memory and exit
    GradDensity_DI = (double *)calloc(3*nRay, sizeof(double)); 
    Density_I = (double *)calloc(nRay, sizeof(double));
    DeltaSNew_I = (double *)calloc(nRay, sizeof(double));
    DistanceToCritSurf_I = (double *)calloc(nRay, sizeof(double));
    GentleRay_I = (short *)calloc(nRay, sizeof(short));
    for (i = 0; i < nRay; i++) {
      DistanceToCritSurf_I[i] = cZero;
      GentleRay_I[i] = TRUE;
    } // for (i = 0; i < nRay; i++) 
    DensityCrInv = cOne/DensityCr;
    Tolerance = min(ToleranceInit,0.1);  // i.e. minimum ten points between a vacuum and a critical surface and
    ToleranceSqr = pow(Tolerance,2);     //  minimum 10 points over 1 rad of the curvature
    AbsoluteMinimumStep = 1e-4*sum(DeltaS_I, nRay)/nRay; // One ten-thousandth of average step
    //write(*,*) 'AbsoluteMinimumStep = ', AbsoluteMinimumStep
    if (deb) { 
      printf("ToleranceInit = %g, Tolerance = %g\n", ToleranceInit, Tolerance);
      printf("DensityCr = %g\n", DensityCr);
      printf("DensityCrInv = %g\n", DensityCrInv);
      //printf("POINTER GradDensity_DI = %i\n", GradDensity_DI);
      //printf("POINTER Density_I = %i\n", Density_I);
      //printf("POINTER DeltaSNew_I = %i\n", DeltaSNew_I);
      printf("Slope_DI = \n");
      print2d(3, nRay, Slope_DI);
    }

  } // if (NewEntry) 

  for (iRay=0; iRay < nRay; iRay++) {
    //Pos_D[0] = Position_DI[0][iRay];
    //Pos_D[1] = Position_DI[1][iRay];
    //Pos_D[2] = Position_DI[2][iRay];
    //if (sqrt(dot_product(Pos_D, Pos_D)) <= SolarRadius) ExcludeRay_I[iRay] = TRUE;
    if (ExcludeRay_I[iRay]) continue;        // Do not process the rays that are done; save time
    HalfDeltaS = cHalf*DeltaS_I[iRay];
    Position_DI[0][iRay] = Position_DI[0][iRay] + Slope_DI[0][iRay]*HalfDeltaS; // Now Position_DI moved by 1/2 DeltaS !!!
    Position_DI[1][iRay] = Position_DI[1][iRay] + Slope_DI[1][iRay]*HalfDeltaS; // Now Position_DI moved by 1/2 DeltaS !!!
    Position_DI[2][iRay] = Position_DI[2][iRay] + Slope_DI[2][iRay]*HalfDeltaS; // Now Position_DI moved by 1/2 DeltaS !!!
  }

  //subroutine plasma_density(Position_DI, nRay, Density_I, GradDensity_DI, DeltaS_I, RayFlag_I)
  //Get_Plasma_Density(
  //    int nRay, double Position_DI[3][nRay], double Density_I[nRay], double GradDensity_DI[3][nRay], 
  //              double DeltaSNew_I[nRay], short RayFlag_I[nRay]);
  //void plasma_density(
  //    int nRay, double Position_DI[3][nRay], double Density_I[nRay], double GradDensity_DI[3][nRay],
  //	          double DeltaS_I[nRay], short RayFlag_I[nRay]) {

  Get_Plasma_Density(nRay, Position_DI, Density_I, GradDensity_DI, DeltaSNew_I, RayFlag_I);

  if (deb) {
    printf("POINTER Density_I = %i\n", Density_I);
    printf("Density_I = "); print1d(nRay, Density_I); //printf("\n");
    printf("Density_I = "); print1d(nRay, Density_I);
    //printf("Density_I = %g %g %g %g\n", Density_I[0], Density_I[1], Density_I[2], Density_I[3]);
    //printf("GradDensity_DI = \n"); 
    //for (i = 0; i < nRay; i++) {
    //  printf("%g %g %g\n", GradDensity_DI[i], GradDensity_DI[nRay+i], GradDensity_DI[n2Ry+i]); 
    //}
  }

  //
  // Flag rays that got into the region with critical density
  //   where (Density_I .ge. DensityCr) RayFlag_I = .true.
  //
  for (iRay = 0; iRay < nRay; iRay++) {
    if (Density_I[iRay] >= DensityCr) RayFlag_I[iRay] = TRUE;   // TRUE indicates "bad ray"
  } // for (iRay=1; iRay < nRay; iRay++) 

  //if (deb) { 
  //printf("Dens2DensCr = %g, Density_I = ", Dens2DensCr); print1d(nRay, Density_I);
  //}


  for (iRay = 0; iRay < nRay; iRay++) {

    if (ExcludeRay_I[iRay]) continue;        // Do not process the rays that are done; save time

    HalfDeltaS = cHalf*DeltaS_I[iRay];
    //PositionHalfBack_D = Position_DI(:,iRay) - Slope_DI(:,iRay)*HalfDeltaS     // Original Position (at an integer point) 
    PositionHalfBack_D[0] = Position_DI[0][iRay] - Slope_DI[0][iRay]*HalfDeltaS; // Original Position (at an integer point)
    PositionHalfBack_D[1] = Position_DI[1][iRay] - Slope_DI[1][iRay]*HalfDeltaS; // Original Position (at an integer point)
    PositionHalfBack_D[2] = Position_DI[2][iRay] - Slope_DI[2][iRay]*HalfDeltaS; // Original Position (at an integer point)
    Dens2DensCr = Density_I[iRay]*DensityCrInv;

    DielPerm = cOne - Dens2DensCr;

    //GradDielPerm_D = -GradDensity_DI(:,iRay)*DensityCrInv
    GradDielPerm_D[0] = -GradDensity_DI[iRay]*DensityCrInv;
    GradDielPerm_D[1] = -GradDensity_DI[nRay+iRay]*DensityCrInv;
    GradDielPerm_D[2] = -GradDensity_DI[n2Ry+iRay]*DensityCrInv;
    //GradDielPerm_D[0] = -GradDensity_DI[0][iRay]*DensityCrInv;
    //GradDielPerm_D[1] = -GradDensity_DI[1][iRay]*DensityCrInv;
    //GradDielPerm_D[2] = -GradDensity_DI[2][iRay]*DensityCrInv;

    //GradDielPermSqr = sum(GradDielPerm_D**2);
    //GradDielPermSqr = pow(GradDielPerm_D[0],2) + pow(GradDielPerm_D[1],2) + pow(GradDielPerm_D[2],2);
    GradDielPermSqr = dot_product(GradDielPerm_D, GradDielPerm_D);

    if (deb) {
      printf("Slope_DI = %g %g %g\n", Slope_DI[0][iRay], Slope_DI[1][iRay], Slope_DI[1][iRay]);
      printf("GradDielPerm_D = "); print1d(3, GradDielPerm_D);
      printf("GradDielPermSqr = %g\n", GradDielPermSqr);
      double GradDielPermSqr1;
      GradDielPermSqr1 = pow(GradDielPerm_D[0],2) + pow(GradDielPerm_D[1],2) + pow(GradDielPerm_D[2],2);
      printf("GradDielPermSqr1 = %g\n", GradDielPermSqr1);
    }


    //GradEpsDotSlope = sum(GradDielPerm_D*Slope_DI(:,iRay))
    GradEpsDotSlope = GradDielPerm_D[0]*Slope_DI[0][iRay] + GradDielPerm_D[1]*Slope_DI[1][iRay]
      + GradDielPerm_D[2]*Slope_DI[2][iRay];

    DielPermHalfBack = DielPerm  - GradEpsDotSlope*HalfDeltaS;

    //Curv = (cHalf*HalfDeltaS/DielPermHalfBack)**2* &
    //       (GradDielPermSqr - GradEpsDotSlope**2)
    // (a x b)^2 = a^2*b^2 - (a dot b)^2
    // if |v|=1, (a x v)^2 = a^2 - (a dot v)^2
    Curv = pow(cHalf*HalfDeltaS/DielPermHalfBack,2) * (GradDielPermSqr - pow(GradEpsDotSlope,2));

    if (GentleRay_I[iRay]) {

      //
      // Check if the trajectory curvature is too sharp to meet the Tolerance
      // If so, reduce the DeltaS step  for the iRay-th ray and leave it until
      // the next call.
      //
      //

      if (Curv >=  ToleranceSqr) {
	DeltaS_I[iRay] = DeltaS_I[iRay]/(cTwo*sqrt(Curv/ToleranceSqr));
	//Position_DI(:,iRay) = PositionHalfBack_D
	Position_DI[0][iRay] = PositionHalfBack_D[0];
	Position_DI[1][iRay] = PositionHalfBack_D[1];
	Position_DI[2][iRay] = PositionHalfBack_D[2];
	if (deb) {
	  printf("Curv = %g, ToleranceSqr = %g\n", Curv, ToleranceSqr);
	}
	continue;    //-------------------------------------------------------------->>
      } // if (Curv >=  ToleranceSqr) 

      //
      // Test if some of the next points can get into the prohibited part of space with "negative"
      // dielectric permittivity
      //

      if (GradEpsDotSlope*HalfDeltaS <= -cThird*DielPermHalfBack) {
	//
	// Mark the ray as steep; memorize the distance to the critical surface;
	// reduce step
	//
	GentleRay_I[iRay] = FALSE;
	DistanceToCritSurf_I[iRay] = DielPermHalfBack/sqrt(GradDielPermSqr);
	DeltaS_I[iRay] =  cHalf*Tolerance*DistanceToCritSurf_I[iRay];
	//Position_DI(:,iRay) = PositionHalfBack_D
	Position_DI[0][iRay] = PositionHalfBack_D[0];
	Position_DI[1][iRay] = PositionHalfBack_D[1];
	Position_DI[2][iRay] = PositionHalfBack_D[2];
	continue;    //-------------------------------------------------------------->>
      } // if (GradEpsDotSlope*HalfDeltaS <= -cThird*DielPermHalfBack)

    } // if (GentleRay_I[iRay])

    //
    // Either switch to opposite branch of parabola
    // or make a Boris step
    //

    if ((GradEpsDotSlope*HalfDeltaS <= -cThird*DielPermHalfBack) || (DeltaS_I[iRay] < AbsoluteMinimumStep)) {

      // Switch to the opposite branch of parabolic trajectory
      //
      // When a next step can drive the ray into the area with the
      // plasma density greater then its critical value, then a special 
      // technique of "parabolic ray reflection" is employed. 
      // It can be shown that a ray trajectory in the medium with constant
      // gradient of dielectric permittivity is a parabola. If the step DeltaS_I
      // is small enough we can assume the grad \epsilon constant and hence
      // assume that the ray approaches the critical surface in a parabolic path.
      // We calculate the parameters of the parabola, namely:
      // StepX_D -- vector along the -grad \epsilon, with the length equal 
      //     to the distance from the PositionHalfBack_D to the parabola extremum;  
      // StepY_D -- vector perpendicular to StepX_D, lying inside of the
      //     parabola plane, pointing at the opposite parabola branch, and with
      //     the length equal the distance from PositionHalfBack_D to the 
      //     opposite branch of parabola.
      // The parabolic reflection is just replacement of the Position_DI with
      // the symmetric point at the opposite branch of parabola, and changing
      // the "incident" direction Slope_DI for the "departing" ray direction
      // according to Snell law.
      //

      if (deb) {
	printf("PARABOLA SWITCHING, iRay = %i, DeltaS = %g\n", iRay, DeltaS_I[iRay]);
	printf("PARABOLA SWITCHING, Slope_DI = %g %g %g\n", Slope_DI[0][iRay], Slope_DI[1][iRay], Slope_DI[2][iRay]);
      }

      LCosAl = -GradEpsDotSlope/GradDielPermSqr;

      if (deb) {printf("LCosAl = %g, GradEpsDotSlope = %g, GradDielPermSqr = %g;\n", LCosAl, GradEpsDotSlope, GradDielPermSqr); }

      ProjSlopeOnMinusGradEps_D[0] = -LCosAl*GradDielPerm_D[0];            // Here v_proj
      ProjSlopeOnMinusGradEps_D[1] = -LCosAl*GradDielPerm_D[1];            // Here v_proj
      ProjSlopeOnMinusGradEps_D[2] = -LCosAl*GradDielPerm_D[2];            // Here v_proj

      //StepY_D = Slope_DI(:,iRay) - ProjSlopeOnMinusGradEps_D        // Here c; |c| = sin \alpha
      StepY_D[0] = Slope_DI[0][iRay] - ProjSlopeOnMinusGradEps_D[0];        // Here c; |c| = sin \alpha
      StepY_D[1] = Slope_DI[1][iRay] - ProjSlopeOnMinusGradEps_D[1];        // Here c; |c| = sin \alpha
      StepY_D[2] = Slope_DI[2][iRay] - ProjSlopeOnMinusGradEps_D[2];        // Here c; |c| = sin \alpha

      //Slope_DI(:,iRay) = Slope_DI(:,iRay) - cTwo*ProjSlopeOnMinusGradEps_D       
      Slope_DI[0][iRay] = Slope_DI[0][iRay] - cTwo*ProjSlopeOnMinusGradEps_D[0];
      Slope_DI[1][iRay] = Slope_DI[1][iRay] - cTwo*ProjSlopeOnMinusGradEps_D[1];
      Slope_DI[2][iRay] = Slope_DI[2][iRay] - cTwo*ProjSlopeOnMinusGradEps_D[2];

      //
      //   We need to have |Step_Y| = 4 * sin(\alpha)*cos(\alpha)*DielPermHalfBack*L
      //   Now the direction of Step_Y is right, the length of it is equal to  sin(\alpha)
      //   Multiply it by L*Cos(\alpha)*DielPermHalfBack
      //

      StepY_D[0] = cFour*StepY_D[0]*LCosAl*DielPermHalfBack; 
      StepY_D[1] = cFour*StepY_D[1]*LCosAl*DielPermHalfBack; 
      StepY_D[2] = cFour*StepY_D[2]*LCosAl*DielPermHalfBack; 

      Position_DI[0][iRay] = PositionHalfBack_D[0] + StepY_D[0];
      Position_DI[1][iRay] = PositionHalfBack_D[1] + StepY_D[1];
      Position_DI[2][iRay] = PositionHalfBack_D[2] + StepY_D[2];
	
      //
      //   Step_X is in the direction of - grad \epsilon, whose vector is of the length of 1/L
      //   The length of Step_X is cos^2(\alpha)*L*DielPermHalfBack 
      //   Thus,
      //

      LCosAl2 = LCosAl*LCosAl;
      StepX_D[0] = (DielPermHalfBack*LCosAl2)*GradDielPerm_D[0];
      StepX_D[1] = (DielPermHalfBack*LCosAl2)*GradDielPerm_D[1];
      StepX_D[2] = (DielPermHalfBack*LCosAl2)*GradDielPerm_D[2];

      if (deb) {
	printf("Slope_DI = %g %g %g\n", Slope_DI[0][iRay], Slope_DI[1][iRay], Slope_DI[2][iRay]);
	printf("StepX_D = "); print1d(3, StepX_D);
	printf("StepY_D = "); print1d(3, StepY_D);
      }

      ////ParabLen = cTwo*sqrt(sum(StepX_D**2))
      //ParabLen = sqrt(sum((cTwo*StepX_D)**2) + sum(StepY_D**2))
      //StepXSqr = pow(cTwo*StepX_D[0],2) + pow(cTwo*StepX_D[1],2) + pow(cTwo*StepX_D[2],2);
      //StepYSqr = pow(StepY_D[0],2) + pow(StepY_D[1],2) + pow(StepY_D[2],2);
      StepXdbl_D[0] = cTwo*StepX_D[0];
      StepXdbl_D[1] = cTwo*StepX_D[1];
      StepXdbl_D[2] = cTwo*StepX_D[2];
      StepXSqr = dot_product(StepXdbl_D, StepXdbl_D); 
      StepYSqr = dot_product(StepY_D, StepY_D); 
      ParabLen = sqrt(StepXSqr + StepYSqr);
      if (deb) { printf("ParabLen = %g\n", ParabLen);}
      //Intensity_I(iRay) = Intensity_I(iRay) + ParabLen*(Dens2DensCr**2)*(cHalf - Dens2DensCr)**2
      Intensity_I[iRay] = Intensity_I[iRay] + ParabLen*pow(Dens2DensCr,2)*pow(cHalf - Dens2DensCr,2);
    } // if ((GradEpsDotSlope*HalfDeltaS <= -cThird*DielPermHalfBack) || (DeltaS_I[iRay] < AbsoluteMinimumStep))
    else {
      //
      // Make a step using Boris' algorithm
      //
	
      if (deb) {printf("// Make a step using Boris' algorithm...;\n");}

      PosBeforeBorisStep_D[0] = Position_DI[0][iRay];
      PosBeforeBorisStep_D[1] = Position_DI[1][iRay];
      PosBeforeBorisStep_D[2] = Position_DI[2][iRay];

      //Coef = cHalf*HalfDeltaS/(cOne - Dens2DensCr);
      Coef = cHalf*HalfDeltaS/(cOne - Dens2DensCr);
      RelGradRefrInx_D[0] = Coef*GradDielPerm_D[0];                       // grad(n) = grad(eps(i+1/2))/(2*eps(i+1/2))
      RelGradRefrInx_D[1] = Coef*GradDielPerm_D[1];                       // grad(n) = grad(eps(i+1/2))/(2*eps(i+1/2))
      RelGradRefrInx_D[2] = Coef*GradDielPerm_D[2];                       // grad(n) = grad(eps(i+1/2))/(2*eps(i+1/2))

      Slope_D[0] = Slope_DI[0][iRay];
      Slope_D[1] = Slope_DI[1][iRay];
      Slope_D[2] = Slope_DI[2][iRay];

      //Omega_D = cross_product(RelGradRefrInx_D, Slope_DI(:,iRay))
      cross_product(RelGradRefrInx_D, Slope_D, Omega_D);

      //Slope1_D = Slope_DI(:,iRay) + cross_product(Slope_DI(:,iRay), Omega_D)
      cross_product(Slope_D, Omega_D, Xprod_D);
      Slope1_D[0] = Slope_D[0] + Xprod_D[0];
      Slope1_D[1] = Slope_D[1] + Xprod_D[1];
      Slope1_D[2] = Slope_D[2] + Xprod_D[2];
	
     //Omega_D = cross_product(RelGradRefrInx_D, Slope1_D)
      cross_product(RelGradRefrInx_D, Slope1_D, Omega_D);

      //Slope1_D = Slope_DI(:,iRay) + cross_product(Slope_DI(:,iRay), Omega_D)
      cross_product(Slope_D, Omega_D, Xprod_D);
      Slope1_D[0] = Slope_D[0] + Xprod_D[0];
      Slope1_D[1] = Slope_D[1] + Xprod_D[1];
      Slope1_D[2] = Slope_D[2] + Xprod_D[2];

      //Curv1 = pow(Omega_D[0],2) + pow(Omega_D[1],2) + pow(Omega_D[2],2);
      Curv1 = dot_product(Omega_D, Omega_D);
      Coef = cTwo/(cOne + Curv1);
      //Slope_DI(:,iRay) = Slope_DI(:,iRay) + Coef*cross_product(Slope1_D, Omega_D)
      cross_product(Slope1_D, Omega_D, Xprod_D);
      Slope_DI[0][iRay] = Slope_DI[0][iRay] + Coef*Xprod_D[0];
      Slope_DI[1][iRay] = Slope_DI[1][iRay] + Coef*Xprod_D[1];
      Slope_DI[2][iRay] = Slope_DI[2][iRay] + Coef*Xprod_D[2];

      //Position_DI(:,iRay) = Position_DI(:,iRay) + Slope_DI(:,iRay)*HalfDeltaS
      Position_DI[0][iRay] = Position_DI[0][iRay] + Slope_DI[0][iRay]*HalfDeltaS;
      Position_DI[1][iRay] = Position_DI[1][iRay] + Slope_DI[1][iRay]*HalfDeltaS;
      Position_DI[2][iRay] = Position_DI[2][iRay] + Slope_DI[2][iRay]*HalfDeltaS;

      //Intensity_I(iRay) = Intensity_I(iRay) + DeltaS_I(iRay)*(Dens2DensCr**2)*(cHalf - Dens2DensCr)**2
      Intensity_I[iRay] = Intensity_I[iRay] + DeltaS_I[iRay]*pow(Dens2DensCr,2)*pow(cHalf - Dens2DensCr, 2);
      if (deb) { printf("DeltaS_I = %g, IntensIncr = %g\n", DeltaS_I[iRay], pow(Dens2DensCr,2)*pow(cHalf - Dens2DensCr, 2));}
    } // else
    //
    // Unsuccessful attempt to stop calculations for rays that reached the solar surface.
    // These can happen for frequencies above ~127 MHz, when the critical surface is 
    // below the solar surface, which is unphysical
    //
    //Pos_D[0] = Position_DI[0][iRay];
    //Pos_D[1] = Position_DI[1][iRay];
    //Pos_D[2] = Position_DI[2][iRay];
    //if (sqrt(dot_product(Pos_D, Pos_D)) <= SolarRadius) {
    //  ExcludeRay_I[iRay] = TRUE;
    //  Position_DI[0][iRay] = PosBeforeBorisStep_D[0];
    //  Position_DI[1][iRay] = PosBeforeBorisStep_D[1];
    //  Position_DI[2][iRay] = PosBeforeBorisStep_D[2];
    //  continue;
    //}

    //
    //   The code below makes gradual increases of the DeltaS up to the value
    // specified in DeltaSNew. The smooth step increase is required so as not to
    // get into the space behind the critical surface, stepping with DeltaS that
    // instantly changes from a very little size to the normal DeltaSNew length.
    // DeltaS is usually reduced in a close vicinity of the critical surface,
    // where the ray is travelling along a very sharp curve with high curvature.
    // For many rays it means fractioning of the DeltaS down several orders of 
    // magnitude, therefore the new step trial should start from a bigger step
    // of the same order of magnitude.
    //   This problem is solved using a non-linear difference equation:
    //           Y(i+1) = [2 - Y(i)/X(i)]*Y(i),
    // where X(i) is the desired final DeltaS length from DeltaSNew, and
    // Y(i) is the actual DeltaS length. A simple analysis of the equation shows
    // that, when Y is small compared to X, the next value of Y will be almost 
    // 2*X, so the DeltaS would grow in a geometrical progression. However, as
    // Y approaches X, its growth rate becomes slower. However, Y always reaches
    // X in several steps. One can check that for Y = X the next value of Y is 
    // always that of X.
    // 
    if (GentleRay_I[iRay]) {
      //
      // For shallow rays the DeltaS is increased unconditionally
      //
      DeltaS_I[iRay] = (cTwo - DeltaS_I[iRay]/DeltaSNew_I[iRay])*DeltaS_I[iRay];
    } // if (GentleRay_I[iRay])
    else {
      //
      // If the iRay-th ray is marked as steep (i.e. "not gentle" or "not shallow")
      // then the code below increases its DeltaS only if the current distance
      // to the critical surface, calculated as \epsilon / grad \epsilon, is greater 
      // than this distance value saved along with marking the ray as steep in
      // the DistanceToCritSurf_I. 
      //   This can happen in two cases: 
      // (1) either the ray was "parabola reflected" and
      // after several steps it went away from the surface by the distance where 
      // the parabola switching occurred; 
      // (2) or the ray is not steep any more because 
      // the current DeltaS is so small, that the next step does not penetrate the
      // critical surface.
      //   The ray is reverted to "gentle" or "shallow"
      //
      if (DielPermHalfBack > DistanceToCritSurf_I[iRay]*sqrt(GradDielPermSqr)) {
	GentleRay_I[iRay] = TRUE;
	DeltaS_I[iRay] = (cTwo - DeltaS_I[iRay]/DeltaSNew_I[iRay])*DeltaS_I[iRay];
      } // if (DielPermHalfBack > DistanceToCritSurf_I[iRay]*sqrt(GradDielPermSqr)) 
    } // else
  } // for (iRay = 1; iRay < nRay; iRay++) 
} // beam_path()
