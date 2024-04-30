subroutine beam_path(Get_Plasma_Density, nRay, ExcludeRay_I, Position_DI, Slope_DI, DeltaS_I, &
     ToleranceInit, DensityCr, Intensity_I, RayFlag_I, NewEntry)

  !
  !   The subroutine beam_path() makes raytracing and emissivity integration along ray paths.
  ! It works on a group of rays (a beam). Each ray is specified by its Cartesian Position_DI
  ! and its direction cosines Slope_DI. The subroutine calculates new Position_DI, which is 
  ! DeltaS_I away from the initial position. It calculates the intensity along the step
  ! as emissivity by DeltaS and adds this to the Intensity_I. Thus, after series of calls
  ! to beam_path(), the Intensity_I contains the result of integration along the paths of 
  ! every ray. 
  !   The electromagnetic rays are refracted in plasma due to its non-zero refractive index,
  ! which is the square root of the diellectric permittivity \epsilon. The \epsilon, in turn,
  ! is a function of the plasma density. The external subprogram, Get_Plasma_Density(),
  ! provides the plasma Density_I along with its gradient, GradDensity_DI at the Position_DI.
  ! The \epsilon can be calculated as 1 - Density_I/DensityCr, where DensityCr is the "critical"
  ! plasma density at which the dielectric permittivity \epsilon falls down to zero. 
  ! The value of DensityCr is proportional to the square of the wave frequency. For example,
  ! for the 42.3 MHz radiowaves the critical density is ~3.71x10^(-17) g/cm^3. 
  ! A surface where the plasma density achieves the critical value acts lika a mirror. 
  ! No radiation penetrates the critical surface. The radiowaves can only travel in the regions
  ! with lower density.
  !   The emissivity w of plasma at a selected plasma frequency is calculated as a polynomial
  ! function w = (Dens2DensCr)^2*(0.5 - Dens2DensCr)^2, where Dens2DensCr is the quotient
  ! Density_I/DensityCr. So, here the plasma density is used for calculation of both emissivity
  ! and dielectric permittyvity.
  !   The parameters of the beam_path() are briefly described below.
  !
  ! Get_Plasma_Density():  external subroutine that returns the plasma Density_I and its gradient
  !     GradDensity_DI. It also provides the recommended step size DeltaSNew_I and asserts
  !     the RayFlag_I (assigns .true.) for a ray should it (illegally) penetrate into the region
  !     with "negative" dielectric permittivity.
  ! nRay:           number of rays being processed.
  ! ExcludeRay_I:  the caller program can use this logical array to stop processing of any individual 
  !     ray when it already finished travelling inside of a specified part of space. Set the 
  !     corresponding element of ExcludeRay_I to .true. to leave it unprocessed during the subsequent
  !     calls to beam_path(). Before the first call to beam_path() all the elements of ExcludeRay_I
  !     should be set to .false.; then all the rays will be processed.
  ! Position_DI:    Cartesian position vectors of the rays.
  ! Slope_DI:       Direction cosines of the rays.
  ! DeltaS_I:       Current step values for the rays. Set the elements of Slope_DI to some reasonable 
  !     value, say, 1.0. The DeltaS_I elements are individually modified by 
  !     beam_path() to satisfy the precision requirements set by ToleranceInit.
  ! ToleranceInit:   determines the precision of ray paths calculation. ToleranceInit is the inverse
  !     of the minimum number of ray trajectory points per one radian of the ray curvature. If this
  !     requirement is not satisfied, the corresponding element of DeltaS_I is decreased. Do not set
  !     the ToleranceInit to any value greater than 0.1 (which means 10 points per curvature radian): 
  !     it will be internally set to 0.1 anyway.   
  ! DensityCr: the plasma density at which its dielectric permittivity becomes zero for chosen wave
  !     frequency.
  ! Intensity_I:    the intensities of each ray calculated as the integral of emissivity along the ray 
  !     path. During each call to ray_path(), each element of the Intensity_I is incremented by the
  !     integral along the step DeltaS_I. Set all the elements of Intensity_I to 0.0 before the 
  !     first call to the beam_path().
  ! RayFlag_I:      the .true. elements of this logical array indicate that the corresponding rays
  !     penetrated into the "prohibited" region of space with the plasma density above its critical 
  !     value. Normally, it should never happen. However, in case the algorithm made such an error,
  !     the flagged rays should be considered as "bad" and thrown away from the resultant Intensity_I.
  !     Set all the elements of RayFlag_I to .false. before calling ray_path() for the first time.
  ! NewEntry:       Set this logical variable to .true. before the first call to beam_path(). This 
  !     value forces the beam_path() to allocate internal dynamic arrays and take several initial
  !     actions. During subsequent calls to beam_path() the NewEntry will keep the value .false..
  !     Setting the NewEntry to .true. during the same run will cause another allocation of the internal
  !     allocatables, with possibly different dimensions. This variable was introduced out of the 
  !     speed considerations: the automatic arrays work slower than the allocatable ones, because
  !     the allocation occurs only once at the first call.
  ! 
  !

  implicit none

  interface
     subroutine Get_Plasma_Density(Position_DI, nRay, Density_I, GradDensity_DI, DeltaSNew_I, RayFlag_I)
       !
       ! Calculates plasma density, Density_I, and its gradient, 
       ! GradDensity_DI(3,nRay), at specified locations Position_DI(3,nRay)
       ! Also, it provides appropriate step, DeltaSNew_I, conceivably dependent
       ! on the numeric grid size
       !
       implicit none
       integer, intent(in) :: nRay
       real, intent(in) :: Position_DI(3,nRay)
       logical, intent(inout) :: RayFlag_I(nRay)
       real, intent(out) :: Density_I(nRay), GradDensity_DI(3,nRay), DeltaSNew_I(nRay)
     end subroutine Get_Plasma_Density

     function cross_product(a, b) result(c)
       implicit none
       real, dimension(3), intent(in) :: a, b
       real, dimension(3) :: c
     end function cross_product
  end interface


  !    external Get_Plasma_Density
  integer, intent(in) :: nRay                                  ! # of pixels in the raster
  real, intent(inout), dimension(3,nRay) :: Position_DI, Slope_DI
  real, intent(inout), dimension(nRay) :: Intensity_I, DeltaS_I
  real, intent(in) ::  ToleranceInit, DensityCr
  logical, intent(inout), dimension(nRay) :: RayFlag_I      ! .true. if a ray is OK, .false. otherwise
  logical, intent(inout) :: NewEntry                        ! Must be set to .true. before a call with new value of nRay
  logical, intent(inout), dimension(nRay) :: ExcludeRay_I   ! A ray is excluded from processing if it is .true.
  integer, parameter :: nSplitDeltaS = 2

  real, parameter :: cZero = 0.0, cOne = 1.0, cTwo = 2.0, cThree = 3.0, cFour = 4.0
  real, parameter :: cHalf = 0.5, cThird = 0.33333333333333333333

  real, save, dimension(3)       :: Slope1_D, Omega_D
  real, save, dimension(3)       :: ProjSlopeOnMinusGradEps_D
  real, save, dimension(3)       :: StepX_D, StepY_D, RelGradRefrInx_D 
  real, save, dimension(3)       :: GradDielPerm_D, PositionHalfBack_D

  real, save, pointer, dimension(:,:) :: GradDensity_DI
  real, save, pointer, dimension(:)   :: Density_I
  real, save, pointer, dimension(:)   :: DeltaSNew_I
  real, save, pointer, dimension(:)   :: DistanceToCritSurf_I
  logical, save, pointer, dimension(:)   :: GentleRay_I       ! .true. for shallow rays; is set to .false. for steep rays.

  real :: HalfDeltaS                        ! DeltaS halved
  real :: DielPerm, DielPermHalfBack, Dens2DensCr, Dens2DensCr1, Dens2DensCr2
  real :: Coef, Curv, Curv1
  real :: LCosAl    !Where L is inverse grad of \epsilon, Alpha is the incidence angle
  real :: GradDielPermSqr, GradEpsDotSlope
  real :: ParabLen, GradDielPerm
  real, save :: Tolerance, ToleranceSqr, DensityCrInv, AbsoluteMinimumStep
  integer :: i, j, iRay 

  if (NewEntry) then
     NewEntry = .false.
     allocate(GradDensity_DI(3,nRay))
     allocate(Density_I(nRay))
     allocate(DeltaSNew_I(nRay))
     allocate(DistanceToCritSurf_I(nRay))
     DistanceToCritSurf_I = cZero
     allocate(GentleRay_I(nRay))
     GentleRay_I = .true.
     DensityCrInv = cOne/DensityCr
     Tolerance = min(ToleranceInit,0.1)  ! i.e. minimum ten points between a vacuum and a critical surface and
     ToleranceSqr = Tolerance**2         !  minimum 10 points over 1 rad of the curvature
     AbsoluteMinimumStep = 1e-4*sum(DeltaS_I)/nRay ! One ten-thousandth of average step
     !write(*,*) 'AbsoluteMinimumStep = ', AbsoluteMinimumStep
  end if


  do iRay=1,nRay
     if (ExcludeRay_I(iRay)) CYCLE        ! Do not process the rays that are done; save time
     HalfDeltaS = cHalf*DeltaS_I(iRay)
     Position_DI(:,iRay) = Position_DI(:,iRay) + Slope_DI(:,iRay)*HalfDeltaS ! Now Position_DI moved by 1/2 DeltaS !!!
  end do

  call Get_Plasma_Density(Position_DI, nRay, Density_I, GradDensity_DI, DeltaSNew_I, RayFlag_I)

  where (Density_I .ge. DensityCr) RayFlag_I = .true.   ! .true. indicates "bad ray"

  do iRay = 1, nRay

     if (ExcludeRay_I(iRay)) CYCLE        ! Do not process the rays that are done; save time

     !if (iRay .eq. 4950) then
     !   write(*,'(a,i5,a,f9.6)') 'iRay = ', iRay, ', DeltaS = ', DeltaS_I(iRay)
     !end if

     !if (Density_I(iRay) .le. cZero)  write(*,*) '++++++++++++++++ iRay = ', iRay, ', Density_I = ', Density_I(iRay), ' < 0 +++'

     HalfDeltaS = cHalf*DeltaS_I(iRay)
     PositionHalfBack_D = Position_DI(:,iRay) - Slope_DI(:,iRay)*HalfDeltaS        ! Original Position (at an integer point) 
     Dens2DensCr = Density_I(iRay)*DensityCrInv

     DielPerm = cOne - Dens2DensCr
     GradDielPerm_D = -GradDensity_DI(:,iRay)*DensityCrInv
     GradDielPermSqr = sum(GradDielPerm_D**2)

     GradEpsDotSlope = sum(GradDielPerm_D*Slope_DI(:,iRay))

     DielPermHalfBack = DielPerm  - GradEpsDotSlope*HalfDeltaS

     Curv = (cHalf*HalfDeltaS/DielPermHalfBack)**2* &  ! (a x b)^2 = a^2*b^2 - (a dot b)^2
          (GradDielPermSqr - GradEpsDotSlope**2)       ! if |v|=1, (a x v)^2 = a^2 - (a dot v)^2

     if (GentleRay_I(iRay)) then

        !
        ! Check if the trajectory curvature is too sharp to meet the Tolerance
        ! If so, reduce the DeltaS step  for the iRay-th ray and leave it until
        ! the next call.
        !
        !

        if (Curv .ge. ToleranceSqr) then
           DeltaS_I(iRay) = DeltaS_I(iRay)/(cTwo*sqrt(Curv/ToleranceSqr))
           Position_DI(:,iRay) = PositionHalfBack_D
           !if (iRay .eq. 4950) then
           !   write(*,*) '<= CYCLED due to Curv test'
           !end if
           CYCLE
        end if

        !
        ! Test if some of the next points can get into the prohibited part of space with "negative"
        ! dielectric permittivity
        !

        if (GradEpsDotSlope*HalfDeltaS .le. -cThird*DielPermHalfBack) then  ! Too close to critical surface?
           !
           ! Mark the ray as steep; memorize the distance to the critical surface;
           ! reduce step
           !
           GentleRay_I(iRay) = .false.
           DistanceToCritSurf_I(iRay) = DielPermHalfBack/sqrt(GradDielPermSqr)
           DeltaS_I(iRay) =  cHalf*Tolerance*DistanceToCritSurf_I(iRay)
           Position_DI(:,iRay) = PositionHalfBack_D
           !if (iRay .eq. 4950) then
           !   write(*,*) '<= CYCLED due to Steep test'
           !end if
           CYCLE
        end if

     end if ! GentleRay_I

     !
     ! Either switch to opposite branch of parabola
     ! or make a Boris step
     !

     if ((GradEpsDotSlope*HalfDeltaS .le. -cThird*DielPermHalfBack) .or. (DeltaS_I(iRay) .lt. AbsoluteMinimumStep)) then

        ! Switch to the opposite branch of parabolic trajectory
        !
        ! When a next step can drive the ray into the area with the
        ! plasma density greater than its critical value, then a special 
        ! technique of "parabolic ray reflection" is employed. 
        ! It can be shown that a ray trajectory in the medium with constant
        ! gradient of dielectric permittivity is a parabola. If the step DeltaS_I
        ! is small enough we can assume the grad \epsilon constant and hence
        ! assume that the ray approaches the critical surface in a parabolic path.
        ! We calculate the parameters of the parabola, namely:
        ! StepX_D -- vector along the -grad \epsilon, with the length equal 
        !     to the distance from the PositionHalfBack_D to the parabola extremum;  
        ! StepY_D -- vector perpendicular to StepX_D, lying inside of the
        !     parabola plane, pointing at the opposite parabola branch, and with
        !     the length equal the distance from PositionHalfBack_D to the 
        !     opposite branch of parabola.
        ! The parabolic reflection is just replacement of the Position_DI with
        ! the symmetric point at the opposite branch of parabola, and changing
        ! the "incident" direction Slope_DI for the "departing" ray direction
        ! according to Snell law.
        !


        !write(*,*) 'PARABOLA SWITCHING, iRay = ', iRay, ', DeltaS = ', DeltaS_I(iRay)

        LCosAl = -GradEpsDotSlope/GradDielPermSqr
        ProjSlopeOnMinusGradEps_D = -LCosAl*GradDielPerm_D            ! Here v_proj

        StepY_D = Slope_DI(:,iRay) - ProjSlopeOnMinusGradEps_D        ! Here c; |c| = sin \alpha

        Slope_DI(:,iRay) = Slope_DI(:,iRay) - cTwo*ProjSlopeOnMinusGradEps_D

        !
        !   We need to have |Step_Y| = 4 * sin(\alpha)*cos(\alpha)*DielPermHalfBack*L
        !   Now the direction of Step_Y is right, the length of it is equal to  sin(\alpha)
        !   Multiply it by L*Cos(\alpha)*DielPermHalfBack
        !

        StepY_D = cFour*StepY_D*LCosAl*DielPermHalfBack 

        Position_DI(:,iRay) = PositionHalfBack_D + StepY_D

        !
        !   Step_X is in the direction of - grad \epsilon, whose vector is of the length of 1/L
        !   The length of Step_X is cos^2(\alpha)*L*DielPermHalfBack 
        !   Thus,
        !

        StepX_D = (DielPermHalfBack*LCosAl**2)*GradDielPerm_D

        !ParabLen = cTwo*sqrt(sum(StepX_D**2))
        ParabLen = sqrt(sum((cTwo*StepX_D)**2) + sum(StepY_D**2))

        Intensity_I(iRay) = Intensity_I(iRay) + ParabLen*(Dens2DensCr**2)*(cHalf - Dens2DensCr)**2

     else 

        ! Make a step using Boris' algorithm

        Coef=cHalf*HalfDeltaS/(cOne - Dens2DensCr)
        RelGradRefrInx_D = Coef*GradDielPerm_D                       ! grad(n) = grad(eps(i+1/2))/(2*eps(i+1/2))

        Omega_D = cross_product(RelGradRefrInx_D, Slope_DI(:,iRay))
        Slope1_D = Slope_DI(:,iRay) + cross_product(Slope_DI(:,iRay), Omega_D)

        Omega_D = cross_product(RelGradRefrInx_D, Slope1_D)
        Slope1_D = Slope_DI(:,iRay) + cross_product(Slope_DI(:,iRay), Omega_D)

        Curv1 = sum(Omega_D**2)
        Coef = cTwo/(cOne + Curv1)
        Slope_DI(:,iRay) = Slope_DI(:,iRay) + Coef*cross_product(Slope1_D, Omega_D)

        Position_DI(:,iRay) = Position_DI(:,iRay) + Slope_DI(:,iRay)*HalfDeltaS

        Intensity_I(iRay) = Intensity_I(iRay) + DeltaS_I(iRay)*(Dens2DensCr**2)*(cHalf - Dens2DensCr)**2

     end if

     !
     !   The code below makes gradual increases of the DeltaS up to the value
     ! specified in DeltaSNew. The smooth step increase is required so as not to
     ! get into the space behind the critical surface, stepping with DeltaS that
     ! instantly changes from a very little size to the normal DeltaSNew length.
     ! DeltaS is usually reduced in a close vicinity of the critical surface,
     ! where the ray is travelling along a very sharp curve with high curvature.
     ! For many rays it means fractioning of the DeltaS down several orders of 
     ! magnitude, therefore the new step trial should start from a bigger step
     ! of the same order of magnitude.
     !   This problem is solved using a non-linear difference equation:
     !           Y(i+1) = [2 - Y(i)/X(i)]*Y(i),
     ! where X(i) is the desired final DeltaS length from DeltaSNew, and
     ! Y(i) is the actual DeltaS length. A simple analysis of the equation shows
     ! that, when Y is small compared to X, the next value of Y will be almost 
     ! 2*X, so the DeltaS would grow in a geometrical progression. However, as
     ! Y approaches X, its growth rate becomes slower. However, Y always reaches
     ! X in several steps. One can check that for Y = X the nexy value of Y is 
     ! always that of X.
     ! 
     if (GentleRay_I(iRay)) then
        !
        ! For shallow rays the DeltaS is increased unconditionally
        !
        DeltaS_I(iRay) = (cTwo - DeltaS_I(iRay)/DeltaSNew_I(iRay))*DeltaS_I(iRay)
     else 
        !
        ! If the iRay-th ray is marked as steep (i.e. "not gentle" or "not shallow")
        ! then the code below increases its DeltaS only if the current distance
        ! to the critical surface, calculated as \epsilon / grad \epsilon, is greater 
        ! than this distance value saved along with marking the ray as steep in
        ! the DistanceToCritSurf_I. 
        !   This can happen in two cases: 
        ! (1) either the ray was "parabola reflected" and
        ! after several steps it went away from the surface by the distance where 
        ! the parabola switching occurred; 
        ! (2)or the ray is not steep any more because 
        ! the current DeltaS is so small, that the next step does not penetrate the
        ! critical surface.
        !   The ray is reverted to "gentle" or "shallow"
        !
        if (DielPermHalfBack .gt. DistanceToCritSurf_I(iRay)*sqrt(GradDielPermSqr)) then
           GentleRay_I(iRay) = .true.
           DeltaS_I(iRay) = (cTwo - DeltaS_I(iRay)/DeltaSNew_I(iRay))*DeltaS_I(iRay)
           !if (iRay .eq. 4950) then
           !   write(*,*) '=> REVERTED to GENTLE: iRay = ', iRay
           !end if
        end if
     end if
  end do

end subroutine beam_path


!
! Vector cross-product: c = a x b
!
function cross_product(a, b) result(c)
  implicit none
  real, dimension(3), intent(in) :: a, b
  real, dimension(3) :: c
  c(1) = a(2)*b(3) - a(3)*b(2)
  c(2) = -a(1)*b(3) + a(3)*b(1)
  c(3) = a(1)*b(2) - a(2)*b(1)
end function cross_product

