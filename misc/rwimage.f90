subroutine beam_intensity(XyzEarth_D, RadioFrequency, ImageRange_I, rIntegration, nXPixel, nYPixel, Intensity_II)

  implicit none

  interface
     subroutine beam_path(Get_Plasma_Density, nRay, ExcludeRay_I, Position_DI, Slope_DI, DeltaS_I, &
          ToleranceInit, DensityCr, Intensity_I, RayFlag_I, NewEntry)
       external Get_Plasma_Density
       integer, intent(in) :: nRay                                  ! # of pixels in the raster
       real, intent(inout), dimension(3,nRay) :: Position_DI, Slope_DI
       real, intent(inout), dimension(nRay) :: Intensity_I, DeltaS_I
       real, intent(in) ::  ToleranceInit, DensityCr
       logical, intent(inout), dimension(nRay) :: RayFlag_I      ! .true. if a ray is OK, .false. otherwise
       logical, intent(inout) :: NewEntry                        ! Must be set to .true. before a call with new value of nRay
       logical, intent(inout), dimension(nRay) :: ExcludeRay_I   ! A ray is excluded from processing if it is .true.
     end subroutine beam_path

     function cross_product(a, b) result(c)
       real, dimension(3), intent(in) :: a, b
       real, dimension(3) :: c
     end function cross_product
  end interface

  external plasma_density
  real, intent(in) :: XyzEarth_D(3)          ! Earth position
  real, intent(in) :: RadioFrequency         ! 
  real, intent(in) :: ImageRange_I(4)        ! (x0, y0, x1, y1), i.e. (XLower, YLower, XUpper, YUpper)
  real, intent(in) :: rIntegration           ! Radius of "integration sphere"
  integer, intent(in) ::  nXPixel, nYPixel   ! Dimensions of the raster in pixels
  real, dimension(nYPixel,nXPixel), intent(out) :: Intensity_II          ! The result from emissivity integration
  logical, dimension(nYPixel*nXPixel) :: ExcludeRay_I   ! A ray is excluded from processing if it is .true.

  real, dimension(nYPixel,nXPixel) :: XPixel_II, YPixel_II   ! Pixel coordinates INSIDE of the image plane
  real :: Normal_D(3)                                          ! Unity vector normal to image plane
  real :: ZAxisOrt_D(3) = (/0, 0, 1/)
  real :: Tau_D(3), Xi_D(3)                                   ! Image plane inner Cartesian orts
  real :: XyzEarthLen
  real, dimension(3,nYPixel,nXPixel) :: XyzPixel_DII          ! SGI pixel coordinates
  real, dimension(3,nYPixel,nXPixel) :: Slope_DII  ! SGI unity slope vectors for all the line-of-sights pointing at the pixels 
  real, dimension(nYPixel,nXPixel) :: EarthToIntSphereDist_II  ! Distance from the radiotelescope to the integration sphere
  real, dimension(nYPixel,nXPixel) :: EarthToIntSphere2_II  ! Distance from the radiotelescope to the integration sphere
  real, dimension(3,nYPixel,nXPixel) :: Position_DII        ! 
  real, dimension(nYPixel,nXPixel) :: XPosition_II, YPosition_II, ZPosition_II, SolarDistSqr_II
  real, dimension(nYPixel*nXPixel) :: SolarDistSqr_I
  real :: XPixelDel, YPixelDel    
  real :: SlopeUnscaled_D(3)
  real :: XLower, XUpper, YLower, YUpper
  real, dimension(3,nXPixel*nYPixel) :: Position_DI, Slope_DI
  real, dimension(nXPixel*nYPixel) :: Intensity_I, DeltaS_I, RayPath_I
  logical, dimension(nXPixel*nYPixel) :: RayFlag_I         !  if a ray is OK, 0 otherwise
  logical :: NewEntry           ! Must be set to .true. before a call with new value of nRay
  real :: OneAU = 215.0, Tolerance = 0.01, DeltaS = 1.0
  real :: MaxRayPath = 60.
  real :: DensityCr 
  real :: XPixel, YPixel, SolarDistMin, MinRayPath, PercentRayLeft, rIntegrationSqr
  integer :: nRay, nIteration, i, j, iRay, nRayInsideIntSphere
  integer, dimension(nXPixel*nYPixel) :: RayInsideIntSphere_I
  logical :: deb = .false.
  !real, parameter :: ProtonChargeSGSe = 4.8e-10 !SGSe

  real, parameter :: cZero = 0.0, cOne = 1.0, cTwo = 2.0, cThree = 3.0, cFour = 4.0
  real, parameter :: cHalf = 0.5, cThird = 0.33333333333333333333
  real, parameter :: cPi = 3.1415926535897931;
  real, parameter :: ProtonChargeSGSe = 4.8e-10      ! StatCoulombs, SGSe
  real, parameter :: cProtonMass = 1.672621636E-24   !g
  real, parameter :: cElectronMass = 9.10938215E-28  !g
  integer ii, kk

  !
  ! Calculate the critical density from the frequency
  !
  write(*,*) 'cPi = ', cPi, ', cProtonMass = ',  cProtonMass, ', cElectronMass = ', cElectronMass, &
       ', RadioFrequency = ', RadioFrequency, ', ProtonChargeSGSe = ', ProtonChargeSGSe, &
       ', RadioFrequency/ProtonChargeSGSe =', RadioFrequency/ProtonChargeSGSe  

  DensityCr = cPi*cProtonMass*cElectronMass*(RadioFrequency/ProtonChargeSGSe)**2

   write(*,*) 'DensityCr = ', DensityCr
  !
  ! Determine the image plane inner coordinates of pixel centers
  !
  XLower = ImageRange_I(1)
  YLower = ImageRange_I(2)
  XUpper = ImageRange_I(3)
  YUpper = ImageRange_I(4)
  XPixelDel = (XUpper - XLower)/nXPixel
  YPixelDel = (YUpper - YLower)/nYPixel
  XPixel_II(1,:) = (/ (XLower + (real(j)-0.5)*XPixelDel, j = 1, nXPixel) /)
  do i = 2, nXPixel
     XPixel_II(i,:) = XPixel_II(1,:)
  end do
  YPixel_II(:,1) = (/ (YLower + (real(j)-0.5)*YPixelDel, j = 1, nYPixel) /)
  do i = 2, nYPixel
     YPixel_II(:,i) = YPixel_II(:,1)
  end do

  ! write(*,'(10f6.1)') (XPixel_II(i,:), i = 1, nYPixel)
  ! write(*,*)
  ! write(*,'(10f6.1)') (YPixel_II(i,:), i = 1, nYPixel)
  ! write(*,*)

  !
  ! Determune the orts, Tau and Xi, of the inner coordinate system of the image plane
  !
  XyzEarthLen = sqrt(sum(XyzEarth_D**2))
  Normal_D = XyzEarth_D/XyzEarthLen
  Tau_D = cross_product(ZAxisOrt_D, Normal_D)
  Xi_D = cross_product(Normal_D, Tau_D)

  write(*,*)
  write(*,'(a,3f8.5)') 'Normal_D = ', Normal_D
  write(*,'(a,3f8.5)') 'Tau_D = ', Tau_D
  write(*,'(a,3f8.5)') 'Xi_D = ', Xi_D
  write(*,*)

  !
  ! Calculate coordinates of all the pixels in the SGI
  !
  do i = 1, nYPixel
     do j = 1, nXPixel
        XyzPixel_DII(:,i,j) = XPixel_II(i,j)*Tau_D + YPixel_II(i,j)*Xi_D
        SlopeUnscaled_D = XyzPixel_DII(:,i,j) - XyzEarth_D
        Slope_DII(:,i,j) = SlopeUnscaled_D/sqrt(sum(SlopeUnscaled_D**2))             ! v
     end do
  end do

  !
  ! Find the points on the integration sphere where it intersects with the straight "rays" 
  !
  do i = 1, nYPixel
     do j = 1, nXPixel
        EarthToIntSphereDist_II(i,j) = -sum(XyzEarth_D*Slope_DII(:,i,j)) - &
             sqrt(sum((Slope_DII(:,i,j)*rIntegration)**2) - sum(cross_product(Slope_DII(:,i,j), XyzEarth_D)**2))
        EarthToIntSphere2_II(i,j) = -sum(XyzEarth_D*Slope_DII(:,i,j)) + &
             sqrt(sum((Slope_DII(:,i,j)*rIntegration)**2) - sum(cross_product(Slope_DII(:,i,j), XyzEarth_D)**2))
        Position_DII(:,i,j) = XyzEarth_D + Slope_DII(:,i,j)*EarthToIntSphereDist_II(i,j)
     end do
  end do


  !write(*,*) 
  !write(*,*) "Position_DII(X,:,:)"
  !write(*,'(10f7.2)') (Position_DII(1,:,j), j = 1, nXPixel)
  !write(*,*) 
  !write(*,*) "Position_DII(Y,:,:)"
  !write(*,'(10f7.2)') (Position_DII(2,:,j), j = 1, nXPixel)
  !write(*,*) 
  !write(*,*) "Position_DII(Z,:,:)"
  !write(*,'(10f7.2)') (Position_DII(3,:,j), j = 1, nXPixel)

  !if (deb) then
  !write(*,*) 'EarthToIntSphereDist_II = '
  !write(*,'(10f7.2)') (EarthToIntSphereDist_II(i,:), i = 1, nYPixel)
  !write(*,*)
  !write(*,*) 'EarthToIntSphere2 = '
  !write(*,'(10f7.2)') (EarthToIntSphere2_II(i,:), i = 1, nYPixel)
  !write(*,*)
  !================
  !write(*,*) 'Slope_X = ' 
  !write(*,'(10f10.5)') (Slope_DII(1,i,:), i = 1, nYPixel)
  !write(*,*)
  !write(*,*) 'Slope_Y = ' 
  !write(*,'(10f10.5)') (Slope_DII(2,i,:), i = 1, nYPixel)
  !write(*,*)
  !write(*,*) 'Slope_Z = ' 
  !write(*,'(10f10.5)') (Slope_DII(3,i,:), i = 1, nYPixel)
  !write(*,*)
  !end if

  !stop

  !
  ! Do emissivity integration inside of the integration sphere 
  !
  nRay = nXPixel*nYPixel

  Position_DI = reshape(Position_DII, (/3, nRay/))

  !write(*,*) "Position_DI:"
  !write(*,'(f7.2)') (Position_DI(:,j), j = 1, nRay)

  Slope_DI = reshape(Slope_DII, (/3, nRay/))
  Intensity_I = 0
  RayFlag_I = .false.
  NewEntry = .true.
  ExcludeRay_I = .false.
  RayInsideIntSphere_I = 1
  DeltaS_I = DeltaS
  RayPath_I = 0
  nIteration = 0
  rIntegrationSqr = rIntegration**2 + 0.01

  MinRayPath = 0.0
  nRayInsideIntSphere = nRay

  !write(*,*) 'nRayInsideIntSphere = ', nRayInsideIntSphere

  do while(nRayInsideIntSphere .gt. 0)
  !do kk = 1, 2
     nIteration = nIteration + 1

     SolarDistSqr_I = sum(Position_DI**2,1)
     where(SolarDistSqr_I .gt. rIntegrationSqr) 
        RayInsideIntSphere_I = 0
        ExcludeRay_I = .true.
     end where

     nRayInsideIntSphere = sum(RayInsideIntSphere_I)

     if (mod(nIteration,int(10.0/DeltaS)) .eq. 0) then

     PercentRayLeft = (real(nRayInsideIntSphere))/real(nRay)*100.0

     write(*,*)
     write(*,'(i5,a,f8.4,a,f21.6)') nIteration,': ======= ', PercentRayLeft, &
          ' % ====== Minimum Distance from Sun: ', minval(sqrt(SolarDistSqr_I)) 
     write(*,*)
     end if

     call beam_path(Plasma_Density, nRay, ExcludeRay_I, Position_DI, Slope_DI, DeltaS_I, &
          Tolerance, DensityCr, Intensity_I, RayFlag_I, NewEntry)     

     RayPath_I = RayPath_I + DeltaS_I
     !write(*,*) RayPath_I
     !write(*,*) "Position_DI:"
     !do ii = 1, nRay
     !   write(*,*) Position_DI(:,ii)
     !end do
     MinRayPath = minval(RayPath_I)
  end do

  Intensity_II = reshape(Intensity_I, (/nYPixel,nXPixel/))

end subroutine beam_intensity



!========================================================================

subroutine plasma_density(Position_DI, nRay, Density_I, GradDensity_DI, DeltaS_I, RayFlag_I)
  !
  ! at a specified location, Position_D(3)
  !

  implicit none

  integer, intent(in) :: nRay
  real, intent(in) :: Position_DI(3,nRay)
  real, intent(out) :: Density_I(nRay), GradDensity_DI(3,nRay),DeltaS_I(nRay)
  logical, intent(inout) :: RayFlag_I(nRay)
  real :: OneAU = 215.0, SolarDistSqr_I(nRay), SolarDistQuad_I(nRay)
  !real, parameter :: DensityCr = 1./9.
  ! DensityAtSolarSurface = Ne(@SolarSurface)*ProtonMass = 2x10^8(cm^-3)*1.6726x10^-24(g) = 3.3452e-16(g/cm^3)
  real, parameter :: DensityAtSolarSurface = 3.3452e-16 ! g/cm^3
  integer :: i

  SolarDistSqr_I = sum(Position_DI**2,1)
  Density_I = DensityAtSolarSurface/SolarDistSqr_I
  SolarDistQuad_I = SolarDistSqr_I**2
  do i = 1, 3
     GradDensity_DI(i,:) = -2.*DensityAtSolarSurface*Position_DI(i,:)/SolarDistQuad_I
  end do

  DeltaS_I = 1.0       ! Revert to original
end subroutine Plasma_Density

