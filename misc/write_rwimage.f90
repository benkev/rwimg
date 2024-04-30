program write_rwimage

  !real REarth_D(3) = (/150.0, 145.7, 50.0/)                ! Earth position in the SGI, in solar radii
  real :: REarth_D(3) = (/150.0, 0.0, 0.0/)                    ! Earth position in the SGI, in solar radii
  real :: RadioFrequency =  42.3e6;                            ! Hz
  real :: ImageRange_I(4) = (/-10.0, -10.0, 10.0, 10.0/)    ! (x0, y0, x1, y1)), i.e. XLower, YLower, XUpper, YUpper
  real :: HalfImageRangeX, HalfImageRangeY
  real :: rIntegration = 25.0                               ! AU, Radius of "integration sphere"
  integer, parameter ::  nXPixel = 1000, nYPixel = 1000         ! Dimensions of the raster in pixels
  real :: Intensity_II(nYPixel, nXPixel) 
  integer fh, i, j;


  call beam_intensity(REarth_D, RadioFrequency, ImageRange_I, rIntegration, nXPixel, nYPixel, Intensity_II);

  fh = 1
  open(fh, file="image.txt", status="replace")

  !write(fh,*) ImageRange_II

  do i = 1, nYPixel
      write(fh,'(1000f10.5)') Intensity_II(i,1:nXPixel)
  end do

  close(fh)

end program write_rwimage
