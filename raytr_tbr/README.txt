This directory contains Fortran and C programs for raytracing
by Benkevitch and Sokolov.

This branch (raytr_tbr) is created in order to try to correct 
the raytracing algorithm so that it could work for exponential 
electron number density distribution in the chromosphere. 
Probably, this will require combining the "parabola-switching"
with the "Snell-law-switching" in the regions with exponential
Ne distribution (so that a ray as if is reflected from a solid
surface). This method can only be used if the "parabola" is 
narrow enough and the incidence point is actually the same as 
the reflection point. 
It is supposed that this correction will allow raytracing for 
the sun images in the radio bands of gigaHertz and tens or 
hundreds gigaHertz. For instance, the "limb brightening" effect
is expected to be reproduced.



Current errors:

(1) for radio frequencies above ~127 MHz the critical surface radius 
is less than that of the sun, so the rays penetrate the solar surface.
It is wrong.

(2) The  intensity is calculated as if there is NO absorption.

(3) The raytracing near steeper rays makes "kinks" in intensity.
 This may be due to wrong conditions for "parabola switching" or so.

2009-Jan-15: The archiving of programs before changing the rwimage.c
	and write_rwimage.c so that the trajrctories of several selected 
	rays be stored to be plotted for visual control and illustration.

2009-Jan-20: Added trajectory saving feature. The beam_intensity()
is replaced by beam_intens_traj() in rwimage.c. Also included
Python programs for image for presentation  creation 

2009-Jul-14: The rays marked as bad (RayFlag_I[] = 1) do not stop
processing, which causes infinite loop in case of the exponential
law for electron number density in the chromosphere.