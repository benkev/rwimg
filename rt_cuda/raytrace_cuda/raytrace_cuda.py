###########################################################
#                                                         #
# raytrace_cuda.py                                        #
#                                                         #
# Module for making images  of space objects              #
# through ray tracing                                     #
# Created: 14 March 2011 by Leonid Bevkevitch             #
# Changed name: 23 Feb 2012 from raytrace.py to           #
#          raytrace_cuda.py                               #
###########################################################
#
import pylab as pl
import numpy as np
import os
import glob
from distutils.sysconfig import get_python_lib
import rtcore
import pyfits
from datetime import date

class Fatal_Error(): pass  # Exception for all fatal errors

class rtconstants():
    """
    Physical constants.
    This class is encapsulates many parameters so as not
    to contaminate the namespaces with many names.
    After instantiation, other parameters can be added.
    """
    def __init__(self):
        self.ProtonChargeCGSe = 4.80320427e-10   # StatCoulombs, CGSe
        self.ProtonMass_g = 1.672621638e-24   # g
        self.ProtonMassInv = 1./self.ProtonMass_g   # g^(-1)
        self.ElectronMass_g = 9.10938215e-28  # g
        self.Boltzmannk_ergK = 1.380650424e-16;  # erg/K, Boltzmann const CGS
        self.Boltzmannk_SI = 1.380650424e-23
        self.Rsun_cm = 6.955e10               # cm, Solar radius
        self.Rsun_km = 6.955e5                # km, Solar radius
        self.c_light_cms = 2.9979245800e10    # cm/s, speed of light in vacuum
        self.c_light_SI = 2.9979245800e8      # m/s, speed of light in vacuum
        self.c_light_Rsun = self.c_light_cms/self.Rsun_cm # in Rsun
        self.h_chromo_km = 10000.     # km, Chromosphere thickness
        self.r_chromo = (self.Rsun_km + self.h_chromo_km - 1000.)/self.Rsun_km
        self.r_chro_cor = (self.Rsun_km + self.h_chromo_km)/self.Rsun_km
        self.r_corona = (self.Rsun_km + self.h_chromo_km + 1000.)/self.Rsun_km
        self.Te_corona_K = 1.0e6      # K
        self.Te_chromo_K = 3.0e4      # K
        self.AU_m = 149597870700.     # Sun-earth distance in meters
        self.AU_Rsun = 215.097        # Sun-earth distance in solar radii
        self.Msun_G = 1.0             # Solar dipole, G/Rsun^3
        # Solar dipole direction vector (unity length): 
        self.Mdir = np.array((0., 0., 1.), dtype=np.double)
        self.Coulomb_log = 20.0       # Coulomb logarithm

class implane():
    """
    Image plane representation.
    """
    def __init__(self, grid=(10,10), rect=(-5., -5, 5., 5.),
                 obs=(215., 0., 0.), rsph = 25., freq=100e6,
                 mode='Tbr', cnu=3.0, msun=1.0, units='K',
                 scattering=False):
        """
        Create an image plane as a two-dimensional array of doubles.
        The pixel sizes are RA x Dec = grid[1] x grid[0] = nx * ny.
        The image dimensions are Dec x RA = dim[0] x dim[1].
        """
        rect = np.double(rect)
        obs = np.double(obs)

        #self.traj = bool(traj);  # True: store trajectories, False: no
        
        
        if mode == 'basic' or mode == 1:
            self.mode = np.int32(1)     # No Tbr or TbrIQUV calculation;
        
        elif mode == 'Tbr' or mode == 2:
            self.mode = np.int32(2)     # Tbr calculation;
        
        elif mode == 'TbrIQUV' or mode == 3:
            self.mode = np.int32(3)     # No Tbr or TbrIQUV calculation;
        
        else:
            print 'Error: mode parameter can be only 1, 2, or 3'
            return
        

        # Exceptions
        self.XCP = Fatal_Error()

        # Determine the location of the raytrace package in the site-packages
        # (/usr/lib/python2.6/site-packages/raytrace or
        # /usr/lib64/python2.6/site-packages/raytrace)
        # NOTE: get_python_lib() returns the non arch-specific directory
        # which is under /usr/lib. One needs to use get_python_lib(1) to
        # get the platform specific directory if that's required.
        
        #self.package = get_python_lib(1) + '/raytrace'
        self.package = os.path.dirname(rtcore.__file__)
        

        # Bit flag positions for self.flags
        self.INACTIVE =  0x0001   # The ray is not being processed  
        self.SHARP =     0x0002   # The ray is sharp (close to normal) 
        self.PENETR =    0x0004   # Eps < 0; Critical surface search
        self.WASREFL =   0x0008   # First step after reflection
        self.BADRAY =    0x0010   # The ray is bad
        self.TRACED =    0x0020   # The ray is traced (trajectory stored)

        # Create structure prm containing physical constants
        self.const = rtconstants()

        # Change some constants
        self.const.Msun_G = msun  # Solar dipole strength in Gauss at equator

        # Counter for calls to the algorithm; roughly number of steps
        self.callcount = 1.0 
        
        # Constants
        self.DeltaS = 0.1   # Initial ray path increment in solar radii

        ## # Stokes components
        ## self.I = None
        ## self.Q = None
        ## self.U = None
        ## self.V = None

        # Grid sizes along Dec and RA: grid=[nDec,nRA]
        self.grid = np.array(grid, dtype=np.int32)

        # Grid size along Dec
        self.ny = ny = grid[0]
        
        # Grid size along RA
        self.nx = nx = grid[1]
        
        # Positions of the image plane left bottom and right top
        # corners in SGI system
        # rect = [xleft, ybottom, xright, ytop]
        #      = [RA0, Dec0, RA1, Dec1]
        #      = [x0,  y0,   x1,  y1]
        self.rect = np.empty(4, dtype=np.double)
        self.rect[:] = rect

        # Image plane size along RA
        self.sizera = rect[2] - rect[0]
        
        # Image plane size along Dec
        self.sizedec = rect[3] - rect[1]
        
        # Frequency of the radio waves
        self.freq = np.double(freq)
        self.freq2 = freq**2
        self.omega = 2.*np.pi*freq
        self.omega2 = self.omega**2

        # Wavenumber in a vacuum
        self.k0 = self.omega/self.const.c_light_cms          # rad/cm
        self.k0_Rsun = self.omega/self.const.c_light_Rsun    # rad/Rsun

        # Coefficient at the Ginzburg's formula for effective collision
        # frequency. Ginzburg recommends Cnu = 5.5
        self.cnu = np.double(cnu)
        
        # Total number of the rays
        self.nrays = nx*ny    # = grid[1]*grid[0]
        
        # Coordinates of the rays
        self.pos = np.empty((ny, nx, 3), dtype=np.double)
        
        # Direction vectors of the rays
        self.dir = np.empty((ny, nx, 3), dtype=np.double)
        print 'ny=%d, nx=%d' % (ny, nx)

        # Initial ray increments (steps)
        self.ds = np.empty((ny, nx), dtype=np.double)
        self.ds[:] = self.DeltaS

        # Critical plasma density where ray at freq cannot penetrate
        # Units: g cm^-3
        self.rhocr = np.pi*self.const.ProtonMass_g* \
                     self.const.ElectronMass_g* \
                     pow(self.freq/self.const.ProtonChargeCGSe,2)
        self.rhocr_inv = 1.0/self.rhocr

        # Initial tolerance in radians per ray step
        self.tol = np.double(0.005)
        self.tol2 = self.tol**2

        # Tolerance and maximum number of iterations for finding
        # the critical surface with the Newton method
        self.toleps = np.double(1e-6)
        self.cntmax = np.double(50.)

        # Absolute minimum step for the rays
        #self.AbsMinStep = 1e-4*sum(self.ds)/self.nrays
        self.absminstep = 1e-4*np.average(self.ds)

        ## Critical plasma density where ray at freq cannot penetrate
        ## Units: g cm^-3
        ##self.rhocr = np.pi*self.const.ProtonMass_g*self.const.ElectronMass_g*\
        ##             pow(self.freq/self.const.ProtonChargeCGSe,2)

        # Radius of the integrarion sphere
        self.rsph = np.double(rsph)
        
        # Maximum number of iterations (steps)
        self.niter = 100
        
        # Brightness temperature
        self.flags = np.empty((ny, nx), dtype=np.short)
        self.flags[:] = 0   # All the rays are ACTIVE at start

        # Number of active rays
        self.nactive = self.nrays

        # Index of the ray closest to the sun
        self.irayminsd = np.int32(-1)

        # Minimum solar distance
        self.minsoldist = np.double(1e16)  # Just a big number

        # The brightness temperature, Tb, calculated along the rays
        if self.mode == 2:
            self.tbr = np.empty((ny, nx), dtype=np.double)
            self.tbr[:] = np.double(0.0)
        else: self.tbr = np.empty(0, dtype=np.double)

        #  TbrIQUV: Tb calculated for the Stokes parameters, I, Q, U, and V
        if self.mode == 3:
            self.tbriquv = np.empty((ny, nx, 4), dtype=np.double)
            self.tbriquv[:] = np.double(0.0)
        else:
            self.tbriquv = np.empty(0, dtype=np.double)

        #  Units: sets the units as either Kelvin or Janskeys
        self.units = units

        #  Scattering: sets whether or not the program should scatter rays
        if(scattering == False):
            self.scattering = np.int32(0)
        else:
            self.scattering = np.int32(1)

        #  Magnetic field that causes the polarization, Bfield
        if self.mode == 3:
            self.bfield = np.empty((ny, nx, 3), dtype=np.double)
            self.bfield[:] = 0.0
        else: self.bfield = np.empty(0, dtype=np.double)

        # Optical depths
        self.opdepth = np.empty((ny, nx), dtype=np.double)
        self.opdepth[:] = 0.0
        
        # Plasma densities
        self.rho = np.empty((ny, nx), dtype=np.double)
        self.rho[:] = 0.0
        
        # Gradients of plasma densities
        self.gradrho = np.empty((ny, nx, 3), dtype=np.double)
        self.gradrho[:] = 0.0
        
        # Coordinates of the rays at previous step
        self.pospr = np.empty((ny, nx, 3), dtype=np.double)
        self.pospr[:] = 0.0
        
        # Directions of the rays at previous step
        self.dirpr = np.empty((ny, nx, 3), dtype=np.double)
        self.dirpr[:] = 0.0
        
        # New ray increments (steps)
        self.dsnew = np.empty((ny, nx), dtype=np.double)
        self.dsnew[:] = 0.0
        
        # Estimated ray distances to the critical surface
        self.distcr = np.empty((ny, nx), dtype=np.double)
        self.distcr[:] = 0.0
        
        # Position of the observer (the earth's coordinates) in SGI system
        self.obs = np.empty(3, dtype=np.double)
        self.obs[:] = obs

        # Pixel X and Y increments on the image plane
        self.dx = np.double(rect[2] - rect[0])/self.nx
        self.dy = np.double(rect[3] - rect[1])/self.ny


        # Rulers: arrays of x and y coordinates of all the pixels
        self.xruler = np.linspace(rect[0]+0.5*self.dx, rect[2]-0.5*self.dx, \
                                  self.nx)  #, dtype=np.double)
        self.yruler = np.linspace(rect[1]+0.5*self.dy, rect[3]-0.5*self.dy, \
                                  self.ny)  #, dtype=np.double)


        self.theta = []
        self.phi = []
        self.orientation = []
        self.density = []
        self.baseStrength = []
        self.stalkStrength = []
        self.scale = []

        # Miscellaneous parameters.
        # The parameters are precalculated to lower the computational
        # overhead in the inner loops of the ray tracing algorithm.
        e2 = (self.const.ProtonChargeCGSe)**2   # e^2
        e2w2 = e2/self.omega2                   # (e/w)^2
        e_ovr_mcw2 = e2w2/((self.const.ElectronMass_g*                  \
                           self.const.c_light_cms)**2)   # (e/mcw)^2
        e2_4pi_ovr_m = 4.*np.pi*e2/self.const.ElectronMass_g # 4pi e^2/m
        e2_4pi_ovr_mw2 = 4.*np.pi*e2w2/self.const.ElectronMass_g # 4pi e^2/m w^2
        twokf2c2 = 2.*self.const.Boltzmannk_ergK*self.freq2/            \
                   self.const.c_light_cms # Rayleigh-Jeans factor, 2kf^2/c^2
        lnLambda_13p7 = 13.7*self.const.Coulomb_log
        
        # Calculate position, pos, and direction, dir, of ray vectors
        # at the ray intersections with the integrarion sphere
        # Old method (too slow):
        ##     init_posdir(obs, rect, grid, self.rsph,
        ##                 self.pos, self.dir)

        # self.isec[i] is number of intersections with sphere of i-th ray
        self.isec = np.empty((ny, nx), dtype=np.short)

        ## First, the directions of the rays:
        ##rtcore.raydir(self.obs, self.xruler, self.yruler, self.dir)
        
        ## Then, based on directions, the ray positions on integration sphere
        ## self.nisec is set to total number of intersections with sphere

        ##self.nisec = rtcore.raysph(self.obs, self.dir, self.rsph, self.isec,
        ##                           self.pos)

        # Calculate position, pos, and direction, dir, of ray vectors
        # at the ray intersections with the integrarion sphere
        self.nisec = rtcore.raypd(self.obs, self.xruler, self.yruler,
                                  self.rsph, self.isec, self.dir, self.pos) 

        if self.nisec < 2*self.nrays:
            print 'Warning: Some rays are outside the integration sphere'

        if self.nisec == 0:
            print 'Warning: No ray intersections with the integration sphere'
        
        # Number of points whose trajectories must be saved
        self.ntrj = np.int32(0)
        
        # Array of coordinate triplets, trajectories[ntrj,niter,3].
        # It is created if the integer array tpts[ntrj,2] of
        # image plane coordinates of the rays whose trajectories must be
        # saved is passed in the call to self.trace(). 
        self.traj = np.empty(0, dtype=float)

        # Array of integer coordinates within the image plane,
        # tpts[ntrj,2], of the rays 
        # whose trajectories to be recorded in self.trajectories array.
        # It is created after the namesake parameter, tpts, is
        # passed at the call to self.trace() method.  
        self.tpts = np.empty(0, dtype=np.int32)

        # The 1D array of integer coordinates in the image plane treat
        self.ppts = np.empty(0, dtype=np.int32)
       
        # Maybe, in some time I will understand the wisdom of introducing
        # this array. However, right now it is beyond my comprehension...
        self.lastStep = np.empty(0, dtype=np.int32)
        
        # Plasma parameters dynamically linked function
        cwd = os.getcwd()
        self.plfname = ''
        plf_soname = glob.glob('plasma_*.so')
        if len(plf_soname) > 0: # plasma_*.so file(s) exist in cwd
            plfname = os.path.basename(plf_soname[0])
            self.plfname = plfname.split('.')[0] # Only the name before .so
            print 'Warning: Plasma parameters function set to '+plfname
        else:  # No plasma_*.so file(s) in cwd
            print 'Warning: no plasma_*.so file in current directory'
            plf_cname = glob.glob('plasma_*.cu') 
            if len(plf_cname) > 0: # A plasma_*.cu file exists
                plf_cname = plf_cname[0]
                plfname = os.path.basename(plf_cname)
                print 'Warn - A plasma_*.cu file exists: ', plf_cname
                self.plfname = plfname.split('.')[0] # Only the name before .cu
                print self.plfname
                print 'plf_cname = ', plf_cname
                self.set_plfunc(plf_cname)
                #print 'Warning: Plasma parameters function set to '+plfname
                #print 'Warning: To provide another plasma parameters function'
                #print 'Warning: use .set_plfunc(file) method.'
            else:  # Neither plasma_*.c nor .so file(s) in cwd
                print 'Warning: no plasma_*.c file in current directory'
                print 'Warning: You need to provide plasma parameters function'
                print 'Warning: using .set_plfunc(file) method.'


        # A single array for parameters: phisical constants, algorithmic
        # settings and work variables
        self.arprm = np.array((
            self.DeltaS,     # Initial ray path increment in solar radii 
            #self.rsph+0.01,  # Radius of the integrarion sphere + a bit more
            self.freq,       # Hz, radio wave frequency
            self.omega,      # rad/s, radio wave frequency
            self.freq2,      # Hz^2, radio wave frequency squared
            self.omega2,     # (rad/s)^2, radio wave frequency squared
            self.k0,         # k0 = omega/c_light_cms, wave number in a vacuum
            self.k0_Rsun,    # k0 = omega/c_light_Rsun, wave # in rad/Rsun
            self.rhocr,           # Critical density at given self.freq
            self.rhocr_inv,  # 1/self.RhoCr
            self.tol,        # Maximum radians per one step 
            self.tol2,       # Tolerance self.tol squared (tol2 = tol^2)
            self.absminstep, # Limit to adaptive step decrease
            self.toleps,     # Tolerance for dielectric permittivity
            self.cntmax,      # Upper limit to Newton iterations number
            self.const.h_chromo_km,     # km, Chromosphere thickness
            self.const.r_chromo,   # in Rsun units: "top" of chromosphere
            self.const.r_chro_cor, # in Rsun units: chromosphere-corona "border"
            self.const.r_corona,   # in Rsun units: "bottom" of corona
            self.const.c_light_cms,      # cm/s, Speed of light
            self.const.ProtonChargeCGSe,   # StatCoulombs, CGSe
            self.const.ProtonMass_g,       # g
            self.const.ProtonMassInv,      # 1/g
            self.const.ElectronMass_g,     # g
            self.const.Boltzmannk_ergK,    # erg/K, Boltzmann constant 
            self.const.Rsun_cm,         # cm, Solar radius
            self.const.Rsun_km,         # km, Solar radius
            self.const.Te_corona_K,     # K
            self.const.Te_chromo_K,     # K
            self.const.AU_m,      # Sun-earth distance (1 AU) in meters
            self.const.AU_Rsun,   # Sun-earth distance (1 AU) in solar radii
            self.const.Msun_G,    # G/Rsun^3, solar dipole field at equator
            self.const.Mdir[0],         # Solar dipole x direction, CGI
            self.const.Mdir[1],         # Solar dipole y direction, CGI
            self.const.Mdir[2],         # Solar dipole z direction, CGI
            e_ovr_mcw2,           # (e/mcw)^2
            e2_4pi_ovr_m,
            e2_4pi_ovr_mw2,
            twokf2c2,
            lnLambda_13p7,
            self.cnu,             # Coef. at Ginzburg's nu_eff
            self.callcount        # Number of calls to advance_beam()
            ), dtype=np.double)


    def trace(self, niter=None, tpts=None):
        """
        Trace the beam of rays comprising the nx by ny image plane grid
        for niter steps.

        parameter: niter-The maximum iterations to perform using the algorithm.
        parameter: tpts-A 2D array of integer coordinates [ix,iy] on the
        image plane to trace throughout the algorithm. Note: ix is horizontal,
        and iy is vertical dimension, so the 2D image plane arrays are a
        accessed as A[iy,ix].
        
        ## returns: if tpts was not None, then the trajectories of all
        ## of those points are returned as a 3D array ([x,y,step]). Otherwise
        ## nothing is returned.        
        """
        errmsg = 'ERROR: tpts must be an array-like structure\n'   \
                 '       convertible to a Nx2 integer array, \n' \
                 '       or a string "all"'
        
        if niter is None: niter = self.niter

        if tpts is not None:
            if isinstance(tpts, str) and tpts.lower() == 'all':
                self.ntrj = self.nx*self.ny
                self.tpts = np.zeros((self.ntrj,2), dtype=int)
                k = 0
                for i in xrange(self.ny):
                    for j in xrange(self.nx):
                        self.tpts[k,:] = i, j
                        k = k + 1
            else:
                try:
                    self.tpts = np.array(tpts, dtype=np.int32)
                except ValueError:
                    print errmsg
                    print 'tpts = ', tpts
                    return
                if self.tpts.ndim <> 2:
                    print errmsg
                    print 'tpts = ', tpts
                    return
                elif len(self.tpts[0,:]) <> 2:
                    print errmsg
                    print 'tpts = ', tpts
                    return
                
            self.ntrj = np.int32(np.size(self.tpts,0))
            ntrj = self.ntrj
            self.ppts = np.empty(ntrj, dtype=np.int32)
            for i in range(ntrj):
                iy = self.tpts[i,0]
                ix = self.tpts[i,1]
                self.ppts[i] = iy*self.nx + ix
                self.flags[iy,ix] |= self.TRACED   # Flag the traced rays
            self.traj = np.empty((ntrj,niter,3), dtype=np.double)
            self.traj[:,:,:] = float('nan')

            self.lastStep = np.empty(ntrj, dtype=np.int32)
            self.lastStep[:] = 0
            
        self.theta.append(float('nan'))
        theta = np.empty((len(self.theta)),dtype=np.double)
        theta[:] = self.theta[:]
        self.phi.append(float('nan'))
        phi = np.empty((len(self.phi)),dtype=np.double)
        phi[:] = self.phi[:]
        self.orientation.append(float('nan'))
        orientation = np.empty((len(self.orientation)),dtype=np.double)
        orientation[:] = self.orientation[:]
        self.density.append(float('nan'))
        density = np.empty((len(self.density)),dtype=np.double)
        density[:] = self.density[:]
        self.baseStrength.append(float('nan'))
        baseStrength = np.empty((len(self.baseStrength)),dtype=np.double)
        baseStrength[:] = self.baseStrength[:]
        self.stalkStrength.append(float('nan'))
        stalkStrength = np.empty((len(self.stalkStrength)),dtype=np.double)
        stalkStrength[:] = self.stalkStrength[:]
        self.scale.append(float('nan'))
        scale = np.empty((len(self.scale)),dtype=np.double)
        scale[:] = self.scale[:]






        rtcore.trace_beam(
            self.arprm,
            self.pos,
            self.dir,
            self.ds,
            self.rsph,
            niter,
            self.mode,
            self.scattering,
            self.flags,
            self.tbr,
            self.tbriquv,
            self.opdepth,
            self.rho,
            self.gradrho,
            self.bfield,
            self.pospr,
            self.dirpr,
            self.dsnew,
            self.distcr,
            self.ppts,
            self.traj,
            self.lastStep,
            theta,
            phi,
            orientation,
            density,
            baseStrength,
            stalkStrength,
            scale)
        
        if(self.units == 'Jy'):
            dist = np.sqrt(np.dot(self.obs,self.obs))
            alpha = (self.rect[2]-self.rect[0])/self.nx/dist
            beta = (self.rect[3]-self.rect[1])/self.ny/dist
            solid_angle = alpha*beta
            self.tbr = 2e26*self.const.Boltzmannk_SI/ \
                       (self.const.c_light_SI**2)*    \
            solid_angle*self.freq**2*self.tbr
            self.tbriquv = 2e26*self.const.Boltzmannk_SI/ \
                           (self.const.c_light_SI**2)*    \
                           solid_angle*self.freq**2*self.tbriquv
            
            pl.figure(); pl.imshow(t.tbr, interpolation='none', origin='lower')


        ## if ntrj != 0:
        ##     return self.traj
        

    def soldist(self):
        """
        Calculate the solar distance of all the rays
        """
        dis = np.empty((self.ny,self.nx), dtype=np.double)
        for i in xrange(self.ny):
            for j in xrange(self.nx):
                dis[i,j] = vmagn(self.pos[i,j,:])
        return dis

    def set_plfunc(self, fname):
        # Check if the file exist
        if not os.path.isfile(fname):
            print 'Error: no such file "'+fname+'".'
            return
        # File exists at this point.
        bfname = os.path.basename(fname) # Strip the base name of dirs
        name = bfname.split('.')[0] # Remove .cu or .so from base file name
        # Prepare compilation and linking strings
        nvcc_compile = 'nvcc -arch=sm_20 -g -Xcompiler -fPIC  -I'+self.package+\
                      '/inc -c '+fname+' -o '+name+'.o'
        nvcc_link = 'nvcc  -arch=sm_20 -Xcompiler -shared '+name+'.o -L'+ \
                   self.package+'/lib -L/usr/lib -lm -lmxv -o '+name+'.so'
        # If it is name.cu, recompile
        if bfname[-3:] == '.cu':
            print 'bfname[-2:] = ', bfname[-2:]
            # Compile .cu file into shared library in current directory
            print nvcc_compile
            os.system(nvcc_compile)
            print nvcc_link
            os.system(nvcc_link)
            self.plfname = name
            return
        
        # If it is name.so, check which is newer, fname.cu or name.so
        # If .cu is newer, recompile
        if bfname[-3:] == '.so': # Shared library
            if bfname == fname:  # Name.so is in current directory
                cname = name+'.cu'
                if os.path.isfile(cname): # If name.cu is in current directory
                    # Recompile, if name.so is older than name.cu
                    if os.path.getctime(cname) > os.path.getctime(bfname):
                        print nvcc_compile
                        os.system(nvcc_compile)
                        print nvcc_link
                        os.system(nvcc_link)
            # Do nothing if name.so is not in current directory
            self.plfname = name
            return
                   
        # The file is neither name.cu nor name.so
        print "Error: file is neither name.cu nor name.so"
        return


    def make_streamer(self,theta,phi,
            orientation=0.0,
            density=2.0,
            baseStrength=5.0,
            stalkStrength=2.0,
            scaleX = .75,
            scaleY = .75,
            scaleZ = .75):
        """
    This function adds a streamer definition to the plasma parameters provided.
    Note that this streamer will apply to all implane's created and that
    remove_streamers must be called to remove them.

    parameter: theta-Colatitude of the streamer (degrees)
    parameter: phi-Azimuth of the streamer (degrees)
    parameter: orientation-The orientation of the monopoles of the streamer in
    degrees where at theta=90 phi=0 orientation=0 the axis of the streamer is
    parallel to the y-axis, and at theta=90 phi=0 orientation=90 the axis of
    the streamer is parallel to the z-axis
    parameter: density-The density of the streamer
    parameter: baseStrength-The magnetic field of the monopoles
    parameter: stalkStrength-The magnetic field of the stalks
    parameter: scales-The scales in the respective directions
    parameter: plfname-The name of the plasma parameters file
    """

        degtorad = 3.14159/180.01
        self.theta.append(degtorad*theta)
        self.phi.append(degtorad*phi)
        self.orientation.append(degtorad*orientation)
        self.density.append(density)
        self.baseStrength.append(baseStrength)
        self.stalkStrength.append(stalkStrength)
        self.scale.append(scaleX)


    def remove_streamers(self):
        """
    This function removes all streamers from the plasma parameters provided.

    parameter: plfname-The name of the plasma parameters file
        """
        self.theta = []
        self.phi = []
        self.orientation = []
        self.density = []
        self.baseStrength = []
        self.stalkStrength = []
        self.scale = []









    def plprofile(self, pos):
        """
        Given an array of positions this returns the density, density gradient,
        and B field at those positions.

        parameter: pos-The array of coordinates.
        returns: A tupple of arrays (rho,gradrho,bfield) where rho is a 2D array
        of the density at any [x,y]. And graderho and bfield are both 3D arrays
        indexed by [x,y,component] where component=0 is x, where 1 is y and
        2 is z.
        """

        pos = np.array(pos)

        if(np.size(pos,1)!=3):
            print("Error, array of coordinates required")
            return


        rho = np.empty((np.size(pos,0)), dtype=np.double)

        gradrho = np.empty((np.size(pos,0), 3), dtype=np.double)

        bfield = np.empty((np.size(pos,0), 3), dtype=np.double)

        ds = np.empty(np.size(pos,0) , dtype=np.double)

        flags = np.empty(np.size(pos,0) , dtype=np.short)

        rtcore.plprofile(
                self.arprm,
                pos,
                rho,
                gradrho,
                bfield,
                self.mode,
                ds,
                flags)
        return rho,gradrho,bfield



def init_posdir(obs, rect, grid, rsph, pos, dir):
    """
    For the observer SGI coordinates in obs and the image plane
    rectangle coordinates rect = [xleft, ybottom, xright, ytop],
    (or, maybe clearer, rect = [xlower, ylower, xupper, yupper],
    calculates position, pos, and direction, dir, of ray vectors
    at the ray intersections with the integrarion sphere with
    the center at SGI coordinates origin (the center of the sun).
    This calculation is based on the assumption that the rays
    outside the intehration sphere are not refracted and are
    actually the straight lines drawn between the observer point
    and every pixel in the image plane. The image plane center
    coincides with the SGI coordinates origin (the center of the
    sun). The image plane is normal to the line between the
    observer and the SGI coordinates origin. 
    """
    nx = grid[1]
    ny = grid[0]

    # HERE MUST BE CHECK FOR CORRECT DIMENSIONS of pos, dir etc.

    #  
    # Find the basis, (tau,xi), of the inner coordinate
    # system in the image plane
    #
    norm2impl = obs/vmagn(obs)               # n is unity normal to image plane
    tau = np.cross((0.,0.,1.), norm2impl) # tau = ez x n
    tau = tau/vmagn(tau) # tau is a unity-length ort, perp to SGI Z and n
    xi =  np.cross(norm2impl, tau)        # xi =  n x tau, unity-length ort  
    #print 'vmagn(tau) = ', vmagn(tau), ', vmagn(xi) = ', vmagn(xi)
    print 'rsph = ', rsph

    #
    # Determine the image plane inner coordinates of the pixel centers
    #
    xlw = rect[0]
    ylw = rect[1]
    xup = rect[2]
    yup = rect[3]
    dx = (xup - xlw)/nx
    dy = (yup - ylw)/ny
    xruler = np.linspace(xlw+0.5*dx, xup-0.5*dx, nx)
    yruler = np.linspace(ylw+0.5*dy, yup-0.5*dy, ny)
    xpix, ypix = np.meshgrid(xruler, yruler)

    #
    # Calculate the SGI coordinates of all the image plane pixels
    #
    for i in xrange(ny):
        for j in xrange(nx):
            rpix = xpix[i,j]*tau + ypix[i,j]*xi # A pixel-vector
            slope = rpix - obs
            #print 'rpix = ', rpix, ', slope = ', slope
            dir[i,j,:] = slope/vmagn(slope) # Directions of Obs-to-pixel, SGI

    #
    # Find the points on the integration sphere where it intersects 
    # with the straight "rays" 
    #
    rsph2 = rsph**2
    for i in xrange(ny):
        for j in xrange(nx):
            v = dir[i,j,:]
            vxo = np.cross(v,obs) # dir x obs
            vdo = np.dot(v,obs)   # dir . obs
            discr = rsph2 - np.dot(vxo,vxo)
            if discr < 0:
                print 'Increase integration sphere radius rsph' 
                raise XCP
            discr = np.sqrt(discr)
            obs2sph = -vdo - discr # Observer-to-integration-sphere distance
            #print 'rsph2 - np.dot(vxo,vxo) = ', rsph2 - np.dot(vxo,vxo)
            #print 'discr = ', discr
            pos[i,j,:] = obs + v*obs2sph



def vmagn(vect):
    """
    Magnitude of a vector as sqrt(dot(vect,vect)).
    Is written exclusively for brevity :)
    """
    return np.sqrt(np.dot(vect, vect))



def run_freq_range(file,
        freq_start=80e6,
        freq_end=300e6,
        steps=12,
        grid=(20,20),
        rect=(-2., -2, 2., 2.),
        obs=(215., 0., 0.),
        rsph = 25.,
        mode='TbrIQUV',
        cnu=3.0,
        msun=1.0,
        units = 'K',
        theta = None,
        phi = None,
        orientation = 0,
        density = 2.0,
        baseStrength = 5.0,
        stalkStrength = 2.0,
        scattering=False):
    """
    This will run a batch job for many frequencies and export the results
    to a fits file.

    parameter: freq_start-The begining of the frequency range to run (Hertz)
    parameter: freq_end-The end of the frequency range to run (Hertz)
    parameter: steps-The number of frequencies between start and end to run.
    parameter: grid-The number of rays to run (x-axis rays,y-axis rays)
    parameter: rect-The imageplane rectangle i.e. the box around the sun from
    which the rays are launched (in solar radii)
    parameter: obs-The location of the earth in terms of solar radii
    parameter: rsph-The size of the integration sphere in terms of solar radii
    parameter: mode-The mode in which to run the simulation.
    parameter: units-K for Kelvin Jy for Janskeys
    parameter: theta,phi,orientation,density,baseStrength,stalkStrength-The
    streamer to put on these simulations, see make_streamer for details on
    the parameters
    parameter: scattering-True to include scattering, False to omit it.
    """

    nx = grid[1]
    ny = grid[0]

    dist = np.sqrt(np.dot(obs,obs))

    const = rtconstants()

    alpha = (rect[2]-rect[0])/nx/dist
    beta =  (rect[3]-rect[1])/ny/dist
    solid_angle = alpha*beta
    conversion = 2e26*const.Boltzmannk_SI/(const.c_light_SI**2)*solid_angle

    if(mode=='Tbr'):
        data = np.empty((1,steps,ny,nx))
    if(mode=='TbrIQUV'):
        data = np.empty((1,steps,ny,nx,2))
    index = 0

    for i in np.linspace(freq_start,freq_end,steps):
        a = implane(grid,rect,obs,rsph,i,mode,cnu,msun,
                units = units,scattering=scattering)
        remove_streamers()
        if(theta!=None and phi!=None and orientation!=None):
            make_streamer(theta,phi,orientation,
                    density,baseStrength,stalkStrength)
            a.trace(5000)


        if(mode=='Tbr'):
            # Make unfinished pixels equal to 0
            a.tbr[np.where(a.flags == 0)] = 0.

            endcap = 1
            data[0,index,:,:] = a.tbr
            if(np.isnan(data[0,index,:,:]).any()):
                continue

        if(mode=='TbrIQUV'):
            endcap = 2
            data[0,index,:,:,0] = a.tbriquv[:,:,0]
            data[0,index,:,:,1] = a.tbriquv[:,:,3]
            if(np.isnan(data[0,index,:,:,0]).any() or
                    np.isnan(data[0,index,:,:,1]).any()):
                continue
        index = index + 1




    for j in [0]:

        if(mode=='Tbr'):
            hdu = pyfits.PrimaryHDU(data[:,:,:,:])
        if(mode=='TbrIQUV'):
            hdu = pyfits.PrimaryHDU(data[:,:,:,:,j])


        prihdr = hdu.header
        prihdr.update('BITPIX', -64, 'FITS BITS/PIXEL') #Not sure
        prihdr.update('NAXIS', 4, 'NUMBER OF AXES')
        prihdr.update('NAXIS1', np.size(data,2), 'X pixels')
        prihdr.update('NAXIS2', np.size(data,3) , 'Y pixels')
        prihdr.update('NAXIS3', steps, 'Number of frequencies')
        prihdr.update('NAXIS4', 1, 'Stokes')

        prihdr.update('OBJECT', 'Sun', 'Source Name')

        prihdr.update('TELESCOP', 'Ray Sim','')
        prihdr.update('INSTRUME', 'HART','')
        prihdr.update('OBSERVER', '','')
        prihdr.update('DATE-OBS', str(date.today()),'')

        prihdr.update('BSCALE', 1.0,'')
        prihdr.update('BZERO', 0.0,'')
        if(units=='K'):
            prihdr.update('BUNIT', 'K', 'Temperature Brightness')
        if(units=='Jy'):
            prihdr.update('BUNIT', 'JY', 'Temperature Brightness')


        #prihdr.update('DATAMAX', ''+np.max(data[:,:,:,:,j]),
        #              'Maximum Pixel Value')
        #prihdr.update('DATAMIN', ''+np.min(data[:,:,:,:,j]),
        #              'Minimum Pixel Value')



        prihdr.update('CTYPE1', 'RA---SIN', 'Coordinates along the ecliptic')
        prihdr.update('CRVAL1', 0.0, '') 
        prihdr.update('CDELT1', (rect[2]-rect[0])/float(np.size(data,2))/dist*180/np.pi, 'Degrees/Pixel')
        prihdr.update('CRPIX1', (-rect[0])/(rect[2]-rect[0])*nx, '')


        prihdr.update('CTYPE2', 'DEC--SIN', 'Coordinates perpendicular to the ecliptic')
        prihdr.update('CRVAL2', 0.0, '')
        prihdr.update('CDELT2', (rect[3]-rect[1])/float(np.size(data,3))/dist*180/np.pi, 'Degrees/Pixel')
        prihdr.update('CRPIX2', (-rect[1])/(rect[3]-rect[1])*ny, '')


        prihdr.update('CTYPE3', 'FREQ', '')
        prihdr.update('CRVAL3', freq_start, 'Frequency in Hz')
        prihdr.update('CDELT3', (freq_end-freq_start)/(steps-1), 'Frequency per index')
        prihdr.update('CRPIX3', 1.0, 'Reference Freq index')



        prihdr.update('CTYPE4', 'STOKES', '')
        prihdr.update('CRVAL4', 1.0+j*3.0, '')
        prihdr.update('CDELT4', 1.0, '')
        prihdr.update('CRPIX4', 1.0, '')
        prihdr.update('CROTA4', 0.0, '')



        prihdr.update('ORIGIN', 'MIT Haystack Observatory', 'Organization')



        if(units == 'K'):
            prihdr.update('KTOJY', conversion , 'This times nu^2 gives the conversion to Jy')
        else:
            prihdr.update('JyTOK', (1/conversion) , 'This divided by nu^2 gives the conversion to K')

        if(theta!=None and phi!=None and orientation!=None):
            prihdr.add_comment('Streamer coordinates: (theta,phi) ('+
                    str(theta)+','+str(phi)+')')
            prihdr.add_comment('Streamer Orientation: '+str(orientation))
            prihdr.add_comment('Streamer Density: '+str(density))
            prihdr.add_comment('Streamer Base Strength: '+str(baseStrength))
            prihdr.add_comment('Streamer Stalk Strength: '+str(stalkStrength))

        hdulist = pyfits.HDUList([hdu])
        if(j==0):
            hdulist.writeto(file+'_I.fits')
        else:
            hdulist.writeto(file+'_V.fits')









