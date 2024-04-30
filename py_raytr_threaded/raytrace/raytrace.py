###########################################################
#                                                         #
# raytrace.py                                             #
#                                                         #
# Module for making images  of space objects              #
# through ray tracing                                     #
# Created: 14 March 2011 by Leonid Bevkevitch             #
###########################################################
#
import pylab as pl
import numpy as np
import numpy.linalg as npl
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


class trajstore():
    """
    Selected ray traces and accompanying parameters 
    """
    def __init__(self, rays=None, grid=None, store=None, npmax=3000):
        """
        rays: 2-Dimensional array or array-like sequence, rays[nr,2],
          of the selected ray indices
          Example: rays=[[2,5],[3,3],[3,5],[5,5]]
        grid: 2-number sequence, numbers of pixels [n_RA,n_DEC] or [n_X,n_Y]   
        store: a sequence of strings with the names of accompanying
          parameters to store along with the ray coordinates.
          Example: store=('Dir', 'Rho', 'gRho', 'Bfield', 'IQUV')
        npmax: maximum number of trajectory points to store.
        """
        if rays == None:
            print 'Error: first argument, rays, must be array-like rays[nr,2]'
            return
        if grid == None:
            print 'Error: grid parameter must be present'
            return
        nx = int(grid[0])
        ny = int(grid[1])
        rays = np.array(rays, dtype=np.int32)
        nr = rays.shape[0] # Number of rays to track
        npmax = int(npmax)
        
        # Always store:
        self.npmax = np.int32(npmax) # max num of trajectory points to store.
        self.ntraj  = np.int32(nr)   # Number of rays to track.
        self.pos = np.empty((nr,npmax,3), dtype=np.double) # Positions
        self.pos[:] = np.NaN
        self.last = np.empty(nr, dtype=np.int32) # Last indices for all rays
        self.last[:] = 0
        self.irays = np.empty(nr, dtype=np.int32) # 1D ray indices

        # Find 1D indices into 2D and 3D beam data arrays
        for i in range(nr):
            self.irays[i] = rays[i,1]*nx + rays[i,0]
                 
        if store == None: return
        
        store = [s.lower() for s in store] # Drop all letters to lowercase
        
        if 'arclen' in store:
            self.arclen = np.empty((nr,npmax), dtype=np.double)
            self.arclen[:] = np.NaN           
        if 'angbd' in store:
            self.angbd = np.empty((nr,npmax), dtype=np.double)
            self.angbd[:] = np.NaN           
        if 'dist' in store:
            self.dist = np.empty((nr,npmax), dtype=np.double)
            self.dist[:] = np.NaN           
        if 'dir' in store:
            self.dir = np.empty((nr,npmax,3), dtype=np.double)
            self.dir[:] = np.NaN
        if 'rho' in store:
            self.rho = np.empty((nr,npmax), dtype=np.double)
            self.rho[:] = np.NaN
        if 'grho' in store:
            self.grho = np.empty((nr,npmax,3), dtype=np.double)
            self.grho[:] = np.NaN
        if 'grhomagn' in store:
            self.grhomagn = np.empty((nr,npmax), dtype=np.double)
            self.grhomagn[:] = np.NaN
        if 'bfield' in store:
            self.bfield = np.empty((nr,npmax,3), dtype=np.double)
            self.bfield[:] = np.NaN
        if 'bmagn' in store:
            self.bmagn = np.empty((nr,npmax), dtype=np.double)
            self.bmagn[:] = np.NaN
        if 'opdepth' in store or 'tau' in store:
            self.opdepth = np.empty((nr,npmax), dtype=np.double)
            self.opdepth[:] = np.NaN
        if 'stokes' in store or 'iquv' in store or 'tbriquv' in store:
            self.tbriquv = np.empty((nr,npmax,4), dtype=np.double)
            self.tbriquv[:] = np.NaN
        elif 'tbr' in store or 'tbri' in store or 'i' in store:
            self.tbr = np.empty((nr,npmax), dtype=np.double)
            self.tbr[:] = np.NaN

        


class implane():
    """
    Image plane representation.
    """
    def __init__(self, grid=(10,10), rect=(-5., -5, 5., 5.),
                 obs=(215., 0., 0.), rsph = 25., freq=100e6,
                 mode='Tbr', cnu=3.0, msun=1.0,
                 trkrays=None, trkparms=None, trknpmax=None,
                 units='K', scattering=False):
        """
        Create an image plane as a two-dimensional array of doubles.
        The size in pixels is RA x Dec = grid[0] x grid[1].
        """
        rect = map(np.double, rect)
        obs = map(np.double, obs)
        grid = map(int, grid)

        # Grid sizes along (RA,Dec)
        self.grid = grid

        # Grid size along RA
        self.nx = nx = grid[0]
        
        # Grid size along Dec
        self.ny = ny = grid[1]
        
        # Indicator of the ray tracing
        self.inprocess = False # True: self.trace() has been called

        # Maximum number of iterations (steps)
        self.niter = 1000

        # Maximum number of points to store in trajectories (if traj=True)
        self.npmax = self.niter

        # Number of stored trajectories
        self.ntraj = 0

        # Tracked rays' image plane coordinates
        if trkrays != None:
            self.trkrays = np.array(trkrays, dtype=np.int32)
        else:
            self.trkrays = None
        
        # Tracked rays' image plane coordinates
        self.trkparms = trkparms

        #
        # Parameters allowed for tracking dependent on the mode
        #
        # Basic set:
        self.allowedtotrack = set(('pos', 'dir', 'arclen', 'rho', 'grho',
                                   'grhomagn', 'dist', 'opdepth', 'tau'))
        
        if mode == 'basic' or mode == 1:
            self.mode = 1     # No Tbr or TbrIQUV calculation;
        elif mode == 'Tbr' or mode == 2:
            self.mode = 2     # Tbr calculation;
            self.allowedtotrack.update(('tbr', 'tbri', 'i'))
        elif mode == 'TbrIQUV' or mode == 3:
            self.mode = 3     # No Tbr or TbrIQUV calculation;
            self.allowedtotrack.update(('bfield', 'bmagn', 'angbd', 'tbriquv',
                                        'iquv', 'stokes'))
        else:
            print 'Error: mode parameter can be only "basic" (or 1), ' \
                  '"Tbr" (or 2), or "TbrIQUV" (or 3). You gave mode = ', mode
            return

        
        #
        # If trkrays given, create the self.traj structure to save trkparms
        #
        # Check if trkrays is array-like
        if trkrays != None:
            if trkparms != None:
                trkparms = [s.lower() for s in trkparms] # Drop to lowercase
            else:
                trkparms = ['pos'] # Always store
            trkrays = np.array(trkrays)        
            # Check if trkrays is array-like, with (ntraj,2) shape
            if np.size(trkrays.shape) == 2 and trkrays.shape[1] == 2:
                if trknpmax == None: trknpmax = self.niter
                else: self.npmax = trknpmax
                # Check if legal params requested to track
                for p in trkparms:
                    if not (p in self.allowedtotrack):
                        print 'Error: illegal parameter name "'+p+'".'
                        return
                self.traj = trajstore(rays=trkrays, grid=grid, store=trkparms,
                                  npmax=trknpmax)
            else:
                print 'Error: trkrays parameter must be array-like with ' \
                      '(ntraj,2) shape.'
                return
        else:
            self.traj = None
        

        # Exceptions
        self.XCP = Fatal_Error()

        # Determine the location of the raytrace package in the site-packages
        # (/usr/lib/python2.6/site-packages/raytrace or
        # /usr/lib64/python2.6/site-packages/raytrace)
        # NOTE: get_python_lib() returns the non arch-specific directory
        # which is under /usr/lib. One needs to use get_python_lib(1) to
        # get the platform specific directory if that's required.

        #self.package = get_python_lib(1) + '/raytrace' 

        #self.package = os.path.dirname(rtcore.__file__)+ '/raytrace'
        self.package = os.path.dirname(rtcore.__file__)

        # Bit flag positions for self.flags
        self.INACTIVE =  0x0001   # The ray is not being processed  
        self.SHARP =     0x0002   # The ray is sharp (close to normal) 
        self.PENETR =    0x0004   # Eps < 0; Critical surface search
        self.WASREFL =   0x0008   # First step after reflection
        self.BADRAY =    0x0010   # The ray is bad

        # Create structure prm containing physical constants
        self.const = rtconstants()

        # Change some constants
        self.const.Msun_G = msun  # Solar dipole strength in Gauss at equator

        # Counter for calls to the algorithm; roughly number of steps
        self.callcount = 1.0 
        
        # Constants
        self.DeltaS = 0.1   # Initial ray path increment in solar radii

        # Positions of the image plane left bottom and right top
        # corners in SGI system
        # rect=[xleft, ybottom, xright, ytop]
        self.rect = np.empty(4, dtype=np.double)
        self.rect[:] = rect
        self.xleft =   xleft =   self.rect[0]
        self.ybottom = ybottom = self.rect[1]
        self.xright =  xright =  self.rect[2]
        self.ytop =    ytop =    self.rect[3]

        # Image plane size along RA
        self.sizera = xright - xleft
        
        # Image plane size along Dec
        self.sizedec = ytop - ybottom
        
        # Frequency of the radio waves
        self.freq = freq
        self.freq2 = freq**2
        self.omega = 2.*np.pi*freq
        self.omega2 = self.omega**2

        # Wavenumber in a vacuum
        self.k0 = self.omega/self.const.c_light_cms          # rad/cm
        self.k0_Rsun = self.omega/self.const.c_light_Rsun    # rad/Rsun

        # Coefficient at the Ginzburg's formula for effective collision
        # frequency. Ginzburg recommends Cnu = 5.5
        self.cnu = cnu
        
        # Total number of the rays
        self.nrays = ny*nx    #grid[1]*grid[0]
        
        # Coordinates of the rays
        self.pos = np.empty((ny, nx, 3), dtype=np.double)
        
        # Direction vectors of the rays
        self.dir = np.empty((ny, nx, 3), dtype=np.double)

        # Ray arc lengths
        self.arclen = np.empty((ny, nx), dtype=np.double)
        self.arclen[:] = 0.0

        # Cosines of angles between ray direction and magnetic field
        self.cosbd = np.empty((ny, nx), dtype=np.double)
        self.cosbd[:] = 0.0

        # Distances of each ray from the sun
        self.dist = np.empty((ny, nx), dtype=np.double)
        self.dist[:] = 0.0
        
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
        self.tol = 0.005
        self.tol2 = self.tol**2

        # Tolerance and maximum number of iterations for finding
        # the critical surface with the Newton method
        self.toleps = 1e-6
        self.cntmax = 50.

        # Absolute minimum step for the rays
        #self.AbsMinStep = 1e-4*sum(self.ds)/self.nrays
        self.absminstep = 1e-4*np.average(self.ds)

        ## Critical plasma density where ray at freq cannot penetrate
        ## Units: g cm^-3
        ##self.rhocr = np.pi*self.const.ProtonMass_g*self.const.ElectronMass_g*\
        ##             pow(self.freq/self.const.ProtonChargeCGSe,2)

        # Radius of the integrarion sphere
        self.rsph = rsph
        
        # Brightness temperature
        self.flags = np.empty((ny, nx), dtype=np.short)
        self.flags[:] = 0   # All the rays are ACTIVE at start

        # Number of active rays
        self.nactive = self.nrays

        # Index of the ray closest to the sun
        self.irayminsd = -1

        # Minimum solar distance
        self.minsoldist = 1e16  # Just a big number

        # The brightness temperature, Tb, calculated along the rays
        if self.mode == 2:
            self.tbr = np.empty((ny, nx), dtype=np.double)
            self.tbr[:] = 0.0
        else: self.tbr = np.empty(0)

        #  TbrIQUV: Tb calculated for the Stokes parameters, I, Q, U, and V
        if self.mode == 3:
            self.tbriquv = np.empty((ny, nx, 4), dtype=np.double)
            self.tbriquv[:] = 0.0
        else:
            self.tbriquv = np.empty(0)

        #  Units: sets the units as either Kelvin or Janskeys
        self.units = units

        #  Scattering: sets whether or not the program should scatter rays
        if(scattering==False):
            self.scattering = 0
        else:
            self.scattering = 1

        #  Magnetic field that causes the polarization, Bfield
        if self.mode == 3:
            self.bfield = np.empty((ny, nx, 3), dtype=np.double)
            self.bfield[:] = 0.0
        else: self.bfield = np.empty(0)

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
        self.dx = dx = (xright - xleft)/self.nx
        self.dy = dy = (ytop - ybottom)/self.ny

        # Rulers: arrays of x and y coordinates of all the pixels
        #             rect=[xleft, ybottom, xright, ytop]
        #                     0       1       2      3
        self.xruler = np.linspace(xleft+0.5*dx, xright-0.5*dx, nx)
        self.yruler = np.linspace(ybottom+0.5*dy, ytop-0.5*dy, ny)

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

        # System type: 32 0r 64-bit
        s = os.popen('file /sbin/init').readline()
        if s.find('64-bit') >= 0:
            self.sysbits = 64
        elif s.find('32-bit') >= 0:
            self.sysbits = 32
        else:
            self.sysbits = -1 # Not determined


        # Plasma parameters dynamically linked function
        cwd = os.getcwd()
        self.plfname = ''
        plf_soname = glob.glob('plasma_*.so')
        
        # Find for what type of system, 32 or 64-bit the *.so built
        sobits = 0 # Assume *.so file does not exist
        if len(plf_soname) > 0: # plasma_*.so file(s) exist in cwd
            s = os.popen('file '+plf_soname[0]).readline() # so file type?
            if s.find('64-bit') >= 0:
                sobits = 64   # ELF type is 64-bit
            elif s.find('32-bit') >= 0:
                sobits = 32   # ELF type is 32-bit
            else:
                sobits = -1 # Not determined (0 if does not exist)
                   
        if sobits == self.sysbits: # plasma_*.so file exist and in right format
            plfname = os.path.basename(plf_soname[0])
            self.plfname = plfname.split('.')[0] # Only the name before .so
            print 'Warning: Plasma parameters function set to '+plfname
        else:  # No plasma_*.so file(s) in cwd
            print 'Warning: either no plasma_*.so file in current directory ' \
                  'or it is in wrong ELF format'
            plf_cname = glob.glob('plasma_*.c') 
            if len(plf_cname) > 0: # A plasma_*.c file exists
                plf_cname = plf_cname[0]
                plfname = os.path.basename(plf_cname)
                self.plfname = plfname.split('.')[0] # Only the name before .c
                print 'plf_cname = ', plf_cname
                self.set_plfunc(plf_cname)
                print 'Warning: Plasma parameters function set to '+plfname
                print 'Warning: To provide another plasma parameters function'
                print 'Warning: use .set_plfunc(file) method.'
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


    def trace(self, niter=None, trkrays=None, trkparms=None, npmax=None):
        """
        Trace the beam of rays comprising the nx by ny image plane grid
        for niter steps.
        Parameters:
          niter-The maximum iterations to perform using the algorithm.
          trkrays-A 2D array of coordinates on the implane to
            trace throughout the algorithm.
          trkparms: a list of parameter names to store with the trajectory.
            The names are 'Dir', 'Arclen', 'Dist', 'Rho', 'gRho', 'Bfield',
            'Tbr' or 'TbrI' or 'Stokes' or 'IQUV' or 'TbrIQUV'.
          npmax: maximum number of trajectory points to store.  
        returns: if trkrays was not None, then the trajectories of all
        of those points are returned as a 3D array ([x,y,step]). Otherwise
        nothing is returned.        
        """
        
        if niter == None:
            niter = self.niter
        else:
            self.niter = niter
        if npmax == None: npmax = self.npmax
        
        #
        # self.traj is a structure containing the arrays
        # to hold the ray trajectory points and other
        # parameters along the rays, specified in trkparms.
        # Example: trkparms=['bfield', 'arclen', 'rho', Stokes']
        #
        #
        # If trkrays given, create or replace the self.traj structure
        # to save trkparms
        #
        # Ignore trkrays, trkparms, and npmax if this call is not
        # the first call
        if (trkrays != None) and (not self.inprocess):
            if trkparms != None:
                trkparms = [s.lower() for s in trkparms] # Drop to lowercase
            else:
                trkparms = ['pos'] # Always store
            trkrays = np.array(trkrays)        
            # Check if trkrays is array-like, with (ntraj,2) shape
            if np.size(trkrays.shape) == 2 and trkrays.shape[1] == 2:
                if trknpmax == None: trknpmax = self.niter
                # Check if legal params requested to track
                for p in trkparms:
                    if not (p in self.allowedtotrack):
                        print 'Error: illegal parameter name "'+p+'".'
                        return
                self.traj = trajstore(rays=trkrays, grid=self.grid,
                                      store=trkparms, npmax=trknpmax)
            else:
                print 'Error: trkrays parameter must be array-like with ' \
                      '(ntraj,2) shape.'
                return


        rtcore.trace_beam(
            self.arprm,
            self.pos,
            self.dir,
            self.arclen,
            self.cosbd,
            self.dist,
            self.ds,
            self.plfname,
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
            self.traj
            )
        
        self.nactive = np.size(np.where((~(self.flags)) & 0x0001)[0])

        #
        # Convert units from Kelvin to Jansky, if requested
        #
        #             rect=[xleft, ybottom, xright, ytop]
        #                     0       1       2      3
        if(self.units == 'Jy'):
            dist = np.sqrt(np.dot(self.obs,self.obs))
            alpha = (self.xright - self.xleft)/self.nx/dist
            beta = (self.ytop - self.ybottom)/self.ny/dist
            solid_angle = alpha*beta
            self.tbr = 2e26*self.const.Boltzmannk_SI/ \
                       (self.const.c_light_SI**2)*    \
            solid_angle*self.freq**2*self.tbr
            self.tbriquv = 2e26*self.const.Boltzmannk_SI/ \
                           (self.const.c_light_SI**2)*    \
                           solid_angle*self.freq**2*self.tbriquv
       
        self.inprocess = True # The rays have been traced

        

    def soldist(self):
        """
        Calculate the solar distance of all the rays
        """
        dis = np.empty((self.ny,self.nx), dtype=np.double)
        for i in xrange(self.ny):
            for j in xrange(self.nx):
                dis[i,j] = npl.norm(self.pos[i,j,:])
        return dis

    def set_plfunc(self, fname):
        # Check if the file exist
        if not os.path.isfile(fname):
            print 'Error: no such file "'+fname+'".'
            return
        # File exists at this point.
        rname = os.path.realpath(fname)   # Full path for fname
        dname = os.path.dirname(rname)   # Directory name
        bfname = os.path.basename(fname) # Strip the base name of dirs
        name = bfname.split('.')[0] # Remove .c or .so from base file name
        print 'fname = ', fname
        print 'dname = ', dname
        print 'bfname = ', bfname
        print 'name = ', name
        print 'self.package = ', self.package
        # Prepare compilation and linking strings

        # Temporal patch: compile streamer.c assuming it is in same dir
        gcc_comp_streamer = 'gcc -g -fPIC  -I'+self.package+'/inc ' \
                            '-c '+dname+'/streamer.c -o streamer.o'
        gcc_compile = 'gcc -g -fPIC  -I'+self.package+ \
                      '/inc -c '+fname+' -o '+name+'.o'
        gcc_link = 'gcc -shared streamer.o '+name+'.o -L'+self.package+ \
                   '/lib -L/usr/lib -lm -lmxv -o '+name+'.so'
        # If it is name.c, recompile
        if bfname[-2:] == '.c':
            # Compile .c file into shared library in current directory
            print gcc_comp_streamer
            os.system(gcc_comp_streamer)
            print gcc_compile
            os.system(gcc_compile)
            print gcc_link
            os.system(gcc_link)
            self.plfname = name
            return
        
        # If it is name.so, check which is newer, fname.c or name.so
        # If .c is newer, recompile
        if bfname[-3:] == '.so': # Shared library
            if bfname == fname:  # Name.so is in current directory
                cname = name+'.c'
                if os.path.isfile(cname): # If name.c is in current directory
                    # Recompile, if name.so is older than name.c
                    if os.path.getctime(cname) > os.path.getctime(bfname):
                        print gcc_compile
                        os.system(gcc_compile)
                        print gcc_link
                        os.system(gcc_link)
            # Do nothing if name.so is not in current directory
            self.plfname = name
            return
                   
        # The file is neither name.c nor name.so
        print "Error: file is neither name.c nor name.so"
        return


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
            self.plfname,
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
    ny = grid[1]
    nx = grid[0]

    # HERE MUST BE CHECK FOR CORRECT DIMENSIONS of pos, dir etc.

    #  
    # Find the basis, (tau,xi), of the inner coordinate
    # system in the image plane
    #
    norm2impl = obs/npl.norm(obs)               # n is unity normal to image plane
    tau = np.cross((0.,0.,1.), norm2impl) # tau = ez x n
    tau = tau/npl.norm(tau) # tau is a unity-length ort, perp to SGI Z and n
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
            dir[i,j,:] = slope/npl.norm(slope) # Directions of Obs-to-pixel, SGI

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


#
# instead of vmagn(vec), use npl.norm() (ie numpy.linalg.norm()) 
#
## def vmagn(vect):
##     """
##     Magnitude of a vector as sqrt(dot(vect,vect)).
##     Is written exclusively for brevity :)
##     """
##     return np.sqrt(np.dot(vect, vect))


def make_streamer(theta,phi,
                  orientation=0.0,
                  density=2.0,
                  baseStrength=5.0,
                  stalkStrength=2.0,
                  scaleX = .75,
                  scaleY = .75,
                  scaleZ = .75,
                  plfname='plasma_parameters'):
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
    rtcore.make_streamer(1.0,degtorad*theta,degtorad*phi,
                         degtorad*orientation,density,baseStrength,
                         stalkStrength,scaleX,scaleY,scaleZ,plfname)

def remove_streamers(plfname='plasma_parameters'):
    """
    This function removes all streamers from the plasma parameters provided.

    parameter: plfname-The name of the plasma parameters file
    """
    rtcore.remove_streamers(plfname)


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


    dist = np.sqrt(np.dot(obs,obs))

    const = rtconstants()
    
    #             rect=[xleft, ybottom, xright, ytop]
    #                     0       1       2      3
    xleft =   rect[0]
    ybottom = rect[1]
    xright =  rect[2]
    ytop =    rect[3]

    alpha = (xright - xleft)/grid[0]/dist
    beta = (ytop - ybottom)/grid[1]/dist
    solid_angle = alpha*beta
    conversion = 2e26*const.Boltzmannk_SI/(const.c_light_SI**2)*solid_angle
    
    if(mode=='Tbr'):
        data = np.empty((1,steps,grid[0],grid[1]))
    if(mode=='TbrIQUV'):
        data = np.empty((1,steps,grid[0],grid[1],2))
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

        xleft =   self.rect[0]
        ybottom = self.rect[1]
        xright =  self.rect[2]
        ytop =    self.rect[3]


        prihdr.update('CTYPE1', 'RA---SIN', 'Coordinates along the ecliptic')
        prihdr.update('CRVAL1', 0.0, '') 
        prihdr.update('CDELT1', (xright - xleft)/float(np.size(data,2))/ \
                      dist*180/np.pi, 'Degrees/Pixel')
        prihdr.update('CRPIX1', (-xleft)/(xright - xleft)*grid[0], '')


        prihdr.update('CTYPE2', 'DEC--SIN', 'Coordinates perpendicular ' \
                      'to the ecliptic')
        prihdr.update('CRVAL2', 0.0, '')
        prihdr.update('CDELT2', (ytop - ybottom)/float(np.size(data,3))/ \
                      dist*180/np.pi, 'Degrees/Pixel')
        prihdr.update('CRPIX2', (-ybottom)/(ytop - ybottom)*grid[1], '')


        prihdr.update('CTYPE3', 'FREQ', '')
        prihdr.update('CRVAL3', freq_start, 'Frequency in Hz')
        prihdr.update('CDELT3', (freq_end-freq_start)/(steps-1),
                      'Frequency per index')
        prihdr.update('CRPIX3', 1.0, 'Reference Freq index')



        prihdr.update('CTYPE4', 'STOKES', '')
        prihdr.update('CRVAL4', 1.0+j*3.0, '')
        prihdr.update('CDELT4', 1.0, '')
        prihdr.update('CRPIX4', 1.0, '')
        prihdr.update('CROTA4', 0.0, '')



        prihdr.update('ORIGIN', 'MIT Haystack Observatory', 'Organization')



        if(units == 'K'):
            prihdr.update('KTOJY', conversion ,
                          'This times nu^2 gives the conversion to Jy')
        else:
            prihdr.update('JyTOK', (1/conversion) ,
                          'This divided by nu^2 gives the conversion to K')

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



def vangle(a, b, amagn=None, bmagn=None):
    """
    Calculates angle(s) between two vectors in arrays a and b.
    The shape of both is (n,3), n is arbitrary number of 3d vectors.
    """
    sha = a.shape
    shb = b.shape
    if len(sha) == 1: # a and b are just 3d vectors, not arrays of vectors
        if amagn == None: amagn = np.sqrt(dot(a,a))
        if bmagn == None: bmagn = np.sqrt(dot(b,b))
        return np.arccos(np.dot(a,b)/(amagn*bmagn))
    else:
        n = sha[0]
        theta = np.empty(n, dtype=np.double)
    
    for i in xrange(n):
        if amagn == None: amagn = np.sqrt(np.dot(a[i,:],a[i,:]))
        if bmagn == None: bmagn = np.sqrt(np.dot(b[i,:],b[i,:]))
        theta[i] = np.arccos(np.dot(a[i,:],b[i,:])/(amagn*bmagn))

    return theta


def line(x0=None, y0=None, x1=None, y1=None, n=None):
    """
    Returns array of (x,y) point coordinates lying along a straight
    line on the image plane. The points in the resulting array are
    evenly distributed.
    The line is specified by the couple of its endpoint coordinates,
    (x0,y0) and (x1,y1). The function may be called as:
      line(x0,y0,x1,y1, n)
      line((x0,y0,x1,y1), n)
      line((x0,y0),(x1,y1), n).
    
    """
    # If called as line((x0,y0,x1,y1), n)
    if len(x0) == 4:
        if y0 == None:
            print 'Error: number of points must be given'
            return None
        n = int(y0)
        y0 = np.double(x0[1])
        x1 = np.double(x0[2])
        y1 = np.double(x0[3])
        x0 = np.double(x0[0])
        print x0, y0, x1, y1

    elif len(x0) == 2 and len(y0) == 2:
        if x1 == None:
            print 'Error: number of points must be given'
            return None
        n = int(x1)
        y = np.double(x0[1])
        x0 = np.double(x0[0])
        x1 = np.double(y0[0])
        y1 = np.double(y0[1])
        y0 = y
        print x0, y0, x1, y1

    elif x0 == None or y0 == None or x1 == None or y1 == None:
        print 'Error: must be given four numbers, x0, y0, x1, y1'
        return None
    if n == None:
        print 'Error: number of points must be given'
        return None

    la = np.empty((n,2), dtype=np.double)
    la[:,0] = np.linspace(x0, x1, n)
    la[:,1] = np.linspace(y0, y1, n)

    return la
