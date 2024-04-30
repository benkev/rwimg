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
import os
import glob
from distutils.sysconfig import get_python_lib
import rtcore

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
        self.c_light_cms = 2.9979245800e10    # cm/s, speed of light in vacuum
        self.Rsun_cm = 6.955e10               # cm, Solar radius
        self.Rsun_km = 6.955e5                # km, Solar radius
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
                 mode='Tbr', cnu=3.0, msun=1.0, traj=False):
        """
        Create an image plane as a two-dimensional array of doubles.
        The pixel sizes are RA x Dec = grid[0] x grid[1].
        The image dimensions are RA x Dec = dim[0] x dim[1].
        """
        rect = map(np.double, rect)
        obs = map(np.double, obs)

        self.traj = bool(traj);  # True: store trajectories, False: no
        self.ccdef = 'BASIC'
        
        if mode == 'basic' or mode == 1:
            self.mode = 1     # No Tbr or TbrIV calculation;
            self.ccdef = 'BASIC'
            self.traj = True  # Store ray trajectories
        elif mode == 'Tbr' or mode == 2:
            self.mode = 2     # Tbr calculation;
            self.ccdef = 'TBR'
        elif mode == 'TbrIV' or mode == 3:
            self.mode = 3     # No Tbr or TbrIV calculation;
            self.ccdef = 'TBRIV'
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
        self.package = get_python_lib(1) + '/raytrace'
        

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

        ## # Stokes components
        ## self.I = None
        ## self.Q = None
        ## self.U = None
        ## self.V = None

        # Grid size along RA
        self.grid = grid

        # Grid size along RA
        self.nx = grid[1]
        
        # Grid size along Dec
        self.ny = grid[0]
        
        # Positions of the image plane left bottom and right top
        # corners in SGI system
        # rect=[xleft, ybottom, xright, ytop]
        self.rect = pl.empty(4, dtype=pl.double)
        self.rect[:] = rect

        # Image plane size along RA
        self.sizera = rect[2] - rect[0]
        
        # Image plane size along Dec
        self.sizedec = rect[3] - rect[1]
        
        # Frequency of the radio waves
        self.freq = freq
        self.freq2 = freq**2
        self.omega = 2.*np.pi*freq
        self.omega2 = self.omega**2

        # Wavenumber in a vacuum
        self.k0 = self.omega/self.const.c_light_cms      # rad/cm

        # Coefficient at the Ginzburg's formula for effective collision
        # frequency. Ginzburg recommends Cnu = 5.5
        self.cnu = cnu
        
        # Total number of the rays
        self.nrays = grid[1]*grid[0]
        
        # Coordinates of the rays
        self.pos = pl.empty((grid[1], grid[0], 3), dtype=pl.double)
        
        # Direction vectors of the rays
        self.dir = pl.empty((grid[1], grid[0], 3), dtype=pl.double)

        # Initial ray increments (steps)
        self.ds = pl.empty((grid[1], grid[0]), dtype=pl.double)
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
        
        # Maximum number of iterations (steps)
        self.niter = 100
        
        # Brightness temperature
        self.flags = pl.empty((grid[1], grid[0]), dtype=pl.short)
        self.flags[:] = 0   # All the rays are ACTIVE at start

        # Number of active rays
        self.nactive = self.nrays

        # Index of the ray closest to the sun
        self.irayminsd = -1

        # Minimum solar distance
        self.minsoldist = 1e16  # Just a big number

        # The brightness temperature, Tb, calculated along the rays
        if self.mode == 2:
            self.tbr = pl.empty((grid[1], grid[0]), dtype=pl.double)
            self.tbr[:] = 0.0
        else: self.tbr = pl.empty(0)

        #  TbrIV: Tb calculated for two Stokes parameters, I and V
        if self.mode == 3:
            self.tbriv = pl.empty((grid[1], grid[0], 2), dtype=pl.double)
            self.tbriv[:] = 0.0
            self.tpriv = pl.empty((grid[1], grid[0], 2), dtype=pl.double)
            self.tpriv[:] = 0.0
        else:
            self.tbriv = pl.empty(0)
            self.tpriv = pl.empty(0)

        #  Magnetic field that causes the polarization, Bfield
        if self.mode == 3:
            self.bfield = pl.empty((grid[1], grid[0], 3), dtype=pl.double)
            self.bfield[:] = 0.0
        else: self.bfield = pl.empty(0)

        # Optical depths
        self.opdepth = pl.empty((grid[1], grid[0]), dtype=pl.double)
        self.opdepth[:] = 0.0
        
        # Plasma densities
        self.rho = pl.empty((grid[1], grid[0]), dtype=pl.double)
        self.rho[:] = 0.0
        
        # Gradients of plasma densities
        self.gradrho = pl.empty((grid[1], grid[0], 3), dtype=pl.double)
        self.gradrho[:] = 0.0
        
        # Coordinates of the rays at previous step
        self.pospr = pl.empty((grid[1], grid[0], 3), dtype=pl.double)
        self.pospr[:] = 0.0
        
        # Directions of the rays at previous step
        self.dirpr = pl.empty((grid[1], grid[0], 3), dtype=pl.double)
        self.dirpr[:] = 0.0
        
        # New ray increments (steps)
        self.dsnew = pl.empty((grid[1], grid[0]), dtype=pl.double)
        self.dsnew[:] = 0.0
        
        # Estimated ray distances to the critical surface
        self.distcr = pl.empty((grid[1], grid[0]), dtype=pl.double)
        self.distcr[:] = 0.0
        
        # Position of the observer (the earth's coordinates) in SGI system
        self.obs = pl.empty(3, dtype=pl.double)
        self.obs[:] = obs

        # Pixel X and Y increments on the image plane
        self.dx = (rect[2] - rect[0])/self.nx
        self.dy = (rect[3] - rect[1])/self.ny

        # Rulers: arrays of x and y coordinates of all the pixels
        self.xruler = np.linspace(rect[0]+0.5*self.dx,
                                     rect[2]-0.5*self.dx, self.nx)
        self.yruler = np.linspace(rect[1]+0.5*self.dy,
                                     rect[3]-0.5*self.dy, self.ny)

        # Miscellaneous parameters.
        # The parameters are precalculated to lower the computational
        # overhead in the inner loops of the ray tracing algorithm.
        e2w2 = (self.const.ProtonChargeCGSe/self.omega)**2   # (e/w)^2
        e_ovr_mcw2 = e2w2/((self.const.ElectronMass_g*                  \
                           self.const.c_light_cms)**2)   # (e/mcw)^2
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
        self.isec = pl.empty((grid[1], grid[0]), dtype=pl.short)

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
            plf_cname = glob.glob('plasma_*.c') 
            if len(plf_cname) > 0: # A plasma_*.c file exists
                plf_cname = plf_cname[0]
                plfname = os.path.basename(plf_cname)
                self.plfname = plfname.split('.')[0] # Only the name before .c
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
            e2_4pi_ovr_mw2,
            twokf2c2,
            lnLambda_13p7,
            self.cnu,             # Coef. at Ginzburg's nu_eff
            self.callcount        # Number of calls to advance_beam()
            ), dtype=np.double)


    def trace(self, niter=None):
        """
        Trace the beam of rays comprising the nx by ny image plane grid
        for niter steps.
        """
        if niter == None: niter = self.niter

        rtcore.trace_beam(
            self.arprm,
            self.pos,
            self.dir,
            self.ds,
            self.plfname,
            self.rsph,
            niter,
            self.mode,
            self.flags,
            self.tbr,
            self.tbriv,
            self.tpriv,
            self.opdepth,
            self.rho,
            self.gradrho,
            self.bfield,
            self.pospr,
            self.dirpr,
            self.dsnew,
            self.distcr
            )

        

    def soldist(self):
        """
        Calculate the solar distance of all the rays
        """
        dis = np.empty((self.grid[1],self.grid[0]), dtype=np.double)
        for i in xrange(self.grid[1]):
            for j in xrange(self.grid[0]):
                dis[i,j] = vmagn(self.pos[i,j,:])
        return dis

    def set_plfunc(self, fname):
        # Check if the file exist
        if not os.path.isfile(fname):
            print 'Error: no such file "'+fname+'".'
            return
        # File exists at this point.
        bfname = os.path.basename(fname) # Strip the base name of dirs
        name = bfname.split('.')[0] # Remove .c or .so from base file name
        # Prepare compilation and linking strings
        gcc_compile = 'gcc -g -fPIC -D'+self.ccdef+' -I'+self.package+ \
                      '/inc -c '+fname+' -o '+name+'.o'
        gcc_link = 'gcc -shared '+name+'.o -L'+self.package+ \
                   '/lib -L/usr/lib -lm -lmxv -o '+name+'.so'
        # If it is name.c, recompile
        if bfname[-2:] == '.c':
            # Compile .c file into shared library in current directory
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
                
       

    def plprof(self, pos):
        """
        For each SGI position in pos[n,3] calculate the profiles of
        density, rho[n], density gradient, gradrho[n,3], and magnetic
        field, bfield[n,3], using the plasma parameters function
        whose name base (without extension) is in self.plfname string.
        self.plprof() returns the tuple of three arrays:
        rho, gradrho, and bfield.
        Created 10 July 2011 by Mark Benjamin
        """
        pos = np.array(pos)

        if len(pos.shape) < 2:
            print("Error, [n,3] array of coordinates required")
            return
        elif pos.shape[1] != 3: # even though pos has 2 or more dimensions
            print("Error, [n,3] array of coordinates required")
            return
        
        rho = pl.empty((pl.size(pos,0)), dtype=pl.double)
        
        gradrho = pl.empty((pl.size(pos,0), 3), dtype=pl.double)
        
        bfield = pl.empty((pl.size(pos,0), 3), dtype=pl.double)
        
        ds = pl.empty(pl.size(pos,0) , dtype=pl.double)

        flags = pl.empty(pl.size(pos,0) , dtype=pl.short)
        
        rtcore.plprofile(
            self.arprm,
            pos,
            rho,
            gradrho,
            bfield,
            self.mode,
            self.plfname,
            ds,
            flags,            
            )
        return rho, gradrho, bfield
        

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




