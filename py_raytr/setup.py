###############################################################
#                                                             #
# setup.py                                                    #
#                                                             #
# Setup program for the module raytrace.                      #
#                                                             #
# Created 15 March 2011 by Leonid Benkevitch                  #
#                                                             #
###############################################################

import distutils
from distutils.core import setup, Extension
import os

# Check if numpy and other modules are installed
try:
    import numpy 
except ImportError:
    print '^^^'
    print '^^^ The numpy module is not installed'
    print '^^^'
    sys.exit(1)

#
# Create variants of advance_beam code:
#   advance_beam_basic.c - calculate ray positions and directions only;
#   advance_beam_Tbr.c - calculate brightness temperature Tbr
#   advance_beam_TbrIV.c - calculate brightness temperature as Stokes I and V
# They differ by only the first line with #define TBR or TBRIV statement.
#
bdir = os.getcwd()
os.chdir('rtcore')
print 'cp advance_beam.c advance_beam_basic.c'
os.system('cp advance_beam.c advance_beam_basic.c')
print 'sed 1i"#define TBR" advance_beam.c > advance_beam_Tbr.c'
os.system('sed 1i"#define TBR" advance_beam.c > advance_beam_Tbr.c')
print 'sed 1i"#define TBRIV" advance_beam.c > advance_beam_TbrIV.c'
os.system('sed 1i"#define TBRIV" advance_beam.c > advance_beam_TbrIV.c')

#
# Compile ray tracing programs and create library
# lib/libmxvec.a if it does not exist
#
#bdir = os.getcwd()
#if not os.path.isfile(bdir+'/core/libmxv.a'):
#os.chdir('rtcore')
os.system('make')
os.chdir(bdir)

# Remove old versions of headers and libs, if there are any,
# and copy the freshest versions to ilib and inc
if os.path.isdir(bdir+'/raytrace/lib'):
    print 'rm -fr '+bdir+'/raytrace/lib'
    os.system('rm -fr '+bdir+'/raytrace/lib')
print 'mkdir raytrace/lib'
os.system('mkdir raytrace/lib')
print 'cp rtcore/lib/libmxv.a raytrace/lib/'
os.system('cp rtcore/lib/libmxv.a raytrace/lib/')

if os.path.isdir(bdir+'/raytrace/inc'):
    print 'rm -fr '+bdir+'/raytrace/inc'
    os.system('rm -fr '+bdir+'/raytrace/inc')
print 'mkdir raytrace/inc'
os.system('mkdir raytrace/inc')
print 'cp rtcore/raytrace.h raytrace/inc/'
os.system('cp rtcore/raytrace.h raytrace/inc/')
 

module_rtcore = Extension('rtcore', \
    library_dirs = ['rtcore/lib'], \
    libraries = ['m', 'mxv'], \
    sources = ['rtcore/rtcoremodule.c',       \
               'rtcore/advance_beam_basic.c', \
               'rtcore/advance_beam_Tbr.c',   \
               'rtcore/advance_beam_TbrIV.c'])

setup(name='raytrace', # Distribution name; will be like "raytrace-0.00.tar.gz
      version='0.00',
      author='Leonid Benkevitch',
      author_email='benkev@haystack.mit.edu',
      description = 'This is a ray tracing package. Can be used for' \
      'obtaining simulated pixel images of the sun.',
      url='http://haystack.mit.edu',
      
      packages = ['raytrace'],
      #package_dir = {'' : 'raytrace'},
      #
      # package_data = {'package name' : ['list of relative path names',..]}
      # is mapping 'package name' to a 'list of relative path names' that
      # should be copied into the package (here 'raytrace')
      # The paths are interpreted as relative to the directory
      # containing the package:
      package_data = {'raytrace' : ['lib/libmxv.a', 'inc/raytrace.h']},
      ext_modules = [module_rtcore]
)

# Back to the base directort
os.chdir(bdir)

# Remove libs and headers. If not to remove, these directories
# may perplex -- they may contain obsolete versions
print 'rm -fr '+bdir+'/raytrace/lib'
os.system('rm -fr '+bdir+'/raytrace/lib')
print 'rm -fr '+bdir+'/raytrace/inc'
os.system('rm -fr '+bdir+'/raytrace/inc')

# Remove variants of advance_beam.c code
print 'rm -f '+bdir+'/rtcore/advance_beam_basic.c'
os.system('rm -f '+bdir+'/rtcore/advance_beam_basic.c')
print 'rm -f '+bdir+'/rtcore/advance_beam_Tbr.c'
os.system('rm -f '+bdir+'/rtcore/advance_beam_Tbr.c')
print 'rm -f '+bdir+'/rtcore/advance_beam_TbrIV.c'
os.system('rm -f '+bdir+'/rtcore/advance_beam_TbrIV.c')

