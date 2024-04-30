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
from distutils.sysconfig import get_python_lib
import os

# Check if numpy and other modules are installed
try:
    import numpy 
except ImportError:
    print '^^^'
    print '^^^ The numpy module is not installed'
    print '^^^'
    sys.exit(1)

AOPATH = get_python_lib(1) + '/numpy/core/include/numpy/arrayobject.h' 

#
# Compile ray tracing programs and create library
# lib/libmxvec.a if it does not exist
#
bdir = os.getcwd()
os.chdir('rtcore')
os.system('make')
os.chdir(bdir)

#
# Remove old versions of headers and libs, if there are any.
# 
if os.path.isdir(bdir+'/raytrace/lib'):
    print 'rm -fr '+bdir+'/raytrace/lib'
    os.system('rm -fr '+bdir+'/raytrace/lib')
if os.path.isdir(bdir+'/raytrace/inc'):
    print 'rm -fr '+bdir+'/raytrace/inc'
    os.system('rm -fr '+bdir+'/raytrace/inc')  
#
# Create inc and lib folders for headers and libs
# and copy the freshest versions to lib and inc
print 'mkdir raytrace/lib'
os.system('mkdir raytrace/lib')
print 'cp rtcore/lib/libmxv.a raytrace/lib/'
os.system('cp rtcore/lib/libmxv.a raytrace/lib/')
print 'mkdir raytrace/inc'
os.system('mkdir raytrace/inc')
print 'cp rtcore/raytrace.h raytrace/inc/'
os.system('cp rtcore/raytrace.h raytrace/inc/')
print 'cp streamer.h raytrace/inc/'
os.system('cp streamer.h raytrace/inc/')
 

module_rtcore = Extension('raytrace/rtcore', \
    library_dirs = ['rtcore/lib'], \
    libraries = ['m', 'mxv'], \
    sources = ['rtcore/rtcoremodule.c',  \
               'rtcore/advance_beam.c',
               'rtcore/calc_tbr.c',
               'rtcore/calc_tbriquv.c'], \
    extra_compile_args = ['-lpthread'],)

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
      package_data = {'raytrace' : ['lib/libmxv.a', 'inc/raytrace.h',  \
                                    'inc/streamer.h']},
      ext_modules = [module_rtcore]
)


# Back to the base directort
os.chdir(bdir)

