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
# Compile ray tracing programs and create library
# lib/libmxvec.a if it does not exist
#
bdir = os.getcwd()
os.chdir('rtcore')
os.system('make')
os.chdir(bdir)

# Remove old versions of headers and libs, if there are any,
# create new ones, and copy the freshest versions to ilib and inc
if os.path.isdir(bdir+'/raytrace_cuda/lib'):
    print 'rm -fr '+bdir+'/raytrace_cuda/lib'
    os.system('rm -fr '+bdir+'/raytrace_cuda/lib')
print 'mkdir raytrace_cuda/lib'
os.system('mkdir raytrace_cuda/lib')
print 'cp rtcore/lib/libmxv.a raytrace_cuda/lib/'
os.system('cp rtcore/lib/libmxv.a raytrace_cuda/lib/')
print 'cp rtcore/lib/libsimulation_core.a raytrace_cuda/lib/'
os.system('cp rtcore/lib/libsimulation_core.a raytrace_cuda/lib/')

if os.path.isdir(bdir+'/raytrace_cuda/inc'):
    print 'rm -fr '+bdir+'/raytrace_cuda/inc'
    os.system('rm -fr '+bdir+'/raytrace_cuda/inc')
print 'mkdir raytrace_cuda/inc'
os.system('mkdir raytrace_cuda/inc')
print 'cp rtcore/raytrace_cuda.h raytrace_cuda/inc/'
os.system('cp rtcore/raytrace_cuda.h raytrace_cuda/inc/')
print 'cp plasma_parameters.h raytrace_cuda/inc/'
os.system('cp plasma_parameters.h raytrace_cuda/inc/')
print 'cp streamer.h raytrace_cuda/inc/'
os.system('cp streamer.h raytrace_cuda/inc/')
print 'cp simulation.h raytrace_cuda/inc/'
os.system('cp simulation.h raytrace_cuda/inc/')
 

module_rtcore = Extension('raytrace_cuda/rtcore', \
    library_dirs = ['rtcore/lib', '/usr/local/cuda/lib64'],\
    libraries = ['m', 'mxv', 'simulation_core', 'cudart', 'cuda' ], \
    sources = ['rtcore/rtcoremodule.c'],  \
    extra_compile_args = ['-lpthread', '-lstdc++'], \
    verbose=True)
    
# Distribution name; will be like "raytrace_cuda-0.00.tar.gz
setup(name='raytrace_cuda',
      version='0.00',
      author='Leonid Benkevitch',
      author_email='benkev@haystack.mit.edu',
      description = 'This is a ray tracing package. Can be used for' \
      'obtaining simulated pixel images of the sun.',
      url='http://haystack.mit.edu',
      verbose=5, 
      packages = ['raytrace_cuda'],
      #package_dir = {'' : 'raytrace_cuda'},
      #
      # package_data = {'package name' : ['list of relative path names',..]}
      # is mapping 'package name' to a 'list of relative path names' that
      # should be copied into the package (here 'raytrace_cuda')
      # The paths are interpreted as relative to the directory
      # containing the package:
      package_data = {'raytrace_cuda' : ['lib/libmxv.a',              \
                                         'lib/libsimulation_core.a',  \
                                         'inc/raytrace_cuda.h',       \
                                         'inc/plasma_parameters.h',       \
                                         'inc/streamer.h',       \
                                         'inc/simulation.h']},
      ext_modules = [module_rtcore]
)

# Back to the base directory
os.chdir(bdir)

