#
# Makefile to build 
#

# Macros
CXX=nvcc
CXXFLAGS=  -arch=sm_20 -g -Xcompiler -fPIC -Xcompiler -fno-exceptions \
                          -Xcompiler -fno-rtti 
OBJECTS=advance_beam.o simulation.o streamer.o 
#       calc.o matrix_math.o plasma_parameters.o 
LIBS = 
#  -lpgplot -lX11 -lftn -lm
BASEDIR = $(PWD)

# Targets
all:    lib/libsimulation_core.a lib/libmxvec.a

lib/libsimulation_core.a:   $(OBJECTS) 
	   $(CXX) -lib $(LIBS) $(OBJECTS) -o lib/libsimulation_core.a 

lib/libmxvec.a:
	cd $(BASEDIR)/mxvec && $(MAKE) install && cd $(BASEDIR)

advance_beam.o: advance_beam.cu
	$(CXX) $(CXXFLAGS) -c advance_beam.cu -o advance_beam.o $(LIBS) 

simulation.o: simulation.cu
	$(CXX) $(CXXFLAGS) -c simulation.cu -o simulation.o $(LIBS) 

#calc.o: calc.cu
#	$(CXX) $(CXXFLAGS) -c calc.cu -o calc.o $(LIBS) 

#plasma_parameters.o: plasma_parameters.cu
#	$(CXX) $(CXXFLAGS) -c plasma_parameters.cu -o plasma_parameters.o \
				$(LIBS) 

streamer.o: streamer.cu
	$(CXX) $(CXXFLAGS) -c streamer.cu -o streamer.o $(LIBS) 

#matrix_math.o: matrix_math.cu
#	$(CXX) $(CXXFLAGS) -c matrix_math.cu -o matrix_math.o $(LIBS) 

clean: 
	rm -f *.o

purge:  clean
	rm -f lib/lib*.a