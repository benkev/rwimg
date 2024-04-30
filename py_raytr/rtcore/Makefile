#
# Top level Makefile for rtcore, the C-language raytracing engine.
# The rtcore Python module is compiled into the shared library 
# rtcore.so, which is used by the pure Python package raytrace.
#

BASEDIR = $(PWD)
SUBDIRS = mxvec

.PHONY: install, clean, purge

install:
	@for subdir in $(SUBDIRS); do \
	  echo "Installing all in $$subdir"; \
	  cd $$subdir && $(MAKE) install && cd $(BASEDIR); \
	done

clean:
	@echo "Cleaning object files in $(BASEDIR)";
	rm -f *.o
	@for subdir in $(SUBDIRS); do \
	  echo "Cleaning object files in $(BASEDIR)/$$subdir"; \
	  cd $$subdir && $(MAKE) clean && cd $(BASEDIR); \
	done

purge:
	@echo 'Cleaning remporary emacs *~ and vi #*# files'; \
	rm -f `find ./ -name '*~'` && rm -f `find ./ -name '#*#'`; \
	echo 'Cleaning libraries in $(BASEDIR)/lib';		\
	cd $(BASEDIR)/lib && rm -f *.a  && cd $(BASEDIR);
	@echo "Cleaning object files in $(BASEDIR)";
	rm -f *.o
	@for subdir in $(SUBDIRS); do \
	  echo "Cleaning object files in $(BASEDIR)/$$subdir"; \
	  cd $$subdir && $(MAKE) clean && cd $(BASEDIR); \
	done


