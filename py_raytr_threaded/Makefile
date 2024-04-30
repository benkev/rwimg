#
# Makefile for the Python and C package raytrace.
# Current version has only one function: cleanup of 
# the whole directory tree.
#

BASEDIR = $(PWD)
SUBDIRS = rtcore

.PHONY: clean, purge

install:
	rm -fr build
	python setup.py build_ext --inplace
	python setup.py install --home=~

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
	cd $(BASEDIR)/rtcore/lib && rm -f *.a  && cd $(BASEDIR);
	@echo "Cleaning object files in $(BASEDIR)";
	rm -f *.o
	@for subdir in $(SUBDIRS); do \
	  echo "Cleaning object files in $(BASEDIR)/$$subdir"; \
	  cd $$subdir && $(MAKE) purge && cd $(BASEDIR); \
	done


