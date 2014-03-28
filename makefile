# COMPILER
FC = gfortran

# COMPILER OPTIONS
FFLAGS = -g -w -std=legacy -mcmodel=large -mtune=native -march=native -masm=intel -freal-4-real-10 -freal-8-real-10 -ffixed-line-length-none -O3 -floop-parallelize-all -ftree-parallelize-loops=4

# SOURCES
SOURCES = $(wildcard *.f) 

# OBJECTS
OBJECTS = $(SOURCES:.f=.o)

PMFCMC: $(OBJECTS)
	$(FC) $(FFLAGS) $(OBJECTS) -o $@

cleanObjects:	
	rm -f $(OBJECTS)

cleanOutput:
	rm -f $(PWD)/res/*.dat

cleanFigures:
	rm -f $(PWD)/res/*.eps
	rm -f $(PWD)/res/*.fig

# CHECK GFORTRAN VERSION

GFV       := $(shell expr `gfortran -dumpversion`)
GFVMAJL4  := $(shell expr `gfortran -dumpversion | cut -f1 -d.` \<  4)
GFVMAJE4  := $(shell expr `gfortran -dumpversion | cut -f1 -d.` \=  4)
GFVMINGE7 := $(shell expr `gfortran -dumpversion | cut -f2 -d.` \>= 7)
checkGfortranVersion:
	$(info Checking if gfortran version supports 80-bit precision...)
# Is version 4.x.y installed? Issue error is older version is installed.
ifeq "$(GFVMAJL4)" "1" 
	$(info Installed gfortran version: $(GFV))
	$(error At least version 4.7 is required for 80-bit precision)
endif
# If version 4.x.y installed, is x >= 7?
ifeq "$(GFVMAJE4)" "1" 
ifeq "$(GFVMINGE7)" "0" 
	$(info Installed gfortran version: $(GFV))
	$(error At least version 4.7 is required for 80-bit precision)
endif
endif
# Do nothing if major release is greater than 4.7.y
	$(info Installed gfortran version: $(GFV). Version is OK!)

run:	
	clear
	make checkGfortranVersion
#	make cleanOutput
#	make cleanFigures
	make cleanObjects
	make
	./PMFCMC

open:
	emacs makefile&
	emacs COM&
	emacs MAIN.f&
	emacs INIT.f&
	emacs GRID.f&
	emacs MISC.f&
	emacs PMFSML.f&
	emacs PMFDSML.f&
	emacs BETA.f&
	emacs CG.f&
	emacs UNFRM.f&
	emacs OTHERMODELS.f&
	emacs MATHUTIL.f&
	emacs ROOTFINDER.f&
	emacs INTEGRATOR.f&
	emacs DIFFERENTIATOR.f&
	emacs NONLINEQSOLVER.f&
