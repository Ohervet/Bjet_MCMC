# Executable name will be bj_mcmc.
# Modified from Makefile for bjet02 in Feb 2022 for use in MCMC by Sarah Youngquist
# Change: changed file names to match mcmc, added support for MacOS Darwin,
# optimization is more aggressive, using O3 instead of O2.
SHELL=/bin/sh

OperatingSystem:=$(shell uname -s)

ifeq (${OperatingSystem},Linux)
	ARCH=linux
endif

ifeq (${OperatingSystem},OSF1)
	ARCH=alphacxx6
endif

# ADDED SY
ifeq (${OperatingSystem},Darwin)
	ARCH=darwin
endif

Dir=./
oDir=./
iDir=-Wall
lDir=-lm

CXX		=

PROG		= bj_mcmc
NO		= 02

ObjSuf		=o
SrcSuf		=cpp
ExeSuf		=exe
DllSuf		=so
OutPutOpt	=-o

SRC		= processes_supp_mcmc.$(SrcSuf) $(PROG).$(SrcSuf)

.SUFFIXES:	.$(SrcSuf)

OBJ		= $(SRC:.$(SrcSuf)=.$(ObjSuf))

ifeq ($(ARCH),alphacxx6)
#Alpha/OSF with cxx6
$(Attention, je fais n importe quoi...)
CXX	=cxx
CXXFLAGS=-O0 -gall

LD	=cxx
LDFLAGS	=-O
SOFLAGS	=-Wl,-expect_unresolved,*-shared
endif

ifeq ($(ARCH),linux)
CXX	=g++
CXXFLAGS=-O3 -Wall -g -fPIC 
LD	=g++
LDFLAGS	=-O3 -Wall -lm -g -fPIC
SOFLAGS	=-shared
endif

# ADDED SY
ifeq ($(ARCH),darwin)
CXX	=g++
CXXFLAGS=-O3 -Wall -g -fPIC 
LD	=g++
LDFLAGS	=-O3 -Wall -lm -g -fPIC
SOFLAGS	=-shared
endif

ifeq ($(CXX),)
$(error $(ARCH) invalid architecture)
endif

LIBS	=$(SYSLIBS)
GLIBS	=$(SYSLIBS)


#---------------------------------------------
.$(SrcSuf).$(ObjSuf):
	$(CXX) $(CXXFLAGS) -c $<

###
$(PROG):	$(OBJ)
		$(LD) $(OutPutOpt) $@ $^ $(LDFLAGS) $(LIBS)
		@echo ""
		@echo "$@ done"
		@echo ""

clean:
		@rm -f $(OBJ) $(PROG)
		@rm -f *.o $(PROG2)
		@echo ""
		@echo "cleaned"
		@echo ""

cleandata:
		@rm -f ./data/*.dat
		@echo ""
		@echo "data removed"
		@echo ""


#----------host gal-----------#

CXX=g++
Dir=./
oDir=./
iDir=-Wall
lDir=-lm

PROG1=gal

N1=eg02

gal:
    EXOBJS1=\
	$(oDir)$(N1).o \
        
    $(PROG1):	$(EXOBJS1)
	$(CXX) -o $@ $(EXOBJS1) $(lDir)

    $(oDir)/$(N1).o: $(N1).cpp
	$(CXX) -c $< $(iDir) -o $@


