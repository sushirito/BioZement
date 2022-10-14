
VERSION = -std=c++11

PROGRAM = lb_ejah
OS := $(shell uname)
HOST := $(shell hostname)


### check if we are on Abel or not
#ifneq (,$(findstring login-,$(HOST)))
### Abel
#DESTDIR = ../badchimp_run/program_files/bin
#else
### not Abel
DESTDIR = run/program_files/bin
#endif

### turn on MPI by default
MPI = 1

### C compilers
#CC = icc         # Intel C compiler, use on cluster
#CC = gcc        # GNU C compiler, use on PC
#CC = g++         # GNU C++ compiler
#CC = icpc       # Intel C++ compiler, use on cluster
#MPICC = /usr/local/bin/mpicc    # build parallel c code on abel cluster
#MPICC = mpicc    # build parallel c code on linux
#MPICC = /usr/local/bin/mpic++  # build parallel c++ code
#OMPI_CXX=gcc-7
MPICC = mpic++  # build parallel c++ code

### clang (Apple) compiler flags
#CFLAGS = -O3 -Wall
#CFLAGS = -g -O0 -Wall
#CFLAGS = -O3 -pipe -fomit-frame-pointer -Wall
#CFLAGS = -O3 -pipe -fomit-frame-pointer -Wall -fslp-vectorize-aggressive

### GNU gcc compiler flags
CFLAGS = -O1 -pipe -fomit-frame-pointer -Wall
CFLAGS += -g 
LDFLAGS += -lm

## Intel compiler flags
#CFLAGS = -O3 -xAVX -mavx -fomit-frame-pointer -fno-alias -Wall
#CFLAGS = -O3 -Wall 
#LDFLAGS = -openmp

ifeq ($(OS), Darwin)
LIBS = -I/usr/include/malloc
endif

ifeq ($(MPI),1)
CC = $(MPICC)
CFLAGS += -D_MPI_
endif

CFLAGS += $(VERSION)
LDFLAGS += $(VERSION)


### Turn on _2D_RUN_ if 2D given in inp.dat
### ---------------------------------------
### Read inp.dat and strip commented lines (s/\#.*$$//), leading whitespace (s/^[[:space:]]*//), 
### and empty lines (/^[[:space:]]*$$/d) (Makefile specifics: protect $ with extra $, and # with \).
### First line of output is dimension of simulation (sed -n '1p')
#DIM := $(shell sed -e 's/\#.*$$//;s/^[[:space:]]*//;/^[[:space:]]*$$/d' inp.dat | sed -n '1p' )
ifeq ($(MAKECMDGOALS),)
MAKECMDGOALS := flow
endif
ifneq ($(MAKECMDGOALS),clean)
DIM := $(shell grep 'nD' run/$(MAKECMDGOALS)/input.dat | cut -d'\#' -f1 | cut -d'D' -f2 | sed -e 's/^[[:space:]]*//;s/[[:space:]]*$$//')
ifeq ($(DIM),2)
CFLAGS += -D_2D_RUN_
endif
endif

### the case targets require specific pre-processor macros  
grad_P: CFLAGS += -D_FLUID_BDRY_PRESS_
perc: CFLAGS += -D_FIND_PERC_NODES_
source: CFLAGS += -D_FLUID_SOURCE_ -D_SPECIES_SOURCE_
calc_magn_replace: CFLAGS += -D_SPECIES_SOURCE_BULK_ -D_RANDOM_MINERAL_ -D_MAGN_ON_MAGN_ 
two-phase: CFLAGS += -D_USE_COLOR_GRAD_
biozement: #CFLAGS += -DDIFFUSIVE_RUN 

flow grad_P perc source calc_magn_replace two-phase biozement: $(PROGRAM)

SRCS = $(wildcard src/*.cpp)
OBJS = $(SRCS:.cpp=.o)
DEPS = $(SRCS:.cpp=.P)

$(PROGRAM):	$(OBJS)
		$(CC) -o $@ $^ $(LDLIBS) $(LDFLAGS)

%.o : %.cpp
	$(CC) -MMD -c -o $@ $< $(LIBS) $(CFLAGS) 
	@cp $*.d $*.P; \
	  sed -e 's/\#.*//' -e 's/^[^:]*: *//' -e 's/ *\\$$//' \
	      -e '/^$$/ d' -e 's/$$/ :/' < $*.d >> $*.P; \
	  rm -f $*.d

### print variable values by command: make print-VARNAME
print-%  : ; @echo $* = $($*)

### copy binary to case-specific folder
flow grad_P perc source calc_magn_replace two-phase biozement:
	mkdir -p $(DESTDIR)/$@/ && cp $(PROGRAM) $(DESTDIR)/$@/

.PHONY: clean
clean:
	@-rm -f src/*.o src/*.P src/*.d

### avoid including DEPS during clean (why create them when they are just deleted?)
ifneq ($(MAKECMDGOALS),clean)
-include $(DEPS) 
endif

# end of Makefile 

