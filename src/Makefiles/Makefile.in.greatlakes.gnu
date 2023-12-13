# Directories
%.o: %.mod
HOMEDIR = $(shell pwd | sed -e 's/\/src.*//')
LIBDIR  = $(HOMEDIR)/lib
MODDIR  = $(HOMEDIR)/mod
OBJDIR  = $(HOMEDIR)/obj
BINDIR  = $(HOMEDIR)/bin
VPATH   = $(LIBDIR) $(BINDIR) $(OBJDIR)

# Compiler and archiver
CC  = mpicc
CXX = mpicxx
F90 = mpif90
F77 = mpif90
LD  = mpif90
AR  = ar rcv
RL  = ranlib

#  Preprocessing options
OPTIONS =

# Compiler flags
CFLAGS   =
F90FLAGS = -std=gnu -fallow-argument-mismatch -fallow-invalid-boz -ffree-line-length-0 -cpp $(OPTIONS)
F77FLAGS = -std=legacy
LDFLAGS  = 
INCFLAGS = -I$(MODDIR)
MODFLAGS = -J $(MODDIR)
DBGFLAGS = -O0 -fbacktrace -Wall -Wextra -Waliasing -Wno-unused-dummy-argument -Wno-unused-parameter -ffree-line-length-none -fall-intrinsics -fcheck=all -g 
OPTFLAGS = -O3 -ftree-vectorize -ffast-math -funroll-loops -fomit-frame-pointer -pipe

# Installation script
INSTDIR = $(HOME)/bin
INSTSCPT = cp $(BINDIR)/* $(INSTDIR)/.
