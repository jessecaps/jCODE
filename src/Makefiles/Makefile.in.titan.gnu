# Directories
HOMEDIR = $(shell pwd | sed -e 's/\/src.*//')
LIBDIR  = $(HOMEDIR)/lib
MODDIR  = $(HOMEDIR)/mod
OBJDIR  = $(HOMEDIR)/obj
BINDIR  = $(HOMEDIR)/bin
VPATH   = $(LIBDIR) $(BINDIR) $(OBJDIR)

# Compiler and archiver
CC  = cc
CXX = CC
F90 = ftn
F77 = ftn
LD  = ftn
AR  = ar rcv
RL  = ranlib

# Preprocessing options
OPTIONS = -DUSE_GPU

# Compiler flags
CFLAGS   = 
F90FLAGS = -std=gnu -ffree-line-length-0 -cpp $(OPTIONS)
F77FLAGS = -std=legacy
LDFLAGS  = #-lSystem -lSystemStubs
INCFLAGS = -I$(MODDIR)
MODFLAGS = -J $(MODDIR)
DBGFLAGS = -O0 -fbacktrace -Wall -Wextra -Waliasing -Wno-unused-dummy-argument -Wno-unused-parameter -ffree-line-length-none -fall-intrinsics -g
OPTFLAGS = -O3 -ftree-vectorize -ffast-math -funroll-loops -fomit-frame-pointer -pipe

# External libraries
#CUDA_DIR = /Users/jcaps/Research/Codes/Builds/cuda/
#CUDA_INC = -I$(CUDA_DIR)/include
#CUDA_LIB = -L$(CUDA_DIR)/lib -lcuda -lcudart 

# Installation script
INSTDIR = $(HOME)/bin
INSTSCPT = cp $(BINDIR)/* $(INSTDIR)/.
