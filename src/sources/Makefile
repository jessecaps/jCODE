include ../Makefile.in

F90FILES = sources.f90 gravity_source.f90 acoustic_source.f90 ignition_source.f90            \
	fuel_source.f90 momentum_source.f90 pressure_gradient.f90 solenoidal_excitation.f90  \
	linear_forcing.f90
LIBFILE = libsources.a

OFILES = $(F90FILES:.f90=.o)
MODFILES = $(F90FILES:.f90=.mod)

.SUFFIXES: .o .f90

$(LIBFILE):$(OFILES) Makefile	
	cd $(OBJDIR); $(AR) $(LIBDIR)/$(LIBFILE) $(OFILES); $(RL) $(LIBDIR)/$(LIBFILE)

.f90.o: 
	$(F90) $(F90FLAGS) $(FLAGS) $(INCFLAGS) -c $*.f90 -o $(OBJDIR)/$*.o $(MODFLAGS)

clean: 
	cd $(OBJDIR); rm -f $(OFILES)
	cd $(MODDIR); rm -f $(MODFILES)
	cd $(LIBDIR); rm -f $(LIBFILE)
