include ../Makefile.in

F90FILES = boundary.f90 sponge.f90 farfield.f90 impenetrable.f90 isothermal.f90 adiabatic.f90 \
	slip.f90 outflow.f90 inflow.f90 boundary_strong.f90 excitation.f90
LIBFILE = libboundary.a

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
