include ../Makefile.in

F90FILES = ibm.f90 ibm_ghost_point.f90 ibm_solver.f90 ibm_integrator.f90 ibm_io.f90

LIBFILE = libibm.a

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
