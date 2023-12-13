include ../Makefile.in

F90FILES = functional.f90 acoustic_noise.f90 form_drag.f90 drag_force.f90 reynolds_stress.f90\
	flame_temperature.f90 heat_release.f90 reactant_depletion.f90 binary_mixing.f90      \
	mixing_norm.f90 velocity_norm.f90 data_assimilation.f90 shadowgraph.f90
LIBFILE = libfunctional.a

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
