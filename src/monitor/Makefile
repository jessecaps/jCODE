include ../Makefile.in

F90FILES = monitor.f90 timing.f90 monitor_state.f90 monitor_particles.f90                    \
	monitor_cases.f90 monitor_shear.f90 monitor_rayleigh_taylor.f90                      \
	monitor_impulse.f90 monitor_time_info.f90 monitor_adjoint.f90                        \
	monitor_homogeneous.f90 monitor_shock_tube.f90 monitor_drag.f90                      \
	monitor_ibm.f90
LIBFILE = libmonitor.a

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
