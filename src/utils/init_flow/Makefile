include ../../Makefile.in

F90FILES = init_flow.f90 quiescent.f90 particle_box.f90 shear_layer.f90 boundary_layer.f90  \
	premixed.f90 rayleigh_taylor.f90 cylinder.f90 jet_crossflow.f90 impulsive_plate.f90 \
	convect_scalar.f90 fluidized_bed.f90 mixing_layer.f90 taylor_green.f90 channel.f90  \
	shock_tube.f90 lib.f90 sphere.f90 flame.f90 pipe.f90 jet.f90 ibm_init.f90 ibm_part.f90 \
	nozzle.f90 erosion.f90 shu_osher.f90 hit.f90 jet_impingement.f90 el_ibm.f90 el_ibm_scalar.f90

BINFILE = init_flow

OFILES   = $(F90FILES:.f90=.o)
MODFILES = $(F90FILES:.f90=.mod)

.SUFFIXES: .f90 .o

#$(BINFILE): $(OFILES)
#	cd $(OBJDIR);$(LD) $(FLAGS) $(LIBRARIES:%=$(LIBDIR)/lib%.a) -o $(BINDIR)/$(BINFILE) $+ $(LDFLAGS)

$(BINFILE): $(OFILES) $(LIBDIR)/lib*.a $(LIBDIR)/lib*.a
	cd $(OBJDIR);$(LD) $(FLAGS) -o $(BINDIR)/$(BINFILE) $+ $(LDFLAGS) $(FFTW_LIB)

clean: 
	cd $(OBJDIR); rm -f $(OFILES)
	cd $(MODDIR); rm -f $(MODFILES)
	cd $(BINDIR); rm -f $(BINFILE)

.f90.o: 
	$(F90) $(F90FLAGS) $(FLAGS) $(INCFLAGS) $(FFTW_INC) -c $*.f90 -o $(OBJDIR)/$*.o $(MODFLAGS)
