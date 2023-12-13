include ../../Makefile.in

F90FILES = spectrum_1D.f90

BINFILE = spectrum

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
