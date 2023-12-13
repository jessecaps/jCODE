include ../../Makefile.in

F90FILES = stl2ensight.f90

BINFILE = stl2ensight

OFILES   = $(F90FILES:.f90=.o)
MODFILES = $(F90FILES:.f90=.mod)

.SUFFIXES: .f90 .o

#$(BINFILE): $(OFILES)
#	cd $(OBJDIR);$(LD) $(FLAGS) $(LIBRARIES:%=$(LIBDIR)/lib%.a) -o $(BINDIR)/$(BINFILE) $+ $(LDFLAGS)

$(BINFILE): $(OFILES) $(LIBDIR)/lib*.a $(LIBDIR)/lib*.a
	cd $(OBJDIR);$(LD) $(FLAGS) -o $(BINDIR)/$(BINFILE) $+ $(LDFLAGS)

clean: 
	cd $(OBJDIR); rm -f $(OFILES)
	cd $(MODDIR); rm -f $(MODFILES)
	cd $(BINDIR); rm -f $(BINFILE)

.f90.o: 
	$(F90) $(F90FLAGS) $(FLAGS) $(INCFLAGS) -c $*.f90 -o $(OBJDIR)/$*.o $(MODFLAGS)
