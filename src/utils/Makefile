include ../Makefile.in

BINFILES = init_flow data2ensight analyzer modify_particles spectrum_1D stl2ensight \
	empty_particle test_adjoint stl2levelset ibm2part modify_ibm

all:	Makefile
	@for i in $(BINFILES); do make -C $$i LIBRARIES="$(LIBRARIES)"; done

clean: 
	@for i in $(BINFILES); do make clean -C $$i; done
