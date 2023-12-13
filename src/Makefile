include Makefile.in

# Files
CODE = jcode
MODULES = library config operator grid state boundary sources combustion particle ibm controller functional postprocess core gradient monitor
LIBRARIES = core postprocess gradient functional controller ibm particle combustion boundary state sources grid operator config monitor library

# Targets

default:Makefile Makefile.in
	@make libraries
	@make $(CODE)
	@make util

all:    Makefile Makefile.in 
	@make libraries
	@make $(CODE)
	@make util

opt:
	@make "FLAGS = $(OPTFLAGS)"

debug:
	@make "FLAGS = $(DBGFLAGS)"

allopt:
	@make all "FLAGS = $(OPTFLAGS)"

alldebug:
	@make all "FLAGS = $(DBGFLAGS)"

$(CODE):$(MODULES:%=lib%.a)
	 $(LD) $(FLAGS) $(OBJDIR)/main.o $(LIBRARIES:%=$(LIBDIR)/lib%.a) -o $(BINDIR)/$(CODE) $(OPTIONAL_LIBS) $(LDFLAGS)

libraries:
	@for i in $(MODULES); do make -C $$i; done

clean:
	@for i in $(MODULES); do make clean -C $$i; done
	rm -f $(BINDIR)/$(CODE)

util:  $(MODULES:%=lib%.a)
	@make -C utils LIBRARIES="$(LIBRARIES)"

utilclean:
	@make -C utils clean

distclean:
	@make clean
	@make utilclean

install:
	$(INSTSCPT)
