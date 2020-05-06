os = $(shell uname -s)

ifeq ($(os),Linux)
CXX          = g++ 
else
CXX          = c++
endif

# uncomment for debug info in the library
# CXXFLAGS     += -g

INCFLAGS      = -I$(EICDIRECTORY)/include -I$(ROOTSYS)/include -I$(FASTJETDIR)/include
CXXFLAGS      = `root-config --cflags`
# turn on all warnings
CXXFLAGS      += -Wall


ROOTLIBS      = $(shell root-config --libs)
LIBPATH       = $(ROOTLIBS) -L$(FASTJETDIR)/lib -L$(BASEDIR)/tmpsmear/lib
LIBS          = -leicsmear
LIBS         += -lfastjet

# clunky way to see if we can turn on grooming
GROOMING = $(if $(wildcard $(FASTJETDIR)/lib/libRecursiveTools.*),1)
ifeq ($(GROOMING),1)
	CXXFLAGS     += -DGROOMING
	LIBS         += -lfastjettools -lRecursiveTools
endif



# includes
INCS          = 

############################ locations ################################
SDIR          = src
ODIR          = src/obj
BDIR          = bin

################### helper to create subdirs  #########################
dir_guard=@mkdir -p $(@D)


######################### standard rules ##############################
$(ODIR)/%.o : $(SDIR)/%.cxx $(INCS)
	$(dir_guard)
	@echo 
	@echo COMPILING
	$(CXX) $(CXXFLAGS) $(INCFLAGS) -c $< -o $@

$(BDIR)/%  : $(ODIR)/%.o 
	$(dir_guard)
	@echo 
	@echo LINKING
	$(CXX) $(LDFLAGS) $(LIBPATH) $^ $(LIBS) -o $@

###############################################################################
############################# Main Targets ####################################
###############################################################################
all    : $(BDIR)/jetspectra $(BDIR)/extendedjetexample

$(ODIR)/jetspectra.o 	: $(SDIR)/jetspectra.cxx $(INCS)

$(ODIR)/extendedjetexample.o 	: $(SDIR)/extendedjetexample.cxx $(INCS)

# bin
$(BDIR)/jetspectra	: $(ODIR)/jetspectra.o
$(BDIR)/extendedjetexample	: $(ODIR)/extendedjetexample.o

###############################################################################
##################################### MISC ####################################
###############################################################################


clean :
	@echo 
	@echo CLEANING
	rm -vf $(ODIR)/*.o
	rm -rvf $(BDIR)/*dSYM
	rm -rvf lib/*dSYM	
	rm -vf $(BDIR)/*

.PHONY : clean
