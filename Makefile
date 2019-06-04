# Makefile

COMPILER	:=  gfortran
VPATH		:=  Poisson1D/Source
BINDIR		:=  Poisson1D/Bin
OBJECTDIR	:=  Poisson1D/Objects
LIB		:=  Poisson1D/Libraries
LIBOBJECTS	:=  Poisson1D/Libraries/Objects
FFLAGS		:=  -g -O3 -J$(OBJECTDIR)
FFLAGSDebug 	:=  -g -Wall -fcheck=all -J$(OBJECTDIR)

OBJECTS := $(BINDIR)/element1DMod.o $(BINDIR)/GLFEM.o $(BINDIR)/domainMod.o $(BINDIR)/quicksort.o $(BINDIR)/ST_TO_CC.o $(BINDIR)/sparse.o $(BINDIR)/GMRES.o $(BINDIR)/Poisson1D.o $(BINDIR)/GidData.o $(BINDIR)/functions.o $(BINDIR)/main.o
LIBRARIES := $(LIB)/mylib.a

main: $(OBJECTS)
	$(COMPILER) $(FFLAGS) $(OBJECTS) -I$(LIBOBJECTS) -L$(LIB) $(LIBRARIES) -o main

mainDebug: $(OBJECTS)
	$(COMPILER) $(FFLAGSDebug) $(OBJECTS) -o main

$(BINDIR)/%.o : $(VPATH)/%.f90
	$(COMPILER) $(FFLAGS) -c $^ -I$(LIBOBJECTS) -o $@ 
clean:
	rm -f $(BINDIR)/*.o $(OBJECTDIR)/*.mod main
