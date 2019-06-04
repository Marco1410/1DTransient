# Makefile

COMPILER	:=  gfortran
VPATH		:=  Source
BINDIR		:=  Bin
OBJECTDIR	:=  Objects
FFLAGS		:=  -g -O3 -J$(OBJECTDIR)

OBJECTS := $(BINDIR)/utilities.o $(BINDIR)/defines.o $(BINDIR)/GWprint.o $(BINDIR)/linearSys.o

main: $(OBJECTS)
	ar rcv mylib.a $(OBJECTS)

$(BINDIR)/%.o : $(VPATH)/%.f90
	$(COMPILER) $(FFLAGS) -c $^ -o $@ -static

%.a : $(BINDIR)/%.o
	ar rcv $^ $@

clean:
	rm -f $(BINDIR)/*.o $(OBJECTDIR)/*.mod *.a
