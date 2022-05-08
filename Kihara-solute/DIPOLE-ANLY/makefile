#COMP= f95 
COMP= ifort
#OPT= -traceback -check all #-fast
OPT= -O3 #-fast
#OPT= -fPIC -shared-intel -mcmodel=large #-fast
OBJ = modules.o main.o aldeArrays.o dcdread.o psfread.o pdbread.o dipole.o atominfo.o maths.o 
all: $(OBJ)
	$(COMP) $(OPT) -o Dipole-anly.e $(OBJ)
modules.o: modules.f90
	$(COMP) $(OPT) -c modules.f90 -o modules.o
main.o: main.f90
	$(COMP) $(OPT) -c main.f90 -o main.o
aldeArrays.o: aldeArrays.f90
	$(COMP) $(OPT) -c aldeArrays.f90 -o aldeArrays.o
dcdread.o: dcdread.f90
	$(COMP) $(OPT) -c dcdread.f90 -o dcdread.o
psfread.o: psfread.f90
	$(COMP) $(OPT) -c psfread.f90 -o psfread.o
pdbread.o: pdbread.f90
	$(COMP) $(OPT) -c pdbread.f90 -o pdbread.o
dipole.o: dipole.f90
	$(COMP) $(OPT) -c dipole.f90 -o dipole.o
atominfo.o: atominfo.f90
	$(COMP) $(OPT) -c atominfo.f90 -o atominfo.o
maths.o: maths.f90
	$(COMP) $(OPT) -c maths.f90 -o maths.o
clean:
	rm -rf *.o *.mod
