F90_FILES := $(wildcard src/*.f90)
OBJ_FILES := $(addprefix obj/,$(notdir $(F90_FILES:.f90=.o)))

F90=gfortran

#FFLAGS=-g -CB -traceback
FFLAGS=-g -ffree-form

main: $(OBJ_FILES)
	$(F90) $(FFLAGS) -o syminfo $(OBJ_FILES)

obj/syminfo.o: src/syminfo.f90 obj/datatypes.o obj/structure.o obj/symmetry.o obj/iodef.o 
	$(F90) $(FFLAGS) -c -o $@ $<

obj/datatypes.o: src/datatypes.f90 obj/iodef.o

obj/%.o: src/%.f90 obj
	$(F90) $(FFLAGS) -c -o $@ $<

obj:
	mkdir obj

clean: 
	rm -rf obj
	rm *.mod
