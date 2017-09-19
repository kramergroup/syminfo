F90_FILES := $(wildcard src/*.f90)
OBJ_FILES := $(addprefix obj/,$(notdir $(F90_FILES:.f90=.o)))

F90=ifort

FFLAGS=-g -CB -traceback

main: $(OBJ_FILES)
	$(F90) $(FFLAGS) -o syminfo $(OBJ_FILES)

obj/syminfo.o: src/syminfo.f90 obj/datatypes.o obj/structure.o obj/symmetry.o 
	$(F90) $(FFLAGS) -c -o $@ $<

obj/%.o: src/%.f90 obj
	$(F90) $(FFLAGS) -c -o $@ $<

obj:
	mkdir obj

clean: 
	rm -rf obj
	rm *.mod
