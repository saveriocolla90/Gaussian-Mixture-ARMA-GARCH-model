PROGRAM = EM_core
PROGRAM_OBJ = $(PROGRAM).o

OBJECTS = EM_core.o mod_PARAMS.o mod_TSERIES.o mod_MODELS.o mod_DERIVATIVES.o mod_STAT.o mod_FITTING.o mod_MATH.o mod_IO.o
MODULES = params.mod tseries.mod models.mod derivatives.mod stat.mod fitting.mod math.mod io.mod

SUBF = ./datetime-fortran/src/
SUBMOD = datetime_module.o

FC = gfortran
FCFLAGS = -c -Wall -Wno-tabs -Wextra -Wconversion -Og -pedantic -fcheck=all -fmax-errors=5

# Make without parameters will make first target found
default : $(PROGRAM).exe


#--EXECUTABLE
$(PROGRAM).exe: $(OBJECTS)
	$(FC) $(OBJECTS) -o $(PROGRAM).exe $(SUBF)$(SUBMOD)


#--MODULES FILES
params.mod: mod_PARAMS.o mod_PARAMS.f90  
	$(FC) $(FCFLAGS) mod_PARAMS.f90
	
tseries.mod: mod_TSERIES.o mod_TSERIES.f90
	$(FC) $(FCFLAGS) mod_TSERIES.f90 -I$(SUBF)

models.mod: mod_MODELS.o mod_MODELS.f90  
	$(FC) $(FCFLAGS) mod_MODELS.f90 

io.mod: mod_IO.o mod_IO.f90
	$(FC) $(FCFLAGS) mod_IO.f90

derivatives.mod: mod_DERIVATIVES.o mod_DERIVATIVES.f90
	$(FC) $(FCFLAGS) mod_DERIVATIVES.f90 

stat.mod: mod_STAT.o mod_STAT.f90
	$(FC) $(FCFLAGS) mod_STAT.f90 

fitting.mod: mod_FITTING.o mod_FITTING.f90
	$(FC) $(FCFLAGS) mod_FITTING.f90 
	
math.mod: mod_MATH.o mod_MATH.f90
	$(FC) $(FCFLAGS) mod_MATH.f90 



#--OBJECTS FILES
mod_PARAMS.o: mod_PARAMS.f90
	$(FC) $(FCFLAGS) mod_PARAMS.f90

mod_TSERIES.o: params.mod mod_TSERIES.f90
	$(FC) $(FCFLAGS) mod_TSERIES.f90 -I$(SUBF)

mod_MODELS.o: tseries.mod params.mod mod_MODELS.f90
	$(FC) $(FCFLAGS) mod_MODELS.f90 

mod_IO.o: models.mod params.mod mod_IO.f90
	$(FC) $(FCFLAGS) mod_IO.f90

mod_DERIVATIVES.o: tseries.mod models.mod params.mod mod_DERIVATIVES.f90
	$(FC) $(FCFLAGS) mod_DERIVATIVES.f90 

mod_STAT.o: tseries.mod params.mod mod_STAT.f90
	$(FC) $(FCFLAGS) mod_STAT.f90 

mod_FITTING.o: tseries.mod params.mod models.mod derivatives.mod math.mod mod_FITTING.f90
	$(FC) $(FCFLAGS) mod_FITTING.f90 

mod_MATH.o: mod_MATH.f90
	$(FC) $(FCFLAGS) mod_MATH.f90 

$(PROGRAM_OBJ): io.mod tseries.mod models.mod params.mod fitting.mod stat.mod EM_core.f90
	$(FC) $(FCFLAGS) $(PROGRAM).f90 


# Debug to show the value of each variable
debug:
	@echo "OBJS = $(OBJECTS)"
	@echo "MODS = $(MODULES)"
	@echo "SUBMODS = $(SUBMOD)"
	@echo "PROGRAM = $(PROGRAM)"
	@echo "PRG_OBJ = $(PROGRAM_OBJ)"
	

clean:
	rm -rf $(OBJECTS)
	rm -rf $(MODULES)
	rm -rf $(PROGRAM).exe
