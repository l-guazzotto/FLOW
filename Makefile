F95 = gfortran
#F95 = pgfortran

#FFLAGS =  -O3 -fPIC -ffree-line-length-0
#FFLAGS = -O3 -fPIC -ffree-line-length-0 -fcheck=all -Wall -g -fbacktrace -fbounds-check 
FFLAGS =  -Og -g  -ffree-line-length-0
FFLAGS2 =   -ffree-line-length-0
LINK	= $(F95) $(FFLAGS)
LINK2	= $(F95) $(FFLAGS2)

OBJECTS =  \
constant.o array_dimensions.o pseudo_IMSL.o  interpolating_functions.o \
	magnetic.o exp_data.o p_d_profile.o triangularity.o flow.o trans_solve.o  \
	readinput.o data_dump.o mgrid.o inter_grid.o ellipt.o \
	nclass_mod.o rarray_zero.o u_erf.o u_lu_backsub.o u_lu_decomp.o main.o

FLOW:  $(OBJECTS) 
	$(LINK) $(OBJECTS)  -o FLOW_2017_exe_mac

debug:  $(OBJECTS) 
	$(LINK2) $(OBJECTS)  -o FLOW_2017_exe_mac


%.o : %.f90
	$(F95) $(FFLAGS)  -c $<

%.o : %.for
	$(F95) $(FFLAGS)  -c $<


data_dump.o: data_dump.f90 trans_solve.o magnetic.o constant.o triangularity.o exp_data.o \
	pseudo_IMSL.o ellipt.o
	$(F95) $(FFLAGS)  -c data_dump.f90

ellipt.o: ellipt.f90 constant.o
	$(F95) $(FFLAGS)  -c ellipt.f90

exp_data.o: exp_data.f90 constant.o
	$(F95) $(FFLAGS)  -c exp_data.f90


flow.o: flow.f90 constant.o
	$(F95) $(FFLAGS)  -c flow.f90


inter_grid.o: inter_grid.f90 constant.o
	$(F95) $(FFLAGS)  -c inter_grid.f90

interpolating_functions.o: interpolating_functions.f90 constant.o pseudo_IMSL.o
	$(F95) $(FFLAGS)  -c interpolating_functions.f90

magnetic.o: magnetic.f90 constant.o
	$(F95) $(FFLAGS)  -c magnetic.f90

mgrid.o: mgrid.f90 constant.o trans_solve.o triangularity.o magnetic.o array_dimensions.o inter_grid.o
	$(F95) $(FFLAGS)  -c mgrid.f90

pseudo_IMSL.o: pseudo_IMSL.f90 constant.o
	$(F95) $(FFLAGS)  -c pseudo_IMSL.f90

p_d_profile.o: p_d_profile.f90 constant.o
	$(F95) $(FFLAGS)  -c p_d_profile.f90

readinput.o: readinput.f90  array_dimensions.o constant.o flow.o magnetic.o p_d_profile.o \
             exp_data.o triangularity.o trans_solve.o interpolating_functions.o pseudo_IMSL.o
	$(F95) $(FFLAGS)  -c readinput.f90

trans_solve.o: trans_solve.f90 constant.o magnetic.o p_d_profile.o flow.o triangularity.o \
	exp_data.o interpolating_functions.o pseudo_IMSL.o
	$(F95) $(FFLAGS)  -c trans_solve.f90

triangularity.o: triangularity.f90 array_dimensions.o constant.o interpolating_functions.o pseudo_IMSL.o \
	exp_data.o
	$(F95) $(FFLAGS)  -c triangularity.f90

main.o: main.f90 constant.o trans_solve.o  exp_data.o triangularity.o magnetic.o
	$(F95) $(FFLAGS)  -c main.f90


clean:
	rm -rf *.o

realclean:
	rm -rf *.o *.mod


cleanall:
	rm -rf *.o *.mod FLOW_2017_exe_mac
