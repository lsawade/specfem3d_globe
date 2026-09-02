# includes default Makefile from previous configuration
include Makefile

# test target
default: test_gf_stf

## compilation directories
O := ./obj

OBJECTS = \
	$(EMPTY_MACRO)

# The test links against the filter subroutines extracted into a local copy.
# The main green_function_stf.F90 has specfem module dependencies that we
# don't need for pure filter testing.
test_gf_stf:
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -o ./bin/test_gf_stf test_gf_stf.f90 $(LDFLAGS)
