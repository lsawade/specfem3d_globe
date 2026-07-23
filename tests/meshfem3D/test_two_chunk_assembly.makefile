# includes default Makefile from previous configuration
include Makefile

# test target
default: test_two_chunk_assembly

## compilation directories
O := ./obj

OBJECTS = \
	$O/assemble_MPI_scalar.shared.o \
	$O/assemble_MPI_vector.shared.o \
	$O/parallel.sharedmpi.o \
	$O/shared_par.shared_module.o \
	$(EMPTY_MACRO)

test_two_chunk_assembly:
	${MPIFCCOMPILE_CHECK} ${FCFLAGS_f90} -o ./bin/test_two_chunk_assembly \
		test_two_chunk_assembly.f90 test_two_chunk_assembly_stubs.f90 -I./obj $(OBJECTS) $(MPILIBS)
