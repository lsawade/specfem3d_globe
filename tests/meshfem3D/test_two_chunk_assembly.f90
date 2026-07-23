!=====================================================================
!
!                       S p e c f e m 3 D  G l o b e
!                       ----------------------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!                and CNRS / University of Marseille, France
! (c) Princeton University and CNRS / University of Marseille, April 2014
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
! (at your option) any later version.
!
!=====================================================================

  program test_two_chunk_assembly

  use mpi
  use constants, only: CUSTOM_REAL,NDIM
  use my_mpi, only: my_local_mpi_comm_world

  implicit none

  integer, parameter :: NGLOB_TEST = 2
  integer, parameter :: NUM_OWNERS = 4
  integer, parameter :: OWNERS_FOUR(NUM_OWNERS) = (/ 0,1,2,3 /)
  integer, parameter :: OWNERS_EIGHT(NUM_OWNERS) = (/ 0,2,5,7 /)
  integer, parameter :: EXPECTED_SUM = 15

  integer :: ierr,myrank,nproc,member,interface_index,num_interfaces
  integer, dimension(:), allocatable :: my_neighbors,nibool_interfaces
  integer, dimension(:,:), allocatable :: ibool_interfaces
  integer, dimension(NUM_OWNERS) :: owners
  logical :: is_owner
  real(kind=CUSTOM_REAL) :: scalar(NGLOB_TEST)
  real(kind=CUSTOM_REAL) :: vector(NDIM,NGLOB_TEST)
  real(kind=CUSTOM_REAL) :: contribution

  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD,myrank,ierr)
  call MPI_Comm_size(MPI_COMM_WORLD,nproc,ierr)
  my_local_mpi_comm_world = MPI_COMM_WORLD

  if (nproc == 4) then
    owners = OWNERS_FOUR
  else if (nproc == 8) then
    owners = OWNERS_EIGHT
  else
    call fail('expected four or eight MPI ranks')
  endif

  is_owner = any(owners == myrank)
  if (is_owner) then
    num_interfaces = NUM_OWNERS - 1
  else
    num_interfaces = 0
  endif
  allocate(my_neighbors(num_interfaces),nibool_interfaces(num_interfaces), &
           ibool_interfaces(1,num_interfaces),stat=ierr)
  if (ierr /= 0) call fail('allocation failed')

  contribution = 0.0_CUSTOM_REAL
  if (is_owner) then
    do member = 1,NUM_OWNERS
      if (owners(member) == myrank) contribution = real(2 ** (member - 1),CUSTOM_REAL)
    enddo

    interface_index = 0
    do member = 1,NUM_OWNERS
      if (owners(member) /= myrank) then
        interface_index = interface_index + 1
        my_neighbors(interface_index) = owners(member)
      endif
    enddo
    nibool_interfaces(:) = 1
    ibool_interfaces(:,:) = 1
  endif

  scalar(:) = real(1000 + myrank,CUSTOM_REAL)
  if (is_owner) scalar(1) = contribution
  call assemble_MPI_scalar(nproc,NGLOB_TEST,scalar,num_interfaces,1, &
                           nibool_interfaces,ibool_interfaces,my_neighbors)
  if (is_owner) then
    call require(scalar(1) == real(EXPECTED_SUM,CUSTOM_REAL),'scalar assembly')
    call require(scalar(2) == real(1000 + myrank,CUSTOM_REAL),'scalar non-owner value')
  endif

  vector(:,:) = real(2000 + myrank,CUSTOM_REAL)
  if (is_owner) vector(:,1) = (/ contribution,2.0_CUSTOM_REAL * contribution, &
                                  4.0_CUSTOM_REAL * contribution /)
  call assemble_MPI_vector(nproc,NGLOB_TEST,vector,num_interfaces,1, &
                           nibool_interfaces,ibool_interfaces,my_neighbors)
  if (is_owner) then
    call require(all(vector(:,1) == (/ real(EXPECTED_SUM,CUSTOM_REAL), &
                                       real(2 * EXPECTED_SUM,CUSTOM_REAL), &
                                       real(4 * EXPECTED_SUM,CUSTOM_REAL) /)), &
                 'vector assembly')
    call require(all(vector(:,2) == real(2000 + myrank,CUSTOM_REAL)), &
                 'vector non-owner value')
  endif

  if (myrank == 0) then
    if (nproc == 8) then
      write(*,'(a)') 'PASS two-chunk owner set {0,2,5,7}: scalar/vector sum = 15'
    else
      write(*,'(a)') 'PASS four-owner scalar/vector assembly sum = 15'
    endif
  endif

  deallocate(my_neighbors,nibool_interfaces,ibool_interfaces)
  call MPI_Finalize(ierr)

contains

  subroutine require(condition,message)

    logical, intent(in) :: condition
    character(len=*), intent(in) :: message

    if (.not. condition) call fail(message)

  end subroutine require

  subroutine fail(message)

    character(len=*), intent(in) :: message

    write(*,'(a,i0,2a)') 'FAIL rank ',myrank,': ',trim(message)
    call MPI_Abort(MPI_COMM_WORLD,1,ierr)

  end subroutine fail

  end program test_two_chunk_assembly
