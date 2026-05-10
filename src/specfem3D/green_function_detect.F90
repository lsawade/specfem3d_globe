!=====================================================================
!
!                       S p e c f e m 3 D  G l o b e
!                       ----------------------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!                and CNRS / University of Marseille, France
!                 (there are currently many more authors!)
! (c) Princeton University and CNRS / University of Marseille, April 2014
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along
! with this program; if not, write to the Free Software Foundation, Inc.,
! 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
!
!=====================================================================

  subroutine gf_detect_force_component()

! Detects which force component (N/E/Z) is being simulated and extracts
! the station identity from the FORCESOLUTION header.
!
! Called after setup_sources() when GF_DATABASE_ENABLED = .true.
! At that point, comp_dir_vect_source_E/N/Z_UP have already been read
! and broadcast to all ranks.

  use constants, only: IIN,IMAIN,MAX_STRING_LEN,TINYVAL,mygroup

  use shared_parameters, only: NSOURCES,NUMBER_OF_SIMULTANEOUS_RUNS

  use specfem_par, only: myrank, &
    comp_dir_vect_source_E,comp_dir_vect_source_N,comp_dir_vect_source_Z_UP

  use green_function_par, only: gf_force_component,gf_network_name,gf_station_name

  implicit none

  ! local variables
  character(len=MAX_STRING_LEN) :: string,FORCESOLUTION,path_to_add
  integer :: ier,idot
  logical :: has_N,has_E,has_Z
  integer :: ncomp

  ! validate: GF database requires exactly one source
  if (NSOURCES /= 1) then
    stop 'Error: GF_DATABASE_ENABLED requires exactly 1 source (NSOURCES must be 1)'
  endif

  ! detect force component from direction vectors (available on all ranks after broadcast)
  has_N = abs(comp_dir_vect_source_N(1)) > TINYVAL
  has_E = abs(comp_dir_vect_source_E(1)) > TINYVAL
  has_Z = abs(comp_dir_vect_source_Z_UP(1)) > TINYVAL

  ncomp = 0
  if (has_N) ncomp = ncomp + 1
  if (has_E) ncomp = ncomp + 1
  if (has_Z) ncomp = ncomp + 1

  if (ncomp /= 1) then
    print *, 'Error: GF database requires exactly one nonzero force direction component'
    print *, '  comp_dir_vect_source_N = ', comp_dir_vect_source_N(1)
    print *, '  comp_dir_vect_source_E = ', comp_dir_vect_source_E(1)
    print *, '  comp_dir_vect_source_Z_UP = ', comp_dir_vect_source_Z_UP(1)
    stop 'Error: GF_DATABASE_ENABLED requires exactly one nonzero force component (N, E, or Z)'
  endif

  if (has_N) gf_force_component = 1
  if (has_E) gf_force_component = 2
  if (has_Z) gf_force_component = 3

  ! parse station identity from FORCESOLUTION header (rank 0 only, then broadcast)
  if (myrank == 0) then
    FORCESOLUTION = 'DATA/FORCESOLUTION'
    if (NUMBER_OF_SIMULTANEOUS_RUNS > 1 .and. mygroup >= 0) then
      write(path_to_add,"('run',i4.4,'/')") mygroup + 1
      FORCESOLUTION = path_to_add(1:len_trim(path_to_add))//FORCESOLUTION(1:len_trim(FORCESOLUTION))
    endif

    open(unit=IIN,file=trim(FORCESOLUTION),status='old',action='read',iostat=ier)
    if (ier /= 0) then
      print *,'Error opening file: ',trim(FORCESOLUTION)
      stop 'Error opening FORCESOLUTION file for GF station detection'
    endif

    ! read header line (e.g., "FORCE  IU.SJG")
    read(IIN,"(a)") string
    ! skip empty lines
    do while (len_trim(string) == 0)
      read(IIN,"(a)") string
    enddo
    close(IIN)

    ! extract station identity: everything after "FORCE" prefix, trimmed
    ! header format: "FORCE  IU.SJG" or "FORCE  001"
    ! skip the first 5 characters ("FORCE"), then trim leading spaces
    string = adjustl(string(6:))

    ! split on '.' to get network and station
    idot = index(trim(string), '.')
    if (idot > 0) then
      gf_network_name = string(1:idot-1)
      gf_station_name = string(idot+1:len_trim(string))
    else
      ! no dot found — use entire label as station name, empty network
      gf_network_name = ''
      gf_station_name = trim(string)
    endif

    ! user output
    write(IMAIN,*)
    write(IMAIN,*) 'Green function database:'
    write(IMAIN,*) '  network: ', trim(gf_network_name)
    write(IMAIN,*) '  station: ', trim(gf_station_name)
    if (gf_force_component == 1) then
      write(IMAIN,*) '  force component: N (1)'
    else if (gf_force_component == 2) then
      write(IMAIN,*) '  force component: E (2)'
    else
      write(IMAIN,*) '  force component: Z (3)'
    endif
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  ! broadcast station identity to all ranks
  call bcast_all_ch(gf_network_name, 8)
  call bcast_all_ch(gf_station_name, 32)

  end subroutine gf_detect_force_component
