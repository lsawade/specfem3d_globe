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

!----
!---- Green function database: Morton code computation and manifest writing
!----

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gf_compute_morton_codes()

! Computes center coordinates and Morton codes for all tagged local elements.
! Must be called while xstore/ystore/zstore_crust_mantle are still allocated.

  use constants, only: CUSTOM_REAL,MIDX,MIDY,MIDZ,IMAIN

  use specfem_par, only: myrank

  use specfem_par_crustmantle, only: &
    ibool => ibool_crust_mantle, &
    xstore => xstore_crust_mantle, &
    ystore => ystore_crust_mantle, &
    zstore => zstore_crust_mantle

  use green_function_par, only: &
    gf_nelem_local, gf_local_elements, &
    gf_morton_codes, gf_morton_hex, gf_center_xyz

  implicit none

  ! external functions
  integer(8), external :: gf_morton_encode

  ! local parameters
  integer :: i, ispec, iglob_center, ier

  ! NOTE: do not early-return on empty ranks. gf_report_morton_codes() below
  ! contains a collective (gather_all_i) that every rank must reach, otherwise
  ! ranks with no local GF elements would skip it and the run would deadlock.
  if (gf_nelem_local > 0) then
    ! allocate arrays
    allocate(gf_center_xyz(3, gf_nelem_local), stat=ier)
    if (ier /= 0) call exit_MPI(myrank, 'Error allocating gf_center_xyz')

    allocate(gf_morton_codes(gf_nelem_local), stat=ier)
    if (ier /= 0) call exit_MPI(myrank, 'Error allocating gf_morton_codes')

    allocate(gf_morton_hex(gf_nelem_local), stat=ier)
    if (ier /= 0) call exit_MPI(myrank, 'Error allocating gf_morton_hex')

    ! extract center GLL node coordinates and compute Morton codes
    do i = 1, gf_nelem_local
      ispec = gf_local_elements(i)
      iglob_center = ibool(MIDX, MIDY, MIDZ, ispec)

      gf_center_xyz(1, i) = xstore(iglob_center)
      gf_center_xyz(2, i) = ystore(iglob_center)
      gf_center_xyz(3, i) = zstore(iglob_center)

      ! gf_morton_encode is a standalone function defined below
      gf_morton_codes(i) = gf_morton_encode( &
        gf_center_xyz(1, i), gf_center_xyz(2, i), gf_center_xyz(3, i))

      write(gf_morton_hex(i), '(Z16.16)') gf_morton_codes(i)
    enddo
  endif

  ! report on rank 0 (collective: all ranks must call this)
  call gf_report_morton_codes()

  end subroutine gf_compute_morton_codes

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gf_report_morton_codes()

! Gathers element counts from all ranks and prints Morton code info on rank 0.

  use constants, only: IMAIN,CUSTOM_REAL
  use constants_solver, only: NPROCTOT_VAL

  use specfem_par, only: myrank

  use green_function_par, only: &
    gf_nelem_local, gf_morton_hex, gf_center_xyz

  implicit none

  integer :: irank, ier, j, nelem_total
  integer, allocatable :: counts_all(:)
  integer :: sendcount_dummy(1)

  ! gather counts on rank 0
  allocate(counts_all(0:NPROCTOT_VAL-1), stat=ier)
  if (ier /= 0) call exit_MPI(myrank, 'Error allocating counts_all in gf_report_morton_codes')

  sendcount_dummy(1) = gf_nelem_local
  call gather_all_i(sendcount_dummy, 1, counts_all, 1, NPROCTOT_VAL)

  if (myrank == 0) then
    nelem_total = sum(counts_all)

    write(IMAIN,*)
    write(IMAIN,*) '********************'
    write(IMAIN,*) ' Green function database: Morton codes'
    write(IMAIN,*) '********************'
    write(IMAIN,*)
    write(IMAIN,*) 'Total tagged elements across all ranks: ', nelem_total
    write(IMAIN,*)

    do irank = 0, NPROCTOT_VAL - 1
      write(IMAIN,'(a,i6,a,i6,a)') '  Rank ', irank, ': ', counts_all(irank), ' elements'
    enddo

    write(IMAIN,*)
    write(IMAIN,*) 'Rank 0 Morton codes (sample):'
    do j = 1, min(gf_nelem_local, 10)
      write(IMAIN,'(a,i6,a,a,a,3(f12.6))') '  elem ', j, ': ', &
        gf_morton_hex(j), '  center=', &
        gf_center_xyz(1,j), gf_center_xyz(2,j), gf_center_xyz(3,j)
    enddo
    if (gf_nelem_local > 10) then
      write(IMAIN,'(a,i6,a)') '  ... (', gf_nelem_local - 10, ' more)'
    endif
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  deallocate(counts_all)

  end subroutine gf_report_morton_codes

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gf_create_directories()

! Creates the database directory structure:
!   {GF_DATABASE_PATH}/elements/{morton_hex}/
! Each rank creates directories for its own elements.

  use constants, only: IMAIN,MAX_STRING_LEN

  use shared_parameters, only: GF_DATABASE_PATH

  use specfem_par, only: myrank

  use green_function_par, only: &
    gf_nelem_local, gf_morton_hex

  implicit none

  integer :: i
  character(len=MAX_STRING_LEN) :: dirpath
  character(len=MAX_STRING_LEN) :: command

  ! rank 0 creates base directories
  if (myrank == 0) then
    ! note: uses the system_command() wrapper, system() is a non-standard GNU extension
    command = 'mkdir -p ' // trim(GF_DATABASE_PATH) // '/elements'
    call system_command(command)

    write(IMAIN,*) 'Green function database: creating directories'
    write(IMAIN,*) '  base path: ', trim(GF_DATABASE_PATH)
    call flush_IMAIN()
  endif

  ! synchronize so base directory exists before ranks create subdirectories
  call synchronize_all()

  ! each rank creates directories for its local elements
  do i = 1, gf_nelem_local
    dirpath = trim(GF_DATABASE_PATH) // '/elements/' // gf_morton_hex(i)
    command = 'mkdir -p ' // trim(dirpath)
    call system_command(command)
  enddo

  call synchronize_all()

  if (myrank == 0) then
    write(IMAIN,*) 'Green function database: directories created'
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  end subroutine gf_create_directories

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gf_write_manifest()

! Writes manifest files:
!   {GF_DATABASE_PATH}/elements/centroids.bin — binary (morton_code + cx,cy,cz as float64)
!   {GF_DATABASE_PATH}/elements/manifest.csv  — human-readable
!
! Only rank 0 writes. Center coordinates (CUSTOM_REAL) are gathered from all
! ranks via send_cr/recv_cr. Rank 0 recomputes Morton codes from the gathered
! float values to avoid needing int(8) MPI transfers.

  use constants, only: IMAIN,MAX_STRING_LEN,CUSTOM_REAL
  use constants_solver, only: NPROCTOT_VAL

  use shared_parameters, only: GF_DATABASE_PATH

  use specfem_par, only: myrank

  use green_function_par, only: &
    gf_nelem_local, gf_center_xyz

  implicit none

  ! external functions
  integer(8), external :: gf_morton_encode

  integer :: irank, i, j, ier, nelem_total
  integer, allocatable :: counts_all(:), offsets(:)
  integer :: sendcount_dummy(1)

  ! gathered data on rank 0
  real(kind=CUSTOM_REAL), allocatable :: center_all(:,:)
  integer(8) :: morton_val
  character(len=16) :: hex_str

  ! file I/O
  character(len=MAX_STRING_LEN) :: filepath_bin, filepath_csv
  integer, parameter :: IUNIT_BIN = 71
  integer, parameter :: IUNIT_CSV = 72

  ! gather counts
  allocate(counts_all(0:NPROCTOT_VAL-1), stat=ier)
  if (ier /= 0) call exit_MPI(myrank, 'Error allocating counts_all in gf_write_manifest')
  sendcount_dummy(1) = gf_nelem_local
  call gather_all_i(sendcount_dummy, 1, counts_all, 1, NPROCTOT_VAL)

  if (myrank == 0) then
    ! build offsets
    allocate(offsets(0:NPROCTOT_VAL-1), stat=ier)
    if (ier /= 0) call exit_MPI(myrank, 'Error allocating offsets in gf_write_manifest')
    offsets(0) = 0
    do irank = 1, NPROCTOT_VAL - 1
      offsets(irank) = offsets(irank - 1) + counts_all(irank - 1)
    enddo
    nelem_total = offsets(NPROCTOT_VAL - 1) + counts_all(NPROCTOT_VAL - 1)
  else
    nelem_total = 0
  endif

  ! gather center coordinates on rank 0
  if (myrank == 0) then
    allocate(center_all(3, nelem_total), stat=ier)
    if (ier /= 0) call exit_MPI(myrank, 'Error allocating center_all')

    ! copy rank 0 data
    do i = 1, counts_all(0)
      center_all(1:3, i) = gf_center_xyz(1:3, i)
    enddo

    ! receive from other ranks
    do irank = 1, NPROCTOT_VAL - 1
      if (counts_all(irank) > 0) then
        j = offsets(irank) + 1
        call recv_cr(center_all(1, j), 3 * counts_all(irank), irank, 100)
      endif
    enddo
  else
    ! non-zero ranks send their data
    if (gf_nelem_local > 0) then
      call send_cr(gf_center_xyz, 3 * gf_nelem_local, 0, 100)
    endif
  endif

  ! rank 0 writes manifest files
  if (myrank == 0) then

    filepath_bin = trim(GF_DATABASE_PATH) // '/elements/centroids.bin'
    filepath_csv = trim(GF_DATABASE_PATH) // '/elements/manifest.csv'

    ! write binary manifest
    open(unit=IUNIT_BIN, file=trim(filepath_bin), form='unformatted', &
         access='stream', status='replace', iostat=ier)
    if (ier /= 0) call exit_MPI(myrank, 'Error opening centroids.bin')

    ! write CSV manifest
    open(unit=IUNIT_CSV, file=trim(filepath_csv), form='formatted', &
         status='replace', iostat=ier)
    if (ier /= 0) call exit_MPI(myrank, 'Error opening manifest.csv')
    write(IUNIT_CSV, '(a)') 'morton_hex,cx,cy,cz'

    do i = 1, nelem_total
      ! recompute Morton code from gathered float coordinates
      morton_val = gf_morton_encode(center_all(1,i), center_all(2,i), center_all(3,i))
      write(hex_str, '(Z16.16)') morton_val

      ! binary: int64 morton + 3x float64 centroid = 32 bytes per record
      write(IUNIT_BIN) morton_val, dble(center_all(1,i)), &
                        dble(center_all(2,i)), dble(center_all(3,i))

      ! CSV
      write(IUNIT_CSV, '(a,a,ES22.14,a,ES22.14,a,ES22.14)') &
        hex_str, ',', dble(center_all(1,i)), ',', dble(center_all(2,i)), &
        ',', dble(center_all(3,i))
    enddo

    close(IUNIT_BIN)
    close(IUNIT_CSV)

    write(IMAIN,*) 'Green function database: manifest files written'
    write(IMAIN,*) '  centroids.bin: ', nelem_total, ' records (32 bytes each)'
    write(IMAIN,*) '  manifest.csv:  ', nelem_total, ' entries'
    write(IMAIN,*)
    call flush_IMAIN()

    deallocate(center_all, offsets)
  endif

  deallocate(counts_all)

  end subroutine gf_write_manifest

!
!-------------------------------------------------------------------------------------------------
!

  function gf_morton_encode(cx, cy, cz) result(morton)

! Standalone Morton encoder (not inside a contains block) so it can be
! called from gf_write_manifest.

  use constants, only: CUSTOM_REAL

  implicit none

  real(kind=CUSTOM_REAL), intent(in) :: cx, cy, cz
  integer(8) :: morton
  integer(8) :: ix, iy, iz

  ! external functions
  integer(8), external :: gf_spread_bits

  integer, parameter :: MORTON_BITS = 21
  integer, parameter :: MORTON_BINS = 2**MORTON_BITS  ! 2,097,152

  ! quantize from [-1, 1] to [0, MORTON_BINS - 1]
  ix = int((real(cx) + 1.0) * 0.5 * real(MORTON_BINS - 1), 8)
  iy = int((real(cy) + 1.0) * 0.5 * real(MORTON_BINS - 1), 8)
  iz = int((real(cz) + 1.0) * 0.5 * real(MORTON_BINS - 1), 8)

  ! clamp to valid range
  ix = max(int(0, 8), min(ix, int(MORTON_BINS - 1, 8)))
  iy = max(int(0, 8), min(iy, int(MORTON_BINS - 1, 8)))
  iz = max(int(0, 8), min(iz, int(MORTON_BINS - 1, 8)))

  ! bit-interleave: x in bits 0,3,6,...; y in bits 1,4,7,...; z in bits 2,5,8,...
  morton = ior(ior(gf_spread_bits(ix), ishft(gf_spread_bits(iy), 1)), &
               ishft(gf_spread_bits(iz), 2))

  end function gf_morton_encode

!
!-------------------------------------------------------------------------------------------------
!

  function gf_spread_bits(v) result(spread)

! Spread 21 low bits of v into every-third-bit positions.
! Standalone version callable from gf_morton_encode and gf_compute_morton_codes.

  implicit none
  integer(8), intent(in) :: v
  integer(8) :: spread

  spread = v
  spread = iand(ior(spread, ishft(spread, 32)), int(Z'001F00000000FFFF', 8))
  spread = iand(ior(spread, ishft(spread, 16)), int(Z'001F0000FF0000FF', 8))
  spread = iand(ior(spread, ishft(spread,  8)), int(Z'100F00F00F00F00F', 8))
  spread = iand(ior(spread, ishft(spread,  4)), int(Z'10C30C30C30C30C3', 8))
  spread = iand(ior(spread, ishft(spread,  2)), int(Z'1249249249249249', 8))

  end function gf_spread_bits
