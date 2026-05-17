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

  subroutine prepare_green_function_storage()

! Central wrapper for all Green function database preparation steps.
! Called from prepare_timerun() before coordinate arrays are deallocated.

  use shared_parameters, only: GF_DATABASE_ENABLED, TOPOGRAPHY

  use specfem_par, only: ibathy_topo

  implicit none

  if (.not. GF_DATABASE_ENABLED) return

  ! precompute Butterworth-filtered STF
  call gf_compute_stf()

  ! compute Morton codes from element center coordinates
  call gf_compute_morton_codes()

  ! create directory structure
  call gf_create_directories()

  ! write per-element coordinate files (idempotent — skips if exists)
  call gf_write_coordinates()

  ! write shared mesh info file (topography + ellipticity, rank 0, idempotent)
  call gf_write_mesh_info()

  ! ibathy_topo was kept alive for gf_write_mesh_info — deallocate now
  if (TOPOGRAPHY) then
    if (allocated(ibathy_topo)) deallocate(ibathy_topo)
  endif

  ! write per-station metadata file (rank 0, create-or-skip)
  call gf_write_station_metadata()

  ! initialize HDF5 files and allocate write buffer
  call gf_init_hdf5()

  end subroutine prepare_green_function_storage
