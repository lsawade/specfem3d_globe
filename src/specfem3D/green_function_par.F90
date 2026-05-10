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

module green_function_par

  use constants, only: MAX_STRING_LEN,CUSTOM_REAL

  implicit none

  ! Force component detection (Stage 2)
  ! 1=N, 2=E, 3=Z
  integer :: gf_force_component = 0
  character(len=8) :: gf_network_name = ''
  character(len=32) :: gf_station_name = ''

  ! STF (Stage 3)
  ! precomputed source time function: Gaussian filtered with Butterworth lowpass
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: gf_stf
  double precision :: gf_f_cutoff = 0.d0   ! lowpass cutoff frequency (Hz)
  double precision :: gf_hdur = 0.d0       ! half-duration used for Gaussian

  ! Element location (Stage 4)
  integer :: gf_nlocations = 0                                     ! total GF locations read
  double precision, dimension(:), allocatable :: gf_lat, gf_lon, gf_depth  ! location coords
  integer, dimension(:), allocatable :: gf_ispec_selected          ! element index per location
  integer, dimension(:), allocatable :: gf_islice_selected         ! MPI rank per location
  integer :: gf_nelem_local = 0                                    ! unique tagged elements on this rank
  integer, dimension(:), allocatable :: gf_local_elements          ! unique local element indices

  ! Morton codes (Stage 6)
  integer(8), dimension(:), allocatable :: gf_morton_codes          ! Morton code per local element
  character(len=16), dimension(:), allocatable :: gf_morton_hex     ! hex string per local element
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: gf_center_xyz  ! (3, gf_nelem_local) center coords

end module green_function_par
