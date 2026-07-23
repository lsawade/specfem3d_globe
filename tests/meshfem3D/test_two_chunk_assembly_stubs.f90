!=====================================================================
!
!                 S p e c f e m 3 D  G l o b e  V e r s i o n  7 . 0
!                 ---------------------------------------------------
!
!     Main authors: Dimitri Komatitsch and Jeroen Tromp
!                    Princeton University, USA and CNRS, France
!                    (c) August 2020
!
!     This program is free software: you can redistribute it and/or modify
!     it under the terms of the GNU General Public License as published by
!     the Free Software Foundation, either version 3 of the License, or
!     (at your option) any later version.
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!     GNU General Public License for more details.
!
!     You should have received a copy of the GNU General Public License
!     along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
!=====================================================================
!
! Link-only stubs for init_mpi() dependencies in parallel.f90.
! The assembly regression initializes MPI directly and never calls init_mpi().

  subroutine open_parameter_file_from_main_only(ier)

  implicit none

  integer, intent(out) :: ier

  ier = 0

  end subroutine open_parameter_file_from_main_only

  subroutine read_value_integer(value,name,ier)

  implicit none

  integer, intent(out) :: value,ier
  character(len=*), intent(in) :: name

  value = 0
  ier = 0

  end subroutine read_value_integer

  subroutine read_value_logical(value,name,ier)

  implicit none

  logical, intent(out) :: value
  integer, intent(out) :: ier
  character(len=*), intent(in) :: name

  value = .false.
  ier = 0

  end subroutine read_value_logical

  subroutine close_parameter_file()

  implicit none

  end subroutine close_parameter_file
