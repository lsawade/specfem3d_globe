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

!----------------------------------------------------------------------
! Green function STF: Gaussian source time function filtered with
! a Butterworth lowpass at the Nyquist frequency of the subsampled output.
!
! Uses cascaded second-order sections (SOS) with forward+backward
! filtering for zero-phase response (equivalent to scipy filtfilt).
!----------------------------------------------------------------------

!
!----------------------------------------------------------------------
!

  subroutine gf_compute_stf()

! Computes a Gaussian STF filtered with a Butterworth lowpass.
! Stores the result in gf_stf(:) array (allocated here).

  use constants, only: IMAIN,PI,CUSTOM_REAL,MAX_STRING_LEN
  use specfem_par, only: DT,NSTEP,t0,hdur,hdur_Gaussian,myrank,OUTPUT_FILES
  use shared_parameters, only: GF_SUBSAMPLE_STEP,T_min_period
  use green_function_par, only: gf_stf,gf_f_cutoff,gf_hdur

  implicit none

  ! local variables
  integer :: it,ier
  double precision :: timeval
  double precision, dimension(:), allocatable :: stf_raw,stf_unfiltered
  character(len=MAX_STRING_LEN) :: filename

  ! Butterworth SOS filter variables
  ! For order 4: 2 second-order sections, each with 6 coefficients (b0,b1,b2,a0,a1,a2)
  integer, parameter :: GF_FILTER_ORDER = 4
  integer :: nsections
  double precision, dimension(GF_FILTER_ORDER/2, 6) :: sos
  integer, parameter :: IOUT_STF = 71

  ! set half-duration from mesh resolution: T_min_period / 10
  ! This keeps the Gaussian broadband while limiting unresolvable energy
  ! (~14% above 1/T_min). The Butterworth lowpass handles anti-aliasing
  ! for the subsampled output.
  gf_hdur = T_min_period / 10.0d0

  ! compute cutoff frequency: Nyquist of subsampled output
  gf_f_cutoff = 1.0d0 / (2.0d0 * DT * dble(GF_SUBSAMPLE_STEP))

  ! allocate STF arrays
  allocate(gf_stf(NSTEP))
  allocate(stf_raw(NSTEP))
  allocate(stf_unfiltered(NSTEP))

  ! generate Gaussian STF (same formula as comp_source_time_function_gauss)
  do it = 1, NSTEP
    timeval = dble(it-1) * DT - t0
    stf_raw(it) = exp(-(timeval**2) / (gf_hdur**2)) / (sqrt(PI) * gf_hdur)
  enddo

  ! save unfiltered for output
  stf_unfiltered(:) = stf_raw(:)

  ! compute Butterworth lowpass SOS coefficients via bilinear transform
  nsections = GF_FILTER_ORDER / 2
  call gf_butterworth_sos(GF_FILTER_ORDER, gf_f_cutoff, 1.0d0/DT, nsections, sos)

  ! apply zero-phase filter (forward + backward pass through each SOS)
  call gf_sosfiltfilt(stf_raw, NSTEP, sos, nsections)

  ! store as CUSTOM_REAL
  do it = 1, NSTEP
    gf_stf(it) = real(stf_raw(it), kind=CUSTOM_REAL)
  enddo

  ! write STF to OUTPUT_FILES for cross-validation
  if (myrank == 0) then
    filename = trim(OUTPUT_FILES) // '/gf_stf.txt'
    open(unit=IOUT_STF, file=trim(filename), status='unknown', iostat=ier)
    if (ier == 0) then
      write(IOUT_STF, '(a)') '# Green function STF'
      write(IOUT_STF, '(a,ES20.10)') '# hdur = ', gf_hdur
      write(IOUT_STF, '(a,ES20.10)') '# f_cutoff = ', gf_f_cutoff
      write(IOUT_STF, '(a,ES20.10)') '# dt = ', DT
      write(IOUT_STF, '(a,I10)') '# nstep = ', NSTEP
      write(IOUT_STF, '(a,ES20.10)') '# t0 = ', t0
      write(IOUT_STF, '(a,I10)') '# filter_order = ', GF_FILTER_ORDER
      write(IOUT_STF, '(a,I10)') '# subsample_step = ', GF_SUBSAMPLE_STEP
      write(IOUT_STF, '(a)') '# time(s)  stf_unfiltered  stf_filtered'
      do it = 1, NSTEP
        timeval = dble(it-1) * DT - t0
        write(IOUT_STF, '(3ES20.10)') timeval, stf_unfiltered(it), stf_raw(it)
      enddo
      close(IOUT_STF)
    endif
  endif

  deallocate(stf_raw)
  deallocate(stf_unfiltered)

  ! user output
  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) 'Green function STF:'
    write(IMAIN,*) '  half-duration (Gaussian):  ', sngl(gf_hdur), ' s'
    write(IMAIN,*) '  lowpass cutoff frequency:  ', sngl(gf_f_cutoff), ' Hz'
    write(IMAIN,*) '  filter order:              ', GF_FILTER_ORDER
    write(IMAIN,*) '  subsample step:            ', GF_SUBSAMPLE_STEP
    write(IMAIN,*) '  STF peak value:            ', maxval(gf_stf)
    write(IMAIN,*)
  endif

  end subroutine gf_compute_stf

!
!----------------------------------------------------------------------
!

  subroutine gf_butterworth_sos(order, fc, fs, nsections, sos)

! Computes second-order section (SOS) coefficients for a Butterworth
! lowpass filter using the bilinear transform.
!
! Each section has 6 coefficients: [b0, b1, b2, a0, a1, a2]
! where a0 is always normalized to 1.0.

  implicit none

  double precision, parameter :: PI = 3.141592653589793238462643383279502884197d0

  integer, intent(in) :: order, nsections
  double precision, intent(in) :: fc, fs
  double precision, intent(out) :: sos(nsections, 6)

  ! local variables
  integer :: k
  double precision :: wc, wc2
  double precision :: pole_real, pole_imag
  double precision :: b0, b1, b2, a0, a1, a2

  ! pre-warp the cutoff frequency for bilinear transform
  ! wc = 2 * fs * tan(pi * fc / fs)
  wc = 2.0d0 * fs * tan(PI * fc / fs)
  wc2 = wc * wc

  ! for each conjugate pole pair, compute a second-order section
  do k = 1, nsections
    ! analog Butterworth poles (unit circle, left half plane)
    ! pole angle for k-th pair: pi * (2*k + order - 1) / (2*order)
    pole_real = cos(PI * dble(2*k + order - 1) / dble(2*order))
    pole_imag = sin(PI * dble(2*k + order - 1) / dble(2*order))

    ! scale to cutoff frequency
    pole_real = wc * pole_real
    pole_imag = wc * pole_imag

    ! bilinear transform: s = 2*fs*(z-1)/(z+1)
    ! For a conjugate pair with analog transfer function:
    !   H(s) = wc^2 / (s^2 - 2*Re(p)*s + |p|^2)
    ! where |p|^2 = pole_real^2 + pole_imag^2 = wc^2 (for Butterworth)
    !
    ! After bilinear transform, the digital filter coefficients are:
    !   b0 = wc^2
    !   b1 = 2*wc^2
    !   b2 = wc^2
    !   a0 = 4*fs^2 - 4*fs*Re(p) + |p|^2
    !   a1 = 2*|p|^2 - 8*fs^2
    !   a2 = 4*fs^2 + 4*fs*Re(p) + |p|^2

    b0 = wc2
    b1 = 2.0d0 * wc2
    b2 = wc2

    a0 = 4.0d0*fs*fs - 4.0d0*fs*pole_real + wc2
    a1 = 2.0d0*wc2 - 8.0d0*fs*fs
    a2 = 4.0d0*fs*fs + 4.0d0*fs*pole_real + wc2

    ! normalize so that a0 = 1
    sos(k, 1) = b0 / a0
    sos(k, 2) = b1 / a0
    sos(k, 3) = b2 / a0
    sos(k, 4) = 1.0d0
    sos(k, 5) = a1 / a0
    sos(k, 6) = a2 / a0
  enddo

  end subroutine gf_butterworth_sos

!
!----------------------------------------------------------------------
!

  subroutine gf_sosfiltfilt(x, n, sos, nsections)

! Applies zero-phase filtering using cascaded second-order sections.
! For each SOS: forward pass then backward pass (like scipy filtfilt).

  implicit none

  integer, intent(in) :: n, nsections
  double precision, intent(inout) :: x(n)
  double precision, intent(in) :: sos(nsections, 6)

  ! local variables
  integer :: isec

  do isec = 1, nsections
    ! forward pass
    call gf_sos_filter_forward(x, n, sos(isec, 1), sos(isec, 2), sos(isec, 3), &
                                      sos(isec, 5), sos(isec, 6))
    ! backward pass (reverse, filter, reverse)
    call gf_reverse_array(x, n)
    call gf_sos_filter_forward(x, n, sos(isec, 1), sos(isec, 2), sos(isec, 3), &
                                      sos(isec, 5), sos(isec, 6))
    call gf_reverse_array(x, n)
  enddo

  end subroutine gf_sosfiltfilt

!
!----------------------------------------------------------------------
!

  subroutine gf_sos_filter_forward(x, n, b0, b1, b2, a1, a2)

! Applies a single second-order section IIR filter (direct form II transposed).
! Assumes a0 = 1.0 (already normalized).
!
! Difference equation:
!   y(i) = b0*x(i) + w1
!   w1   = b1*x(i) - a1*y(i) + w2
!   w2   = b2*x(i) - a2*y(i)

  implicit none

  integer, intent(in) :: n
  double precision, intent(inout) :: x(n)
  double precision, intent(in) :: b0, b1, b2, a1, a2

  ! local variables
  integer :: i
  double precision :: w1, w2, yi

  w1 = 0.0d0
  w2 = 0.0d0

  do i = 1, n
    yi = b0 * x(i) + w1
    w1 = b1 * x(i) - a1 * yi + w2
    w2 = b2 * x(i) - a2 * yi
    x(i) = yi
  enddo

  end subroutine gf_sos_filter_forward

!
!----------------------------------------------------------------------
!

  subroutine gf_reverse_array(x, n)

! Reverses an array in place.

  implicit none

  integer, intent(in) :: n
  double precision, intent(inout) :: x(n)

  ! local variables
  integer :: i
  double precision :: tmp

  do i = 1, n/2
    tmp = x(i)
    x(i) = x(n - i + 1)
    x(n - i + 1) = tmp
  enddo

  end subroutine gf_reverse_array
