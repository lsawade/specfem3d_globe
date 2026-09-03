program test_gf_stf

! Unit test for Green function STF: Butterworth filter + Gaussian STF.
!
! Tests:
! 1. Butterworth SOS coefficient computation
! 2. Filter application (lowpass removes high frequencies)
! 3. Full STF generation: Gaussian + Butterworth lowpass
! 4. STF properties: peak location, area, frequency content
! 5. Write output for Python cross-validation

  implicit none

  double precision, parameter :: PI = 3.141592653589793d0

  ! test parameters
  double precision, parameter :: dt = 0.05d0        ! sampling interval (s)
  integer, parameter :: nstep = 2000                 ! number of timesteps
  double precision, parameter :: hdur = 0.25d0       ! Gaussian half-duration (s)
  double precision, parameter :: f_cutoff = 2.5d0    ! lowpass cutoff (Hz)
  double precision, parameter :: fs = 1.0d0 / dt     ! sampling rate (Hz)
  integer, parameter :: filter_order = 4
  integer, parameter :: nsections = filter_order / 2

  ! arrays
  double precision :: sos(nsections, 6)
  double precision :: stf(nstep), stf_filtered(nstep)
  double precision :: t(nstep)

  ! test variables
  double precision :: peak_val, peak_time, area
  double precision :: timeval
  integer :: i, ipeak
  integer :: nfail

  ! for frequency content check
  integer :: nfft, k
  double precision :: freq, amp_sum_above, amp_sum_below
  double precision, allocatable :: fft_real(:), fft_imag(:)

  ! output file
  integer, parameter :: IOUT = 20

  nfail = 0

  print *, 'program: test_gf_stf'
  print *

  !--------------------------------------------------------------------
  ! Test 1: Butterworth SOS coefficients
  !--------------------------------------------------------------------
  print *, 'Test 1: Butterworth SOS coefficients'

  call gf_butterworth_sos(filter_order, f_cutoff, fs, nsections, sos)

  ! verify basic properties
  do i = 1, nsections
    print *, '  Section ', i
    print *, '    b = ', sos(i,1), sos(i,2), sos(i,3)
    print *, '    a = ', sos(i,4), sos(i,5), sos(i,6)

    ! a0 must be 1.0
    if (abs(sos(i,4) - 1.0d0) > 1.0d-12) then
      print *, '  FAIL: a0 is not 1.0'
      nfail = nfail + 1
    endif

    ! b0, b2 must be positive (lowpass)
    if (sos(i,1) <= 0.0d0 .or. sos(i,3) <= 0.0d0) then
      print *, '  FAIL: b0 or b2 is not positive'
      nfail = nfail + 1
    endif

    ! b1 = 2*b0 for lowpass Butterworth
    if (abs(sos(i,2) - 2.0d0*sos(i,1)) > 1.0d-12) then
      print *, '  FAIL: b1 != 2*b0'
      nfail = nfail + 1
    endif
  enddo

  ! DC gain check: product of section DC gains should be 1.0
  call check_dc_gain(sos, nsections, nfail)

  print *, '  done.'
  print *

  !--------------------------------------------------------------------
  ! Test 2: Lowpass filter removes high frequencies
  !--------------------------------------------------------------------
  print *, 'Test 2: Lowpass filter effect on impulse'

  ! create an impulse signal
  stf(:) = 0.0d0
  stf(nstep/2) = 1.0d0

  ! apply filter
  stf_filtered(:) = stf(:)
  call gf_sosfiltfilt(stf_filtered, nstep, sos, nsections)

  ! the filtered impulse should be spread out and have lower peak
  if (maxval(abs(stf_filtered)) >= 1.0d0) then
    print *, '  FAIL: filtered impulse peak should be less than 1.0, got ', maxval(abs(stf_filtered))
    nfail = nfail + 1
  else
    print *, '  filtered impulse peak: ', maxval(abs(stf_filtered)), ' (< 1.0, OK)'
  endif

  ! the filtered signal should be symmetric (zero-phase)
  call check_symmetry(stf_filtered, nstep, nstep/2, nfail)

  print *, '  done.'
  print *

  !--------------------------------------------------------------------
  ! Test 3: Full Gaussian + Butterworth STF
  !--------------------------------------------------------------------
  print *, 'Test 3: Gaussian + Butterworth STF'

  ! generate time array and Gaussian STF
  ! center the Gaussian at t = 5*hdur
  do i = 1, nstep
    t(i) = dble(i-1) * dt
    timeval = t(i) - 5.0d0 * hdur
    stf(i) = exp(-(timeval**2) / (hdur**2)) / (sqrt(PI) * hdur)
  enddo

  ! make a copy and filter it
  stf_filtered(:) = stf(:)
  call gf_sosfiltfilt(stf_filtered, nstep, sos, nsections)

  ! check peak location (should be near 5*hdur)
  ipeak = 1
  peak_val = stf_filtered(1)
  do i = 2, nstep
    if (stf_filtered(i) > peak_val) then
      peak_val = stf_filtered(i)
      ipeak = i
    endif
  enddo
  peak_time = t(ipeak)

  print *, '  peak value:    ', peak_val
  print *, '  peak time:     ', peak_time, ' s (expected ~', 5.0d0*hdur, ' s)'

  if (abs(peak_time - 5.0d0*hdur) > 2.0d0*dt) then
    print *, '  FAIL: peak time shifted too much'
    nfail = nfail + 1
  endif

  ! check area (trapezoidal integration, should be ~1.0 for normalized Gaussian)
  area = 0.0d0
  do i = 1, nstep - 1
    area = area + 0.5d0 * (stf_filtered(i) + stf_filtered(i+1)) * dt
  enddo
  print *, '  area (integral): ', area, ' (expected ~1.0)'

  if (abs(area - 1.0d0) > 0.05d0) then
    print *, '  FAIL: area deviates from 1.0 by more than 5%'
    nfail = nfail + 1
  endif

  ! check peak is positive
  if (peak_val <= 0.0d0) then
    print *, '  FAIL: peak value should be positive'
    nfail = nfail + 1
  endif

  print *, '  done.'
  print *

  !--------------------------------------------------------------------
  ! Test 4: Frequency content check (DFT)
  !--------------------------------------------------------------------
  print *, 'Test 4: Frequency content of filtered STF'

  nfft = nstep
  allocate(fft_real(nfft), fft_imag(nfft))

  ! compute DFT of filtered STF
  call compute_dft(stf_filtered, nstep, fft_real, fft_imag, nfft)

  ! compute total energy below and above cutoff
  amp_sum_below = 0.0d0
  amp_sum_above = 0.0d0
  do k = 1, nfft/2
    freq = dble(k-1) * fs / dble(nfft)
    if (freq > 0.0d0 .and. freq <= f_cutoff) then
      amp_sum_below = amp_sum_below + fft_real(k)**2 + fft_imag(k)**2
    else if (freq > f_cutoff) then
      amp_sum_above = amp_sum_above + fft_real(k)**2 + fft_imag(k)**2
    endif
  enddo

  print *, '  energy below cutoff: ', amp_sum_below
  print *, '  energy above cutoff: ', amp_sum_above

  ! energy above cutoff should be negligible compared to below
  if (amp_sum_below > 0.0d0) then
    if (amp_sum_above / amp_sum_below > 1.0d-4) then
      print *, '  FAIL: too much energy above cutoff (ratio = ', amp_sum_above/amp_sum_below, ')'
      nfail = nfail + 1
    else
      print *, '  energy ratio above/below cutoff: ', amp_sum_above/amp_sum_below, ' (< 1e-4, OK)'
    endif
  endif

  deallocate(fft_real, fft_imag)

  print *, '  done.'
  print *

  !--------------------------------------------------------------------
  ! Test 5: Write output for Python cross-validation
  !--------------------------------------------------------------------
  print *, 'Writing STF to test_stf_output.txt'

  ! regenerate unfiltered for comparison
  do i = 1, nstep
    timeval = t(i) - 5.0d0 * hdur
    stf(i) = exp(-(timeval**2) / (hdur**2)) / (sqrt(PI) * hdur)
  enddo

  open(unit=IOUT, file='test_stf_output.txt', status='unknown')
  write(IOUT, '(a)') '# Green function STF unit test output'
  write(IOUT, '(a,ES20.10)') '# hdur = ', hdur
  write(IOUT, '(a,ES20.10)') '# f_cutoff = ', f_cutoff
  write(IOUT, '(a,ES20.10)') '# dt = ', dt
  write(IOUT, '(a,I10)') '# nstep = ', nstep
  write(IOUT, '(a,ES20.10)') '# t_center = ', 5.0d0*hdur
  write(IOUT, '(a,I10)') '# filter_order = ', filter_order
  write(IOUT, '(a)') '# time(s)  stf_unfiltered  stf_filtered'
  do i = 1, nstep
    write(IOUT, '(3ES20.10)') t(i), stf(i), stf_filtered(i)
  enddo
  close(IOUT)

  print *, '  done.'
  print *

  !--------------------------------------------------------------------
  ! Summary
  !--------------------------------------------------------------------
  if (nfail == 0) then
    print *, '--- all tests passed ---'
  else
    print *, '--- FAILED: ', nfail, ' test(s) ---'
    stop 1
  endif

end program test_gf_stf

!======================================================================
! Helper subroutines for the test
!======================================================================

subroutine check_dc_gain(sos, nsections, nfail)

  implicit none
  integer, intent(in) :: nsections
  double precision, intent(in) :: sos(nsections, 6)
  integer, intent(inout) :: nfail

  double precision :: dc_gain, section_gain
  integer :: i

  dc_gain = 1.0d0
  do i = 1, nsections
    section_gain = (sos(i,1) + sos(i,2) + sos(i,3)) / &
                   (sos(i,4) + sos(i,5) + sos(i,6))
    dc_gain = dc_gain * section_gain
  enddo

  print *, '  total DC gain: ', dc_gain, ' (expected 1.0)'

  if (abs(dc_gain - 1.0d0) > 1.0d-10) then
    print *, '  FAIL: DC gain is not 1.0'
    nfail = nfail + 1
  endif

end subroutine check_dc_gain

!----------------------------------------------------------------------

subroutine check_symmetry(x, n, icenter, nfail)

  implicit none
  integer, intent(in) :: n, icenter
  double precision, intent(in) :: x(n)
  integer, intent(inout) :: nfail

  integer :: i, ncheck
  double precision :: max_asym

  ! check symmetry around icenter
  ncheck = min(icenter - 1, n - icenter)
  ncheck = min(ncheck, 100)

  max_asym = 0.0d0
  do i = 1, ncheck
    max_asym = max(max_asym, abs(x(icenter + i) - x(icenter - i)))
  enddo

  print *, '  max asymmetry around center: ', max_asym

  if (max_asym > 1.0d-12) then
    print *, '  FAIL: filtered impulse is not symmetric (zero-phase violated)'
    nfail = nfail + 1
  endif

end subroutine check_symmetry

!----------------------------------------------------------------------

subroutine compute_dft(x, n, dft_real, dft_imag, nfft)

! Simple DFT computation (not FFT, but sufficient for testing)

  implicit none

  double precision, parameter :: PI = 3.141592653589793d0

  integer, intent(in) :: n, nfft
  double precision, intent(in) :: x(n)
  double precision, intent(out) :: dft_real(nfft), dft_imag(nfft)

  integer :: k, j
  double precision :: angle

  do k = 1, nfft
    dft_real(k) = 0.0d0
    dft_imag(k) = 0.0d0
    do j = 1, n
      angle = 2.0d0 * PI * dble(k-1) * dble(j-1) / dble(nfft)
      dft_real(k) = dft_real(k) + x(j) * cos(angle)
      dft_imag(k) = dft_imag(k) - x(j) * sin(angle)
    enddo
  enddo

end subroutine compute_dft

!======================================================================
! Filter subroutines (identical to src/specfem3D/green_function_stf.F90)
!======================================================================

subroutine gf_butterworth_sos(order, fc, fs, nsections, sos)

! Computes second-order section (SOS) coefficients for a Butterworth
! lowpass filter using the bilinear transform.

  implicit none

  double precision, parameter :: PI = 3.141592653589793238462643383279502884197d0

  integer, intent(in) :: order, nsections
  double precision, intent(in) :: fc, fs
  double precision, intent(out) :: sos(nsections, 6)

  integer :: k
  double precision :: wc, wc2
  double precision :: pole_real, pole_imag
  double precision :: b0, b1, b2, a0, a1, a2

  wc = 2.0d0 * fs * tan(PI * fc / fs)
  wc2 = wc * wc

  do k = 1, nsections
    pole_real = cos(PI * dble(2*k + order - 1) / dble(2*order))
    pole_imag = sin(PI * dble(2*k + order - 1) / dble(2*order))

    pole_real = wc * pole_real
    pole_imag = wc * pole_imag

    b0 = wc2
    b1 = 2.0d0 * wc2
    b2 = wc2

    a0 = 4.0d0*fs*fs - 4.0d0*fs*pole_real + wc2
    a1 = 2.0d0*wc2 - 8.0d0*fs*fs
    a2 = 4.0d0*fs*fs + 4.0d0*fs*pole_real + wc2

    sos(k, 1) = b0 / a0
    sos(k, 2) = b1 / a0
    sos(k, 3) = b2 / a0
    sos(k, 4) = 1.0d0
    sos(k, 5) = a1 / a0
    sos(k, 6) = a2 / a0
  enddo

end subroutine gf_butterworth_sos

!----------------------------------------------------------------------

subroutine gf_sosfiltfilt(x, n, sos, nsections)

! Applies zero-phase filtering using cascaded second-order sections.

  implicit none

  integer, intent(in) :: n, nsections
  double precision, intent(inout) :: x(n)
  double precision, intent(in) :: sos(nsections, 6)

  integer :: isec

  do isec = 1, nsections
    call gf_sos_filter_forward(x, n, sos(isec, 1), sos(isec, 2), sos(isec, 3), &
                                      sos(isec, 5), sos(isec, 6))
    call gf_reverse_array(x, n)
    call gf_sos_filter_forward(x, n, sos(isec, 1), sos(isec, 2), sos(isec, 3), &
                                      sos(isec, 5), sos(isec, 6))
    call gf_reverse_array(x, n)
  enddo

end subroutine gf_sosfiltfilt

!----------------------------------------------------------------------

subroutine gf_sos_filter_forward(x, n, b0, b1, b2, a1, a2)

! Single second-order section IIR filter (direct form II transposed).

  implicit none

  integer, intent(in) :: n
  double precision, intent(inout) :: x(n)
  double precision, intent(in) :: b0, b1, b2, a1, a2

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

!----------------------------------------------------------------------

subroutine gf_reverse_array(x, n)

! Reverses an array in place.

  implicit none

  integer, intent(in) :: n
  double precision, intent(inout) :: x(n)

  integer :: i
  double precision :: tmp

  do i = 1, n/2
    tmp = x(i)
    x(i) = x(n - i + 1)
    x(n - i + 1) = tmp
  enddo

end subroutine gf_reverse_array
