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
!---- Green function database: read GF_LOCATIONS and locate containing elements
!----

  subroutine gf_read_locations()

! Reads DATA/GF_LOCATIONS file on rank 0 and broadcasts to all ranks.
! Each line: latitude longitude depth_km
! Lines starting with '#' are comments, empty lines are skipped.

  use constants, only: IMAIN,MAX_STRING_LEN,IIN
  use specfem_par, only: myrank

  use green_function_par, only: gf_nlocations,gf_lat,gf_lon,gf_depth

  implicit none

  ! local parameters
  integer :: ier,iline,igf
  character(len=MAX_STRING_LEN) :: line
  character(len=MAX_STRING_LEN) :: gf_locations_file
  double precision :: lat,lon,depth_km

  ! counts valid lines and reads coordinates on main process
  gf_locations_file = 'DATA/GF_LOCATIONS'

  if (myrank == 0) then
    ! first pass: count valid lines
    gf_nlocations = 0
    open(unit=IIN,file=trim(gf_locations_file),status='old',action='read',iostat=ier)
    if (ier /= 0) then
      print *, 'Error: GF_DATABASE_ENABLED is .true. but could not open ',trim(gf_locations_file)
      call exit_MPI(myrank,'GF_DATABASE_ENABLED requires a DATA/GF_LOCATIONS file')
    else
      do
        read(IIN,'(a)',iostat=ier) line
        if (ier /= 0) exit
        ! skip comments and blank lines
        line = adjustl(line)
        if (len_trim(line) == 0) cycle
        if (line(1:1) == '#') cycle
        gf_nlocations = gf_nlocations + 1
      enddo
      close(IIN)
    endif

    write(IMAIN,*)
    write(IMAIN,*) '********************'
    write(IMAIN,*) ' Green function database: reading GF_LOCATIONS'
    write(IMAIN,*) '********************'
    write(IMAIN,*)
    write(IMAIN,*) 'Number of GF locations: ',gf_nlocations
    call flush_IMAIN()
  endif

  ! broadcast count
  call bcast_all_singlei(gf_nlocations)

  ! error if no locations found
  if (gf_nlocations == 0) then
    if (myrank == 0) print *, 'Error: GF_LOCATIONS file has no valid entries'
    call exit_MPI(myrank,'GF_DATABASE_ENABLED requires at least one location in DATA/GF_LOCATIONS')
  endif

  ! allocate coordinate arrays on all ranks
  allocate(gf_lat(gf_nlocations), &
           gf_lon(gf_nlocations), &
           gf_depth(gf_nlocations),stat=ier)
  if (ier /= 0) call exit_MPI(myrank,'Error allocating GF location arrays')

  ! second pass: read coordinates on main
  if (myrank == 0) then
    open(unit=IIN,file=trim(gf_locations_file),status='old',action='read',iostat=ier)
    if (ier /= 0) call exit_MPI(myrank,'Error re-opening GF_LOCATIONS file')

    igf = 0
    do
      read(IIN,'(a)',iostat=ier) line
      if (ier /= 0) exit
      line = adjustl(line)
      if (len_trim(line) == 0) cycle
      if (line(1:1) == '#') cycle

      igf = igf + 1
      read(line,*,iostat=ier) lat, lon, depth_km
      if (ier /= 0) then
        write(IMAIN,*) 'Error parsing GF_LOCATIONS line ',igf,': ',trim(line)
        call exit_MPI(myrank,'Error reading GF_LOCATIONS coordinates')
      endif

      gf_lat(igf) = lat
      gf_lon(igf) = lon
      gf_depth(igf) = depth_km
    enddo
    close(IIN)

    ! print locations
    do igf = 1,gf_nlocations
      write(IMAIN,'(a,i6,a,f10.4,a,f10.4,a,f10.2,a)') &
        '  Location #',igf, &
        ':  lat=',gf_lat(igf), &
        '  lon=',gf_lon(igf), &
        '  depth=',gf_depth(igf),' km'
    enddo
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  ! broadcast coordinates to all ranks
  call bcast_all_dp(gf_lat,gf_nlocations)
  call bcast_all_dp(gf_lon,gf_nlocations)
  call bcast_all_dp(gf_depth,gf_nlocations)

  end subroutine gf_read_locations

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gf_locate_elements()

! Locates which spectral elements contain the GF locations.
! Follows the pattern of locate_receivers.f90 but only needs ispec/islice,
! not interpolation coefficients.
! After locating, builds a unique list of tagged elements per rank using
! O(n) boolean mask deduplication.

  use constants_solver, only: &
    ELLIPTICITY_VAL,NDIM,HUGEVAL,IMAIN, &
    DEGREES_TO_RADIANS,R_UNIT_SPHERE, &
    nrec_SUBSET_MAX

  use shared_parameters, only: R_PLANET

  use specfem_par, only: &
    myrank, &
    nspec => NSPEC_CRUST_MANTLE, &
    ibathy_topo,TOPOGRAPHY

  use specfem_par, only: rspl_ellip,ellipicity_spline,ellipicity_spline2,nspl_ellip

  use green_function_par, only: &
    gf_nlocations,gf_lat,gf_lon,gf_depth, &
    gf_ispec_selected,gf_islice_selected, &
    gf_nelem_local,gf_local_elements

  implicit none

  ! local parameters
  integer :: igf,igf_in_this_subset,igf_already_done,ier
  integer :: ngf_SUBSET_current_size
  double precision :: lat,lon,depth_km,depth_m
  double precision :: theta,phi
  double precision :: sint,cost,sinp,cosp
  double precision :: r0,r_target,elevation
  double precision :: x_target,y_target,z_target
  double precision :: x,y,z
  double precision :: xi,eta,gamma
  double precision :: distmin_not_squared
  integer :: ispec_selected,ispec,j

  ! subset arrays for locate_MPI_slice
  integer, allocatable, dimension(:) :: ispec_selected_subset
  double precision, allocatable, dimension(:,:) :: xyz_found_subset
  double precision, allocatable, dimension(:) :: xi_subset,eta_subset,gamma_subset
  double precision, allocatable, dimension(:) :: final_distance_subset

  ! full arrays for MPI gathering
  double precision, allocatable, dimension(:,:) :: xyz_found
  double precision, allocatable, dimension(:) :: xi_all,eta_all,gamma_all
  double precision, allocatable, dimension(:) :: final_distance

  ! deduplication mask
  logical, allocatable, dimension(:) :: ispec_tagged

  logical, parameter :: POINT_CAN_BE_BURIED = .true.

  ! nothing to do
  if (gf_nlocations == 0) then
    gf_nelem_local = 0
    return
  endif

  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) '********************'
    write(IMAIN,*) ' Green function database: locating elements'
    write(IMAIN,*) '********************'
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  ! allocate result arrays
  allocate(gf_ispec_selected(gf_nlocations), &
           gf_islice_selected(gf_nlocations), &
           xyz_found(NDIM,gf_nlocations), &
           xi_all(gf_nlocations), &
           eta_all(gf_nlocations), &
           gamma_all(gf_nlocations), &
           final_distance(gf_nlocations),stat=ier)
  if (ier /= 0) call exit_MPI(myrank,'Error allocating GF element location arrays')

  gf_ispec_selected(:) = 0
  gf_islice_selected(:) = -1
  xyz_found(:,:) = 0.d0
  xi_all(:) = 0.d0
  eta_all(:) = 0.d0
  gamma_all(:) = 0.d0
  final_distance(:) = HUGEVAL

  ! loop over subsets (same pattern as locate_receivers)
  do igf_already_done = 0, gf_nlocations, nrec_SUBSET_MAX

    ngf_SUBSET_current_size = min(nrec_SUBSET_MAX, gf_nlocations - igf_already_done)
    if (ngf_SUBSET_current_size <= 0) exit

    ! allocate subset arrays
    allocate(ispec_selected_subset(ngf_SUBSET_current_size), &
             xi_subset(ngf_SUBSET_current_size), &
             eta_subset(ngf_SUBSET_current_size), &
             gamma_subset(ngf_SUBSET_current_size), &
             xyz_found_subset(NDIM,ngf_SUBSET_current_size), &
             final_distance_subset(ngf_SUBSET_current_size),stat=ier)
    if (ier /= 0) call exit_MPI(myrank,'Error allocating GF subset arrays')

    ispec_selected_subset(:) = 0
    final_distance_subset(:) = HUGEVAL

    ! locate points in this subset
    do igf_in_this_subset = 1, ngf_SUBSET_current_size
      igf = igf_in_this_subset + igf_already_done

      lat = gf_lat(igf)
      lon = gf_lon(igf)
      depth_km = gf_depth(igf)

      ! depth in meters
      depth_m = depth_km * 1000.d0

      ! normalize longitude to [0, 360]
      if (lon < 0.d0) lon = lon + 360.d0
      if (lon > 360.d0) lon = lon - 360.d0

      ! geographic latitude -> geocentric colatitude (radians)
      call lat_2_geocentric_colat_dble(lat,theta,ELLIPTICITY_VAL)

      ! longitude to radians
      phi = lon * DEGREES_TO_RADIANS

      ! reduce to valid ranges
      call reduce(theta,phi)

      sint = sin(theta)
      cost = cos(theta)
      sinp = sin(phi)
      cosp = cos(phi)

      ! start with unit sphere radius
      r0 = R_UNIT_SPHERE

      ! topography
      if (TOPOGRAPHY) then
        call get_topo_bathy(lat,lon,elevation,ibathy_topo)
        r0 = r0 + elevation / R_PLANET
      endif

      ! ellipticity
      if (ELLIPTICITY_VAL) then
        call add_ellipticity_rtheta(r0,theta,nspl_ellip,rspl_ellip, &
                                    ellipicity_spline,ellipicity_spline2)
      endif

      ! subtract depth (in meters)
      r0 = r0 - depth_m / R_PLANET
      r_target = r0

      ! Cartesian coordinates
      x_target = r_target * sint * cosp
      y_target = r_target * sint * sinp
      z_target = r_target * cost

      ! locate element
      call locate_point(x_target,y_target,z_target,lat,lon, &
                        ispec_selected,xi,eta,gamma, &
                        x,y,z,distmin_not_squared,POINT_CAN_BE_BURIED)

      ! store results
      ispec_selected_subset(igf_in_this_subset) = ispec_selected
      xi_subset(igf_in_this_subset) = xi
      eta_subset(igf_in_this_subset) = eta
      gamma_subset(igf_in_this_subset) = gamma
      xyz_found_subset(1,igf_in_this_subset) = x
      xyz_found_subset(2,igf_in_this_subset) = y
      xyz_found_subset(3,igf_in_this_subset) = z
      final_distance_subset(igf_in_this_subset) = distmin_not_squared

    enddo

    ! gather across MPI ranks to find best location
    call locate_MPI_slice(ngf_SUBSET_current_size,igf_already_done, &
                          ispec_selected_subset, &
                          xyz_found_subset, &
                          xi_subset,eta_subset,gamma_subset, &
                          final_distance_subset, &
                          gf_nlocations,gf_ispec_selected,gf_islice_selected, &
                          xyz_found, &
                          xi_all,eta_all,gamma_all, &
                          final_distance)

    deallocate(ispec_selected_subset)
    deallocate(xi_subset,eta_subset,gamma_subset)
    deallocate(final_distance_subset)
    deallocate(xyz_found_subset)

  enddo

  ! user output on main
  if (myrank == 0) then
    do igf = 1,gf_nlocations
      write(IMAIN,'(a,i6,a,f10.4,a,f10.4,a,f10.2,a,i6,a,i8,a,f8.3,a)') &
        '  GF location #',igf, &
        ':  lat=',gf_lat(igf), &
        '  lon=',gf_lon(igf), &
        '  depth=',gf_depth(igf), &
        ' km  -> slice ',gf_islice_selected(igf), &
        '  element ',gf_ispec_selected(igf), &
        '  (dist=',sngl(final_distance(igf)),' km)'
    enddo
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  ! broadcast results to all ranks
  call bcast_all_i(gf_ispec_selected,gf_nlocations)
  call bcast_all_i(gf_islice_selected,gf_nlocations)

  ! free temporary arrays
  deallocate(xyz_found)
  deallocate(xi_all,eta_all,gamma_all)
  deallocate(final_distance)

  ! O(n) deduplication: build unique element list per rank using boolean mask
  allocate(ispec_tagged(nspec),stat=ier)
  if (ier /= 0) call exit_MPI(myrank,'Error allocating ispec_tagged')
  ispec_tagged(:) = .false.

  do igf = 1, gf_nlocations
    if (gf_islice_selected(igf) == myrank) then
      ispec_tagged(gf_ispec_selected(igf)) = .true.
    endif
  enddo

  ! count unique elements
  gf_nelem_local = count(ispec_tagged)

  ! pack into compact array
  if (gf_nelem_local > 0) then
    allocate(gf_local_elements(gf_nelem_local),stat=ier)
    if (ier /= 0) call exit_MPI(myrank,'Error allocating gf_local_elements')

    j = 0
    do ispec = 1, nspec
      if (ispec_tagged(ispec)) then
        j = j + 1
        gf_local_elements(j) = ispec
      endif
    enddo
  endif

  deallocate(ispec_tagged)

  ! report per-rank element counts
  if (myrank == 0) then
    write(IMAIN,*) 'Green function database: element tagging summary'
    write(IMAIN,*) '  Rank 0 has ',gf_nelem_local,' unique tagged elements'
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  end subroutine gf_locate_elements
