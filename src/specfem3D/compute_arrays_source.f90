!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  7 . 0
!          --------------------------------------------------
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

  subroutine compute_arrays_source(sourcearray, &
                                   xi_source,eta_source,gamma_source, &
                                   Mxx,Myy,Mzz,Mxy,Mxz,Myz, &
                                   xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz, &
                                   xigll,yigll,zigll)

  use constants

  implicit none

  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLY,NGLLZ) :: sourcearray

  double precision :: xi_source,eta_source,gamma_source
  double precision :: Mxx,Myy,Mzz,Mxy,Mxz,Myz

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: xix,xiy,xiz,etax,etay,etaz, &
        gammax,gammay,gammaz

  ! Gauss-Lobatto-Legendre points of integration and weights
  double precision, dimension(NGLLX) :: xigll
  double precision, dimension(NGLLY) :: yigll
  double precision, dimension(NGLLZ) :: zigll

  ! local parameters
  double precision :: xixd,xiyd,xizd,etaxd,etayd,etazd,gammaxd,gammayd,gammazd

  ! source arrays
  double precision, dimension(NDIM,NGLLX,NGLLY,NGLLZ) :: sourcearrayd
  double precision, dimension(NGLLX) :: hxis,hpxis
  double precision, dimension(NGLLY) :: hetas,hpetas
  double precision, dimension(NGLLZ) :: hgammas,hpgammas

  double precision :: hlagrange
  double precision :: dsrc_dx, dsrc_dy, dsrc_dz
  double precision :: dxis_dx, detas_dx, dgammas_dx
  double precision :: dxis_dy, detas_dy, dgammas_dy
  double precision :: dxis_dz, detas_dz, dgammas_dz

  integer :: k,l,m

! compute Lagrange polynomials at the source location
! the source does not necessarily correspond to a Gauss-Lobatto point
  call lagrange_any(xi_source,NGLLX,xigll,hxis,hpxis)
  call lagrange_any(eta_source,NGLLY,yigll,hetas,hpetas)
  call lagrange_any(gamma_source,NGLLZ,zigll,hgammas,hpgammas)

  dxis_dx = ZERO
  dxis_dy = ZERO
  dxis_dz = ZERO
  detas_dx = ZERO
  detas_dy = ZERO
  detas_dz = ZERO
  dgammas_dx = ZERO
  dgammas_dy = ZERO
  dgammas_dz = ZERO

  do m = 1,NGLLZ
     do l = 1,NGLLY
        do k = 1,NGLLX

           xixd    = dble(xix(k,l,m))
           xiyd    = dble(xiy(k,l,m))
           xizd    = dble(xiz(k,l,m))
           etaxd   = dble(etax(k,l,m))
           etayd   = dble(etay(k,l,m))
           etazd   = dble(etaz(k,l,m))
           gammaxd = dble(gammax(k,l,m))
           gammayd = dble(gammay(k,l,m))
           gammazd = dble(gammaz(k,l,m))

           hlagrange = hxis(k) * hetas(l) * hgammas(m)

           dxis_dx = dxis_dx + hlagrange * xixd
           dxis_dy = dxis_dy + hlagrange * xiyd
           dxis_dz = dxis_dz + hlagrange * xizd

           detas_dx = detas_dx + hlagrange * etaxd
           detas_dy = detas_dy + hlagrange * etayd
           detas_dz = detas_dz + hlagrange * etazd

           dgammas_dx = dgammas_dx + hlagrange * gammaxd
           dgammas_dy = dgammas_dy + hlagrange * gammayd
           dgammas_dz = dgammas_dz + hlagrange * gammazd

       enddo
     enddo
  enddo

! calculate source array
  sourcearrayd(:,:,:,:) = ZERO
  do m = 1,NGLLZ
     do l = 1,NGLLY
        do k = 1,NGLLX

           dsrc_dx = (hpxis(k)*dxis_dx)*hetas(l)*hgammas(m) + hxis(k)*(hpetas(l)*detas_dx)*hgammas(m) + &
                                                                hxis(k)*hetas(l)*(hpgammas(m)*dgammas_dx)
           dsrc_dy = (hpxis(k)*dxis_dy)*hetas(l)*hgammas(m) + hxis(k)*(hpetas(l)*detas_dy)*hgammas(m) + &
                                                                hxis(k)*hetas(l)*(hpgammas(m)*dgammas_dy)
           dsrc_dz = (hpxis(k)*dxis_dz)*hetas(l)*hgammas(m) + hxis(k)*(hpetas(l)*detas_dz)*hgammas(m) + &
                                                                hxis(k)*hetas(l)*(hpgammas(m)*dgammas_dz)

           sourcearrayd(1,k,l,m) = sourcearrayd(1,k,l,m) + (Mxx*dsrc_dx + Mxy*dsrc_dy + Mxz*dsrc_dz)
           sourcearrayd(2,k,l,m) = sourcearrayd(2,k,l,m) + (Mxy*dsrc_dx + Myy*dsrc_dy + Myz*dsrc_dz)
           sourcearrayd(3,k,l,m) = sourcearrayd(3,k,l,m) + (Mxz*dsrc_dx + Myz*dsrc_dy + Mzz*dsrc_dz)

       enddo
     enddo
  enddo

  ! distinguish between single and double precision for reals
  sourcearray(:,:,:,:) = real(sourcearrayd(:,:,:,:), kind=CUSTOM_REAL)

  end subroutine compute_arrays_source

!================================================================

  subroutine compute_arrays_source_derivative(sourcearray, &
                                   xi_source,eta_source,gamma_source, &
                                   Mxx,Myy,Mzz,Mxy,Mxz,Myz, &
                                   xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz, &
                                   xigll,yigll,zigll, &
                                   direction, theta, phi, depth)

  use constants
  implicit none

  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLY,NGLLZ) :: sourcearray

  double precision :: xi_source,eta_source,gamma_source
  double precision :: Mxx,Myy,Mzz,Mxy,Mxz,Myz

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: xix,xiy,xiz,etax,etay,etaz, &
        gammax,gammay,gammaz

  ! Gauss-Lobatto-Legendre points of integration and weights
  double precision, dimension(NGLLX) :: xigll
  double precision, dimension(NGLLY) :: yigll
  double precision, dimension(NGLLZ) :: zigll

  ! local parameters
  double precision :: xixd,xiyd,xizd,etaxd,etayd,etazd,gammaxd,gammayd,gammazd

  ! source arrays
  double precision, dimension(NDIM,NGLLX,NGLLY,NGLLZ) :: sourcearrayd
  double precision, dimension(NGLLX) :: hxis,hpxis,hppxis
  double precision, dimension(NGLLY) :: hetas,hpetas,hppetas
  double precision, dimension(NGLLZ) :: hgammas,hpgammas,hppgammas

  double precision :: hlagrange
  double precision :: dxis_dx, detas_dx, dgammas_dx
  double precision :: dxis_dy, detas_dy, dgammas_dy
  double precision :: dxis_dz, detas_dz, dgammas_dz
  double precision :: d2src_dxi2, d2src_deta2, d2src_dgamma2
  double precision :: d2src_dxideta
  double precision :: d2src_dxidgamma
  double precision :: d2src_detadgamma
  double precision :: d2src_dx2, d2src_dy2, d2src_dz2 
  double precision :: d2src_dxy, d2src_dxz, d2src_dyz
  double precision :: fx, fxx, fy, fyy, fz, fzz
  double precision :: fyx, fzx, fxy, fzy, fxz, fyz
  double precision :: fac_x, fac_y, fac_z
  double precision :: theta, phi, depth
  double precision :: sint, cost, sinp, cosp
  double precision :: grr_inv, gtt_inv, gpp_inv

  integer :: k,l,m
  integer :: direction

! compute Lagrange polynomials at the source location
! the source does not necessarily correspond to a Gauss-Lobatto point
  call lagrange_any(xi_source,NGLLX,xigll,hxis,hpxis)
  call lagrange_any(eta_source,NGLLY,yigll,hetas,hpetas)
  call lagrange_any(gamma_source,NGLLZ,zigll,hgammas,hpgammas)
  call lagrange_2prime(xi_source,NGLLX,xigll,hppxis)
  call lagrange_2prime(eta_source,NGLLY,yigll,hppetas)
  call lagrange_2prime(gamma_source,NGLLZ,zigll,hppgammas)

  sourcearrayd(:,:,:,:) = sourcearray(:,:,:,:)

  dxis_dx = ZERO
  dxis_dy = ZERO
  dxis_dz = ZERO
  detas_dx = ZERO
  detas_dy = ZERO
  detas_dz = ZERO
  dgammas_dx = ZERO
  dgammas_dy = ZERO
  dgammas_dz = ZERO

  do m = 1,NGLLZ
     do l = 1,NGLLY
        do k = 1,NGLLX

           xixd    = dble(xix(k,l,m))
           xiyd    = dble(xiy(k,l,m))
           xizd    = dble(xiz(k,l,m))
           etaxd   = dble(etax(k,l,m))
           etayd   = dble(etay(k,l,m))
           etazd   = dble(etaz(k,l,m))
           gammaxd = dble(gammax(k,l,m))
           gammayd = dble(gammay(k,l,m))
           gammazd = dble(gammaz(k,l,m))

           hlagrange = hxis(k) * hetas(l) * hgammas(m)

           dxis_dx = dxis_dx + hlagrange * xixd
           dxis_dy = dxis_dy + hlagrange * xiyd
           dxis_dz = dxis_dz + hlagrange * xizd

           detas_dx = detas_dx + hlagrange * etaxd
           detas_dy = detas_dy + hlagrange * etayd
           detas_dz = detas_dz + hlagrange * etazd

           dgammas_dx = dgammas_dx + hlagrange * gammaxd
           dgammas_dy = dgammas_dy + hlagrange * gammayd
           dgammas_dz = dgammas_dz + hlagrange * gammazd

       enddo
     enddo
  enddo

  do m = 1,NGLLZ
    do l = 1,NGLLY
      do k = 1,NGLLX

        ! Natural derivative of source location
        d2src_dxi2 = hppxis(k) * hetas(l) * hgammas(m) 
        d2src_deta2 = hxis(k) * hppetas(l) * hgammas(m) 
        d2src_dgamma2 = hxis(k) * hetas(l) * hppgammas(m) 
        d2src_dxideta = hpxis(k) * hpetas(l) * hgammas(m)
        d2src_dxidgamma = hpxis(k) * hetas(l) * hpgammas(m) 
        d2src_detadgamma = hxis(k) * hpetas(l) * hpgammas(m)

        ! Real derivatives
        d2src_dx2 = d2src_dxi2 * dxis_dx * dxis_dx + &
                    d2src_deta2 * detas_dx * detas_dx + &
                    d2src_dgamma2 * dgammas_dx * dgammas_dx + &
                    2.d0 * d2src_dxideta * detas_dx * dxis_dx + &
                    2.d0 * d2src_dxidgamma * dgammas_dx * dxis_dx + &
                    2.d0 * d2src_detadgamma * dgammas_dx * detas_dx
        d2src_dy2 = d2src_dxi2 * dxis_dy * dxis_dy + &
                    d2src_deta2 * detas_dy * detas_dy + &
                    d2src_dgamma2 * dgammas_dy * dgammas_dy + &
                    2.d0 * d2src_dxideta * detas_dy * dxis_dy + &
                    2.d0 * d2src_dxidgamma * dgammas_dy * dxis_dy + &
                    2.d0 * d2src_detadgamma * dgammas_dy * detas_dy
        d2src_dz2 = d2src_dxi2 * dxis_dz * dxis_dz + &
                    d2src_deta2 * detas_dz * detas_dz + &
                    d2src_dgamma2 * dgammas_dz * dgammas_dz + &
                    2.d0 * d2src_dxideta * detas_dz * dxis_dz + &
                    2.d0 * d2src_dxidgamma * dgammas_dz * dxis_dz + &
                    2.d0 * d2src_detadgamma * dgammas_dz * detas_dz
        d2src_dxy = d2src_dxi2 * dxis_dx * dxis_dy + &
                    d2src_deta2 * detas_dx * detas_dy + &
                    d2src_dgamma2 * dgammas_dx * dgammas_dy + &
                    d2src_dxideta * detas_dy * dxis_dx + &
                    d2src_dxideta * detas_dx * dxis_dy + &
                    d2src_dxidgamma * dgammas_dy * dxis_dx + &
                    d2src_dxidgamma * dgammas_dx * dxis_dy + &
                    d2src_detadgamma * dgammas_dy * detas_dx + &
                    d2src_detadgamma * dgammas_dx * detas_dy
        d2src_dxz = d2src_dxi2 * dxis_dx * dxis_dz + &
                    d2src_deta2 * detas_dx * detas_dz + &
                    d2src_dgamma2 * dgammas_dx * dgammas_dz + &
                    d2src_dxideta * detas_dz * dxis_dx + &
                    d2src_dxideta * detas_dx * dxis_dz + &
                    d2src_dxidgamma * dgammas_dz * dxis_dx + &
                    d2src_dxidgamma * dgammas_dx * dxis_dz + &
                    d2src_detadgamma * dgammas_dz * detas_dx + &
                    d2src_detadgamma * dgammas_dx * detas_dz
        d2src_dyz = d2src_dxi2 * dxis_dy * dxis_dz + &
                    d2src_deta2 * detas_dy * detas_dz + &
                    d2src_dgamma2 * dgammas_dy * dgammas_dz + &
                    d2src_dxideta * detas_dz * dxis_dy + &
                    d2src_dxideta * detas_dy * dxis_dz + &
                    d2src_dxidgamma * dgammas_dz * dxis_dy + &
                    d2src_dxidgamma * dgammas_dy * dxis_dz + &
                    d2src_detadgamma * dgammas_dz * detas_dy + &
                    d2src_detadgamma * dgammas_dy * detas_dz
        
        if (m == 3 .AND. l == 3 .AND. k == 3) then
             print *, "Hessian"
             print *, d2src_dx2, d2src_dy2, d2src_dz2, d2src_dxy, d2src_dxz, d2src_dyz
        endif

        
        ! With respect to x
        fxx = (Mxx * d2src_dx2 + Mxy * d2src_dxy + Mxz * d2src_dxz)
        fyx = (Mxy * d2src_dx2 + Myy * d2src_dxy + Myz * d2src_dxz)
        fzx = (Mxz * d2src_dx2 + Myz * d2src_dxy + Mzz * d2src_dxz)
          
        ! With respect to y
        fxy = (Mxx * d2src_dxy + Mxy * d2src_dy2 + Mxz * d2src_dyz)
        fyy = (Mxy * d2src_dxy + Myy * d2src_dy2 + Myz * d2src_dyz)
        fzy = (Mxz * d2src_dxy + Myz * d2src_dy2 + Mzz * d2src_dyz)

        ! With respect to z
        fxz = (Mxx * d2src_dxz + Mxy * d2src_dyz + Mxz * d2src_dz2)
        fyz = (Mxy * d2src_dxz + Myy * d2src_dyz + Myz * d2src_dz2)
        fzz = (Mxz * d2src_dxz + Myz * d2src_dyz + Mzz * d2src_dz2)

        ! Compute rotation factors
        ! note cos(theta) = sin(90 - theta)
        ! This is important since we want the derivative with respect to 
        ! the latitude
        if (m == 1 .AND. l == 1 .AND. k == 1) then
             print *, 'depth: ', depth
             print *, 'phi:   ', phi
             print *, 'theta: ', theta
        endif
        !dthetadlambda = 
        ! thetaprime = 
        sint = sin(theta)
        cost = cos(theta)
        sinp = sin(phi)
        cosp = cos(phi)
        grr_inv = ONE
        gtt_inv = ONE   
        gpp_inv = ONE / sint
        if (direction == 1) then
            fac_x = -1.d0 * sint * cosp / EARTH_R_KM
            fac_y = -1.d0 * sint * sinp / EARTH_R_KM
            fac_z = -1.d0 * cost / EARTH_R_KM
        else if (direction == 2) then
            fac_x = -1.d0 * cost * cosp * depth * DEGREES_TO_RADIANS
            fac_y = -1.d0 * cost * sinp * depth * DEGREES_TO_RADIANS
            fac_z =  1.d0 * sint * depth * DEGREES_TO_RADIANS
        else if (direction == 3) then
            fac_x = -1.d0 * sint * sinp * depth * DEGREES_TO_RADIANS
            fac_y = sint * cosp  * depth * DEGREES_TO_RADIANS
            fac_z = 0.d0 * depth * DEGREES_TO_RADIANS
        else 
            stop "Wrong direction. Should 1 for depth, 2 for lat, 3 for lon."
        endif
        
        ! Rotate
        fx = (fxx * fac_x + fxy * fac_y + fxz * fac_z)
        fy = (fyx * fac_x + fyy * fac_y + fyz * fac_z)
        fz = (fzx * fac_x + fzy * fac_y + fzz * fac_z)

        ! Add to sourcearray
        sourcearrayd(1,k,l,m) = sourcearrayd(1,k,l,m) + fx
        sourcearrayd(2,k,l,m) = sourcearrayd(2,k,l,m) + fy
        sourcearrayd(3,k,l,m) = sourcearrayd(3,k,l,m) + fz

      enddo
    enddo
  enddo

  ! distinguish between single and double precision for reals
  sourcearray(:,:,:,:) = real(sourcearrayd(:,:,:,:), kind=CUSTOM_REAL)

  end subroutine compute_arrays_source_derivative

! =============================================================================

  subroutine compute_arrays_source_derivative2(sourcearray, &
    xi_source,eta_source,gamma_source, &
    Mxx,Myy,Mzz,Mxy,Mxz,Myz, &
    xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz, &
    xigll,yigll,zigll, &
    direction, theta, phi, depth)

  use constants

  implicit none

  double precision, external :: lagrange_deriv_GLL
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLY,NGLLZ) :: sourcearray

  double precision :: xi_source,eta_source,gamma_source
  double precision :: Mxx,Myy,Mzz,Mxy,Mxz,Myz

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: xix,xiy,xiz,etax,etay,etaz, &
  gammax,gammay,gammaz

  ! Gauss-Lobatto-Legendre points of integration and weights
  double precision, dimension(NGLLX) :: xigll
  double precision, dimension(NGLLY) :: yigll
  double precision, dimension(NGLLZ) :: zigll

  ! local parameters
  double precision :: xixd,xiyd,xizd,etaxd,etayd,etazd,gammaxd,gammayd,gammazd

  ! source arrays
  double precision, dimension(NDIM,NGLLX,NGLLY,NGLLZ) :: sourcearrayd
  double precision, dimension(NGLLX) :: hxis,hpxis
  double precision, dimension(NGLLY) :: hetas,hpetas
  double precision, dimension(NGLLZ) :: hgammas,hpgammas

  ! GLL derivative arrays
  double precision, dimension(NGLLX, NGLLX) :: lagx, dlagx
  double precision, dimension(NGLLY, NGLLY) :: lagy, dlagy
  double precision, dimension(NGLLZ, NGLLZ) :: lagz, dlagz

  ! grad array 
  double precision, dimension(3, NGLLX, NGLLY, NGLLZ) :: grad

  ! single values 
  double precision :: hlagrange
  double precision :: dsrc_dx, dsrc_dy, dsrc_dz
  double precision :: dxis_dx, detas_dx, dgammas_dx
  double precision :: dxis_dy, detas_dy, dgammas_dy
  double precision :: dxis_dz, detas_dz, dgammas_dz
  double precision :: dxs_dxsi, dxs_deta, dxs_dgamma
  double precision :: dys_dxsi, dys_deta, dys_dgamma
  double precision :: dzs_dxsi, dzs_deta, dzs_dgamma
  double precision :: d2src_dx2, d2src_dy2, d2src_dz2 
  double precision :: d2src_dxy, d2src_dxz, d2src_dyz  
  double precision :: fx, fxx, fy, fyy, fz, fzz
  double precision :: fyx, fzx, fxy, fzy, fxz, fyz
  double precision :: fac_x, fac_y, fac_z
  double precision :: theta, phi, depth
  double precision :: sint, cost, sinp, cosp
  double precision :: grr_inv, gtt_inv, gpp_inv

  integer :: k,l,m, n, o, p, j
  integer :: i1, i2, k1, k2, j1, j2
  integer :: direction

  ! compute Lagrange polynomials at the source location
  ! the source does not necessarily correspond to a Gauss-Lobatto point
  call lagrange_any(xi_source,NGLLX,xigll,hxis,hpxis)
  call lagrange_any(eta_source,NGLLY,yigll,hetas,hpetas)
  call lagrange_any(gamma_source,NGLLZ,zigll,hgammas,hpgammas)

  ! calculate derivatives of the Lagrange polynomials
  ! and precalculate some products in double precision
  ! hprime(i,j) = h'_j(xigll_i) by definition of the derivation matrix
  do i1 = 1,NGLLX
    do i2 = 1,NGLLX
      dlagx(i2,i1) = real(lagrange_deriv_GLL(i1-1,i2-1,xigll,NGLLX), kind=CUSTOM_REAL)
    enddo
  enddo

  do j1 = 1,NGLLY
    do j2 = 1,NGLLY
      dlagy(j2,j1) = real(lagrange_deriv_GLL(j1-1,j2-1,yigll,NGLLY), kind=CUSTOM_REAL)
    enddo
  enddo

  do k1 = 1,NGLLZ
    do k2 = 1,NGLLZ
      dlagz(k2,k1) = real(lagrange_deriv_GLL(k1-1,k2-1,zigll,NGLLZ), kind=CUSTOM_REAL)
    enddo
  enddo

  
  dxis_dx = ZERO
  dxis_dy = ZERO
  dxis_dz = ZERO
  detas_dx = ZERO
  detas_dy = ZERO
  detas_dz = ZERO
  dgammas_dx = ZERO
  dgammas_dy = ZERO
  dgammas_dz = ZERO

  do m = 1,NGLLZ
    do l = 1,NGLLY
      do k = 1,NGLLX

        xixd    = dble(xix(k,l,m))
        xiyd    = dble(xiy(k,l,m))
        xizd    = dble(xiz(k,l,m))
        etaxd   = dble(etax(k,l,m))
        etayd   = dble(etay(k,l,m))
        etazd   = dble(etaz(k,l,m))
        gammaxd = dble(gammax(k,l,m))
        gammayd = dble(gammay(k,l,m))
        gammazd = dble(gammaz(k,l,m))

        hlagrange = hxis(k) * hetas(l) * hgammas(m)

        dxis_dx = dxis_dx + hlagrange * xixd
        dxis_dy = dxis_dy + hlagrange * xiyd
        dxis_dz = dxis_dz + hlagrange * xizd

        detas_dx = detas_dx + hlagrange * etaxd
        detas_dy = detas_dy + hlagrange * etayd
        detas_dz = detas_dz + hlagrange * etazd

        dgammas_dx = dgammas_dx + hlagrange * gammaxd
        dgammas_dy = dgammas_dy + hlagrange * gammayd
        dgammas_dz = dgammas_dz + hlagrange * gammazd

      enddo
    enddo
  enddo

  ! Differentiate with respect to source location
  sourcearrayd(:,:,:,:) = ZERO
  grad(:,:,:,:) = ZERO

  do m = 1,NGLLZ
    do l = 1,NGLLY
      do k = 1,NGLLX

        dsrc_dx = (hpxis(k)*dxis_dx)*hetas(l)*hgammas(m) + hxis(k)*(hpetas(l)*detas_dx)*hgammas(m) + &
                                        hxis(k)*hetas(l)*(hpgammas(m)*dgammas_dx)
        dsrc_dy = (hpxis(k)*dxis_dy)*hetas(l)*hgammas(m) + hxis(k)*(hpetas(l)*detas_dy)*hgammas(m) + &
                                        hxis(k)*hetas(l)*(hpgammas(m)*dgammas_dy)
        dsrc_dz = (hpxis(k)*dxis_dz)*hetas(l)*hgammas(m) + hxis(k)*(hpetas(l)*detas_dz)*hgammas(m) + &
                                        hxis(k)*hetas(l)*(hpgammas(m)*dgammas_dz)

        ! Up until here these are the same steps as computing the moment 
        ! source. But instead of using the moment tensor here now, we save  
        ! the gradient with respect to the source at all gll locations, so that 
        ! this gradient can be numerically differentiated again.
        grad(1,k,l,m) = dsrc_dx
        grad(2,k,l,m) = dsrc_dy
        grad(3,k,l,m) = dsrc_dz

        ! print *, "Gradient", k, l, m
        ! print *, grad(1,k,l,m)
        ! print *, grad(2,k,l,m)
        ! print *, grad(3,k,l,m)
        
      enddo
    enddo
  enddo


  do m = 1,NGLLZ
    do l = 1,NGLLY
      do k = 1,NGLLX


        ! Differentiate the save gradient again at GLL locations
        ! This computation is very similar to expression A2 in the Appendix
        ! of Komatitsch 1999. We differentiate a multivariate vecotr valued
        ! function at the GLL nodes.
        dxs_dxsi = ZERO
        dxs_deta = ZERO
        dxs_dgamma = ZERO
        dys_dxsi = ZERO
        dys_deta = ZERO
        dys_dgamma = ZERO
        dzs_dxsi = ZERO
        dzs_deta = ZERO
        dzs_dgamma = ZERO

        do n = 1,NGLLX
          dxs_dxsi = dxs_dxsi + grad(1, n, l, m) * dlagx(n,k)
          dys_dxsi = dys_dxsi + grad(2, n, l, m) * dlagx(n,k)
          dzs_dxsi = dzs_dxsi + grad(3, n, l, m) * dlagx(n,k)
        enddo

        do o = 1,NGLLY
          dxs_deta = dxs_deta + grad(1, k, o, m) * dlagy(o,l)
          dys_deta = dys_deta + grad(2, k, o, m) * dlagy(o,l)
          dzs_deta = dzs_deta + grad(3, k, o, m) * dlagy(o,l)
        enddo

        do p = 1,NGLLZ
          dxs_dgamma = dxs_dgamma + grad(1, k, l, p) * dlagz(p,m)
          dys_dgamma = dys_dgamma + grad(2, k, l, p) * dlagz(p,m)
          dzs_dgamma = dzs_dgamma + grad(3, k, l, p) * dlagz(p,m)
        enddo

        if (m == 1 .AND. l == 1 .AND. k == 1) then
           ! print *, 'Local derivative'
           ! print *, 'x: ', dxs_dxsi, dxs_deta, dxs_dgamma
           ! print *, 'y: ', dys_dxsi, dys_deta, dys_dgamma
           ! print *, 'z: ', dzs_dxsi, dzs_deta, dzs_dgamma
        endif
        
        ! Compute full expressions at gll node (multiply with Jacobian)
        d2src_dx2 = dxs_dxsi * dble(xix(k,l,m)) + dxs_deta * dble(etax(k,l,m)) + dxs_dgamma * dble(gammax(k,l,m))
        d2src_dy2 = dys_dxsi * dble(xiy(k,l,m)) + dys_deta * dble(etay(k,l,m)) + dys_dgamma * dble(gammay(k,l,m))
        d2src_dz2 = dzs_dxsi * dble(xiz(k,l,m)) + dzs_deta * dble(etaz(k,l,m)) + dzs_dgamma * dble(gammaz(k,l,m))
        d2src_dxy = dxs_dxsi * dble(xiy(k,l,m)) + dxs_deta * dble(etay(k,l,m)) + dxs_dgamma * dble(gammay(k,l,m))
        d2src_dxz = dxs_dxsi * dble(xiz(k,l,m)) + dxs_deta * dble(etaz(k,l,m)) + dxs_dgamma * dble(gammaz(k,l,m))
        d2src_dyz = dys_dxsi * dble(xiz(k,l,m)) + dys_deta * dble(etaz(k,l,m)) + dys_dgamma * dble(gammaz(k,l,m))

        if (m == 3 .AND. l == 3 .AND. k == 3) then
           print *, "Hessian"
           print *, d2src_dx2, d2src_dy2, d2src_dz2, d2src_dxy, d2src_dxz, d2src_dyz
        endif
        
        ! With respect to x
        fxx = (Mxx * d2src_dx2 + Mxy * d2src_dxy + Mxz * d2src_dxz)
        fyx = (Mxy * d2src_dx2 + Myy * d2src_dxy + Myz * d2src_dxz)
        fzx = (Mxz * d2src_dx2 + Myz * d2src_dxy + Mzz * d2src_dxz)
          
        ! With respect to y
        fxy = (Mxx * d2src_dxy + Mxy * d2src_dy2 + Mxz * d2src_dyz)
        fyy = (Mxy * d2src_dxy + Myy * d2src_dy2 + Myz * d2src_dyz)
        fzy = (Mxz * d2src_dxy + Myz * d2src_dy2 + Mzz * d2src_dyz)

        ! With respect to z
        fxz = (Mxx * d2src_dxz + Mxy * d2src_dyz + Mxz * d2src_dz2)
        fyz = (Mxy * d2src_dxz + Myy * d2src_dyz + Myz * d2src_dz2)
        fzz = (Mxz * d2src_dxz + Myz * d2src_dyz + Mzz * d2src_dz2)

        ! Compute rotation factors
        ! note cos(theta) = sin(90 - theta)
        ! This is important since we want the derivative with respect to 
        ! the latitude
        if (m == 1 .AND. l == 1 .AND. k == 1) then
             print *, 'depth: ', depth
             print *, 'phi:   ', phi
             print *, 'theta: ', theta
        endif
        !dthetadlambda = 
        ! thetaprime = 
        sint = sin(theta)
        cost = cos(theta)
        sinp = sin(phi)
        cosp = cos(phi)
        grr_inv = ONE
        gtt_inv = ONE   
        gpp_inv = ONE / sint
        if (direction == 1) then
            fac_x = -1.d0 * sint * cosp / EARTH_R_KM
            fac_y = -1.d0 * sint * sinp / EARTH_R_KM
            fac_z = -1.d0 * cost / EARTH_R_KM
        else if (direction == 2) then
            fac_x = -1.d0 * cost * cosp * depth * DEGREES_TO_RADIANS
            fac_y = -1.d0 * cost * sinp * depth * DEGREES_TO_RADIANS
            fac_z =  1.d0 * sint * depth * DEGREES_TO_RADIANS
        else if (direction == 3) then
            fac_x = -1.d0 * sint * sinp * depth * DEGREES_TO_RADIANS
            fac_y = sint * cosp  * depth * DEGREES_TO_RADIANS
            fac_z = 0.d0 * depth * DEGREES_TO_RADIANS
        else 
            stop "Wrong direction. Should 1 for depth, 2 for lat, 3 for lon."
        endif
        
        ! Rotate
        fx = (fxx * fac_x + fxy * fac_y + fxz * fac_z)
        fy = (fyx * fac_x + fyy * fac_y + fyz * fac_z)
        fz = (fzx * fac_x + fzy * fac_y + fzz * fac_z)

        ! Add to sourcearray
        sourcearrayd(1,k,l,m) = sourcearrayd(1,k,l,m) + fx
        sourcearrayd(2,k,l,m) = sourcearrayd(2,k,l,m) + fy
        sourcearrayd(3,k,l,m) = sourcearrayd(3,k,l,m) + fz

      enddo
    enddo
  enddo

  ! distinguish between single and double precision for reals
  sourcearray(:,:,:,:) = real(sourcearrayd(:,:,:,:), kind=CUSTOM_REAL)

end subroutine compute_arrays_source_derivative2

!================================================================

!================================================================

  subroutine compute_arrays_source_adjoint(myrank, adj_source_file, &
                                           nu,source_adjoint, &
                                           NSTEP_BLOCK,iadjsrc,it_sub_adj,NSTEP_SUB_ADJ, &
                                           NTSTEP_BETWEEN_READ_ADJSRC,DT)

  use constants, only: CUSTOM_REAL,NDIM,IIN_ADJ,MAX_STRING_LEN

  use specfem_par, only: scale_displ_inv, NUMBER_OF_SIMULTANEOUS_RUNS, READ_ADJSRC_ASDF, mygroup

!  use iso_c_binding, only: C_NULL_CHAR

  implicit none

! input -- notice here NSTEP_BLOCK is different from the NSTEP in the main program
! instead NSTEP_BLOCK = iadjsrc_len(it_sub_adj), the length of this specific block

  integer,intent(in) :: myrank
  character(len=*),intent(in) :: adj_source_file

  double precision, dimension(NDIM,NDIM),intent(in) :: nu


  ! output
  integer,intent(in) :: NTSTEP_BETWEEN_READ_ADJSRC
  real(kind=CUSTOM_REAL),intent(out) :: source_adjoint(NDIM,NTSTEP_BETWEEN_READ_ADJSRC)

  integer,intent(in) :: NSTEP_SUB_ADJ
  integer, dimension(NSTEP_SUB_ADJ,2),intent(in) :: iadjsrc

  integer,intent(in) :: NSTEP_BLOCK
  integer,intent(in) :: it_sub_adj

  double precision,intent(in) :: DT

  ! local parameters
  double precision, dimension(NDIM,NSTEP_BLOCK) :: adj_src_u

  real(kind=CUSTOM_REAL), dimension(NDIM,NSTEP_BLOCK) :: adj_src
  real(kind=CUSTOM_REAL), dimension(NSTEP_BLOCK) :: adj_source_asdf
  real(kind=CUSTOM_REAL) :: junk

  integer :: icomp, itime, ios
  integer :: index_start,index_end,index_i
  character(len=3),dimension(NDIM) :: comp
  character(len=MAX_STRING_LEN) :: filename, path_to_add
  character(len=80) :: adj_source_name
  character(len=2) :: bic

  call band_instrument_code(DT,bic)
  comp(1) = bic(1:2)//'N'
  comp(2) = bic(1:2)//'E'
  comp(3) = bic(1:2)//'Z'

  ! safety check
  if (NSTEP_BLOCK > NTSTEP_BETWEEN_READ_ADJSRC) then
    print *,'Error invalid NSTEP_BLOCK ',NSTEP_BLOCK,' compared to NTSTEP_BETWEEN_READ_ADJSRC ',NTSTEP_BETWEEN_READ_ADJSRC
    call exit_MPI(myrank,'Error invalid NSTEP_BLOCK size in compute_array_source_adjoint')
  endif

  ! (sub)trace start and end
  ! reading starts in chunks of NSTEP_BLOCK from the end of the trace,
  ! i.e. as an example: total length NSTEP = 3000, chunk length NSTEP_BLOCK= 1000
  !                                then it will read in first index_start=2001 to index_end=3000,
  !                                second time, it will be index_start=1001 to index_end=2000 and so on...
  index_start = iadjsrc(it_sub_adj,1)
  index_end = iadjsrc(it_sub_adj,1)+NSTEP_BLOCK-1

  ! unfortunately, things become more tricky because of the Newmark time scheme at
  ! the very beginning of the time loop. however, when we read in the backward/reconstructed
  ! wavefields at the end of the first time loop, we can use the adjoint source index from 3000 down to 1.
  !
  ! see the comment on where we add the adjoint source (compute_add_sources_adjoint()).
  !
  ! otherwise,
  ! we would have to shift this indices by minus 1, to read in the adjoint source trace between 0 to 2999.
  ! since 0 index is out of bounds, we would have to put that adjoint source displacement artificially to zero
  !
  ! here now, index_start is now 2001 and index_end = 3000, then 1001 to 2000, then 1 to 1000.
  index_start = index_start
  index_end = index_end

  itime = 0
  adj_src(:,:) = 0._CUSTOM_REAL

  if (READ_ADJSRC_ASDF) then
    ! ASDF format
    do icomp = 1, NDIM ! 3 components

      ! print *, "READING ADJOINT SOURCES USING ASDF"

      adj_source_name = trim(adj_source_file) // '_' // comp(icomp)

      ! would skip read and set source artificially to zero if out of bounds, see comments above
      if (index_start == 0 .and. itime == 0) then
        adj_src(icomp,1) = 0._CUSTOM_REAL
        cycle
      endif

      call read_adjoint_sources_ASDF(adj_source_name, adj_source_asdf, index_start, index_end)

      adj_src(icomp,:) = real(adj_source_asdf(1:NSTEP_BLOCK))

    enddo

  else
    ! ASCII format
    do icomp = 1, NDIM

      ! opens adjoint component file
      filename = 'SEM/'//trim(adj_source_file) // '.'// comp(icomp) // '.adj'

      if (NUMBER_OF_SIMULTANEOUS_RUNS > 1 .and. mygroup >= 0) then
        write(path_to_add,"('run',i4.4,'/')") mygroup + 1
        filename = path_to_add(1:len_trim(path_to_add))//filename(1:len_trim(filename))
      endif

      open(unit=IIN_ADJ,file=trim(filename),status='old',action='read',iostat=ios)

      ! note: adjoint source files must be available for all three components E/N/Z, even
      !          if a component is just zeroed out
      if (ios /= 0) then
        ! adjoint source file not found
        ! stops simulation
        call exit_MPI(myrank, &
            'file '//trim(filename)//' not found, please check with your STATIONS_ADJOINT file')
      endif
      !if (ios /= 0) cycle ! cycles to next file - this is too error prone and users might easily end up with wrong results

      ! jumps over unused trace length
      do itime  = 1,index_start-1
        read(IIN_ADJ,*,iostat=ios) junk,junk
        if (ios /= 0) &
          call exit_MPI(myrank, &
            'file '//trim(filename)//' has wrong length, please check with your simulation duration')
      enddo

      ! reads in (sub)trace
      do itime = index_start,index_end

        ! index will run from 1 to NSTEP_BLOCK
        index_i = itime - index_start + 1

        ! would skip read and set source artificially to zero if out of bounds, see comments above
        if (index_start == 0 .and. itime == 0) then
          adj_src(icomp,1) = 0._CUSTOM_REAL
          cycle
        endif

        ! reads in adjoint source trace
        !read(IIN_ADJ,*,iostat=ios) junk, adj_src(icomp,itime-index_start+1)
        read(IIN_ADJ,*,iostat=ios) junk, adj_src(icomp,index_i)

        if (ios /= 0) then
          print *,'Error reading adjoint source: ',trim(filename)
          print *,'rank ',myrank,' - time step: ',itime,' index_start: ',index_start,' index_end: ',index_end
          print *,'  ',trim(filename)//'has wrong length, please check with your simulation duration'
          call exit_MPI(myrank,'file '//trim(filename)//' has wrong length, please check with your simulation duration')
        endif
      enddo

      close(IIN_ADJ)

    enddo
  endif

  ! non-dimensionalize
  adj_src(:,:) = adj_src(:,:) * scale_displ_inv

  ! rotates to Cartesian
  do itime = 1, NSTEP_BLOCK
    adj_src_u(:,itime) = nu(1,:) * adj_src(1,itime) &
                       + nu(2,:) * adj_src(2,itime) &
                       + nu(3,:) * adj_src(3,itime)
  enddo

  do icomp = 1, NDIM
    source_adjoint(icomp,1:NSTEP_BLOCK) = adj_src_u(icomp,1:NSTEP_BLOCK)
  enddo

! not used, but for reference in case lagrange interpolators will be added here again ...
!
!  contains
!
!    subroutine multiply_arrays_adjoint(sourcearrayd,hxir,hetar,hgammar,adj_src_ud)
!
!    use constants
!
!    implicit none
!
!    double precision, dimension(NDIM,NGLLX,NGLLY,NGLLZ) :: sourcearrayd
!    double precision, dimension(NGLLX) :: hxir
!    double precision, dimension(NGLLY) :: hetar
!    double precision, dimension(NGLLZ) :: hgammar
!    double precision, dimension(NDIM) :: adj_src_ud
!
!    integer :: i,j,k
!
!    ! adds interpolated source contribution to all GLL points within this element
!    do k = 1, NGLLZ
!      do j = 1, NGLLY
!        do i = 1, NGLLX
!          sourcearrayd(:,i,j,k) = hxir(i) * hetar(j) * hgammar(k) * adj_src_ud(:)
!        enddo
!      enddo
!    enddo
!
!    end subroutine multiply_arrays_adjoint

  end subroutine compute_arrays_source_adjoint
