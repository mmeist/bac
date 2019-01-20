module axis_mod
implicit none

private :: o_deriv, x_search_func, psi_func, calc_o_point, calc_x_point
public :: axis, default_theta_scaling

contains

subroutine axis(nt_core, theta_scaling_func, points)
  use period_mod
  use field_eq_mod, only : psif
  use psi4root, only : psi_root, R_o,  Z_o, theta, R_x, Z_x, R0, Z0
  use const, only : twopi

  implicit none

  interface
    function theta_scaling_func(x)
      implicit none
      integer, INTENT(IN) :: x
      double precision, dimension(x) :: theta_scaling_func
    end function theta_scaling_func
  END interface

  integer, dimension(:), intent(in) :: nt_core
  double precision, dimension(:, :), intent(out) :: points !

  integer, parameter :: iterat=200

  double precision :: hn1, hn1min, relerr, rrr, zzz, ppp
  double precision :: Br,Bp,Bz,dBrdR,dBrdp,dBrdZ,dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ
  integer i, j

  double precision :: theta_o_x
  double precision :: delta_xo, b_root, rho, tlr_root, rtbis
  double precision, dimension(:), allocatable :: R_o_x, Z_o_x, psif_o_x, delta_theta
  double precision, dimension(:), allocatable :: R_coord, Z_coord

  double precision :: tau_rad
  integer :: nr_core, n_points, point_idx

  nr_core = size(nt_core)
  n_points = sum(nt_core) + 1
  
  delta_xo  = 2.d0/dfloat(30-1)
  allocate(R_o_x(nr_core), Z_o_x(nr_core), psif_o_x(nr_core), &
  R_coord(n_points), Z_coord(n_points), delta_theta(maxval(nt_core)))
  
  
  ! parameters for odeint in calc_o_point
  hn1 = 0.1d0
  hn1min = 0.d0
  relerr = 1.d-12

  call calc_o_point(hn1, hn1min, relerr, R_o, Z_o)
  
  call calc_x_point(-341.d0, iterat, 1.d-14, R_x, Z_x)
  
  ! calculate Coords and psi values on o-x line
  do i=1, nr_core
    tau_rad = dfloat(i) / dfloat(nr_core)
    if(tau_rad >= 1.d0) then
      R_o_x(i) = R_x
      Z_o_x(i) = Z_x
    else
      R_o_x(i) = R_o + tau_rad*(R_x - R_o)
      Z_o_x(i) = Z_o + tau_rad*(Z_x - Z_o)
    end if

    ppp = 0.d0
    rrr = R_o_x(i)
    zzz = Z_o_x(i)
    call field(rrr, ppp, zzz, Br, Bp, Bz, dBrdR, dBrdp, dBrdZ   &
          ,dBpdR, dBpdp, dBpdZ, dBzdR, dBzdp, dBzdZ)
    psif_o_x(i) = psif

    if(tau_rad >= 1.d0) then
        !nr_core = i
        exit
    end if
  end do

  theta_o_x = atan2(Z_x-Z_o,R_x-R_o)
  print*, 'theta_o_x=', theta_o_x

  ! --- core region, calculate coords for curves of equal psi value for each psi in psif_o_x ---

  ! parameters of psi_func
  R0 = R_o
  Z0 = Z_o
  tlr_root = 1.d-7

  ! first point is O-Point
  R_coord(1) = R_o
  Z_coord(1) = Z_o

  point_idx = 2
  do i = 1, nr_core
    ! first point for every psi value is the point on the R_o_x line
    R_coord(point_idx) = R_o_x(i)
    Z_coord(point_idx) = Z_o_x(i)
    
    point_idx = point_idx + 1
    
    psi_root = psif_o_x(i)
    if (i < nr_core) then
      b_root = sqrt((R_o_x(i+1)-R_o)**2 + (Z_o_x(i+1)-Z_o)**2)
    else 
      b_root = sqrt((R_o_x(i)-R_o)**2 + (Z_o_x(i)-Z_o)**2)
    end if
    
    theta = theta_o_x
    delta_theta(:nt_core(i)) = theta_scaling_func(nt_core(i))
    do j = 2, nt_core(i)
      theta = theta_o_x + delta_theta(j - 1)
      if (theta >= twopi) theta = theta - (int(theta/twopi))*twopi
      if (theta < 0.) theta = theta + (int(abs(theta)/twopi) +1)*twopi
      rho = rtbis(psi_func, 0.d0, b_root, tlr_root)
      R_coord(point_idx) = rho*cos(theta) + R_o
      Z_coord(point_idx) = rho*sin(theta) + Z_o

      point_idx = point_idx + 1

    end do
  end do

  points(1, :) = R_coord
  points(2, :) = Z_coord
  points(3, :) = 0.d0

  !open(212,file='points.fmt')
  !do point_idx = 1, n_points
  !  write(212,101) R_coord(point_idx), Z_coord(point_idx)
  !end do
  !close(212)

  deallocate(Z_o_x, R_o_x, psif_o_x, R_coord, Z_coord, delta_theta)

  101 format(1000(e21.14,x))

end subroutine axis

subroutine calc_o_point(hn1, hn1min, relerr, raxis, zaxis)
  use const, only : twopi

  external rkqs

  double precision, intent(in) :: hn1, hn1min, relerr
  double precision, intent(out) :: raxis, zaxis

  double precision, dimension(2) :: y
  double precision :: phi0, dist_max, phi, phiout, per_phi, r, z, rmin, rmax, zmin, zmax
  integer :: nplot, nturns, npoint, idir, i, j, nok, nbad

  open(11, file='flt.inp')
  read(11,*)
  read(11,*)
  read(11,*)
  read(11,*) rmin
  read(11,*) rmax
  read(11,*) zmin
  read(11,*) zmax
  read(11,*) raxis
  read(11,*) zaxis
  read(11,*) per_phi
  read(11,*)
  read(11,*)
  read(11,*) idir
  read(11,*) nplot, nturns
  read(11,*)
  read(11,*)
  read(11,*) phi0
  close(11)
  
  phi0 = phi0 * twopi
  per_phi = per_phi * twopi

  do i=1, nplot
    npoint = 0
    phi = phi0

    y(1) = raxis
    y(2) = zaxis

    r = y(1)
    z = y(2)
    
    write(*,*) i, raxis, zaxis

    dist_max = 0.d0
    do j=1,nturns

        phiout = phi + per_phi * idir
        
        call odeint(y, 2, phi, phiout, relerr, hn1 * idir, hn1min, nok, nbad, o_deriv, rkqs) 
        
        npoint = npoint + 1
        phi = phiout
        r = y(1)  ! cylindic R
        z = y(2)  ! cylindic Z
        if (phi >= per_phi) phi = phi - (int(phi/per_phi)) * per_phi
        if (phi < 0.) phi = phi + (int(abs(phi)/per_phi) +1) * per_phi
        
        if( y(1) < rmin .or. y(1) > rmax .or.            &
            y(2) < zmin .or. y(2) > zmax ) exit

        raxis = raxis + r
        zaxis = zaxis + z
    end do

    raxis = raxis/dfloat(nturns+1)
    zaxis = zaxis/dfloat(nturns+1)
  end do

  print *, "Axis: "
  write(*,*) raxis, zaxis
end subroutine calc_o_point

subroutine calc_x_point(x_min_start, iterat, tlr_min, x1_min, x2_min)
  use fixed_coord, only : fixed, num_fixed

  double precision, intent(in) :: x_min_start, tlr_min
  integer, intent(in) :: iterat
  double precision, intent(out) :: x1_min, x2_min

  double precision :: a_min, b_min, c_min, z_min, x_min, z1_min, brent
  integer :: i, j

  x_min = x_min_start
  do i=1, iterat
     do j=1,2
!        num_fixed = modulo(j,2) + 1
        if(j .eq. 1) then
           num_fixed = 2
           a_min = 138.d0 
           b_min = 143.d0 
           c_min = 147.d0 
        else
           num_fixed = 1
           a_min = -93.d0 !50.d0
           b_min = -95.d0 !60.d0
           c_min = -97.d0 !70.d0
        endif
           fixed = x_min
           z_min = brent(a_min, b_min, c_min, x_search_func, tlr_min, x_min)
        if(j.eq.1) then
           x1_min = x_min
           z1_min = -z_min
        else
           x2_min = x_min
        endif
     end do
  end do

  print *, "X-point: "
  print *, x1_min, x2_min, z_min

end subroutine calc_x_point

subroutine o_deriv(phi,y,yp)
  use period_mod

  implicit none

  integer, parameter :: nequat = 2
  double precision y(nequat),yp(nequat)
  double precision :: phi, rrr, zzz, ppp
  double precision :: Br,Bp,Bz,dBrdR,dBrdp,dBrdZ,dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ
  rrr = y(1)
  zzz = y(2)
  ppp = phi
  if (ppp .ge. per_phi) ppp = ppp - (int(ppp/per_phi))*per_phi
  if (ppp .lt. 0.) ppp = ppp + (int(abs(ppp)/per_phi) +1)*per_phi

  call field(rrr,ppp,zzz,Br,Bp,Bz,dBrdR,dBrdp,dBrdZ   &
       ,dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ)
  yp(1) = Br*rrr/Bp
  yp(2) = Bz*rrr/Bp

  return
end subroutine o_deriv
! -----------------------------------------------------------------
double precision function x_search_func(s)
  use fixed_coord, only : fixed, num_fixed
  use field_eq_mod, only : psif
  implicit none
  double precision :: s
  double precision :: rrr, zzz, ppp
  double precision :: Br,Bp,Bz,dBrdR,dBrdp,dBrdZ,dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ

  ppp = 0.d0
  if(num_fixed .eq. 1) then
     rrr = fixed
     zzz = s
    ! call field to calculate psif
     call field(rrr, ppp, zzz, Br, Bp, Bz, dBrdR, dBrdp, dBrdZ   &
          ,dBpdR, dBpdp, dBpdZ, dBzdR, dBzdp, dBzdZ)
     x_search_func = psif
  else
     rrr = s
     zzz = fixed   
    ! call field to calculate psif
     call field(rrr, ppp, zzz, Br, Bp, Bz, dBrdR, dBrdp, dBrdZ   &
          ,dBpdR, dBpdp, dBpdZ, dBzdR, dBzdp, dBzdZ)
     x_search_func = -psif
  endif
  return
end function x_search_func
! --------------------------------------------
double precision function psi_func(rho)
  use field_eq_mod, only : psif
  use psi4root, only : psi_root, R0,  Z0, theta
  implicit none
  double precision :: rho
  double precision :: rrr, zzz, ppp
  double precision :: Br,Bp,Bz,dBrdR,dBrdp,dBrdZ,dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ

  ppp = 0.d0
  rrr = rho*cos(theta) + R0
  zzz = rho*sin(theta) + Z0

  ! call field to calculate psif
  call field(rrr, ppp, zzz, Br, Bp, Bz, dBrdR, dBrdp, dBrdZ   &
       ,dBpdR, dBpdp, dBpdZ, dBzdR, dBzdp, dBzdZ)
  psi_func = psif - psi_root

  return
end function psi_func

function default_theta_scaling(n_theta) result(delta_theta)
  use const, only : pi
  double precision, parameter :: th_mesh_co=-1.5d0
  integer, intent(in) :: n_theta

  double precision, dimension(n_theta) :: delta_theta
  integer :: i, j

  do j=1, n_theta/2
    delta_theta(j) = pi*(exp(-th_mesh_co*dfloat(j)/dfloat(n_theta/2)) - 1.d0) &
                                              /(exp(-th_mesh_co) - 1.d0)
  end do

  i=0
  do j=n_theta/2+1,n_theta-1
    i = n_theta - j + 1
    delta_theta(j) =  delta_theta(j-1) + (delta_theta(i)-delta_theta(i-1)) 
  end do
end function default_theta_scaling

end module axis_mod

!program axis_test
!  use axis_mod, only : axis, default_theta_scaling
!
!  double precision, allocatable, dimension(:, :) :: points
!  integer :: n_points
!  
!  integer, dimension(15) :: vpr
!  !vpr = (/(i, i=6, 6 + size(vpr) * 2 - 1, 2)/)
!  vpr = 50
!
!  n_points = sum(vpr) + 1
!
!  allocate(points(3, n_points))
!
!  call axis(vpr, default_theta_scaling, points)
!
!  deallocate(points)
!end program axis_test