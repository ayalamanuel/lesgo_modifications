!!
!!  Module: wave_file_io
!!
!!  Reads a single snapshot of eta(x) and detadx(x) from a text file,
!!  interpolates them onto the simulation grid using cubic splines,
!!  and stores the base profiles for use in time evolution (wave_type=3).
!!
!!  File format expected:
!!    Line 1:  N  (number of points in the file)
!!    Lines 2 to N+1:  x_i   eta_i   detadx_i
!!
!!  The x-coordinates in the file should be DIMENSIONAL (meters).
!!  The eta and detadx values should also be dimensional.
!!  They will be non-dimensionalized internally by z_i and z_i respectively
!!  (detadx is already dimensionless).
!!
!!  The profile is assumed to be 1D (x-direction), uniform in y,
!!  consistent with mono_wave and stokes_wave.
!!
!*******************************************************************************
module read_wave
!*******************************************************************************
use types, only : rprec
implicit none

private
public :: read_wave_file, eta_0, detadx_0

! Base profiles on the simulation grid (nx, ny), non-dimensional
real(rprec), dimension(:,:), allocatable :: eta_0, detadx_0

contains

!*******************************************************************************
subroutine read_wave_file()
!*******************************************************************************
! Reads the wave file, interpolates to the simulation grid, and stores
! the base profiles eta_0 and detadx_0.

use param, only : nx, ny, ld, dx, z_i, L_x, path, coord
implicit none

integer :: n_file, i, jy
real(rprec), allocatable :: x_file(:), eta_file(:), detadx_file(:)
real(rprec), dimension(nx) :: x_sim, eta_interp, detadx_interp
character(256) :: fname
integer :: io_stat

! Only coord 0 needs the wave surface, but allocate everywhere for safety
if (.not. allocated(eta_0)) allocate(eta_0(nx, ny))
if (.not. allocated(detadx_0)) allocate(detadx_0(nx, ny))
eta_0 = 0.0_rprec
detadx_0 = 0.0_rprec

! Build simulation x-grid (dimensional, in meters)
do i = 1, nx
    x_sim(i) = (i - 1) * dx * z_i
end do

! Read the file
fname = trim(path) // 'wave_input.dat'
if (coord == 0) write(*,*) 'Reading wave input file: ', trim(fname)

open(unit=99, file=trim(fname), status='old', action='read', iostat=io_stat)
if (io_stat /= 0) then
    write(*,*) 'ERROR: Cannot open wave input file: ', trim(fname)
    stop
end if

read(99, *) n_file
if (coord == 0) write(*,*) '  Number of points in wave file: ', n_file

allocate(x_file(n_file), eta_file(n_file), detadx_file(n_file))

do i = 1, n_file
    read(99, *) x_file(i), eta_file(i), detadx_file(i)
end do
close(99)

! Interpolate eta and detadx from file grid to simulation grid
call cubic_spline_interp(x_file, eta_file, n_file, x_sim, eta_interp, nx)
call cubic_spline_interp(x_file, detadx_file, n_file, x_sim, detadx_interp, nx)

! Non-dimensionalize: eta by z_i, detadx is already dimensionless
! (detadx = d(eta)/d(x), both in meters, so the ratio is dimensionless)
eta_interp = eta_interp / z_i

! Fill the 2D arrays (uniform in y)
do jy = 1, ny
    eta_0(:, jy) = eta_interp(:)
    detadx_0(:, jy) = detadx_interp(:)
end do

if (coord == 0) then
    write(*,*) '  Wave file read and interpolated successfully'
    write(*,*) '  max|eta_0| (non-dim) = ', maxval(abs(eta_0))
    write(*,*) '  max|detadx_0|        = ', maxval(abs(detadx_0))
end if

deallocate(x_file, eta_file, detadx_file)

end subroutine read_wave_file

!*******************************************************************************
subroutine cubic_spline_interp(x_in, y_in, n_in, x_out, y_out, n_out)
!*******************************************************************************
! Natural cubic spline interpolation.
! Input:  x_in(n_in), y_in(n_in) - data points (must be sorted in ascending x)
! Output: y_out(n_out) - interpolated values at x_out(n_out)
!
! For points outside the input range, linear extrapolation from the
! nearest interval is used.

implicit none
integer, intent(in) :: n_in, n_out
real(rprec), intent(in) :: x_in(n_in), y_in(n_in), x_out(n_out)
real(rprec), intent(out) :: y_out(n_out)

real(rprec), allocatable :: h(:), alpha(:), l(:), mu(:), z(:), b(:), c(:), d(:)
integer :: i, j, lo, hi, mid
real(rprec) :: dx_loc, a_loc

allocate(h(n_in-1), alpha(n_in), l(n_in), mu(n_in), z(n_in))
allocate(b(n_in), c(n_in), d(n_in-1))

! Step 1: compute h_i = x_{i+1} - x_i
do i = 1, n_in - 1
    h(i) = x_in(i+1) - x_in(i)
end do

! Step 2: compute alpha (RHS for tridiagonal system)
alpha(1) = 0.0_rprec
alpha(n_in) = 0.0_rprec
do i = 2, n_in - 1
    alpha(i) = 3.0_rprec / h(i) * (y_in(i+1) - y_in(i)) &
             - 3.0_rprec / h(i-1) * (y_in(i) - y_in(i-1))
end do

! Step 3: solve tridiagonal system for c (natural spline: c(1)=c(n)=0)
l(1) = 1.0_rprec
mu(1) = 0.0_rprec
z(1) = 0.0_rprec

do i = 2, n_in - 1
    l(i) = 2.0_rprec * (x_in(i+1) - x_in(i-1)) - h(i-1) * mu(i-1)
    mu(i) = h(i) / l(i)
    z(i) = (alpha(i) - h(i-1) * z(i-1)) / l(i)
end do

l(n_in) = 1.0_rprec
z(n_in) = 0.0_rprec
c(n_in) = 0.0_rprec

! Back substitution
do i = n_in - 1, 1, -1
    c(i) = z(i) - mu(i) * c(i+1)
end do

! Step 4: compute b and d coefficients
do i = 1, n_in - 1
    b(i) = (y_in(i+1) - y_in(i)) / h(i) - h(i) * (c(i+1) + 2.0_rprec * c(i)) / 3.0_rprec
    d(i) = (c(i+1) - c(i)) / (3.0_rprec * h(i))
end do

! Step 5: evaluate spline at output points
do j = 1, n_out
    ! Binary search for the interval containing x_out(j)
    if (x_out(j) <= x_in(1)) then
        ! Extrapolate from first interval
        i = 1
    else if (x_out(j) >= x_in(n_in)) then
        ! Extrapolate from last interval
        i = n_in - 1
    else
        ! Binary search
        lo = 1
        hi = n_in
        do while (hi - lo > 1)
            mid = (lo + hi) / 2
            if (x_in(mid) <= x_out(j)) then
                lo = mid
            else
                hi = mid
            end if
        end do
        i = lo
    end if

    dx_loc = x_out(j) - x_in(i)
    y_out(j) = y_in(i) + b(i) * dx_loc + c(i) * dx_loc**2 + d(i) * dx_loc**3
end do

deallocate(h, alpha, l, mu, z, b, c, d)

end subroutine cubic_spline_interp

end module read_wave
