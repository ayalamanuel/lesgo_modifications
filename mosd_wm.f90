!!
!!  Copyright (C) 2011-2017  Johns Hopkins University
!!
!!  This file is part of lesgo.
!!
!!  lesgo is free software: you can redistribute it and/or modify
!!  it under the terms of the GNU General Public License as published by
!!  the Free Software Foundation, either version 3 of the License, or
!!  (at your option) any later version.
!!
!!  lesgo is distributed in the hope that it will be useful,
!!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!  GNU General Public License for more details.
!!
!!  You should have received a copy of the GNU General Public License
!!  along with lesgo.  If not, see <http://www.gnu.org/licenses/>.
!!
!! Manuel trial
!*******************************************************************************
module mosd_wm
!*******************************************************************************
! This module has all the calculations of the MOSD model, the equilibirum 
! wall model (EQM) is called from eqmfit_wm.f90 and the windward potential
! model (WPM) is calculated here.

use types, only : rprec
use param, only : wave_type, nx, ny, ld, dx, dz, dt, total_time_dim, total_time,&
                  amp, wave_n, wave_freq, u_star, z_i, pi, kp_spec, zgrid_match,&
                  nu_molec, vonk, zo, smooth_eqm, dy, coord, nz
use sim_param, only : eqmxz, eqmyz, s_wpmxz, s_wpmyz, eta, &
                      detadx, detady, detadt, u_orb, w_orb, deta2dx2, deta2dy2, &
                      detadxdy, detadydx, grad_eta_mag, dgrad_etadt,&
                      Cx_wave,Cy_wave, ur_mag_wpm, ur_mag_eqm, &
                      u_delta_m, u, v, detadx_dt, detady_dt
use sim_param, only : eta_hat_o, eta_hat_o_filter, omega_wave
use sim_param, only : delta_wpm
use grid_m
use test_filtermodule
use functions, only : cell_indx
use fft
use messages
use derivatives, only : ddx_fd, ddy_fd
use wave_spectrum, only : frequency_shift
use eqmfit_wm, only : eqmfit_calc
use read_wave, only : eta_0, detadx_0
implicit none
integer :: fid

private  
public :: wpm_calc,  mosd_finalize, mono_wave, spectrum_wave, stokes_wave, file_wave

contains
!*******************************************************************************
subroutine mosd_finalize()
!*******************************************************************************
use param, only : nx, ny
use sim_param, only : txz, tyz
implicit none

call wpm_calc ()
call eqmfit_calc ()

txz(1:nx,1:ny,1) = eqmxz(:,:) +  s_wpmxz(:,:) 
tyz(1:nx,1:ny,1) = eqmyz(:,:) +  s_wpmyz(:,:) 

end subroutine mosd_finalize

!*******************************************************************************
subroutine wave_selec() 
!*******************************************************************************
! This subroutine checks what type of wave will be modeled and builds the 
! surface distribution, gradients, etc.

implicit none
character(*), parameter :: sub_name = 'wave type'

select case(wave_type)

        ! Monochromatic wave
        case (0)
        call mono_wave()
        
        ! Spectrum wave
        case (1)
        call spectrum_wave()

        ! Stokes wave
        case (2)
        call stokes_wave()

        ! File-based wave (read from wave_input.dat)
        case (3)
        call file_wave()

        ! Otherwise, invalid
        case default
        call error (sub_name, 'invalid')

end select

end subroutine wave_selec

!*******************************************************************************
subroutine mono_wave() 
!*******************************************************************************
implicit none
integer :: jx, jy
real(rprec), dimension(ld, ny) :: x_grid
real(rprec), dimension(nx, ny) :: detadx_new

do jy = 1,ny
   do jx= 1,ld
        
      x_grid(jx,jy) = (jx-1)*dx
        
   end do
end do


eta = amp*cos(wave_n*x_grid(1:nx,1:ny) - wave_freq*total_time) 
detadx_new = amp*wave_n*sin(wave_freq*total_time - wave_n*x_grid(1:nx,1:ny))
detady = 0.0_rprec
detadt = amp*wave_freq*sin(wave_n*x_grid(1:nx,1:ny)   - wave_freq*total_time)
u_orb = amp*wave_freq*cos(wave_n*x_grid(1:nx,1:ny)  - wave_freq*total_time)
w_orb = amp*wave_freq*sin(wave_n*x_grid(1:nx,1:ny)  - wave_freq*total_time)
deta2dx2 = -amp*wave_n**2*cos(-wave_freq*total_time + wave_n*x_grid(1:nx,1:ny))
deta2dy2 = 0.0_rprec
detadxdy = 0.0_rprec
detadydx = 0.0_rprec

detadx_dt = (detadx_new(:,:) - detadx(:,:))/dt
detadx = detadx_new(:,:)

end subroutine mono_wave

!*******************************************************************************
subroutine stokes_wave()
!*******************************************************************************
implicit none
integer :: jx, jy
real(rprec), dimension(ld, ny) :: x_grid
real(rprec), dimension(nx, ny) :: detadx_new
real(rprec), dimension(nx, ny) :: theta
real(rprec) :: eps, k, w, a

! build x_grid
do jy = 1, ny
   do jx = 1, ld
      x_grid(jx, jy) = (jx-1)*dx
   end do
end do

! shorthand
a   = amp
k   = wave_n
w   = wave_freq
eps = k*a

theta = k*x_grid(1:nx,1:ny) - w*total_time

! ---------------------------
! 3rd-order Stokes surface eta
! eta = a cosθ + 1/2 (ka) a cos2θ + 3/8 (ka)^2 a cos3θ
! ---------------------------
eta = a*cos(theta) + 0.5_rprec*eps*a*cos(2.0_rprec*theta) + (3.0_rprec/8.0_rprec)*eps**2*a*cos(3.0_rprec*theta)

! slope in x: dη/dx = -a k [ sinθ + eps sin2θ + (9/8) eps^2 sin3θ ]
detadx_new = -a*k*( sin(theta) + eps*sin(2.0_rprec*theta) + (9.0_rprec/8.0_rprec)*eps**2*sin(3.0_rprec*theta) )

detady = 0.0_rprec

! time derivative: dη/dt = a w [ sinθ + eps sin2θ + (9/8) eps^2 sin3θ ]
detadt =  a*w*( sin(theta) + eps*sin(2.0_rprec*theta) + (9.0_rprec/8.0_rprec)*eps**2*sin(3.0_rprec*theta) )

! orbital velocities at z=0 (deep-water Stokes 3rd order)
! u = a w [ cosθ + eps cos2θ + (9/8) eps^2 cos3θ ]
! w = a w [ sinθ + eps sin2θ + (9/8) eps^2 sin3θ ]
u_orb = a*w*( cos(theta) + eps*cos(2.0_rprec*theta) + (9.0_rprec/8.0_rprec)*eps**2*cos(3.0_rprec*theta) )
w_orb = a*w*( sin(theta) + eps*sin(2.0_rprec*theta) + (9.0_rprec/8.0_rprec)*eps**2*sin(3.0_rprec*theta) )

! curvature in x: d2η/dx2 = -a k^2[ cosθ + 2 eps cos2θ + (27/8) eps^2 cos3θ ]
deta2dx2 = -a*k**2*( cos(theta) + 2.0_rprec*eps*cos(2.0_rprec*theta) + (27.0_rprec/8.0_rprec)*eps**2*cos(3.0_rprec*theta) )

deta2dy2  = 0.0_rprec
detadxdy  = 0.0_rprec
detadydx  = 0.0_rprec

! time derivative of slope (as you had it)
detadx_dt = (detadx_new(:,:) - detadx(:,:))/dt
detadx    = detadx_new(:,:)

end subroutine stokes_wave


!*******************************************************************************
subroutine file_wave()
!*******************************************************************************
! Frozen-shape translation: the wave profile read from file propagates
! at constant phase speed c without changing shape.
!
!   eta(x, t) = eta_0(x - c*t)
!   detadx(x, t) = detadx_0(x - c*t)
!
! This is equivalent to Sullivan et al. (2018) Eq. (8a)-(8b), where each
! Fourier mode of the stored profile is advanced by exp[j*k_m*(x - c*dt)].
!
! Implementation: the base profile eta_0 is stored on the simulation grid
! at initialization. At each timestep, we compute the shift distance
! c*t, convert to grid-cell units, and do a circular lookup with linear
! interpolation for the fractional cell offset. Since we always read
! from the pristine base profile (not from the previous timestep), there
! is no accumulated numerical diffusion.
!
! Parameters needed from param module:
!   wave_n    -> k (wavenumber), used to compute c = omega/k
!   wave_freq -> omega (angular frequency)
!
implicit none
integer :: jx, jy, jx_src, jx_src_p1
real(rprec), dimension(nx, ny) :: eta_new, detadx_new
real(rprec) :: c_phase, shift, shift_cells, frac
integer :: shift_int

! Phase speed (non-dimensional): c = omega / k
c_phase = wave_freq / wave_n

! Total shift in non-dimensional x-units
shift = c_phase * total_time

! Convert to grid-cell units
shift_cells = shift / dx
shift_int = floor(shift_cells)
frac = shift_cells - real(shift_int, rprec)

! Circular shift with linear interpolation from base profile
do jy = 1, ny
   do jx = 1, nx
      ! Source indices: position (jx - shift), wrapped periodically
      jx_src    = modulo(jx - 1 - shift_int,     nx) + 1
      jx_src_p1 = modulo(jx - 1 - shift_int - 1, nx) + 1

      ! Linear interpolation for fractional grid-cell offset
      eta_new(jx, jy)    = (1.0_rprec - frac) * eta_0(jx_src, jy) &
                         +  frac * eta_0(jx_src_p1, jy)
      detadx_new(jx, jy) = (1.0_rprec - frac) * detadx_0(jx_src, jy) &
                         +  frac * detadx_0(jx_src_p1, jy)
   end do
end do

! No y-component for 2D (x-z) waves
detady = 0.0_rprec

! Time derivatives via finite differences (same approach as spectrum_wave)
detadt    = (eta_new(:,:) - eta(:,:)) / dt
detadx_dt = (detadx_new(:,:) - detadx(:,:)) / dt

! Update stored values
eta    = eta_new(:,:)
detadx = detadx_new(:,:)

! Orbital velocities (deep water, linear theory)
! u_orb from Sullivan Eq. (10): u_o = (h_p * pi * c / lambda) * cos(2*pi*x/lambda)
! Simplified: u_orb ~ omega * eta,  w_orb ~ detadt
u_orb = wave_freq * eta(:,:)
w_orb = detadt(:,:)

! Second spatial derivative via finite differences
call ddx_fd(detadx, deta2dx2)
deta2dx2(nx,:) = -detadx(nx,:) / dx

deta2dy2  = 0.0_rprec
detadxdy  = 0.0_rprec
detadydx  = 0.0_rprec
detady_dt = 0.0_rprec

end subroutine file_wave


!*******************************************************************************
subroutine spectrum_wave() 
!*******************************************************************************
implicit none
integer :: jx, jy
real(rprec), dimension(nx, ny) :: u_orb_1, eta_new, eta_1, detadx_new, detady_new
complex(rprec), dimension(nx,ny) :: eta_hat, eta_shift, u_orb_hat, u_orb_shift
complex(rprec), dimension(nx/2+1,ny) :: eta_hat_2, u_orb_hat_2, eta_hat_filter_2
complex(kind=8) :: ii

ii = cmplx(0.0d0,1.0d0)

do jy = 1,ny
   do jx = 1,nx
    
      eta_hat(jx,jy) = eta_hat_o(jx,jy)*exp(-ii*omega_wave(jx,jy)*total_time_dim)
      u_orb_hat(jx,jy) = eta_hat(jx,jy)*omega_wave(jx,jy)

   end do
end do

! The spectrum needs to be shifted before going through fft and then only grabing
! half of the spectrum for this type of fft
call frequency_shift(eta_hat,eta_shift)
call frequency_shift(u_orb_hat,u_orb_shift)

eta_hat_2 = eta_shift(:nx/2+1,1:ny)
u_orb_hat_2 = u_orb_shift(:nx/2+1,1:ny)

call dfftw_execute_dft_c2r(back_wave, eta_hat_2, eta_1)
call dfftw_execute_dft_c2r(back_wave, u_orb_hat_2, u_orb_1)

! Since the surface spectrum and orbital velocities are calculated using dimensions,
! a normalization must be applied. Because of the FFT, the surface spectrum has to be 
! devided by 2
eta_new = eta_1(:,:)/2.0_rprec/z_i
u_orb = u_orb_1(:,:)/2.0_rprec/u_star
w_orb = 0.0_rprec !This needs to be changed at soem point

! Calculating the time derivative and spatial derivatives of the surface
detadt = (eta_new(:,:) - eta(:,:))/dt
eta = eta_new(:,:)

call ddx_fd(eta, detadx_new)
call ddy_fd(eta, detady_new)
detadx_new(nx,:) = -eta(nx,:)/dx
detady_new(:,ny) = -eta(:,ny)/dy

call ddx_fd(detadx,deta2dx2)
call ddy_fd(detady,deta2dy2)
deta2dx2(nx,:) = -detadx(nx,:)/dx
deta2dy2(:,ny) = -detady(:,ny)/dy

call ddy_fd(detadx,detadxdy)
call ddx_fd(detady,detadydx)
detadxdy(:,ny) = -detadx(:,ny)/dy
detadydx(nx,:) = -detady(nx,:)/dx

detadx_dt = (detadx_new(:,:) - detadx(:,:))/dt
detady_dt = (detady_new(:,:) - detady(:,:))/dt
detadx = detadx_new(:,:)
detady = detady_new(:,:)


end subroutine spectrum_wave

!*******************************************************************************
subroutine wpm_calc ()
!*******************************************************************************
implicit none
integer :: jx, jy
real(rprec), pointer, dimension(:) :: z
real(rprec) :: threshold_wave_speed, k_min, L_x_wave, L_y_wave,   &
               z_diff, H_p, eta_mean
real(rprec), dimension(ld, ny) :: u_delta, v_delta
real(rprec), dimension(nx, ny) :: grad_eta_mag_new, ur_wpm, vr_wpm, nx_wpm,    &
                                  ny_wpm, H_arg, H_wpm, s_wpm, alpha
real(rprec), dimension(nx,ny) :: eta_p, ramp_eta_p
integer :: z_delta_low, z_delta_up
character(*), parameter :: sub_name = 'delta outside coord=0'
nullify(z)
z => grid % z

! Getting the wave profile, gradients, phase velocities, etc
call wave_selec()

! Calculating the velocity input for WPM at a height 3H_p. Where H_s is the significant
! wave height Not using the linear_interp, because it does not work. The velocity is also 
! filtered at delta scale. NOTE: ensure that height 3H_p lies in the coord = 0 since this 
! subroutine is only called at first processor.
eta_mean = sum(eta)/(nx*ny)
eta_p = eta - eta_mean
ramp_eta_p = max(eta_p, 0.0_rprec)**8_rprec  
H_p = ( sum(ramp_eta_p) / (nx*ny) )**0.125_rprec
delta_wpm = 3.0_rprec*H_p

! Check if delta_wpm is within the first processor domain (coord=0).
! z(nz-1) might be to cautious but thats ok!
if (delta_wpm > z(nz-1)) then
    write(*,*) '***********************************************'
    write(*,*) 'ERROR in wpm_calc: delta_wpm exceeds first processor domain!'
    write(*,*) 'delta_wpm = ', delta_wpm
    write(*,*) 'z(nz-1) for coord=0 = ', z(nz-1)
    write(*,*) '***********************************************'
call error (sub_name, 'invalid')
end if

! Taking the velocity for wpm at a height Delta = 3*Hp. If for some reason this height location is 
! below the first uv grid point 0.5*dz, then it will take the velocity at that first uv grid point.
! In the scenario where delta<0.5dz, one just has to make sure that the first grid point is roughly
! 3*hp-2*Hp
if (delta_wpm > z(1)) then 

z_delta_low = cell_indx('k',dz,delta_wpm)
z_delta_up = z_delta_low + 1

z_diff = delta_wpm - z(z_delta_low)
u_delta = u(1:ld,1:ny,z_delta_low) + (u(1:ld,1:ny,z_delta_up) - u(1:ld,1:ny,z_delta_low))*z_diff/dz
v_delta = v(1:ld,1:ny,z_delta_low) + (v(1:ld,1:ny,z_delta_up) - v(1:ld,1:ny,z_delta_low))*z_diff/dz

else
write(*,*) 'using the first grid point since 0.5*dz>delta_wpm'

u_delta = u(1:ld,1:ny,1)
v_delta = v(1:ld,1:ny,1)
end if

call test_filter(u_delta)
call test_filter(v_delta)

! Calculating the magnitude of the gradient, the temporal derivative of the gradient
! and the velocity components of the wave. The phase velocity is clipped to the
! velocity of the wave of size 0.25*kp
grad_eta_mag_new = sqrt((detadx(:,:))**2 + (detady(:,:))**2)
dgrad_etadt = (grad_eta_mag_new(:,:) - grad_eta_mag(:,:))/dt
grad_eta_mag = grad_eta_mag_new(:,:)

Cx_wave = -detadt(:,:)*detadx(:,:)*(1/grad_eta_mag(:,:)**2)
Cy_wave = -detadt(:,:)*detady(:,:)*(1/grad_eta_mag(:,:)**2)

if (wave_type == 1) then        
k_min = 0.25_rprec*(kp_spec)
threshold_wave_speed = sqrt(9.81_rprec/k_min)/u_star
Cx_wave = min(max(Cx_wave(:,:),-threshold_wave_speed),threshold_wave_speed)
Cy_wave = min(max(Cy_wave(:,:),-threshold_wave_speed),threshold_wave_speed)
end if

if (wave_type == 3) then
Cx_wave = wave_freq/wave_n
Cy_wave = 0.0_rprec
end if

! Calcualting the realtive velocities. u_delta and v_delta are already filtered 
! and we CANNOT filter Cx_wave, Cy_wave because they won't cancel out with 
! the surface gradients specially for spectrum wave. We also calculate the 
! magnitude of u,v flow velocity.
ur_wpm = u_delta(1:nx,1:ny) - Cx_wave(:,:)
vr_wpm = v_delta(1:nx,1:ny) - Cy_wave(:,:)

! Calculating the remainder of the steady WPM component
nx_wpm = detadx(:,:)/grad_eta_mag(:,:)
ny_wpm = detady(:,:)/grad_eta_mag(:,:)
ur_mag_wpm = sqrt((ur_wpm(:,:)*nx_wpm(:,:))**2 + (vr_wpm(:,:)*ny_wpm(:,:))**2)
alpha = atan(grad_eta_mag(:,:))

H_arg = (ur_wpm(:,:)*detadx(:,:) + vr_wpm(:,:)*detady(:,:))
H_wpm = (H_arg(:,:) + abs(H_arg(:,:)))/(2*H_arg(:,:))

s_wpm = alpha(:,:)/(pi+alpha(:,:))*(ur_mag_wpm(:,:)**2)*grad_eta_mag(:,:)

s_wpmxz  = -s_wpm(:,:)*nx_wpm(:,:)*H_wpm(:,:)
s_wpmyz  = -s_wpm(:,:)*ny_wpm(:,:)*H_wpm(:,:)


end subroutine wpm_calc

end module mosd_wm
