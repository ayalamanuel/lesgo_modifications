!!
!!  Copyright (C) 2009-2017  Johns Hopkins University
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

!*******************************************************************************
module sim_param
!*******************************************************************************
use types, only : rprec
use param, only : lh, ld, nx, ny, nz, lbz, u_star
implicit none

save
public

logical :: sim_param_initialized = .false.

real(rprec), dimension(:,:,:), allocatable :: u, v, w,                         &
    dudx, dudy, dudz, dvdx, dvdy, dvdz,  dwdx, dwdy, dwdz,                     &
    RHSx, RHSy, RHSz, RHSx_f, RHSy_f, RHSz_f,                                  &
    dpdx, dpdy, dpdz, txx, txy, tyy,                                           &
    txz, tyz, tzz, divtx, divty, divtz,                                        &
    fx, fy, fz, fxa, fya, fza, cc_z, v_omegaz, u_omegay,                       &
    dtxdx, dtydy, dtzdz

real(rprec), target, dimension(:,:,:), allocatable :: p

real(rprec), dimension(:,:), allocatable :: eqmxz, eqmyz, s_wpmxz, s_wpmyz,    &
    eta, detadx, detady, detadt, u_orb, w_orb, deta2dx2, &
    deta2dy2, detadxdy, detadydx, grad_eta_mag, dgrad_etadt,       &
    Cx_wave,Cy_wave, ur_mag_wpm, ur_mag_eqm, u_delta_m, ustar_lbc,&
    detadx_dt, detady_dt

real(rprec), dimension(:,:), allocatable :: K_wave, D_spread, phi, omega_wave ,&
    S_kx_ky

complex(rprec), dimension(:,:), allocatable :: eta_hat_o, eta_hat_o_filter

real(rprec) :: delta_wpm, delta_m
contains

!*******************************************************************************
subroutine sim_param_init ()
!*******************************************************************************
!
! This subroutine initilizes all global arrays defined in the sim_param
! module. Here they are allocated and initialized to zero.
!
implicit none

allocate ( u(ld, ny, lbz:nz) ); u = 0.0_rprec
allocate ( v(ld, ny, lbz:nz) ); v = 0.0_rprec
allocate ( w(ld, ny, lbz:nz) ); w = 0.0_rprec
allocate( dudx(ld, ny, lbz:nz) ); dudx = 0.0_rprec
allocate( dudy(ld, ny, lbz:nz) ); dudy = 0.0_rprec
allocate( dudz(ld, ny, lbz:nz) ); dudz = 0.0_rprec
allocate( dvdx(ld, ny, lbz:nz) ); dvdx = 0.0_rprec
allocate( dvdy(ld, ny, lbz:nz) ); dvdy = 0.0_rprec
allocate( dvdz(ld, ny, lbz:nz) ); dvdz = 0.0_rprec
allocate( dwdx(ld, ny, lbz:nz) ); dwdx = 0.0_rprec
allocate( dwdy(ld, ny, lbz:nz) ); dwdy = 0.0_rprec
allocate( dwdz(ld, ny, lbz:nz) ); dwdz = 0.0_rprec
allocate( RHSx(ld, ny, lbz:nz) ); RHSx = 0.0_rprec
allocate( RHSy(ld, ny, lbz:nz) ); RHSy = 0.0_rprec
allocate( RHSz(ld, ny, lbz:nz) ); RHSz = 0.0_rprec
allocate( cc_z(ld, ny, lbz:nz) ); cc_z = 0.0_rprec
allocate( v_omegaz(ld, ny, lbz:nz) ); v_omegaz = 0.0_rprec
allocate( u_omegay(ld, ny, lbz:nz) ); u_omegay = 0.0_rprec
allocate( RHSx_f(ld, ny, lbz:nz) ); RHSx_f = 0.0_rprec
allocate( RHSy_f(ld, ny, lbz:nz) ); RHSy_f = 0.0_rprec
allocate( RHSz_f(ld, ny, lbz:nz) ); RHSz_f = 0.0_rprec
allocate ( dpdx(ld, ny, nz) ); dpdx = 0.0_rprec
allocate ( dpdy(ld, ny, nz) ); dpdy = 0.0_rprec
allocate ( dpdz(ld, ny, nz) ); dpdz = 0.0_rprec
allocate ( txx(ld, ny, lbz:nz) ); txx = 0.0_rprec
allocate ( txy(ld, ny, lbz:nz) ); txy = 0.0_rprec
allocate ( tyy(ld, ny, lbz:nz) ); tyy = 0.0_rprec
allocate ( txz(ld, ny, lbz:nz) ); txz = 0.0_rprec
allocate ( tyz(ld, ny, lbz:nz) ); tyz = 0.0_rprec
allocate ( tzz(ld, ny, lbz:nz) ); tzz = 0.0_rprec
allocate ( p(ld, ny, 0:nz) ); p = 0.0_rprec
allocate ( divtx(ld, ny, lbz:nz) ); divtx = 0.0_rprec
allocate ( divty(ld, ny, lbz:nz) ); divty = 0.0_rprec
allocate ( divtz(ld, ny, lbz:nz) ); divtz = 0.0_rprec
allocate ( dtxdx(ld, ny, lbz:nz) ); dtxdx = 0.0_rprec
allocate ( dtydy(ld, ny, lbz:nz) ); dtydy = 0.0_rprec
allocate ( dtzdz(ld, ny, lbz:nz) ); dtzdz = 0.0_rprec

allocate ( eqmxz(nx, ny) ); eqmxz = 0.0_rprec
allocate ( eqmyz(nx, ny) ); eqmyz = 0.0_rprec
allocate ( s_wpmxz(nx, ny) ); s_wpmxz = 0.0_rprec
allocate ( s_wpmyz(nx, ny) ); s_wpmyz = 0.0_rprec
allocate ( eta(nx, ny) ); eta = 0.0_rprec
allocate ( detadx(nx, ny) ); detadx = 0.0_rprec
allocate ( detady(nx, ny) ); detady = 0.0_rprec
allocate ( detadt(nx, ny) ); detadt = 0.0_rprec
allocate ( u_orb(nx, ny) ); u_orb = 0.0_rprec
allocate ( w_orb(nx, ny) ); w_orb = 0.0_rprec
allocate ( grad_eta_mag(nx,ny) ); grad_eta_mag = 0.0_rprec
allocate ( Cx_wave(nx,ny) ); Cx_wave = 0.0_rprec
allocate ( Cy_wave(nx,ny) ); Cy_wave = 0.0_rprec
allocate ( ur_mag_wpm(nx,ny) ); ur_mag_wpm = 0.0_rprec
allocate ( ur_mag_eqm(nx,ny) ); ur_mag_eqm = 0.0_rprec
allocate ( u_delta_m(nx,ny) ); u_delta_m= 0.0_rprec

allocate ( K_wave(nx,ny) ); K_wave = 0.0_rprec
allocate ( D_spread(nx,ny) ); D_spread = 0.0_rprec
allocate ( S_kx_ky(nx,ny) ); S_kx_ky = 0.0_rprec
allocate ( phi(nx,ny) ); phi = 0.0_rprec
allocate ( omega_wave(nx,ny) ); omega_wave = 0.0_rprec
allocate ( eta_hat_o(nx,ny) ); eta_hat_o = 0.0_rprec

#if defined(PPTURBINES) || defined(PPATM) || defined(PPLVLSET)
allocate ( fxa(ld, ny, lbz:nz) ); fxa = 0.0_rprec
allocate ( fya(ld, ny, lbz:nz) ); fya = 0.0_rprec
allocate ( fza(ld, ny, lbz:nz) ); fza = 0.0_rprec
#endif

#if defined(PPLVLSET) || defined(PPATM)
allocate ( fx(ld, ny, nz) ); fx = 0.0_rprec
allocate ( fy(ld, ny, nz) ); fy = 0.0_rprec
allocate ( fz(ld, ny, nz) ); fz = 0.0_rprec
#endif

allocate( ustar_lbc(nx, ny) ); ustar_lbc = u_star

delta_wpm = 0.0_rprec
delta_m = 0.0_rprec

sim_param_initialized = .true.

end subroutine sim_param_init

end module sim_param
