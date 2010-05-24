
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  glimmer_paramets.f90 - part of the GLIMMER ice model     + 
! +                                                           +
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! 
! Copyright (C) 2004 GLIMMER contributors - see COPYRIGHT file 
! for list of contributors.
!
! This program is free software; you can redistribute it and/or 
! modify it under the terms of the GNU General Public License as 
! published by the Free Software Foundation; either version 2 of 
! the License, or (at your option) any later version.
!
! This program is distributed in the hope that it will be useful, 
! but WITHOUT ANY WARRANTY; without even the implied warranty of 
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the 
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License 
! along with this program; if not, write to the Free Software 
! Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 
! 02111-1307 USA
!
! GLIMMER is maintained by:
!
! Ian Rutt
! School of Geographical Sciences
! University of Bristol
! University Road
! Bristol
! BS8 1SS
! UK
!
! email: <i.c.rutt@bristol.ac.uk> or <ian.rutt@physics.org>
!
! GLIMMER is hosted on berliOS.de:
!
! https://developer.berlios.de/projects/glimmer-cism/
!
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#ifdef HAVE_CONFIG_H
#include "config.inc"
#endif

!> model scaling constants
module glimmer_paramets

  use glimmer_global, only : sp, dp
  use glimmer_physcon, only : scyr, rhoi, grav, gn

  implicit none; save

#ifdef NO_RESCALE
  real(dp), parameter :: thk0 = 1          ! m
  real(dp), parameter :: len0 = 1        ! m
  real(dp), parameter :: vel0 = 1 / scyr    ! m yr^{-1} converted to S.I. units
  real(dp), parameter :: vis0 = 1 / scyr 
#else
  real(dp), parameter :: thk0 = 2000.0d0          ! m
  real(dp), parameter :: len0 = 200.0d3        ! m
  real(dp), parameter :: vel0 = 500.0 / scyr    ! m yr^{-1} converted to S.I. units
  !real(dp), parameter :: vis0 = 5.70d-18 / scyr  ! yr^{-1} Pa^{-3} converted to S.I. units
  real(dp), parameter :: vis0 = 1d-16 / scyr 
#endif

  ! *sfp* defined these to convert scales to values used by GLAM
  real(dp), parameter :: tau0_glam = rhoi*grav*thk0                   ! stress scale in GLAM ( Pa )  
  real(dp), parameter :: vis0_glam = tau0_glam**(-gn) * (vel0/len0)   ! rate factor scale in GLAM ( Pa^-3 s^-1 )
  real(dp), parameter :: evs0 = tau0_glam * (vel0/len0) ! eff. visc. scale in GLAM ( Pa s )


  real(dp), parameter :: acc0 = thk0 * vel0 / len0  ! m s^{-1} 
  ! ** for zero order model real(dp), parameter :: tim0 = thk0 / acc0      ! s
  real(dp), parameter :: tim0 = len0 / vel0      ! s
  real(dp) :: tau0                        ! Pa note cannot define here as f90 wont allow
                                          ! parameters with noninteger powers in - look
                                          ! in initial in blah.f90 (not sure this applies now...)

!whl - to do - remove lambda0
  real(dp), parameter :: lambda0 = 80.0d3 / len0    ! basal topo/friction parameter for ismip-hom tests 

  real(sp), parameter :: conv = tim0 / scyr

! some parameters for debugging
   integer, parameter ::   &
      itest = 133, jtest = 84,  &          ! in Greenland (FV2), lat 67.3 N, lon 330 E
                  jjtest = 97 - jtest,  &  ! reversed for N to S indexing (FV2, ny = 96)
      itest_local = 60, jtest_local = 54   ! Greenland 20 deg grid, initial usrf = 491 m

   integer, parameter :: idiag = 30, jdiag = 50  ! point for diagnostic output
   integer, parameter :: stdout = 6

!MH!  !MAKE_RESTART
!MH!#ifdef RESTARTS
!MH!#define RST_PARAMETS
!MH!#include "glimmer_rst_head.inc"
!MH!#undef RST_PARAMETS
!MH!#endif
!MH!
!MH!#ifdef RESTARTS
!MH!contains
!MH!
!MH!#define RST_PARAMETS
!MH!#include "glimmer_rst_body.inc"
!MH!#undef RST_PARAMETS
!MH!#endif

end module glimmer_paramets