! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  isostasy_el.f90 - part of the GLIMMER ice model          + 
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
! GLIMMER is hosted on NeSCForge:
!
! http://forge.nesc.ac.uk/projects/glimmer/
!
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

module isostasy_el

  !*FD handle elastic lithosphere

  use glimmer_global, only : dp
  use isostasy_types

  real, private, parameter :: r_lr=6.0        ! influence of disk load at (0,0) is felt within a radius of rbel_r_lr*rbel_r
  
  private :: init_rbel,rbel_ow, rbel_iw

contains
  
  subroutine init_elastic(rbel, deltax)
    !*FD initialise elastic lithosphere calculations
    use physcon, only : pi
    implicit none
    type(isostasy_elastic) :: rbel     !*FD structure holding elastic litho data    
    real(kind=dp), intent(in) :: deltax        !*FD grid spacing

    ! local variables
    real :: a     ! radius of disk
    real :: r     ! distance from centre
    integer :: i,j

    ! calculate a so that a circle of radius a is equivalent to a square with size deltax
    a = deltax/sqrt(pi)
    ! initialise w
    call init_rbel(rbel, a)
    ! calculate size of operator
    rbel%wsize = int(r_lr*rbel%lr/deltax)
    ! allocate memory for operator
    allocate(rbel%w(0:rbel%wsize,0:rbel%wsize))

    ! calculating points within disk
    rbel%w(0,0) = rbel_iw(rbel,0.)
    r = deltax/rbel%lr
    rbel%w(0,1) = rbel_iw(rbel,r)
    rbel%w(1,0) = rbel%w(0,1)
    ! calculating points outside disk
    do j=0,rbel%wsize
       do i=2,rbel%wsize
          r = deltax*sqrt(real(i)**2+real(j)**2)/rbel%lr
          rbel%w(i,j) = rbel_ow(rbel,r)
       end do
    end do
    do j=2,rbel%wsize
       do i=0,1
          r = deltax*sqrt(real(i)**2+real(j)**2)/rbel%lr
          rbel%w(i,j) = rbel_ow(rbel,r)
       end do
    end do
    i=1
    j=1
    r = deltax*sqrt(real(i)**2+real(j)**2)/rbel%lr
    rbel%w(i,j) = rbel_ow(rbel,r)
#ifdef DEB_REBOUND
    open(1,file='w.dat',status='UNKNOWN')
    do j=0,rbel%wsize
       do i=0,rbel%wsize
          write(1,*) i,j,rbel%w(i,j)
       end do
    end do
    close(1)
#endif
    !rbel%w=rbel%w/len0
  end subroutine init_elastic

  subroutine calc_elastic(rbel,load,load_factors)
    !*FD calculate surface loading effect using elastic lithosphere approximation
    implicit none
    type(isostasy_elastic) :: rbel                     !*FD structure holding elastic litho data
    real(kind=dp), dimension(:,:), intent(out) :: load !*FD loading effect due to load_factors
    real(kind=dp), dimension(:,:), intent(in)  :: load_factors !*FD load mass divided by mantle density

    integer ewn,nsn
    integer i,j,n,m

    ewn = size(load,1)
    nsn = size(load,2)

    load = 0.
    do j=1,nsn
       do i=1,ewn
          do n=max(1,j-rbel%wsize),min(nsn,j+rbel%wsize)
             do m=max(1,i-rbel%wsize),min(ewn,i+rbel%wsize)
                load(i,j) = load(i,j) + load_factors(m,n) * rbel%w(abs(m-i),abs(n-j))
             end do
          end do
       end do
    end do
  end subroutine calc_elastic

  subroutine finialise_elastic(rbel)
    !*FD clean-up data structure
    implicit none
    type(isostasy_elastic) :: rbel     !*FD structure holding elastic litho data    

    deallocate(rbel%w)
  end subroutine finialise_elastic

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! private subroutines
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine init_rbel(rbel, a)
    !*FD initialise elastic lithosphere calculations
    use paramets, only: len0
    use physcon, only: rhom,grav
    use kelvin
    implicit none
    type(isostasy_elastic) :: rbel     !*FD structure holding elastic litho data
    real, intent(in) :: a             !*FD radius of disk

    real dummy_a

    call set_kelvin(1.d-10,40)

    rbel%lr = ((rbel%d/(rhom*grav))**0.25)/len0
    rbel%a  = a

    dummy_a = rbel%a/rbel%lr
    
    rbel%c1  =  dummy_a * dker(dummy_a)
    rbel%c2  = -dummy_a * dkei(dummy_a)
    rbel%cd3 =  dummy_a * dber(dummy_a)
    rbel%cd4 = -dummy_a * dbei(dummy_a)
    
  end subroutine init_rbel

  function rbel_ow(rbel,r)
    use kelvin
    !*FD calculating deflection outside disk
    implicit none
    real :: rbel_ow
    real, intent(in) :: r             !*FD radius, r should be scaled with lr
    type(isostasy_elastic) :: rbel     !*FD structure holding elastic litho data
    
    rbel_ow = rbel%cd3*ker(r) + rbel%cd4*kei(r)
  end function rbel_ow

  function rbel_iw(rbel,r)
    use kelvin
    !*FD calculating deflection inside disk
    implicit none
    real :: rbel_iw
    real, intent(in) :: r             !*FD radius, r should be scaled with lr
    type(isostasy_elastic) :: rbel     !*FD structure holding elastic litho data
    
    rbel_iw = 1.0 + rbel%c1*ber(r) + rbel%c2*bei(r)
  end function rbel_iw
end module isostasy_el
