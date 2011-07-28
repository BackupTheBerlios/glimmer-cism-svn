module ssflux_h2o_spc
	use glimmer_global, only : dp,sp 
    use glimmer_paramets, only : thk0,len0
    use glide_thck
    use glide_types
	use glimmer_physcon, only : muw, rhow, grav, scyr
	use glimmer_scales, only : scale2d_f1

	contains
	
!-------------------------

	subroutine ssflux_h2o_SPC_init
	!use glimmer_global, only : dp 
	!initialize water routing subroutine
	implicit none
	end subroutine ssflux_h2o_spc_init
	
!-------------------------

	subroutine ssflux_h2o_spc_final
!	use glimmer_global, only : dp 
	!initialize water routing subroutine
	implicit none
	end subroutine ssflux_h2o_spc_final
	
!-------------------------

	subroutine ssflux_h2o_spc_driver(model,bmlt,bwat,thck,topg,btem,floater)
    !driver routine for the steady state H2o flux model (this is stable)

 !   use glide_temp
    implicit none
    type(glide_global_type),intent(inout) :: model
    real(dp), dimension(:,:), intent(out) :: bwat
    real(dp), dimension(:,:), intent(in) :: bmlt, thck, topg, btem
    logical, dimension(:,:), intent(in) :: floater
	call ssflux_h2o_spc_main(bmlt,bwat,thck,topg,btem,floater,model%numerics%dew,model%numerics%dns,&
	model%general%ewn,model%general%nsn)
    !print *, model%numerics%dew,model%numerics%dns
	end subroutine ssflux_h2o_spc_driver
	

!-------------------------

	subroutine ssflux_h2o_spc_main(bmlt,bwat,thck,topg,btem,floater,dx,dy,ewn,nsn)
!	!initialize water routing subroutine
!	use glimmer_global, only : dp
!	use glimmer_paramets, only : thk0, len0
	!use glimmer_routing
	!use glint_routing
	implicit none
	real(dp), dimension(:,:), intent(inout) :: bwat
    real(dp), dimension(:,:), intent(in) :: bmlt,thck,topg,btem
    real(dp), intent(in) :: dx, dy
    integer, intent(in) :: ewn, nsn
    logical, dimension(:,:), intent(in) :: floater
    real(dp),dimension(:,:),allocatable :: hpot,ssflux
    integer :: ns, ew
   ! print *, ewn, nsn
   !fldp_router(surface,input,output,mask,dx,dy)
    allocate(hpot(ewn,nsn))
    allocate(ssflux(ewn,nsn))
    !print *, 666
    do ns = 1,nsn
             do ew = 1,ewn
             !scale issue
             hpot(ew,ns)=(thck(ew,ns)*0.917+topg(ew,ns))*thk0
            ! print *, topg(31,ns), thck(ew,ns), hpot(ew,ns)
             end do
    end do 
    !print *, topg(50,31), thck(50,31), hpot(50,31)
    !print *, topg(11,31), thck(11,31), hpot(11,31)
    call fldp_router(hpot,bmlt,ssflux,floater,bwat,dx,dy)
    deallocate(hpot,ssflux)
	end subroutine ssflux_h2o_spc_main
	
!-------------------------

!  subroutine fldp_router(surface,input,output,mask,dx,dy)
   subroutine fldp_router(surface,input,output,mask,hsheet,dx,dy)
!    use glimmer_global, only : dp, sp
!	use glimmer_physcon, only : muw, rhow, grav, scyr
!	use glimmer_paramets, only : thk0, len0
!	use glimmer_scales, only : scale2d_f1
    !*FD Routes water from input field to its destination, 
    !*FD according to a surface elevation field. The method used 
    !*FD is by Quinn et. al. (1991)
    !designed for double precision numbers
    implicit none
    real(dp),dimension(:,:),intent(in)  :: surface !*FD hydopotential
    real(dp),dimension(:,:),intent(in)  :: input   !*FD Input water field
    real(dp),dimension(:,:),intent(inout) :: output, hsheet  !*FD Output water field
    logical, dimension(:,:),intent(in)  :: mask    !*FD Masked points!=.false. if matter, .true. ifthey dont
    real(dp),               intent(in)  :: dx      !*FD $x$ grid-length
    real(dp),               intent(in)  :: dy      !*FD $y$ grid-length
      ! Internal variables --------------------------------------
	
    integer :: nx,ny,k,nn,cx,cy,px,py,x,y,ned
    integer, dimension(:,:),allocatable :: sorted
    real(dp),dimension(:,:),allocatable :: flats,surfcopy
    real(dp),dimension(-1:1,-1:1) :: slopes
    real(dp),dimension(-1:1,-1:1) :: dists
    real(dp) :: sqtwo,sloptot,numdn,slpav, hspre3,hsh4now,turd!total slope for dh/dl calculation
    logical :: flag
	sqtwo=sqrt(dx**2+dy**2)*len0
	turd=0.333333333333333333333333333333333333333333333333333333333333333
	!print *, sqtwo
    ! Set up grid dimensions ----------------------------------
	!print *, dx, dy
    nx=size(surface,1) ; ny=size(surface,2)
    nn=nx*ny
	!scale issue
	!dists(-1,:)=(/4d0,2d0*dx/dy,4d0/)
    !dists(0,:)=(/2d0*dy/dx,0d0,2d0*dy/dx/)
    !this is a weird way to quantify distance
    dists(-1,:)=(/sqtwo,dy*len0,sqtwo/)
    dists(0,:)=(/dx*len0,0d0,dx*len0/)
    dists(1,:)=dists(-1,:)
    
    !print *, 4d0

    ! Allocate internal arrays and copy data ------------------

    allocate(sorted(nn,2),flats(nx,ny),surfcopy(nx,ny))
    surfcopy=surface

    ! Fill holes in data, and sort heights --------------------

    call phillholes(surfcopy,flats,mask)
    call heights_sort(surfcopy,sorted)

    ! Initialise output with input, which will then be --------
    ! redistributed -------------------------------------------

    output=input

    ! Begin loop over points, highest first -------------------

    do k=nn,1,-1
    
      ! Get location of current point -------------------------

      x=sorted(k,1)
      y=sorted(k,2)
      
      ! Reset flags and slope arrays --------------------------

      flag=.true.
      slopes=0.0

      ! Loop over adjacent points, and calculate slopes -------
      sloptot=0
      numdn=0
     ! print *, x,y,surfcopy(x,y)

      do cx=-1,1
        do cy=-1,1
          ! If this is the centre point, ignore
          if (cx==0.and.cy==0) continue
          ! Otherwise do slope calculation 
          px=x+cx ; py=y+cy
          if (px>0.and.px<=nx.and.py>0.and.py<=ny) then
              if (surfcopy(px,py)<surfcopy(x,y)) then
                 slopes(cx,cy)=(surfcopy(x,y)-surfcopy(px,py))/dists(cx,cy)
                sloptot=sloptot+slopes(cx,cy)**2
                numdn=numdn+slopes(cx,cy);
              endif
          endif
        enddo
      enddo
      slpav=sloptot/numdn
    !  if (x==11.and.y==31) then
!		print *, surfcopy(x,y)
!		print *, input(x,y)*scale2d_f1,output(x,y),input(x,y)*scale2d_f1*dx*dy*len0*len0/scyr+output(x,y)
!		do ned=-1,1
!		print *, slopes(ned,:)
!		enddo
!	endif
		!output=flux in + melt*dx*dy/secyr
      output(x,y)=output(x,y)+input(x,y)*scale2d_f1*dx*len0*dy*len0/scyr
      
      ! If there are places for the water to drain to, --------
      ! distribute it accordingly -----------------------------
      if (sum(slopes)/=0.0) then

        slopes=slopes/sum(slopes)
        do cx=-1,1
          do cy=-1,1
            px=x+cx ;py=y+cy
            if (slopes(cx,cy)/=0.0) then
              output(px,py)=output(px,py)+output(x,y)*slopes(cx,cy)
            endif
          enddo
        enddo
		!print *, slopes(cx,cy)
        ! Having distributed the water, zero the source -------

        !output(x,y)=0.0

      endif
      slpav=sloptot/numdn
      hspre3=(output(x,y)*12*muw/(((dx+dy)*len0/2)*rhow*grav*slpav))
      hsh4now=(abs(hspre3))**turd
      hsheet(x,y)=hsh4now/thk0
   	  !credit David g. simpson http://www.davidgsimpson.com/index.html
      if (x==11.and.y==31) then
		!print *, hspre3, hsh4now, hsheet(x,y)
		
	endif
      !hsheet(x,y)=output(x,y)
      !(nthroot(ssflxgrd(iy,ix)*12*Muw/(xstep*Rau_w*g*dhdl_grd(iy,ix)),3))


      ! End of main loop --------------------------------------
	
    enddo
    !print *, hsheet(11,31)

    ! Tidy up -------------------------------------------------

    deallocate(sorted,flats)

  end subroutine fldp_router

!==============================================================
! Internal subroutines
!==============================================================


  subroutine phillholes(phi,flats,mask)
!deals with the logical designation of floater versus the 1, 0 designation of mask
   ! use glimmer_global, only : dp
    implicit none

    real(dp),dimension(:,:),intent(inout) :: phi
    real(dp),dimension(:,:),intent(inout) :: flats
    logical, dimension(:,:),intent(in)    :: mask !=1 if points matter, .true. ifthey dont

    ! Internal variables --------------------------------------

    real(dp),allocatable,dimension(:,:) :: old_phi
    integer, allocatable,dimension(:,:) :: pool !identifies Pools!

    real(dp) :: pvs(9), max_val
    real(dp), parameter :: null = 1e+20
    integer :: flag,nx,ny,i,j

    ! ---------------------------------------------------------

    nx=size(phi,1) ; ny=size(phi,2)

    allocate(pool(nx,ny),old_phi(nx,ny))

    flag = 1

    ! ---------------------------------------------------------

    do while (flag .eq. 1)

       flag = 0

       old_phi = phi

       do i=2,nx-1
          do j=2,ny-1

             flats(i,j) = 0

             if (mask(i,j) .eqv. .false.) then
!mask = 1 & at least 1 surrounding point is lower 

                if (any(old_phi(i-1:i+1,j-1:j+1) < old_phi(i,j))) then
                   pool(i,j) = 0
                else
                   pool(i,j) = 1
                end if

                if (pool(i,j) .eq. 1) then

                   flag = 1

                   pvs = (/ old_phi(i-1:i+1,j-1), old_phi(i-1:i+1,j+1), old_phi(i-1:i+1,j) /)

                   where (pvs == old_phi(i,j))
                      pvs = null
                   end where

                   max_val = minval(pvs)

                   if (max_val .ne. null) then
                      phi(i,j) = max_val
                   else
                      flag = 0
                      flats(i,j) = 1
                   end if

                end if

             end if
          end do
       end do

    end do

    deallocate(pool,old_phi)

  end subroutine phillholes
!==============================================================

  subroutine heights_sort(surface,sorted)
 !   use glimmer_global, only : dp
    
    real(dp),dimension(:,:) :: surface
    integer,dimension(:,:) :: sorted

    integer :: nx,ny,nn,i,j,k
    real(dp),dimension(:),allocatable :: vect
    integer,dimension(:),allocatable :: ind

    nx=size(surface,1) ; ny=size(surface,2)
    nn=size(sorted,1)

    allocate(vect(nn),ind(nn)) 

    if (nn/=nx*ny.or.size(sorted,2).ne.2) then
      print*,'Wrong dimensions'
      stop
    endif

    k=1

    do i=1,nx
      do j=1,ny
        vect(k)=surface(i,j)
        k=k+1
      enddo
    enddo

    call indexx(vect,ind)

    do k=1,nn
      sorted(k,1)=floor(real(ind(k)-1)/real(ny))+1
      sorted(k,2)=mod(ind(k)-1,ny)+1
    enddo

    do k=1,nn
      vect(k)=surface(sorted(k,1),sorted(k,2))
    enddo
    
  end subroutine heights_sort

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !
  ! The following two subroutines perform an index-sort of an array. 
  ! They are a GPL-licenced replacement for the Numerical Recipes routine indexx. 
  ! They are not derived from any NR code, but are based on a quicksort routine by
  ! Michael Lamont (http://linux.wku.edu/~lamonml/kb.html), originally written
  ! in C, and issued under the GNU General Public License. The conversion to 
  ! Fortran 90, and modification to do an index sort was done by Ian Rutt.
  !
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine indexx(array,index)

    use glimmer_log
    use glimmer_global, only : dp 

    !*FD Performs an index sort of \texttt{array} and returns the result in
    !*FD \texttt{index}. The order of elements in \texttt{array} is unchanged.
    !*FD
    !*FD This is a GPL-licenced replacement for the Numerical Recipes routine indexx. 
    !*FD It is not derived from any NR code, but are based on a quicksort routine by
    !*FD Michael Lamont (http://linux.wku.edu/~lamonml/kb.html), originally written
    !*FD in C, and issued under the GNU General Public License. The conversion to 
    !*FD Fortran 90, and modification to do an index sort was done by Ian Rutt.

    real(dp),dimension(:) :: array !*FD Array to be indexed.
    integer, dimension(:) :: index !*FD Index of elements of \texttt{array}.
    integer :: i

    if (size(array).ne.size(index)) then
      call write_log('ERROR: INDEXX size mismatch.',GM_FATAL,__FILE__,__LINE__)
    endif

    do i=1,size(index)
       index(i)=i
    enddo

    call q_sort_index(array,index,1,size(array))

  end subroutine indexx

!==============================================================

  recursive subroutine q_sort_index(numbers,index,left,right)

    !*FD This is the recursive subroutine actually used by \texttt{indexx}. 
    !*FD
    !*FD This is a GPL-licenced replacement for the Numerical Recipes routine indexx. 
    !*FD It is not derived from any NR code, but are based on a quicksort routine by
    !*FD Michael Lamont (http://linux.wku.edu/~lamonml/kb.html), originally written
    !*FD in C, and issued under the GNU General Public License. The conversion to 
    !*FD Fortran 90, and modification to do an index sort was done by Ian Rutt.
    use glimmer_global, only : dp 

    implicit none

    real(dp),dimension(:) :: numbers !*FD Numbers being sorted
    integer, dimension(:) :: index   !*FD Returned index
    integer :: left, right           !*FD Limit of sort region

    integer :: ll,rr
    integer :: pv_int,l_hold, r_hold,pivpos
    real(dp) :: pivot

    ll=left
    rr=right

    l_hold = ll
    r_hold = rr
    pivot = numbers(index(ll))
    pivpos=index(ll)

    do
       if (.not.(ll < rr)) exit

       do 
          if  (.not.((numbers(index(rr)) >= pivot) .and. (ll < rr))) exit
          rr=rr-1
       enddo

       if (ll.ne.rr) then
          index(ll) = index(rr)
          ll=ll+1
       endif

       do
          if (.not.((numbers(index(ll)) <= pivot) .and. (ll < rr))) exit
          ll=ll+1
       enddo

       if (ll.ne.rr) then
          index(rr) = index(ll)
          rr=rr-1
       endif
    enddo

    index(ll) = pivpos
    pv_int = ll
    ll = l_hold
    rr = r_hold
    if (ll < pv_int)  call q_sort_index(numbers, index,ll, pv_int-1)
    if (rr > pv_int)  call q_sort_index(numbers, index,pv_int+1, rr)

  end subroutine q_sort_index

end module ssflux_h2o_spc