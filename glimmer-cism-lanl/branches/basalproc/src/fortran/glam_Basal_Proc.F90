module glam_Basal_Proc


use glide_types
use glimmer_paramets, only : dp,sp
use glimmer_physcon,  only : grav, rhow, rhos, scyr
use glimmer_paramets, only : vel0
use glimmer_log,      only : write_log

  implicit none;save

!!Variables
  real (kind = dp), dimension(:,:,:), allocatable:: u,N,etill,dy
  integer, parameter :: unin = 90
  
contains
  
  subroutine Basal_Proc_init(ewn,nsn,basalproc,hotstart,ntem,tillfile)
    implicit none
    
    !Arguments
    integer,intent(in) :: ewn,nsn,hotstart
	real (kind = sp),intent (in) :: ntem
    type(glide_basalproc),intent(inout) :: basalproc
    character (len=30), intent(in), optional :: tillfile
    
    !Variables
    real (kind = dp), dimension (ewn,nsn,basalproc%tnodes) :: por,Ztot
    real (kind = dp), dimension (ewn-1,nsn-1) :: stagHwater
    real (kind = sp), dimension (ewn,nsn) :: TillMap
    integer :: x,y,i
    character(len=512) :: message

!allocate basal processes variables
	allocate (dy(ewn-1,nsn-1,basalproc%tnodes-1));		dy=0.6d0
	allocate (u(ewn-1,nsn-1,basalproc%tnodes));			u=41.0d3
	allocate (N(ewn-1,nsn-1,basalproc%tnodes));			N=11.0d3
	allocate (etill(ewn-1,nsn-1,basalproc%tnodes));		etill=0.5d0




    if (hotstart.eq.1) then
       !From restart file, the following variables are known: dy, u, N and minTauf and etill
  !     basalproc%minTauf=5000  !Dummy value for now
       por=etill/(1+etill)
       stagHwater=0.0
       do i=2,basalproc%tnodes
          stagHwater(:,:)=stagHwater+dy(:,:,i-1)*(por(:,:,i-1)+por(:,:,i))/2
       enddo

       write(message,*) 'Till layer has been initialized from restart file'
       call write_log(message)
       
       
    else if (hotstart.eq.0) then

          open (unit=unin,file=tillfile,access='direct',form='unformatted',recl=ewn*nsn*4)
          read(unin,rec=1) TillMap
          close(unin)
          basalproc%minTauf=dble(TillMap)
          
     
       
       N(:,:,1)=basalproc%minTauf/basalproc%fric
       do i=2,basalproc%tnodes
          N(:,:,i)=N(:,:,1)
       end do
       
       etill=basalproc%etillo-basalproc%Comp*log10(N/basalproc%No)
       por=etill/(1+etill)
       dy=basalproc%Zs*(1+etill(:,:,1:basalproc%tnodes-1))/(basalproc%tnodes-1) 
       Ztot(:,:,1)=0.0
       stagHwater=0.0
       do i=2,basalproc%tnodes
          Ztot(:,:,i)=Ztot(:,:,i-1)+dy(:,:,i-1)
          stagHwater(:,:)=stagHwater+dy(:,:,i-1)*(por(:,:,i-1)+por(:,:,i))/2
       enddo
     
       u=(rhos-rhow)*(1-por)*grav*Ztot-N

         write(message,*) 'Till layer has been initialized using namelist file: ',tillfile
       	 call write_log(message)       
    end if
    
    	!Calculate Hwater on normal grid - using zero gradient as BC
		call stag2norm(ewn,nsn,stagHwater,basalproc%Hwater)	    
    
  end subroutine Basal_Proc_init
  
  
  subroutine Basal_Proc_driver (ewn,      nsn,      upn,  				&
								dt,ubas,vbas,  what, bmlt, basalproc)
								
	use glide_grids, only: stagvarb								
  implicit none
    
    !Arguments
    integer, intent (in) ::ewn, nsn, upn, what
    real (kind = sp), intent(in) :: dt
    real (kind = dp), dimension(:,:), intent (in) :: ubas,vbas
    real (kind = dp), dimension(:,:), intent (in) :: bmlt
    type(glide_basalproc),intent(inout) :: basalproc
    
    !Variables
    real (kind = dp), dimension (ewn-1,nsn-1) :: Ub,stagHwater,stagbmlt
    real (kind = dp) :: f1
    integer :: i

    !Calculate basal melt rate on staggered grid
  	call stagvarb(bmlt, stagbmlt,ewn,nsn)
      
!    !Calculate the magnitude of basal velocity, in m/yr
	Ub=scyr*vel0* (sqrt(ubas(:,:)**2+vbas(:,:)**2))

!	print*,'begin basal_proc_driver'
!	print*,'mean Tauf=',sum(basalproc%minTauf)/((ewn-1)*(nsn-1))
!	print*,'mean bmlt=',sum(stagbmlt)/((ewn-1)*(nsn-1))
!	print*,'mean Ub=',sum(Ub)/((ewn-1)*(nsn-1))
	
    select case(what)

    case(1)
    call Till_FullRes  (ewn,nsn,dt,stagbmlt,Ub,basalproc%tnodes,			&
    					basalproc%Kh,basalproc%Cv, basalproc%etillo,        &
    					basalproc%fric, basalproc%Comp, basalproc%No,  		&
    					basalproc%Zs, stagHwater,basalproc%minTauf)

    case(2)
    call Till_FastCalc (dt,stagbmlt,basalproc%aconst,basalproc%bconst,      &
    					basalproc%Zs,basalproc%minTauf,stagHwater) 
    
    end select

	!Calculate Hwater on normal grid - using zero gradient as BC
	call stag2norm(ewn,nsn,stagHwater,basalproc%Hwater)
	

 
    
  end subroutine Basal_Proc_driver
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  

  
  subroutine Basal_Proc_final

  ! Deallocate till variables
    deallocate(u)
    deallocate(N)
    deallocate(etill)
    deallocate(dy)
    
  end subroutine Basal_Proc_final
  
  
  
  subroutine stag2norm(ewn,nsn,stagvar,normvar)
  use glimmer_paramets, only : dp	
  implicit none
  
  integer, intent(in) :: ewn,nsn
  real (kind = dp), intent(in), dimension(:,:) :: stagvar
  real (kind = dp), intent(out), dimension (:,:) :: normvar

    normvar(1:ewn,1:nsn)=0.0
    normvar(2:ewn-1,2:nsn-1)=(stagvar(1:ewn-2,2:nsn-1)+stagvar(2:ewn-1,2:nsn-1)+ & 
    						stagvar(1:ewn-2,1:nsn-2)+stagvar(2:ewn-1,1:nsn-2)) / 4.0d0

    !Apply zero-gradient to the velocity field on normal grid - using swapbnmelt
    call swapbnmelt(0,size(normvar(1,:)),normvar(1,:),normvar(2,:),normvar(ewn,:),normvar(ewn-1,:))
    call swapbnmelt(0,size(normvar(:,1)),normvar(:,1),normvar(:,2),normvar(:,nsn),normvar(:,nsn-1))

  end subroutine stag2norm	

  subroutine swapbnmelt(bc,lgth,a,b,c,d)
	use glimmer_paramets, only : dp        
    implicit none
    
    integer, intent(in) :: bc,lgth
    real (kind = dp), intent(in), dimension(lgth) :: b, d
    real (kind = dp), intent(out), dimension (lgth):: a, c
    
    if (bc==0) then
       a = b
       c = d
    end if
    return
    
  end subroutine swapbnmelt


subroutine Till_FullRes(ewn,nsn,dt,bdot,Ub,						&
						tnodes,Kh,Cv,etillo,fric,Comp,No,Zs,		&
						Hwater,minTauf)
 !Largely follows Bougamont et al. 2003 (JGR), with explicit till layer properties calculations 
 !Includes vertical mixing as in Christoffersen et al. 2003
	
	use glimmer_paramets, only : dp,sp
	use glimmer_physcon,  only : grav, rhow, rhos, scyr

  implicit none
  
  integer, intent(in) :: ewn,nsn,tnodes
  real (kind=sp), intent(in):: dt
  real (kind=dp),intent(in) ::Kh,Cv,etillo,fric,Comp,No,Zs
  real (kind = dp), dimension(:,:),intent(in) :: bdot ! basal melt rate m.yr
  real (kind = dp), dimension(:,:),intent(in) :: Ub   ! velocity m/yr
 real (kind = dp), dimension(:,:),intent(out) :: minTauf   ! 
 real (kind = dp), dimension(:,:),intent(out) :: Hwater   ! 

  !Local variables
  real, parameter :: f=1e-3
  real (kind = dp), dimension(ewn-1,nsn-1,tnodes) :: uold,Tauf,deltaU,por
  real (kind = dp), dimension(ewn-1,nsn-1,tnodes-1) :: vw   !vertical water flow, m/s
  real (kind = dp), dimension(ewn-1,nsn-1) :: du  
  integer :: i

  !Boundary conditions at the ice/till interface
  
  uold=u 
  du=dy(:,:,1)*(bdot/scyr)*rhow*grav/Kh  ! Removed *dt so that du has indeed unit Pa
  u(:,:,1)=uold(:,:,1) + (Cv*scyr*dt/dy(:,:,1)**2)*(uold(:,:,2)-uold(:,:,1)+du) + &
       f*Ub*dt/(dy(:,:,1))*du
  
  i=2
  do while (i.lt.tnodes)
     u(:,:,i)=uold(:,:,i) + (Cv*scyr*dt/(dy(:,:,i)**2))*(uold(:,:,i+1)-2*uold(:,:,i)+uold(:,:,i-1)) + &
          (f*Ub(:,:)*dt/dy(:,:,i))*(uold(:,:,i-1)-uold(:,:,i))
     i=i+1
  end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Set up two possible lower boundary conditions:

  ! 1. ZERO-FLUX assuming an imaginary nodes at (tnodes+1) 
  u(:,:,tnodes)=uold(:,:,tnodes)+(Cv*scyr*dt/(dy(:,:,tnodes-1)**2))*(uold(:,:,tnodes-1)-uold(:,:,tnodes)) + &
          (f*Ub*dt/dy(:,:,tnodes-1))*(uold(:,:,tnodes-1)-uold(:,:,tnodes))

  ! 2. CONSTANT FLUX assuming an imaginary nodes at (tnodes+1)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!Calculate vertical velocity from darcy equation:
  vw(:,:,1:tnodes-1)=Kh/(rhow*grav)*(1/dy(:,:,1:tnodes-1))*(u(:,:,1:tnodes-1)-u(:,:,2:tnodes))

!From there, new void ratio distribution
!if vw > 0, then water flowing downward, and added to till layer (increases etill)
  etill(:,:,1:tnodes-1)=etill(:,:,1:tnodes-1)+(vw*scyr)*dt/dy

  where (etill.lt.0.15)
     etill=0.15
  end where
  
  Tauf=No*fric*10**(-(etill-etillo)/Comp)
  minTauf=minval(Tauf(:,:,1:tnodes-1),3)  
  dy(:,:,1:tnodes-1)=(Zs/(tnodes-1))*(1+(etill(:,:,1:tnodes-1)+etill(:,:,2:tnodes))/2)
  por=etill/(1+etill)

  Hwater=0.
  do i=2,tnodes
     Hwater=Hwater+dy(:,:,i-1)*((por(:,:,i-1)+por(:,:,i))/2)
  enddo


end subroutine Till_FullRes

subroutine Till_FastCalc(dt,bdot,aconst,bconst,Zs,minTauf,Hwater) 
 !Largely follows Bougamont et al. 2003 (JGR)
 !>>>>>>>Till layer is represented with only one node!!  

  use glimmer_paramets, only : dp,sp
  
  implicit none
  
  !Arguments
  real (kind=sp), intent(in):: dt
  real (kind=dp), intent(in):: aconst,bconst,Zs
  real (kind = dp), dimension(:,:),intent(in) :: bdot ! basal melt rate m.yr
  real (kind = dp), dimension(:,:),intent(out) :: minTauf ! basal melt rate m.yr
  real (kind = dp), dimension(:,:),intent(out) :: Hwater ! basal melt rate m.yr
 
  etill(:,:,1)=etill(:,:,1)+dt*bdot/Zs
  where (etill.lt.0.15)
     etill=0.15
  end where

  minTauf=aconst*exp(-bconst*etill(:,:,1))
  Hwater=Zs*(etill(:,:,1)/(1+etill(:,:,1)))


end subroutine Till_FastCalc

 
end module glam_Basal_Proc

  
