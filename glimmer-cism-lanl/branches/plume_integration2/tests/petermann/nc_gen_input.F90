program nc_gen_input

  use netcdf

  implicit none

  real,parameter :: pi = 3.1415926
  real,parameter :: g = 9.81
  !real,parameter :: rhoi_rhow = 1030.d0 / 920.d0

  character(len=512) :: gen_usage = "nc_gen_input <type_code> &
       &nc_filename [params] &
       &<type_code> = [cs|ts|gs|ss]"

  character(len=512) :: conf_shelf_params = "<fname> <nx> <ny> <hx> <hy> &
       &<slope_start_pos> <grounded_ice_thk> &
       &<ice_front_pos> <ice_front_thick> &
       &<ocean_topg> <land_topg> <k_x> & 
       &<chan_amp> <chan_init_length> & 
       &<kinbc_width> <nsswitch>"

  character(len=512) :: ghost_shelf_params = "<fname> <nx> <ny> <kinbcw> <n_level> <hx> <hy> &
       &<upstream_thk> <upstream_vel> <front_thk> &
       &<k_x> <chan_amp> <ramp_width> <landw> <ifpos>"


  character (len=2) :: type_code

  !variables to be used in all cases and in netcdf writing
  character(len=1024) :: argstr
  character (len=2056) :: fname
  integer :: nx,ny,n_level

  !data arrays
  real,dimension(:),allocatable :: xs,ys,xstag,ystag,level
  real,dimension(:,:),allocatable :: thck,topog
  real,dimension(:,:),allocatable :: kinbcmask
  real,dimension(:,:,:),allocatable :: uvelhom,vvelhom

  !get arguments from command line

  if (command_argument_count() < 1) then

     write(*,*) "Usage: ", trim(gen_usage)
     stop 1

  end if

  call get_command_argument(1,argstr)

  if (trim(argstr) == '-h' .or. trim(argstr) == '--help') then
     write(*,*) "Usage: ", trim(gen_usage)
     stop
  end if

  type_code = argstr

  if (type_code == 'cs') then

     call make_confined_shelf()

  else if (type_code =='ts') then

     call make_two_sided_shelf()

  else if (type_code == 'ls') then
  
     call make_linear_shelf()

  else if (type_code == 'ss') then

     call make_steady_shelf()

  else if (type_code == 'gs') then

     call make_ghost_shelf()

  else

     write(*,*) "Unrecognized type_code: ", type_code

  end if

contains

  subroutine make_ghost_shelf()

    ! <fname> <nx> <ny> <kinbcw> <n_level> <hx> <hy> 
    ! <upstream_thk> <front_thk> 
    ! <k_x> <chan_amp> <ramp_width> <landw> <ifpos>

    ! local variables
    integer :: landw, landw_north 
    integer :: kinbcw 
    integer :: ifpos

    integer :: i,j
    integer :: firstarg = 2
    integer :: ramp_width
    real :: hx,hy,up_thk,if_thk
    real :: k_x,chan_amp,chan_depth
    real :: upstream_vel
    real :: otopg,ltopg
    real :: cutoff_usrf,rhoi_rhow,central_amp

    if (command_argument_count() /= 16) then
       write(*,*)"Incorrect number of parameters. Ghost shelf requires:  ",trim(ghost_shelf_params)
       stop 1
    end if

    call get_command_argument(firstarg,argstr)
    read(argstr,'(a512)') fname
    write(*,*) 'fname ',trim(fname)

    call get_command_argument(firstarg+1, argstr)
    read(argstr,'(i5)') nx
    write(*,*) 'nx:',nx

    call get_command_argument(firstarg + 2,argstr)
    read(argstr,'(i5)') ny
    write(*,*) 'ny:',ny

    call get_command_argument(firstarg + 3,argstr)
    read(argstr,'(i5)') kinbcw
    write(*,*) 'kinbcw:', kinbcw

    call get_command_argument(firstarg + 4,argstr)
    read(argstr,'(i5)') n_level
    write(*,*) 'n_level:',n_level

    call get_command_argument(firstarg + 5,argstr)
    read(argstr,'(f18.12)') hx
    write(*,*) 'hx',hx

    call get_command_argument(firstarg + 6,argstr)
    read(argstr,'(f18.12)') hy
    write(*,*) 'hy',hy

    call get_command_argument(firstarg + 7,argstr)
    read(argstr,'(f18.12)') up_thk
    write(*,*) 'up_thk',up_thk

    call get_command_argument(firstarg + 8,argstr)
    read(argstr,'(f18.12)') upstream_vel
    write(*,*) 'upstream_vel',upstream_vel

    call get_command_argument(firstarg + 9,argstr)
    read(argstr,'(f18.12)') if_thk
    write(*,*) 'if_thk',if_thk

    call get_command_argument(firstarg + 10,argstr)
    read(argstr,'(f18.12)') k_x
    write(*,*) 'k_x',k_x

    call get_command_argument(firstarg + 11,argstr)
    read(argstr,'(f18.12)') chan_amp
    write(*,*) 'chan_amp',chan_amp

    call get_command_argument(firstarg + 12,argstr)
    read(argstr,'(i5)') ramp_width
    write(*,*) 'ramp_width',ramp_width

    call get_command_argument(firstarg+13,argstr)
    read(argstr,'(i5)') landw
    write(*,*) 'landw', landw
    
    call get_command_argument(firstarg+14,argstr)
    read(argstr,'(i5)') ifpos
    write(*,*) 'ifpos', ifpos

    allocate(xs(nx),ys(ny),xstag(nx-1),ystag(ny-1))
    allocate(level(n_level))
    allocate(topog(nx,ny),thck(nx,ny))

    allocate(kinbcmask(nx,ny-1))

    allocate(uvelhom(nx-1,ny-1,n_level), vvelhom(nx-1,ny-1,n_level))

    !now populate the dimension variables
    xs = (/ ( (i-1)*hx,i=1,nx ) /)
    ys = (/ ( (j-1)*hy,j=1,ny ) /)
    level = (/ ( real(i)/real((n_level-1)), i=0,(n_level-1) ) /)
    xstag = (/ ( ((real(i)-0.5)*hx),i=1,nx-1 ) /)
    ystag = (/ ( ((real(j)-0.5)*hy),j=1,ny-1 ) /)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! define the topograph, thickness, kinbcmask       !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! safe values so the ice never touches the bottom or breaches the land

    landw_north = 0

    otopg = -up_thk * 2.0
    ltopg = up_thk * 1.0


    rhoi_rhow = 910.0 / 1028.0
    central_amp = 0.0

    !define topg
    topog = -abs(otopg)
    topog(:,(ny+1-landw_north):ny) = abs(ltopg)   ! North edge is flat and high

    if (landw > 0) then
        topog(1:landw,:) = abs(ltopg)
        topog((nx+1-landw):nx,:) = abs(ltopg)
    
        do j=1,(ny-landw_north)
            topog(landw:(ramp_width-1+landw),j) = (/ ( ltopg +(otopg-ltopg)*i/(ramp_width-1),i=0,(ramp_width-1) ) /)
            topog((nx-landw+2-ramp_width):(nx-landw+1),j) = (/ ( otopg +(ltopg-otopg)*i/(ramp_width-1),i=0,(ramp_width-1) ) /)
        end do
    end if

    !define kinbcmask
    kinbcmask = 0
    kinbcmask(:,(ny-kinbcw):(ny-1)) = 1  !north edge has imposed velocity
    
    !define thickness

    do j=ifpos,ny
       thck(:,j) = if_thk + &
            (up_thk - if_thk)*(real(j-ifpos)/real(ny-ifpos))
    end do

    do i=1,nx
       do j=ifpos,ny
          chan_depth = 	0.5 * chan_amp &
               * (1-cos((real(i-1)/real(nx-1))*2*pi*k_x))
          thck(i,j) = thck(i,j) - chan_depth
       end do
    end do

    thck(:,1:(ifpos-1)) = 0.d0 ! zero ahead of ice front position

    do j= ifpos,(ny - landw_north)

       cutoff_usrf = maxval(thck((ramp_width+1):(nx+1-ramp_width),j)*(1 - rhoi_rhow))

       do i=1,nx
          thck(i,j) = max( thck(i,j) - &
                           max( thck(i,j) + topog(i,j) - cutoff_usrf, 0.0 ), &
                           0.0 )
       end do
    end do

    uvelhom = 0.0
    vvelhom = 0.0

    do i = 1,(nx-1)
        vvelhom(i,(ny-kinbcw):(ny-1),:) = upstream_vel* &
                   (1.0 + central_amp*4.0*real(i-landw)/real(nx-1-2*landw) * real(nx-landw-i)/real(nx-1-2*landw))
    end do
!    vvelhom = upstream_vel

    call write_nc_file(.true.)

    deallocate(xs,ys,level,xstag,ystag,thck,topog,kinbcmask,uvelhom,vvelhom)

  end subroutine make_ghost_shelf

  subroutine make_confined_shelf()

    ! local variables
    integer :: landw = 2
    integer :: i,j
    integer :: firstarg = 3
    integer :: startpos,ifpos,kinbcw
    real :: hx,hy,githk,ifthk,kx,chan_amp,chan_init_length
    real :: chan_depth,otopg,ltopg
    logical :: ns_switch

    if (command_argument_count() /= 17) then
       write(*,*)"Incorrect number of parameters. Confined shelf requires:  ",trim(conf_shelf_params)
       stop 1
    end if

    call get_command_argument(2,argstr)
    read(argstr,'(a512)') fname
    write(*,*) 'fname ',trim(fname)

    call get_command_argument(firstarg, argstr)
    read(argstr,'(i5)') nx
    write(*,*) 'nx:',nx

    call get_command_argument(firstarg + 1,argstr)
    read(argstr,'(i5)') ny
    write(*,*) 'ny:',ny

    call get_command_argument(firstarg + 2,argstr)
    read(argstr,'(f18.12)') hx
    write(*,*) 'hx',hx

    call get_command_argument(firstarg + 3,argstr)
    read(argstr,'(f18.12)') hy
    write(*,*) 'hy',hy

    call get_command_argument(firstarg + 4,argstr)
    read(argstr,'(i5)') startpos
    write(*,*) 'startpos',startpos

    call get_command_argument(firstarg + 5,argstr)
    read(argstr,'(f18.12)') githk
    write(*,*) 'githk',githk

    call get_command_argument(firstarg + 6,argstr)
    read(argstr,'(i5)') ifpos
    write(*,*) 'ifpos',ifpos

    call get_command_argument(firstarg + 7,argstr)
    read(argstr,'(f18.12)') ifthk
    write(*,*) 'ifthk',ifthk

    call get_command_argument(firstarg + 8,argstr)
    read(argstr,'(f18.12)') otopg
    write(*,*) 'otopg', otopg

    call get_command_argument(firstarg + 9,argstr)
    read(argstr,'(f18.12)') ltopg
    write(*,*) 'ltopg', ltopg

    call get_command_argument(firstarg + 10,argstr)
    read(argstr,'(f18.12)') kx
    write(*,*) 'kx',kx

    call get_command_argument(firstarg + 11,argstr)
    read(argstr,'(f18.12)') chan_amp
    write(*,*) 'chan_amp',chan_amp

    call get_command_argument(firstarg + 12,argstr)
    read(argstr,'(f18.12)') chan_init_length
    write(*,*) 'chan_init_length',chan_init_length  

    call get_command_argument(firstarg + 13,argstr)
    read(argstr,'(i5)') kinbcw
    write(*,*) 'kinbc_width',kinbcw

    call get_command_argument(firstarg + 14,argstr)
    read(argstr,'(l1)') ns_switch
    write(*,*) 'ns_switch',ns_switch

    allocate(xs(nx),ys(ny),xstag(nx-1),ystag(ny-1))
    allocate(topog(nx,ny),thck(nx,ny),kinbcmask(nx-1,ny-1))

    !now populate the dimension variables
    xs = (/ ( (i-1)*hx,i=1,nx ) /)
    ys = (/ ( (j-1)*hy,j=1,ny ) /)
    xstag = (/ ( ((real(i)-0.5)*hx),i=1,nx-1 ) /)
    ystag = (/ ( ((real(j)-0.5)*hy),j=1,ny-1 ) /)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! define the topograph, thickness, kinbcmask       !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (ns_switch) then
       !define topg
       topog = -abs(otopg)
       topog(1:landw,:) = abs(ltopg)
       topog((nx+1-landw):nx,:) = abs(ltopg)
       topog(:,1:landw) = abs(ltopg)

       !define thickness
       thck(:,1:startpos) = githk 
       do j=startpos+1,ifpos
          if (ifpos /= startpos) then
             thck(:,j) = githk + &
                  (ifthk -githk)*(real(j-startpos)/real(ifpos-startpos))
          end if
       end do
       do i=1,nx
          do j=startpos+1,ifpos
             chan_depth = 	0.5 * chan_amp &
                  * (1-exp(- real(j-startpos)*hy/chan_init_length))&
                  * (1-cos((real(i-1)/real(nx-1))*2*pi*kx))
             thck(i,j) = thck(i,j) - chan_depth
          end do
       end do

       thck(1:landw,:) = 0.d0
       thck((nx+1-landw):nx,:) = 0.d0
       thck(:,1:landw) = 0.d0
       thck(:,(ifpos+1):ny) = 0.d0 ! zero ahead of ice front position

       ! define kinbcmask
       kinbcmask = 0
       kinbcmask(1:kinbcw,:) = 1            !west
       kinbcmask((nx-kinbcw):(nx-1),:) = 1  !east
       kinbcmask(:,1:kinbcw) = 1 !south

    else

       !define topg
       topog = -abs(otopg)
       topog(:,(ny+1-landw):ny) = abs(ltopg)   !north edge
       topog((nx+1-landw):nx,:) = abs(ltopg)   !east
       topog(1:landw,:) = abs(ltopg)           !west

       !define thickness
       thck(:,startpos:ny) = githk 
       do j=ifpos,(startpos-1)
          thck(:,j) = githk + &
               (ifthk -githk)*(real(startpos-j)/real(startpos-ifpos))
       end do
       do i=1,nx
          do j=ifpos,startpos-1
             chan_depth = 	0.5 * chan_amp &
                  * (1-exp(- real(startpos-j)*hy/chan_init_length))&
                  * (1-cos((real(i-1)/real(nx-1))*2*pi*kx))
             thck(i,j) = thck(i,j) - chan_depth
          end do
       end do

       thck(1:landw,:) = 0.d0            !west edge
       thck((nx+1-landw):nx,:) = 0.d0    !east
       thck(:,(ny+1-landw):ny) = 0.d0    !north
       thck(:,1:(ifpos-1)) = 0.d0 ! zero ahead of ice front position

       ! define kinbcmask
       kinbcmask = 0
       kinbcmask(1:kinbcw,:) = 1           !west
       kinbcmask((nx-kinbcw):(nx-1),:) = 1 !east
       kinbcmask(:,(ny-kinbcw):(ny-1)) = 1   !north

    end if
    call write_nc_file(.false.)

    deallocate(xs,ys,xstag,ystag,thck,topog,kinbcmask)

  end subroutine make_confined_shelf

  subroutine make_steady_shelf()

    ! local variables
    integer :: i,j
    integer :: ifpos,kinbcw
    integer,parameter :: zero_buf = 1
    real :: hx,hy
    real :: rhoi, rhoo, glen_A
    real :: upstream_vel, upstream_thk
    real :: kx,chan_amp,chan_depth,chan_init_length
    real :: otopg,ltopg,acab_per_year
    real,dimension(:),allocatable :: rand_row

    real :: k,tmp
    real :: rand_amp

    character (len=512) :: steady_shelf_params = "<fname> <nx> <ny> <n_level> <hx> <hy> <upstream_thk> &
       & <upstream_vel> <ifpos> <ocean_depth> <kinbcw> <acab_per_year> &
       & <rhoi> <rhoo> <glen_A> <rand_amp> <kx> <chan_amp> <chan_init_length>"

    if (command_argument_count() /= 20) then
       write(*,*)"Incorrect number of parameters. Steady shelf requires: &
            &  ",trim(steady_shelf_params)
       stop 1
    end if

    call get_command_argument(2,argstr)
    read(argstr,'(a512)') fname
    write(*,*) 'fname ',trim(fname)

    call get_command_argument(3, argstr)
    read(argstr,'(i5)') nx
    write(*,*) 'nx:',nx

    call get_command_argument(4,argstr)
    read(argstr,'(i5)') ny
    write(*,*) 'ny:',ny

    call get_command_argument(5,argstr)
    read(argstr,'(i5)') n_level
    write(*,*) 'n_level:',n_level

    call get_command_argument(6,argstr)
    read(argstr,'(f18.12)') hx
    write(*,*) 'hx',hx

    call get_command_argument(7,argstr)
    read(argstr,'(f18.12)') hy
    write(*,*) 'hy',hy

    call get_command_argument(8,argstr)
    read(argstr,'(f18.12)') upstream_thk
    write(*,*) 'upstream_thk',upstream_thk

    call get_command_argument(9,argstr)
    read(argstr,'(f18.12)') upstream_vel
    write(*,*) 'upstream_vel',upstream_vel

    call get_command_argument(10,argstr)
    read(argstr,'(i5)') ifpos
    write(*,*) 'ifpos',ifpos

    call get_command_argument(11,argstr)
    read(argstr,'(f18.12)') otopg
    write(*,*) 'otopg', otopg

    call get_command_argument(12,argstr)
    read(argstr,'(i5)') kinbcw
    write(*,*) 'kinbc_width',kinbcw

    call get_command_argument(13, argstr)
    read(argstr,'(f18.12)') acab_per_year
    write(*,*) 'acab_per_year', acab_per_year

    call get_command_argument(14, argstr)
    read(argstr,'(f18.12)') rhoi
    write(*,*) 'rhoi', rhoi

    call get_command_argument(15, argstr)
    read(argstr,'(f18.12)') rhoo
    write(*,*) 'rhoo', rhoo

    call get_command_argument(16, argstr)
    read(argstr,'(f18.12)') glen_A
    write(*,*) 'glen_A', glen_A

    call get_command_argument(17, argstr)
    read(argstr,'(f18.12)') rand_amp
    write(*,*) 'rand_amp', rand_amp

    call get_command_argument(18, argstr)
    read(argstr,'(f18.12)') kx
    write(*,*) 'kx', kx

    call get_command_argument(19, argstr)
    read(argstr,'(f18.12)') chan_amp
    write(*,*) 'chan_amp', chan_amp

    call get_command_argument(20, argstr)
    read(argstr,'(f18.12)') chan_init_length
    write(*,*) 'chan_init_length', chan_init_length

    allocate(xs(nx),ys(ny),level(n_level),xstag(nx-1),ystag(ny-1))
    allocate(topog(nx,ny),thck(nx,ny),kinbcmask(nx-1,ny-1))
    allocate(uvelhom(nx-1,ny-1,n_level), vvelhom(nx-1,ny-1,n_level))
    allocate(rand_row(nx-2*zero_buf))

    !now populate the dimension variables
    xs = (/ ( (i-1)*hx,i=1,nx ) /)
    ys = (/ ( (j-1)*hy,j=1,ny ) /)
    level = (/ ( real(i)/real((n_level-1)), i=0,(n_level-1) ) /)
    xstag = (/ ( ((real(i)-0.5)*hx),i=1,nx-1 ) /)
    ystag = (/ ( ((real(j)-0.5)*hy),j=1,ny-1 ) /)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! define the topography, thickness, kinbcmask       !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !define topg
    topog = -abs(otopg)
    
    !define thickness

    thck = 0.0
    thck(1+zero_buf:nx-zero_buf,(ny-kinbcw):ny) = upstream_thk 
    
    k = (rhoi*g*0.25*(1.0-rhoi/rhoo))**3.0

    if (acab_per_year .eq. 0.0) then

       do j=ifpos,(ny-kinbcw + 1)       
          tmp = upstream_vel ** 4.0 - 4.0 * glen_A * k * &
               (upstream_thk*upstream_vel)**3.0 * &
               (ny - kinbcw + 1 - j)*hy
          tmp = tmp ** 0.25  
          thck(1+zero_buf:nx-zero_buf,j) = abs(upstream_vel*upstream_thk) / tmp
          call random_number(rand_row)
          thck(1+zero_buf:nx-zero_buf,j) = thck(1+zero_buf:nx-zero_buf,j) * &
	                                   (1+rand_amp*(0.5-rand_row))
       end do

    else

      do j=ifpos,(ny-kinbcw)
         tmp = upstream_vel ** 4.0 - (glen_A/acab_per_year)*k* &
              ((upstream_vel*upstream_thk)**4.0 - (upstream_vel*upstream_thk - &
                                                   acab_per_year*(ny-kinbcw+1-j)*hy)**4.0)
         tmp = tmp ** 0.25
         thck(1+zero_buf:nx-zero_buf,j) = abs(upstream_vel*upstream_thk - &
                                          acab_per_year*(ny-kinbcw+1-j)*hy)/tmp
         call random_number(rand_row)
         thck(1+zero_buf:nx-zero_buf,j) = thck(1+zero_buf:nx-zero_buf,j)* &
                                          (1+rand_amp*(0.5-rand_row))

      end do

    end if

    do i=1,nx
    	do j=ifpos,(ny-kinbcw-1)
		chan_depth = 0.5 * chan_amp &
	                  * (1-exp( real(j-(ny-kinbcw))*hy/chan_init_length))&
        	          * (1-cos((real(i-1)/real(nx-1))*2*pi*kx))
             	thck(i,j) = thck(i,j) - chan_depth
          end do
       end do

    ! define kinbcmask
    kinbcmask = 0
    kinbcmask(:,(ny-kinbcw):(ny-1)) = 1 !north edge
    
    uvelhom = 0.0
    vvelhom = 0.0

    vvelhom(:,(ny-kinbcw):(ny-1),:) = upstream_vel

    call write_nc_file(.true.)

    deallocate(level,xs,ys,xstag,ystag)
    deallocate(thck,topog,kinbcmask,uvelhom,vvelhom)

  end subroutine make_steady_shelf

  subroutine make_linear_shelf()

    ! local variables
    integer :: i,j,k
    integer :: ifpos,kinbcw
    integer,parameter :: zero_buf = 1
    real :: hx,hy
    real :: rhoi, rhoo
    real :: upstream_thk, if_thk, upstream_vel,inflow_a
    real :: chan_depth,otopg,ltopg,acab_per_year
    real,dimension(:),allocatable :: rand_row
    logical :: noslip

    real :: tmp
    real :: rand_amp

    character (len=512) :: linear_shelf_params = "<fname> <nx> <ny> <n_level> <hx> <hy> <upstream_thk> &
       & <upstream_vel> <inflow_a> <ifpos> <ifthk> <ocean_depth> <kinbcw> [<noslip>]"

    if (command_argument_count() < 14) then
       write(*,*)"Incorrect number of parameters. Linear shelf requires: &
            &  ",trim(linear_shelf_params)
       stop 1
    end if

    call get_command_argument(2,argstr)
    read(argstr,'(a512)') fname
    write(*,*) 'fname ',trim(fname)

    call get_command_argument(3, argstr)
    read(argstr,'(i5)') nx
    write(*,*) 'nx:',nx

    call get_command_argument(4,argstr)
    read(argstr,'(i5)') ny
    write(*,*) 'ny:',ny

    call get_command_argument(5,argstr)
    read(argstr,'(i5)') n_level
    write(*,*) 'n_level:',n_level

    call get_command_argument(6,argstr)
    read(argstr,'(f18.12)') hx
    write(*,*) 'hx',hx

    call get_command_argument(7,argstr)
    read(argstr,'(f18.12)') hy
    write(*,*) 'hy',hy

    call get_command_argument(8,argstr)
    read(argstr,'(f18.12)') upstream_thk
    write(*,*) 'upstream_thk',upstream_thk

    call get_command_argument(9,argstr)
    read(argstr,'(f18.12)') upstream_vel
    write(*,*) 'upstream_vel',upstream_vel

    call get_command_argument(10,argstr)
    read(argstr,'(f18.12)') inflow_a
    write(*,*) 'inflow_a', inflow_a

    call get_command_argument(11,argstr)
    read(argstr,'(i5)') ifpos
    write(*,*) 'ifpos',ifpos

    call get_command_argument(12,argstr)
    read(argstr,'(f18.12)') if_thk
    write(*,*) 'if_thk', if_thk

    call get_command_argument(13,argstr)
    read(argstr,'(f18.12)') otopg
    write(*,*) 'otopg', otopg

    call get_command_argument(14,argstr)
    read(argstr,'(i5)') kinbcw
    write(*,*) 'kinbc_width',kinbcw

    if (command_argument_count() == 15) then
       call get_command_argument(15,argstr)
       read(argstr, '(l1)') noslip
    else
	noslip = .false.
    end if
    write(*,*) 'noslip', noslip

    allocate(xs(nx),ys(ny),level(n_level),xstag(nx-1),ystag(ny-1))
    allocate(topog(nx,ny),thck(nx,ny),kinbcmask(nx-1,ny-1))
    allocate(uvelhom(nx-1,ny-1,n_level), vvelhom(nx-1,ny-1,n_level))
    allocate(rand_row(nx-2*zero_buf))

    !now populate the dimension variables
    xs = (/ ( (i-1)*hx,i=1,nx ) /)
    ys = (/ ( (j-1)*hy,j=1,ny ) /)
    level = (/ ( real(i)/real((n_level-1)), i=0,(n_level-1) ) /)
    xstag = (/ ( ((real(i)-0.5)*hx),i=1,nx-1 ) /)
    ystag = (/ ( ((real(j)-0.5)*hy),j=1,ny-1 ) /)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! define the topography, thickness, kinbcmask       !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !define topg
    topog = -abs(otopg)
    
    !define thickness

    thck = 0.0
    thck(1+zero_buf:(nx-zero_buf),(ny-kinbcw):ny) = upstream_thk 
    
    do j=ifpos,(ny-kinbcw)
       thck(1+zero_buf:(nx-zero_buf),j) = if_thk + real(j-ifpos)/(ny-kinbcw+1-ifpos) * (upstream_thk-if_thk)
    end do

    ! define kinbcmask
    kinbcmask = 0
    kinbcmask(:,(ny-kinbcw):(ny-1)) = 1 !north edge
    if (noslip) then		    
       kinbcmask(1:(1+zero_buf),:) = 1
       kinbcmask((nx-1-zero_buf):(nx-1),:) = 1
    end if

    uvelhom = 0.d0
    vvelhom = 0.d0

    if (noslip) then
       vvelhom(2:(nx-2),(ny-kinbcw):ny,:) = upstream_vel
    else
        vvelhom(:,ny-kinbcw:ny-1,:) = upstream_vel
    end if

    do j=ny-kinbcw-10,ny
       do i=1+zero_buf,nx-zero_buf
          thck(i,j) = thck(i,j) + real(j-(ny-kinbcw-10))/real(ny-(ny-kinbcw-10)) & 
                                  *inflow_a * (xs(i)-xs((nx-1)/2+1))*(xs(i)-xs((nx-1)/2+1))
       end do
    end do

    call write_nc_file(.true.)

    deallocate(level,xs,ys,xstag,ystag)
    deallocate(thck,topog,kinbcmask,uvelhom,vvelhom)

  end subroutine make_linear_shelf

  subroutine make_two_sided_shelf()

    ! local variables
    integer :: landw = 2
    integer :: i,j
    integer :: firstarg = 3
    integer :: divpos,ifpos,kinbcw
    real :: hx,hy,githk,ifthk,kx,chan_amp,chan_init_length
    real :: chan_depth,otopg,ltopg

    character(len=512) :: two_sided_shelf_params = "<fname> <nx> <ny> <hx> <hy> &
         &<center_thk> <icefront_pos> &
         &<icefront_thk> <otopg> <ltopg> &
         &<kx> <chan_amp> <chan_init_l> <kinbcw>"

    if (command_argument_count() /= 15) then
       write(*,*)"Incorrect number of parameters. Two-sided shelf requires: ",&
            trim(two_sided_shelf_params)
       stop 1
    end if

    call get_command_argument(2,argstr)
    read(argstr,'(a512)') fname
    write(*,*) 'fname:',trim(fname)

    call get_command_argument(firstarg, argstr)
    read(argstr,'(i5)') nx
    write(*,*) 'nx:',nx

    call get_command_argument(firstarg + 1,argstr)
    read(argstr,'(i5)') ny
    write(*,*) 'ny:',ny

    call get_command_argument(firstarg + 2,argstr)
    read(argstr,'(f18.12)') hx
    write(*,*) 'hx',hx

    call get_command_argument(firstarg + 3,argstr)
    read(argstr,'(f18.12)') hy
    write(*,*) 'hy',hy

    call get_command_argument(firstarg + 4,argstr)
    read(argstr,'(f18.12)') githk
    write(*,*) 'githk',githk

    call get_command_argument(firstarg + 5,argstr)
    read(argstr,'(i5)') ifpos
    write(*,*) 'ifpos',ifpos

    call get_command_argument(firstarg + 6,argstr)
    read(argstr,'(f18.12)') ifthk
    write(*,*) 'ifthk',ifthk

    call get_command_argument(firstarg + 7,argstr)
    read(argstr,'(f18.12)') otopg
    write(*,*) 'otopg', otopg

    call get_command_argument(firstarg + 8,argstr)
    read(argstr,'(f18.12)') ltopg
    write(*,*) 'ltopg', ltopg

    call get_command_argument(firstarg + 9,argstr)
    read(argstr,'(f18.12)') kx
    write(*,*) 'kx',kx

    call get_command_argument(firstarg + 10,argstr)
    read(argstr,'(f18.12)') chan_amp
    write(*,*) 'chan_amp',chan_amp

    call get_command_argument(firstarg + 11,argstr)
    read(argstr,'(f18.12)') chan_init_length
    write(*,*) 'chan_init_length',chan_init_length  

    call get_command_argument(firstarg + 12,argstr)
    read(argstr,'(i5)') kinbcw
    write(*,*) 'kinbc_width',kinbcw

    allocate(xs(nx),ys(ny),xstag(nx-1),ystag(ny-1))
    allocate(topog(nx,ny),thck(nx,ny),kinbcmask(nx-1,ny-1))

    !now populate the dimension variables
    xs = (/ ( (i-1)*hx,i=1,nx ) /)
    ys = (/ ( (j-1)*hy,j=1,ny ) /)
    xstag = (/ ( ((real(i)-0.5)*hx),i=1,nx-1 ) /)
    ystag = (/ ( ((real(j)-0.5)*hy),j=1,ny-1 ) /)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! define the topograph, thickness, kinbcmask       !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !define topg
    topog = -abs(otopg)
    topog(1:landw,:) = abs(ltopg)         !west wall
    topog((nx+1-landw):nx,:) = abs(ltopg) !east wall

    !define thickness
    thck = 0.d0
    divpos = ny/2

    do j=divpos,ifpos
       do i=1,nx
          thck(i,j) = githk + &
               (ifthk -githk)*(real(j-divpos)/real(ifpos-divpos))
          chan_depth = 0.5 * chan_amp &
               * (1-exp(- real(j-divpos)*hy/chan_init_length))&
               * (1-cos((real(i-1)/real(nx-1))*2*pi*kx))
          thck(i,j) = thck(i,j) - chan_depth
       end do
    end do

    do j = (ny-ifpos+1),(divpos-1)
       do i=1,nx
          thck(i,j) = githk + &
               (ifthk -githk)*(real(divpos-j)/real(divpos-(ny-ifpos+1)))
          chan_depth = 0.5 * chan_amp &
               * (1-exp(- real(divpos-j)*hy/chan_init_length))&
               * (1-cos((real(i-1)/real(nx-1))*2*pi*kx))

          thck(i,j) = thck(i,j) - chan_depth
       end do
    end do

    thck(1:landw,:) = 0.d0
    thck((nx+1-landw):nx,:) = 0.d0
    thck(:,1:landw) = 0.d0
    thck(:,(ifpos+1):ny) = 0.d0 ! zero ahead of ice front position
    thck(:,1:(ny-ifpos)) = 0.d0 !zero before south ice front

    ! define kinbcmask
    kinbcmask = 0
    kinbcmask(1:kinbcw,:) = 1
    kinbcmask((nx-kinbcw):(nx-1),:) = 1

    call write_nc_file(.false.)

    deallocate(xs,ys,xstag,ystag,thck,topog,kinbcmask)


  end subroutine make_two_sided_shelf

  subroutine write_nc_file(write_vel)

    logical,intent(in) :: write_vel

    !local variables
    integer :: nc_id
    integer :: time_dimid,x_dimid,y_dimid,xstag_dimid,ystag_dimid,level_dimid
    integer :: x_varid,y_varid,time_varid,level_varid
    integer :: thck_varid,topog_varid,kinbcmask_varid,uvel_varid,vvel_varid
    integer :: xstag_varid,ystag_varid

    call check( nf90_create(fname, NF90_CLOBBER, nc_id) )

    ! defining dimensions
    call check( nf90_def_dim(nc_id,'time',1,time_dimid) )
    call check( nf90_def_dim(nc_id,'x1',nx,x_dimid) )
    call check( nf90_def_dim(nc_id,'y1',ny,y_dimid) )
    call check( nf90_def_dim(nc_id,'x0',nx-1,xstag_dimid) )
    call check( nf90_def_dim(nc_id,'y0',ny-1,ystag_dimid) )


    ! define variables
    call check( nf90_def_var(nc_id,'time',NF90_DOUBLE,(/time_dimid/),time_varid) )
    call check( nf90_put_att(nc_id, time_varid, 'long_name', 'time') )
    call check( nf90_put_att(nc_id, time_varid, 'units', 'seconds') )

    call check( nf90_def_var(nc_id,'x1',NF90_DOUBLE,(/x_dimid/),x_varid) )
    call check( nf90_put_att(nc_id, x_varid, 'long_name', 'Cartisian x-coordinate') )
    call check( nf90_put_att(nc_id, x_varid, 'units', 'meter') )

    call check( nf90_def_var(nc_id,'y1',NF90_DOUBLE,(/y_dimid/),y_varid) )
    call check( nf90_put_att(nc_id, y_varid, 'long_name', 'Cartisian y-coordinate') )
    call check( nf90_put_att(nc_id, y_varid, 'units', 'meter') )

    call check( nf90_def_var(nc_id,'x0',NF90_DOUBLE,(/xstag_dimid/),xstag_varid) )
    call check( nf90_put_att(nc_id, xstag_varid, 'long_name', 'Cartisian y-coordinate velocity grid') )
    call check( nf90_put_att(nc_id, xstag_varid, 'units', 'meter') )

    call check( nf90_def_var(nc_id,'y0',NF90_DOUBLE,(/ystag_dimid/),ystag_varid) )
    call check( nf90_put_att(nc_id, ystag_varid, 'long_name', 'Cartisian y-coordinate velocity grid') )
    call check( nf90_put_att(nc_id, ystag_varid, 'units', 'meter') )

    call check( nf90_def_var(nc_id,'thk',NF90_DOUBLE,(/x_dimid,y_dimid,time_dimid/),thck_varid) )
    call check( nf90_put_att(nc_id, thck_varid, 'long_name', 'ice thickness') )
    call check( nf90_put_att(nc_id, thck_varid, 'units', 'meter') )

    call check( nf90_def_var(nc_id,'topg',NF90_DOUBLE,(/x_dimid,y_dimid/),topog_varid) )
    call check( nf90_put_att(nc_id, topog_varid, 'long_name', 'topography') )
    call check( nf90_put_att(nc_id, topog_varid, 'units', 'meter') )

    call check( nf90_def_var(nc_id,'kinbcmask',NF90_DOUBLE,(/xstag_dimid,ystag_dimid,time_dimid/),kinbcmask_varid) )
    call check( nf90_put_att(nc_id, kinbcmask_varid, 'long_name', 'kinematic boundary condition mask') )

    if (write_vel) then
       call check( nf90_def_dim(nc_id,'level',n_level,level_dimid) )
       call check( nf90_def_var(nc_id,'level',NF90_DOUBLE,(/level_dimid/),level_varid) )
       call check( nf90_put_att(nc_id, level_varid, 'long_name', 'level') )
       call check( nf90_put_att(nc_id, level_varid, 'units', '1') )

       call check( nf90_def_var(nc_id,'uvelhom',NF90_DOUBLE,(/xstag_dimid,ystag_dimid,level_dimid,time_dimid/),uvel_varid) )
       call check( nf90_put_att(nc_id, uvel_varid, 'long_name', 'ice velocity in x direction') )
       call check( nf90_put_att(nc_id, uvel_varid, 'units', 'meter/year') )

       call check( nf90_def_var(nc_id,'vvelhom',NF90_DOUBLE,(/xstag_dimid,ystag_dimid,level_dimid,time_dimid/),vvel_varid) )
       call check( nf90_put_att(nc_id, vvel_varid, 'long_name', 'ice velocity in y direction') )
       call check( nf90_put_att(nc_id, vvel_varid, 'units', 'meter/year') )

    end if

    call check( nf90_enddef(nc_id) )

    ! populate variables

    call check( nf90_put_var(nc_id,x_varid,xs) )
    call check( nf90_put_var(nc_id,y_varid,ys) )
    call check( nf90_put_var(nc_id,xstag_varid,xstag) )
    call check( nf90_put_var(nc_id,ystag_varid,ystag) )
    call check( nf90_put_var(nc_id,time_varid, (/ 1 /) ) )


    ! write the arrays out to the netcdf file and close it
    call check( nf90_put_var(nc_id, thck_varid, thck) )
    call check( nf90_put_var(nc_id, topog_varid, topog) )
    call check( nf90_put_var(nc_id, kinbcmask_varid, kinbcmask) )
    
    if (write_vel) then
       call check( nf90_put_var(nc_id,level_varid, level) )
       call check( nf90_put_var(nc_id, uvel_varid, uvelhom) )
       call check( nf90_put_var(nc_id, vvel_varid, vvelhom) )
    end if

    call check( nf90_close(nc_id) )

  end subroutine write_nc_file

  subroutine check(status_code)

    integer,intent(in) :: status_code

    if (status_code /= 0) then
       write(*,*) 'fatal netcdf error:',nf90_strerror(status_code)
       stop 1
    end if

  end subroutine check

end program nc_gen_input
