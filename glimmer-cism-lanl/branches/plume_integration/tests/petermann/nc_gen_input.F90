program nc_gen_input

  use netcdf

  implicit none

  real,parameter :: pi = 3.1415926
  real,parameter :: rhoi_rhow = 1030.d0 / 920.d0

  integer :: iarg_count
  character(len=512) :: usage = "nc_gen_input <nx> <ny> <hx> <hy> <glpos>  &
                               & <gldep> <ifpos> <ifdep> <kx> <chan_amp> <fname>"
  character(len=128) :: argstr

  integer :: nx,ny,glpos,ifpos
  real :: hx,hy,gldep,ifdep,kx,chan_amp
  character (len=512) :: fname

  !get arguments from command line
  iarg_count = command_argument_count()

  if (iarg_count < 11) then
     write(*,*) "Usage: ", trim(usage)
     stop
  end if

  call get_command_argument(1,argstr)
  read(argstr,'(i5)') nx
  write(*,*) 'nx:',nx

  call get_command_argument(2,argstr)
  read(argstr,'(i5)') ny
  write(*,*) 'ny:',ny

  call get_command_argument(3,argstr)
  read(argstr,'(f8.2)') hx
  write(*,*) 'hx',hx

  call get_command_argument(4,argstr)
  read(argstr,'(f8.2)') hy
  write(*,*) 'hy',hy

  call get_command_argument(5,argstr)
  read(argstr,'(i5)') glpos
  write(*,*) 'glpos',glpos

  call get_command_argument(6,argstr)
  read(argstr,'(f8.2)') gldep
  write(*,*) 'gldep',gldep

  call get_command_argument(7,argstr)
  read(argstr,'(i5)') ifpos
  write(*,*) 'ifpos',ifpos

  call get_command_argument(8,argstr)
  read(argstr,'(f8.2)') ifdep
  write(*,*) 'ifdep',ifdep

  call get_command_argument(9,argstr)
  read(argstr,'(f8.2)') kx
  write(*,*) 'kx',kx
  
  call get_command_argument(10,argstr)
  read(argstr,'(f8.2)') chan_amp
  write(*,*) 'chan_amp',chan_amp
  
  call get_command_argument(11,argstr)
  read(argstr,'(a512)') fname
  write(*,*) 'fname ',trim(fname)

  call make_nc_file(nx,ny,hx,hy,glpos,gldep,ifpos,ifdep,kx,chan_amp,trim(fname))


contains

  subroutine make_nc_file(nx,ny,hx,hy,glpos,gldep,ifpos,ifdep,kx,chan_amp,fname)

    character(len=*),intent(in) :: fname
    integer,intent(in) :: nx,ny,glpos,ifpos
    real,intent(in) :: hx,hy,gldep,ifdep,kx,chan_amp
    
    ! local variables
    integer :: i,j,nc_id
    integer :: time_dimid,x_dimid,y_dimid
    integer :: x_varid,y_varid,time_varid,thck_varid,topog_varid,kinbcmask_varid

    !data arrays
    real,dimension(nx) :: xs
    real,dimension(ny) :: ys
    real,dimension(nx,ny) :: thck,topog,kinbcmask
    
    call check( nf90_create(fname, NF90_CLOBBER, nc_id) )

    ! defining dimensions
    call check( nf90_def_dim(nc_id,'time',1,time_dimid) )
    call check( nf90_def_dim(nc_id,'x1',nx,x_dimid) )
    call check( nf90_def_dim(nc_id,'y1',ny,y_dimid) )

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

    call check( nf90_def_var(nc_id,'thk',NF90_DOUBLE,(/x_dimid,y_dimid,time_dimid/),thck_varid) )
    call check( nf90_put_att(nc_id, thck_varid, 'long_name', 'ice thickness') )
    call check( nf90_put_att(nc_id, thck_varid, 'units', 'meter') )

    call check( nf90_def_var(nc_id,'topg',NF90_DOUBLE,(/x_dimid,y_dimid/),topog_varid) )
    call check( nf90_put_att(nc_id, topog_varid, 'long_name', 'topography') )
    call check( nf90_put_att(nc_id, topog_varid, 'units', 'meter') )
    
    call check( nf90_def_var(nc_id,'kinbcmask',NF90_DOUBLE,(/x_dimid,y_dimid,time_dimid/),kinbcmask_varid) )
    call check( nf90_put_att(nc_id, kinbcmask_varid, 'long_name', 'kinematic boundary condition mask') )

    call check( nf90_enddef(nc_id) )

    !now populate the dimension variables
    xs = (/ ( (i-1)*hx,i=1,nx ) /)
    ys = (/ ( (j-1)*hy,j=1,ny ) /)
    
    call check( nf90_put_var(nc_id,x_varid,xs) )
    call check( nf90_put_var(nc_id,y_varid,ys) )
    call check( nf90_put_var(nc_id,time_varid, (/ 1 /) ) )
   
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! define the topograph, thickness, kinbcmask       !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    !define topg
    topog = -gldep
    topog(1:2,:) = 0.d0 !at sea level
    topog((nx-1):nx,:) = 0.d0
    topog(:,1:2) = 0.d0
    
    !define thickness
    thck(:,1:glpos) = gldep 
    do j=glpos+1,ifpos
	if (ifpos /= glpos) then
           thck(:,j) = gldep + (ifdep - gldep)*(real(j-glpos)/real(ifpos-glpos))
        end if
    end do
    do i=1,nx
       thck(i,(glpos+1):ifpos) = thck(i,(glpos+1):ifpos) - &
                         0.5*chan_amp*(1-cos((real(i-1)/real(nx-1))*2*pi*kx))
    end do
    thck(1:2,:) = 0.d0
    thck((nx-1):nx,:) = 0.d0
    thck(:,1:2) = 0.d0
    thck(:,(ifpos+1):ny) = 0.d0 ! zero ahead of ice front position
    thck = thck * rhoi_rhow
    
    ! define kinbcmask
    kinbcmask(1:3,:) = 1
    kinbcmask(:,1:3) = 1
    kinbcmask((nx-2):nx,:) = 1

    ! write the arrays out to the netcdf file and close it
    call check( nf90_put_var(nc_id, thck_varid, thck) )
    call check( nf90_put_var(nc_id, topog_varid, topog) )
    call check( nf90_put_var(nc_id, kinbcmask_varid, kinbcmask) )

    call check( nf90_close(nc_id) )

  end subroutine make_nc_file

  subroutine check(status_code)

    integer,intent(in) :: status_code

    if (status_code /= 0) then
       write(*,*) 'fatal netcdf error:',nf90_strerror(status_code)
       stop
    end if

  end subroutine check

end program nc_gen_input
