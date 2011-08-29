program nc_uniform_end

  use netcdf

  implicit none

  ! module variables
  character (len=512) :: nc_file_in, nc_file_out

  integer :: thk_n_margin, kin_n_margin

  integer :: t_read
  integer,parameter :: dp = kind(1.0d0)

  integer :: ny, nx, nt
  integer :: nlevel
  real(kind=dp) :: hx, hy

  real(kind=dp) :: domain_xmin, domain_xmax, domain_ymin, domain_ymax

  !data arrays
  real(kind=dp),dimension(:),allocatable :: levels
  real(kind=dp),dimension(:),allocatable :: xs, ys
  real(kind=dp),dimension(:),allocatable :: xstag,ystag
  real(kind=dp),dimension(:,:),allocatable :: thck,topog
  real(kind=dp),dimension(:,:,:),allocatable :: uvelhom,vvelhom
  real(kind=dp),dimension(:,:,:),allocatable :: temp
  integer,dimension(:,:),allocatable :: kinbcmask

  integer :: template_row_index
  
  call main()

contains

subroutine main()

  character(len=512) :: gen_usage = "nc_uniform_end &
                                    &<nc_file_in> <nc_file_out> &
                                    &<t_read>&
                                    &<thk_n_margin> <kin_n_margin> &
                                    &<template_row_index>"

  character(len=512) :: argstr
  integer :: n,i,j
  integer,parameter :: n_baseargs = 6

  if (command_argument_count() < n_baseargs) then
     write(*,*) "Not enough arguments received (",command_argument_count(),").  Usage: ", trim(gen_usage)
     stop 1
  end if

  call get_command_argument(1,argstr)
  read(argstr,'(a)') nc_file_in
  write(*,*) 'nc_file_in:', trim(nc_file_in)

  call get_command_argument(2,argstr)
  read(argstr,'(a)') nc_file_out
  write(*,*) 'nc_file_out:', trim(nc_file_out)
		 
  call get_command_argument(3, argstr)
  read(argstr,'(i5)') t_read
  write(*,*) 't_read:', t_read

  call get_command_argument(4,argstr)
  read(argstr,'(i5)') thk_n_margin
  write(*,*)'thk_n_marg', thk_n_margin

  call get_command_argument(5,argstr)
  read(argstr,'(i5)') kin_n_margin
  write(*,*)'kin_n_marg', kin_n_margin

  call get_command_argument(6,argstr)
  read(argstr,'(i5)') template_row_index
  write(*,*)'template_row_index', template_row_index

  call read_old_nc_file() 

  allocate(xstag(nx-1),ystag(ny-1))
  xstag = (/ ( 0.5d0*(xs(i)+xs(i+1))  ,i=1,nx-1 ) /)
  ystag = (/ ( 0.5d0*(ys(j)+ys(j+1))  ,j=1,ny-1 ) /)


  call modify_thickness_velocity()
  
  call write_nc_file()

  deallocate(xs,ys)
  deallocate(xstag,ystag)

end subroutine main

subroutine read_old_nc_file()

  !local variables  
  integer :: nc_id
  integer :: time_dimid,time_varid,x_dimid,y_dimid,level_dimid
  integer :: x_varid,y_varid,level_varid
  integer :: thck_varid,topog_varid,kinbcmask_varid
  integer :: uvelhom_varid,vvelhom_varid,temp_varid

  character(len=NF90_MAX_NAME) :: x1_name,y1_name,t_name,level_name
  
  call check(nf90_open(trim(nc_file_in), NF90_NOWRITE, nc_id))  

  call check(nf90_inq_dimid(nc_id, 'x1', x_dimid))
  call check(nf90_inq_dimid(nc_id, 'y1', y_dimid))
  call check(nf90_inq_dimid(nc_id, 'time', time_dimid))
  call check(nf90_inq_dimid(nc_id, 'level', level_dimid))

  call check(nf90_inq_varid(nc_id, 'x1', x_varid))
  call check(nf90_inq_varid(nc_id, 'y1', y_varid))
  call check(nf90_inq_varid(nc_id, 'time', time_varid))
  call check(nf90_inq_varid(nc_id, 'level', level_varid))
  
  call check(nf90_inq_varid(nc_id, 'thk', thck_varid))
  call check(nf90_inq_varid(nc_id, 'topg', topog_varid))
  call check(nf90_inq_varid(nc_id, 'kinbcmask', kinbcmask_varid))
  call check(nf90_inq_varid(nc_id, 'uvelhom', uvelhom_varid))
  call check(nf90_inq_varid(nc_id, 'vvelhom', vvelhom_varid))
  call check(nf90_inq_varid(nc_id, 'temp', temp_varid))

  call check(nf90_inquire_dimension(nc_id, x_dimid, x1_name, nx))
  call check(nf90_inquire_dimension(nc_id, y_dimid, y1_name, ny))
  call check(nf90_inquire_dimension(nc_id, time_dimid, t_name, nt))
  call check(nf90_inquire_dimension(nc_id, level_dimid, level_name, nlevel))

  allocate(xs(nx))
  allocate(ys(ny))

  call check(nf90_get_var(nc_id, x_varid, xs, start= (/ 1 /), count=(/ nx /)))
  call check(nf90_get_var(nc_id, y_varid, ys, start= (/ 1 /), count=(/ ny /)))
        
  hx = xs(2) - xs(1)
  hy = ys(2) - ys(1)
 
  if (t_read < 1) then
	write(*,*) 'reading last time slice'
	t_read = nt
  end if

  if (nt < t_read) then
	write(*,*) 't_read exceeds available time slices: ', nt
	stop 1
  end if

  allocate(thck(nx,ny),topog(nx,ny))
  allocate(kinbcmask(nx-1, ny-1))
  allocate(uvelhom(nx-1,ny-1,nlevel),vvelhom(nx-1,ny-1,nlevel))
  allocate(temp(nx, ny, nlevel))

  allocate(levels(nlevel))

  call check(nf90_get_var(nc_id, level_varid, levels, &
                          start= (/ 1 /), &
                          count= (/ nlevel /)))

  call check(nf90_get_var(nc_id, thck_varid, thck,&
                          start= (/ 1,1,t_read /), &
 			  count= (/ nx,ny,1 /) ))
  call check(nf90_get_var(nc_id, topog_varid,topog, &
                           start= (/ 1,1 /), &
			   count= (/ nx,ny /) ))
  call check(nf90_get_var(nc_id, kinbcmask_varid, kinbcmask, &
			start= (/ 1,1,t_read /), &
		        count= (/ (nx-1),(ny-1),1 /) ))

  call check(nf90_get_var(nc_id, uvelhom_varid, uvelhom, &
                        start= (/ 1,1,1,t_read /), &
                        count= (/ (nx-1),(ny-1),nlevel,1 /) ))
  call check(nf90_get_var(nc_id, vvelhom_varid, vvelhom, &
                        start= (/ 1,1,1,t_read /), &
                        count= (/ (nx-1),(ny-1),nlevel,1 /) ))
  
  call check(nf90_get_var(nc_id, temp_varid, temp, &
                        start= (/ 1,1,1,t_read /), &
                        count= (/ (nx),(ny), nlevel, 1 /) ))

end subroutine read_old_nc_file

subroutine modify_thickness_velocity()

  integer :: j

  do j=(template_row_index+1),(ny-thk_n_margin)
     thck(:,j) = thck(:,template_row_index)
  end do

end subroutine modify_thickness_velocity

subroutine write_nc_file()
    
    !local variables
    integer :: nc_id
    integer :: time_dimid,x_dimid,y_dimid,xstag_dimid,ystag_dimid,level_dimid
    integer :: x_varid,y_varid,time_varid,level_varid
    integer :: thck_varid,topog_varid,kinbcmask_varid,uvelhom_varid,vvelhom_varid
    integer :: temp_varid
    integer :: xstag_varid,ystag_varid
    
    call check( nf90_create(trim(nc_file_out), NF90_CLOBBER, nc_id) )

    ! defining dimensions
    call check( nf90_def_dim(nc_id,'time',1,time_dimid) )
    call check( nf90_def_dim(nc_id,'x1',nx,x_dimid) )
    call check( nf90_def_dim(nc_id,'y1',ny,y_dimid) )
    call check( nf90_def_dim(nc_id,'x0',nx-1,xstag_dimid) )
    call check( nf90_def_dim(nc_id,'y0',ny-1,ystag_dimid) )
    call check( nf90_def_dim(nc_id,'level',nlevel,level_dimid) )

    ! define variables
    call check( nf90_def_var(nc_id,'level',NF90_DOUBLE,(/level_dimid/),level_varid) )
    call check( nf90_put_att(nc_id, level_varid, 'long_name', 'sigma level') )
  
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
    call check( nf90_put_att(nc_id, xstag_varid, 'long_name', &
				  'Cartisian y-coordinate velocity grid') )
    call check( nf90_put_att(nc_id, xstag_varid, 'units', 'meter') )

    call check( nf90_def_var(nc_id,'y0',NF90_DOUBLE,(/ystag_dimid/),ystag_varid) )
    call check( nf90_put_att(nc_id, ystag_varid, 'long_name',  &
				 'Cartisian y-coordinate velocity grid') )
    call check( nf90_put_att(nc_id, ystag_varid, 'units', 'meter') )

    call check( nf90_def_var(nc_id,'thk',NF90_DOUBLE, &
 			(/x_dimid,y_dimid,time_dimid/),thck_varid) )
    call check( nf90_put_att(nc_id, thck_varid, 'long_name', 'ice thickness') )
    call check( nf90_put_att(nc_id, thck_varid, 'units', 'meter') )

    call check( nf90_def_var(nc_id,'topg',NF90_DOUBLE,(/x_dimid,y_dimid/),topog_varid) )
    call check( nf90_put_att(nc_id, topog_varid, 'long_name', 'topography') )
    call check( nf90_put_att(nc_id, topog_varid, 'units', 'meter') )
    
    call check( nf90_def_var(nc_id,'kinbcmask',NF90_DOUBLE, &
                             (/xstag_dimid,ystag_dimid,time_dimid/),kinbcmask_varid) )
    call check( nf90_put_att(nc_id, kinbcmask_varid, 'long_name',  &
                                 'kinematic boundary condition mask') )

    call check( nf90_def_var(nc_id,'uvelhom', NF90_DOUBLE, &
                             (/xstag_dimid,ystag_dimid,level_dimid,time_dimid/),uvelhom_varid))
    call check( nf90_put_att(nc_id, uvelhom_varid, 'long_name', &
                                 'x velocity') )

    call check( nf90_def_var(nc_id,'vvelhom', NF90_DOUBLE, &
                             (/xstag_dimid,ystag_dimid,level_dimid,time_dimid/),vvelhom_varid))
    call check( nf90_put_att(nc_id, vvelhom_varid, 'long_name', &
                                 'y velocity') )
    
    call check( nf90_def_var(nc_id,'temp', NF90_DOUBLE, &
                             (/ x_dimid, y_dimid,level_dimid, time_dimid/), temp_varid))
    call check( nf90_put_att(nc_id, temp_varid, 'long_name', &
                              'temperature'))

    call check( nf90_enddef(nc_id) )

    call check( nf90_put_var(nc_id,x_varid,xs) )
    call check( nf90_put_var(nc_id,y_varid,ys) )
    call check( nf90_put_var(nc_id,xstag_varid,xstag) )
    call check( nf90_put_var(nc_id,ystag_varid,ystag) )
    call check( nf90_put_var(nc_id,level_varid,levels) )
    call check( nf90_put_var(nc_id,time_varid, (/ 1 /) ) )

    ! write the arrays out to the netcdf file and close it
    call check( nf90_put_var(nc_id, thck_varid, thck) )
    call check( nf90_put_var(nc_id, topog_varid, topog) )
    call check( nf90_put_var(nc_id, kinbcmask_varid, kinbcmask) )

    call check( nf90_put_var(nc_id, uvelhom_varid, uvelhom) )
    call check( nf90_put_var(nc_id, vvelhom_varid, vvelhom) )
    
    call check( nf90_put_var(nc_id, temp_varid, temp) )

    call check( nf90_close(nc_id) )

  end subroutine write_nc_file

  subroutine check(status_code)

    integer,intent(in) :: status_code

    if (status_code /= 0) then
       write(*,*) 'Fatal netcdf error:',nf90_strerror(status_code)
       stop 1
    end if

  end subroutine check

end program nc_uniform_end

