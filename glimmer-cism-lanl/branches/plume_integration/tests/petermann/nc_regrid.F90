program nc_regrid

  use netcdf

  implicit none

  ! module variables

  character (len=512) :: nc_file_in, nc_file_out
  integer :: n_margin, s_margin, w_margin, e_margin
  integer :: kmask_width
  integer :: t_read

  integer :: nx_new, ny_new, ny_old, nx_old, nt_old
  real :: hx_old, hy_old, hx_new, hy_new

  ! first two old x values, and old y values
  real,dimension(2) :: x1_old,y1_old

  !data arrays
  real,dimension(:),allocatable :: xs_new,ys_new
  real,dimension(:),allocatable :: xstag_new,ystag_new
  real,dimension(:,:),allocatable :: thck_old,thck_new,topog_old,topog_new
  integer,dimension(:,:),allocatable :: kinbcmask_old,kinbcmask_new

  call main()

contains

subroutine main()

  character(len=512) :: gen_usage = "nc_regrid &
                                    &nc_file_in nc_file_out [params]"
  character(len=512) :: params = "<t_read> <n_margin> <s_margin> &
                                    &<w_margin> <e_margin> <new_m> <new_n> &
                                    &<kmask_width>"

  character(len=128) :: argstr

  if (command_argument_count() /= 10) then
     write(*,*) "Usage: ", trim(gen_usage)
     write(*,*) "     params: ", trim(params)
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

  call get_command_argument(4, argstr)
  read(argstr,'(i5)') n_margin
  write(*,*) 'n_margin:',n_margin

  call get_command_argument(5,argstr)
  read(argstr,'(i5)') s_margin
  write(*,*) 's_margin:',s_margin

  call get_command_argument(6, argstr)
  read(argstr,'(i5)') w_margin
  write(*,*) 'w_margin:',w_margin

  call get_command_argument(7,argstr)
  read(argstr,'(i5)') e_margin
  write(*,*) 'e_margin:',e_margin

  call get_command_argument(8,argstr)
  read(argstr,'(i5)') kmask_width
  write(*,*)'kmask_width', kmask_width

  call get_command_argument(9, argstr)
  read(argstr,'(i5)') nx_new
  write(*,*) 'new_m:',nx_new

  call get_command_argument(10,argstr)
  read(argstr,'(i5)') ny_new
  write(*,*) 'new_n:', ny_new

  call read_old_nc_file() 

  allocate(xs_new(nx_new),ys_new(ny_new))
  allocate(xstag_new(nx_new-1),ystag_new(ny_new-1))
  allocate(topog_new(nx_new,ny_new),thck_new(nx_new,ny_new))
  allocate(kinbcmask_new(nx_new-1,ny_new-1))

  call define_new_data()
  
  call write_nc_file()

  deallocate(xs_new,ys_new,xstag_new,ystag_new,&
             thck_new,topog_new,kinbcmask_new)

end subroutine main

elemental function is_land(w)

	integer :: is_land
	integer,intent(in) :: w

	if (w > 0.0) then
	   is_land = 1
        else
	   is_land = 0
        end if

end function is_land

subroutine write_margins(data_old,data_new, &
                         test_sign, threshold, &
                         nx0,nx1,ny0,ny1, &
	                 n_marg,s_marg,w_marg,e_marg)

  real,dimension(:,:),intent(in) :: data_old
  real,dimension(:,:),intent(inout) :: data_new
  real,intent(in) :: test_sign, threshold
  integer,intent(in) :: nx0,nx1,ny0,ny1
  integer,intent(in) :: n_marg,s_marg,w_marg,e_marg

  if (test_sign * data_old(1,ny0/2) > test_sign * threshold) then
	!land margin on west
	data_new(1:w_marg,:) = data_old(1,ny0/2)
  else
	!ocean margin on west
	data_new(1:w_marg, (s_marg+1):(ny1-n_marg)) = &
				          data_old(1,ny0/2)
  end if
  
  if (test_sign * data_old(nx0,ny0/2) > test_sign * threshold) then
	!land margin on east side
	data_new((nx1+1-e_marg):nx1,:) = data_old(nx0,ny0/2)
  else
      ! ocean margin
      data_new((nx1+1-e_marg):nx1,(s_marg+1):(ny1-n_marg)) = &
	 					 data_old(nx0,ny0/2)
  end if

  if(test_sign * data_old(nx0/2, 1) > test_sign * threshold) then
	!land margin on south side
	data_new(1:nx1, 1:s_marg) = data_old(nx0/2,1)
  else
	data_new((w_marg+1):(nx1-e_marg),1:s_marg) = &
					  data_old(nx0/2,1)
  end if

  if(test_sign * data_old(nx0/2, ny0) > test_sign * threshold) then
	!land margin on north side
	data_new(1:nx1, (ny1+1-n_marg):ny1) = &
					data_old(nx0/2,ny0)
  else
	data_new((w_marg+1)+1:(nx1-e_marg),&
	         (ny1+1-n_marg):ny1) = &
                               data_old(nx0/2,ny0)
  end if

end subroutine write_margins


subroutine write_interior(data_old, data_new, nx0,nx1,ny0,ny1)

  real,dimension(:,:),intent(in) :: data_old
  real,dimension(:,:),intent(inout) :: data_new
  integer,intent(in) :: nx0,nx1,ny0,ny1

  !local vars
  integer :: i00,i01,i10,i11,j00,j01,j10,j11 
  integer :: i0,j0
  real :: a,b, temp
  integer :: i,j
  i00 = w_margin+1
  i01 = w_margin+1
  i10 = nx0 - e_margin
  i11 = nx1 - e_margin
  
  j00 = s_margin + 1
  j01 = s_margin + 1
  j10 = ny0 - n_margin
  j11 = ny1 - n_margin

  do i=i01,i11
     do j=j01,j11
 	 
         i0 = aint((real(i) - i01)*(hx_new/hx_old)) + real(i00)
         j0 = aint((real(j) - j01)*(hy_new/hx_old)) + real(j01)
         a = (real(i)-i01)-aint((real(i)-i01)*(hx_new/hx_old))*(hx_old/hx_new)
	 b = (real(j)-j01)-aint((real(j)-j01)*(hy_new/hy_old))*(hy_old/hy_new)

	 data_new(i,j) = (1-a)*(1-b)*data_old(i0,j0)     + &
	                 (1-a)*   b *data_old(i0,j0+1)   + &
                            a *(1-b)*data_old(i0+1,j0)     + &
	                    a *   b *data_old(i0+1,j0+1)
     end do
  end do

end subroutine write_interior

subroutine read_old_nc_file()

  !local variables  
  integer :: nc_id
  integer :: time_dimid,time_varid,x_dimid,y_dimid
  integer :: x_varid,y_varid
  integer :: thck_varid,topog_varid,kinbcmask_varid



  character(len=NF90_MAX_NAME) :: x1_name,y1_name,t_name
  
  call check(nf90_open(trim(nc_file_in), NF90_NOWRITE, nc_id))  

  call check(nf90_inq_dimid(nc_id, 'x1', x_dimid))
  call check(nf90_inq_dimid(nc_id, 'y1', y_dimid))
  call check(nf90_inq_varid(nc_id, 'x1', x_varid))
  call check(nf90_inq_varid(nc_id, 'y1', y_varid))
  call check(nf90_inq_dimid(nc_id, 'time', time_dimid))
  call check(nf90_inq_varid(nc_id, 'time', time_varid))

  call check(nf90_inq_varid(nc_id, 'thk', thck_varid))
  call check(nf90_inq_varid(nc_id, 'topg', topog_varid))
  call check(nf90_inq_varid(nc_id, 'kinbcmask', kinbcmask_varid))

  call check(nf90_inquire_dimension(nc_id, x_dimid, x1_name, nx_old))
  call check(nf90_inquire_dimension(nc_id, y_dimid, y1_name, ny_old))
  call check(nf90_inquire_dimension(nc_id, time_dimid, t_name, nt_old))

  call check(nf90_get_var(nc_id, x_varid, x1_old, start= (/ 1 /), &
						  count= (/ 2 /) ))
  hx_old = x1_old(2) - x1_old(1)

  call check(nf90_get_var(nc_id, y_varid, y1_old, start= (/ 1 /), &
      						  count= (/ 2 /) ))
  hy_old = y1_old(2) - y1_old(1)

  if (nt_old < t_read) then
	write(*,*) 't_read exceeds available time slices: ', nt_old
	stop 1
  end if

  allocate(thck_old(nx_old,ny_old),topog_old(nx_old,ny_old))
  allocate(kinbcmask_old(nx_old-1, ny_old-1))

  call check(nf90_get_var(nc_id, thck_varid, thck_old,&
                          start= (/ 1,1,t_read /), &
 			  count= (/ nx_old,ny_old,1 /) ))
  call check(nf90_get_var(nc_id, topog_varid,topog_old, &
                           start= (/ 1,1 /), &
			   count= (/ nx_old,ny_old /) ))
  call check(nf90_get_var(nc_id, kinbcmask_varid, kinbcmask_old, &
			start= (/ 1,1,t_read /), &
		        count= (/ (nx_old-1),(ny_old-1),1 /) ))

end subroutine read_old_nc_file


subroutine define_new_data()

  ! local variables
  integer :: i,j
  real,dimension(nx_new-1,ny_new-1) :: kinbcmask_new_real

  hx_new = hx_old * (nx_old - w_margin - e_margin - 1) / &
                    (nx_new - w_margin - e_margin - 1)
  hy_new = hy_old * (ny_old - s_margin - n_margin - 1) / &
                    (ny_new - s_margin - n_margin - 1)

  !now populate the dimension variables
  xs_new = (/ ( x1_old(1) + (i-1)*hx_new,i=1,nx_new ) /)
  ys_new = (/ ( y1_old(1) + (j-1)*hy_new,j=1,ny_new ) /)
  xstag_new = (/ ( x1_old(1) + ((real(i)-0.5)*hx_new),i=1,nx_new-1 ) /)
  ystag_new = (/ ( y1_old(1) + ((real(j)-0.5)*hy_new),j=1,ny_new-1 ) /)
  
  ! now define the thck, topog, kinbcmask arrays
  call write_margins(topog_old,topog_new, +1.0, 0.0,&
                          nx_old, nx_new,ny_old, ny_new, &
	                  n_margin,s_margin,w_margin,e_margin)

  call write_margins(thck_old,thck_new, -1.0, 0.0, &
                          nx_old, nx_new,ny_old, ny_new, &
                          n_margin,s_margin,w_margin,e_margin)
  call write_margins(real(kinbcmask_old), kinbcmask_new_real, +1.0,0.0,&
                          nx_old-1, nx_new-1, ny_old-1, ny_new-1, &
			  is_land(n_margin)*kmask_width, &
	                  is_land(s_margin)*kmask_width, &
	 		  is_land(w_margin)*kmask_width, &
	                  is_land(e_margin)*kmask_width)
  
  call write_interior(topog_old,topog_new, nx_old,nx_new,ny_old,ny_new)
  call write_interior(thck_old,thck_new,nx_old,nx_new,ny_old,ny_new)
  call write_interior(real(kinbcmask_old),kinbcmask_new_real, &
		 	nx_old,nx_new,ny_old,ny_new)

  kinbcmask_new = int(kinbcmask_new_real)

end subroutine define_new_data

subroutine write_nc_file()
    
    !local variables
    integer :: nc_id
    integer :: time_dimid,x_dimid,y_dimid,xstag_dimid,ystag_dimid
    integer :: x_varid,y_varid,time_varid
    integer :: thck_varid,topog_varid,kinbcmask_varid
    integer :: xstag_varid,ystag_varid
    
    call check( nf90_create(trim(nc_file_out), NF90_CLOBBER, nc_id) )

    ! defining dimensions
    call check( nf90_def_dim(nc_id,'time',1,time_dimid) )
    call check( nf90_def_dim(nc_id,'x1',nx_new,x_dimid) )
    call check( nf90_def_dim(nc_id,'y1',ny_new,y_dimid) )
    call check( nf90_def_dim(nc_id,'x0',nx_new-1,xstag_dimid) )
    call check( nf90_def_dim(nc_id,'y0',ny_new-1,ystag_dimid) )
 
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

    call check( nf90_enddef(nc_id) )

    call check( nf90_put_var(nc_id,x_varid,xs_new) )
    call check( nf90_put_var(nc_id,y_varid,ys_new) )
    call check( nf90_put_var(nc_id,xstag_varid,xstag_new) )
    call check( nf90_put_var(nc_id,ystag_varid,ystag_new) )
    call check( nf90_put_var(nc_id,time_varid, (/ 1 /) ) )

    ! write the arrays out to the netcdf file and close it
    call check( nf90_put_var(nc_id, thck_varid, thck_new) )
    call check( nf90_put_var(nc_id, topog_varid, topog_new) )
    call check( nf90_put_var(nc_id, kinbcmask_varid, kinbcmask_new) )

    call check( nf90_close(nc_id) )

  end subroutine write_nc_file

  subroutine check(status_code)

    integer,intent(in) :: status_code

    if (status_code /= 0) then
       write(*,*) 'fatal netcdf error:',nf90_strerror(status_code)
       stop 1
    end if

  end subroutine check



end program nc_regrid
