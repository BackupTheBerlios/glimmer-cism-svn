program nc_regrid

  use netcdf

  implicit none

  ! module variables
  character (len=512) :: nc_file_in, nc_file_out

  integer :: top_n_margin, top_s_margin, top_w_margin, top_e_margin
  integer :: thk_n_margin, thk_s_margin, thk_w_margin, thk_e_margin
  integer :: kin_n_margin, kin_s_margin, kin_w_margin, kin_e_margin

  integer :: t_read
  integer,parameter :: dp = kind(1.0d0)

  integer :: nx_new, ny_new, ny_old, nx_old, nt_old
  integer :: nlevel                                !only allow new nlevel = old nlevel
  real(kind=dp) :: hx_old, hy_old, hx_new, hy_new

!  real(kind=dp) :: hx_new_kin, hy_new_kin

  ! first two old x values, and old y values
  real(kind=dp),dimension(2) :: x1_old,y1_old


  !data arrays
  real(kind=dp),dimension(:),allocatable :: levels
  real(kind=dp),dimension(:),allocatable :: xs_new,ys_new
  real(kind=dp),dimension(:),allocatable :: xstag_new,ystag_new
  real(kind=dp),dimension(:,:),allocatable :: thck_old,thck_new,topog_old,topog_new
  real(kind=dp),dimension(:,:,:),allocatable :: uvelhom_old,uvelhom_new,vvelhom_old,vvelhom_new
  real(kind=dp),dimension(:,:,:),allocatable :: temp_old, temp_new
  integer,dimension(:,:),allocatable :: kinbcmask_old,kinbcmask_new
  
  call main()

contains

subroutine main()

  character(len=512) :: gen_usage = "nc_regrid &
                                    &<nc_file_in> <nc_file_out> &
                                    &<t_read> <new_m> <new_n> &
                                    &<top_n_margin> <top_s_margin> &
                                    &<top_w_margin> <top_e_margin> &
                                    &<thk_n_margin> <thk_s_margin> &
                                    &<thk_w_margin> <thk_e_margin> &
                                    &<kin_n_margin> <kin_s_margin> &
                                    &<kin_w_margin> <kin_e_margin> "

  character(len=128) :: argstr


  if (command_argument_count() /= 17) then
     write(*,*) "Usage: ", trim(gen_usage)
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
  read(argstr,'(i5)') nx_new
  write(*,*) 'new_m:',nx_new

  call get_command_argument(5,argstr)
  read(argstr,'(i5)') ny_new
  write(*,*) 'new_n:', ny_new

  call get_command_argument(6, argstr)
  read(argstr,'(i5)') top_n_margin
  write(*,*) 'top_n_margin:',top_n_margin

  call get_command_argument(7,argstr)
  read(argstr,'(i5)') top_s_margin
  write(*,*) 'top_s_margin:',top_s_margin

  call get_command_argument(8, argstr)
  read(argstr,'(i5)') top_w_margin
  write(*,*) 'top_w_margin:',top_w_margin

  call get_command_argument(9,argstr)
  read(argstr,'(i5)') top_e_margin
  write(*,*) 'top_e_margin:',top_e_margin

  call get_command_argument(10,argstr)
  read(argstr,'(i5)') thk_n_margin
  write(*,*)'thk_n_marg', thk_n_margin

  call get_command_argument(11,argstr)
  read(argstr,'(i5)') thk_s_margin
  write(*,*)'thk_s_margin', thk_s_margin

  call get_command_argument(12,argstr)
  read(argstr,'(i5)') thk_w_margin
  write(*,*)'thk_w_margin', thk_w_margin

  call get_command_argument(13,argstr)
  read(argstr,'(i5)') thk_e_margin
  write(*,*)'thk_e_margin', thk_e_margin

  call get_command_argument(14,argstr)
  read(argstr,'(i5)') kin_n_margin
  write(*,*)'kin_n_marg', kin_n_margin

  call get_command_argument(15,argstr)
  read(argstr,'(i5)') kin_s_margin
  write(*,*)'kin_s_margin', kin_s_margin

  call get_command_argument(16,argstr)
  read(argstr,'(i5)') kin_w_margin
  write(*,*)'kin_w_margin', kin_w_margin

  call get_command_argument(17,argstr)
  read(argstr,'(i5)') kin_e_margin
  write(*,*)'kin_e_margin', kin_e_margin

  call read_old_nc_file() 

  allocate(xs_new(nx_new),ys_new(ny_new))
  allocate(xstag_new(nx_new-1),ystag_new(ny_new-1))
  allocate(topog_new(nx_new,ny_new),thck_new(nx_new,ny_new))
  allocate(kinbcmask_new(nx_new-1,ny_new-1))
  allocate(uvelhom_new(nx_new-1,ny_new-1,nlevel),vvelhom_new(nx_new-1,ny_new-1,nlevel))
  allocate(temp_new(nx_new,ny_new,nlevel))

  call define_new_data()
  
  call write_nc_file()

  deallocate(xs_new,ys_new,xstag_new,ystag_new,&
             thck_new,topog_new,kinbcmask_new, &
             uvelhom_new,vvelhom_new,temp_new)

end subroutine main

subroutine write_real_margins(data_old,data_new, &
                         test_sign, threshold, &
                         nx0,nx1,ny0,ny1, &
	                 n_marg,s_marg,w_marg,e_marg)

  real(kind=dp),dimension(:,:),intent(in) :: data_old
  real(kind=dp),dimension(:,:),intent(inout) :: data_new
  real(kind=dp),intent(in) :: test_sign, threshold
  integer,intent(in) :: nx0,nx1,ny0,ny1
  integer,intent(in) :: n_marg,s_marg,w_marg,e_marg

  if (test_sign * data_old(1,ny0/2) > test_sign * threshold) then
	!land margin on west
	data_new(1:w_marg,:) = data_old(1,ny0/2)
  else
	!ocean margin on west
	data_new(1:w_marg, (s_marg+1):(ny1-n_marg)) = data_old(1,ny0/2)
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
	data_new((w_marg+1):(nx1-e_marg),&
	         (ny1+1-n_marg):ny1) = &
                               data_old(nx0/2,ny0)
  end if

end subroutine write_real_margins


subroutine write_interior(data_old, data_new, nx_old,nx_new,ny_old,ny_new,&
                          n_margin,s_margin,w_margin,e_margin)

  real(kind=dp),dimension(:,:),intent(in) :: data_old
  real(kind=dp),dimension(:,:),intent(inout) :: data_new

  integer,intent(in) :: nx_old,nx_new,ny_old,ny_new
  integer,intent(in) :: n_margin,s_margin,w_margin,e_margin

  !local vars
  integer :: i_old_left,i_new_left,i_old_right,i_new_right
  integer :: j_old_bot,j_new_bot,j_old_top,j_new_top 
  integer :: i_prev,j_prev
  real(kind=dp) :: a,b
  real(kind=dp),dimension(size(data_new,1),size(data_new,2)) :: data_new_temp
  integer :: i,j,l

  data_new_temp = 0.d0

  i_old_left = w_margin+1
  i_new_left = w_margin+1
  i_old_right = nx_old - e_margin
  i_new_right = nx_new - e_margin
  
  j_old_bot = s_margin + 1
  j_new_bot = s_margin + 1
  j_old_top = ny_old - n_margin
  j_new_top = ny_new - n_margin

  do i=i_new_left,i_new_right
     do j=j_new_bot,j_new_top
 	 
         i_prev = int((real(i) - i_new_left)*(hx_new/hx_old)) + i_old_left
         j_prev = int((real(j) - j_new_bot) *(hy_new/hy_old)) + j_old_bot

         ! a is the fractional number of old cells in the x-direction
         ! from the previous old grid point
         ! to the current new grid point. b is the same, in the y-direction

         a = (real(i)-i_new_left)*(hx_new/hx_old) - &
                aint((real(i)-i_new_left)*(hx_new/hx_old))

	 b = (real(j)-j_new_bot)*(hy_new/hy_old) - &
		aint((real(j)-j_new_bot) *(hy_new/hy_old))

	 data_new(i,j) = (1-a)*(1-b) * data_old(i_prev,j_prev)     + &
	                 (1-a)*   b  * data_old(i_prev,j_prev+1)   + &
                            a *(1-b) * data_old(i_prev+1,j_prev)     + &
	                    a *   b  * data_old(i_prev+1,j_prev+1)
     end do
  end do

  ! smooth data_new twice
  do l=1,3
  do i=(i_new_left+1),(i_new_right-1)
     do j=(j_new_bot+1),(j_new_top-1)
 	 data_new_temp(i,j) = data_new(i,j) + &
              (4.d0*data_new(i,j) &
                  -data_new(i-1,j) &
                  -data_new(i+1,j) &
                  -data_new(i,j-1) &
                  -data_new(i,j+1))*(-0.1d0)

     end do
  end do

  data_new(i_new_left+1:i_new_right-1,j_new_bot+1:j_new_top-1) = &
       data_new_temp(i_new_left+1:i_new_right-1,j_new_bot+1:j_new_top-1)
      
end do 

end subroutine write_interior

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

  call check(nf90_inquire_dimension(nc_id, x_dimid, x1_name, nx_old))
  call check(nf90_inquire_dimension(nc_id, y_dimid, y1_name, ny_old))
  call check(nf90_inquire_dimension(nc_id, time_dimid, t_name, nt_old))
  call check(nf90_inquire_dimension(nc_id, level_dimid, level_name, nlevel))

  call check(nf90_get_var(nc_id, x_varid, x1_old, start= (/ 1 /), &
						  count= (/ 2 /) ))
  hx_old = x1_old(2) - x1_old(1)

  call check(nf90_get_var(nc_id, y_varid, y1_old, start= (/ 1 /), &
      						  count= (/ 2 /) ))
  hy_old = y1_old(2) - y1_old(1)
 
  if (t_read < 1) then
	write(*,*) 'reading last time slice'
	t_read = nt_old
  end if

  if (nt_old < t_read) then
	write(*,*) 't_read exceeds available time slices: ', nt_old
	stop 1
  end if

  allocate(thck_old(nx_old,ny_old),topog_old(nx_old,ny_old))
  allocate(kinbcmask_old(nx_old-1, ny_old-1))
  allocate(uvelhom_old(nx_old-1,ny_old-1,nlevel),vvelhom_old(nx_old-1,ny_old-1,nlevel))
  allocate(temp_old(nx_old, ny_old,nlevel))

  allocate(levels(nlevel))

  call check(nf90_get_var(nc_id, level_varid, levels, &
                          start= (/ 1 /), &
                          count= (/ nlevel /)))

  call check(nf90_get_var(nc_id, thck_varid, thck_old,&
                          start= (/ 1,1,t_read /), &
 			  count= (/ nx_old,ny_old,1 /) ))
  call check(nf90_get_var(nc_id, topog_varid,topog_old, &
                           start= (/ 1,1 /), &
			   count= (/ nx_old,ny_old /) ))
  call check(nf90_get_var(nc_id, kinbcmask_varid, kinbcmask_old, &
			start= (/ 1,1,t_read /), &
		        count= (/ (nx_old-1),(ny_old-1),1 /) ))

  call check(nf90_get_var(nc_id, uvelhom_varid, uvelhom_old, &
                        start= (/ 1,1,1,t_read /), &
                        count= (/ (nx_old-1),(ny_old-1),nlevel,1 /) ))
  call check(nf90_get_var(nc_id, vvelhom_varid, vvelhom_old, &
                        start= (/ 1,1,1,t_read /), &
                        count= (/ (nx_old-1),(ny_old-1),nlevel,1 /) ))
  
  call check(nf90_get_var(nc_id, temp_varid, temp_old, &
                        start= (/ 1,1,1,t_read /), &
                        count= (/ (nx_old),(ny_old), nlevel, 1 /) ))

end subroutine read_old_nc_file


subroutine define_new_data()

  ! local variables
  integer :: i,j
  real(kind=dp),dimension(nx_new-1,ny_new-1) :: kinbcmask_new_real

  hx_new = hx_old * (nx_old - thk_w_margin - thk_e_margin - 1) / &
                    (nx_new - thk_w_margin - thk_e_margin - 1)
  hy_new = hy_old * (ny_old - thk_s_margin - thk_n_margin - 1) / &
                    (ny_new - thk_s_margin - thk_n_margin - 1)

!  hx_new_kin = hx_old * (nx_old - thk_w_margin - thk_e_margin - 1) / &
!                        (nx_new - thk_w_margin - thk_e_margin - 1)
!  hy_new_kin = hy_old * (ny_old - thk_s_margin - thk_n_margin - 1) / &
!                        (ny_new - thk_s_margin - thk_n_margin - 1)

  !now populate the dimension variables
  xs_new = (/ ( x1_old(1) + (i-1)*hx_new,i=1,nx_new ) /)
  ys_new = (/ ( y1_old(1) + (j-1)*hy_new,j=1,ny_new ) /)
  xstag_new = (/ ( x1_old(1) + ((real(i)-0.5)*hx_new),i=1,nx_new-1 ) /)
  ystag_new = (/ ( y1_old(1) + ((real(j)-0.5)*hy_new),j=1,ny_new-1 ) /)
  
  ! now define the thck, topog, kinbcmask arrays
  topog_new = 0.0d0
  thck_new = 0.0d0
  kinbcmask_new_real = 0.0d0
  uvelhom_new = 0.0d0
  vvelhom_new = 0.0d0
  temp_new = 0.0d0

  call write_real_margins(topog_old,topog_new, +1.0d0, 0.0d0,&
                          nx_old, nx_new,ny_old, ny_new, &
	                  top_n_margin,top_s_margin,top_w_margin,top_e_margin)

  call write_interior(topog_old,topog_new, nx_old,nx_new,ny_old,ny_new,&
		      top_n_margin,top_s_margin,top_w_margin,top_e_margin)

  call write_real_margins(thck_old,thck_new, -1.0d0, 0.0d0, &
                          nx_old, nx_new,ny_old, ny_new, &
                          thk_n_margin,thk_s_margin,thk_w_margin,thk_e_margin)

  call write_interior(thck_old,thck_new,nx_old,nx_new,ny_old,ny_new,&
                      thk_n_margin,thk_s_margin,thk_w_margin,thk_e_margin)

  call write_real_margins(kinbcmask_old*1.0d0, kinbcmask_new_real, +1.0d0,0.0d0,&
                          nx_old-1, nx_new-1, ny_old-1, ny_new-1, &
			  kin_n_margin,kin_s_margin,kin_w_margin,kin_e_margin)
  kinbcmask_new = int(kinbcmask_new_real)


  do j=1,nlevel
     call write_real_margins(uvelhom_old(:,:,j),uvelhom_new(:,:,j), 1.0d0, 0.0d0,&
                             nx_old-1, nx_new-1,ny_old-1,ny_new-1, &
                             kin_n_margin,kin_s_margin,kin_w_margin,kin_e_margin)
     call write_real_margins(vvelhom_old(:,:,j),vvelhom_new(:,:,j), 1.0d0, 0.0d0,&
                             nx_old-1, nx_new-1,ny_old-1,ny_new-1, &
                             kin_n_margin,kin_s_margin-1,kin_w_margin,kin_e_margin)
     call write_interior(uvelhom_old(:,:,j),uvelhom_new(:,:,j), nx_old-1,nx_new-1, &
                         ny_old-1,ny_new-1,&
                         kin_n_margin,kin_s_margin,kin_w_margin,kin_e_margin)

     call write_interior(vvelhom_old(:,:,j),vvelhom_new(:,:,j), nx_old-1,nx_new-1, &
                         ny_old-1,ny_new-1, &
                         kin_n_margin,kin_s_margin-1,kin_w_margin,kin_e_margin)

     call write_real_margins(temp_old(:,:,j), temp_new(:,:,j), 1.0d0, 0.0d0, &
                             nx_old, nx_new, ny_old, ny_new, &
                             0,0,0,0)

     call write_interior(temp_old(:,:,j), temp_new(:,:,j), nx_old,nx_new, &
                         ny_old, ny_new, 0,0,0,0)

  end do

end subroutine define_new_data

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
    call check( nf90_def_dim(nc_id,'x1',nx_new,x_dimid) )
    call check( nf90_def_dim(nc_id,'y1',ny_new,y_dimid) )
    call check( nf90_def_dim(nc_id,'x0',nx_new-1,xstag_dimid) )
    call check( nf90_def_dim(nc_id,'y0',ny_new-1,ystag_dimid) )
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

    call check( nf90_put_var(nc_id,x_varid,xs_new) )
    call check( nf90_put_var(nc_id,y_varid,ys_new) )
    call check( nf90_put_var(nc_id,xstag_varid,xstag_new) )
    call check( nf90_put_var(nc_id,ystag_varid,ystag_new) )
    call check( nf90_put_var(nc_id,level_varid,levels) )
    call check( nf90_put_var(nc_id,time_varid, (/ 1 /) ) )

    ! write the arrays out to the netcdf file and close it
    call check( nf90_put_var(nc_id, thck_varid, thck_new) )
    call check( nf90_put_var(nc_id, topog_varid, topog_new) )
    call check( nf90_put_var(nc_id, kinbcmask_varid, kinbcmask_new) )

    call check( nf90_put_var(nc_id, uvelhom_varid, uvelhom_new) )
    call check( nf90_put_var(nc_id, vvelhom_varid, vvelhom_new) )
    
    call check( nf90_put_var(nc_id, temp_varid, temp_new) )

    call check( nf90_close(nc_id) )

  end subroutine write_nc_file

  subroutine check(status_code)

    integer,intent(in) :: status_code

    if (status_code /= 0) then
       write(*,*) 'Fatal netcdf error:',nf90_strerror(status_code)
       stop 1
    end if

  end subroutine check

!!$  subroutine local_linear(raw_data,ocean_mask,short_waves, long_waves, as, bs)
!!$
!!$    implicit none
!!$
!!$    real(kind=kdp),dimension(:,:),intent(in)  :: raw_data
!!$    integer,dimension(:,:),intent(in)         :: ocean_mask
!!$    real(kind=kdp),dimension(:,:),intent(out) ::short_waves, long_waves, as, bs
!!$
!!$    ! raw_data is the field to be filtered
!!$    ! ocean_mask is equal to 1 where the raw_data is valid, 0 otherwise
!!$
!!$    integer :: m,n,i,k
!!$    integer :: lower_ilim,lower_klim, upper_ilim,upper_klim
!!$    integer :: win_imin,win_kmin,win_imax,win_kmax
!!$    integer,parameter :: window_width = 3
!!$
!!$    real(kind=kdp),dimension(window_width,window_width) :: data_window
!!$    integer(kind=kdp),dimension(window_width,window_width) :: ocean_mask_window
!!$
!!$    real(kind=kdp),dimension(window_width,window_width) :: xmat,ymat
!!$    real(kind=kdp),dimension(3,3) :: matrix
!!$    real(kind=kdp),dimension(3)   :: rhs
!!$    real(kind=kdp)                :: multiplier, x0, y0
!!$
!!$    as = 0.d0
!!$    bs = 0.d0
!!$    long_waves = 0.d0
!!$    short_waves = 0.d0
!!$
!!$    m = size(raw_data,1)
!!$    n = size(raw_data,2)
!!$
!!$    x0 = floor( (window_width-1) /2.d0 )
!!$    y0 = floor( (window_width-1) /2.d0 )
!!$
!!$    do i= 1,window_width
!!$       xmat(i,:) = real(i) - (x0+1.d0)
!!$       ymat(:,i) = real(i) - (y0+1.d0)
!!$    end do
!!$
!!$    !$omp parallel default(none) &
!!$    !$omp private(i,k,lower_ilim,lower_klim,upper_ilim,upper_klim, &
!!$    !$omp         win_imin,win_kmin,win_imax,win_kmax, & 
!!$    !$omp         data_window, ocean_mask_window,matrix,rhs, &
!!$    !$omp         multiplier) &
!!$    !$omp shared(raw_data, ocean_mask, short_waves, long_waves, as,bs, &
!!$    !$omp        x0,y0, xmat,ymat, &
!!$    !$omp        domain_imin,domain_imax,domain_kmin,domain_kmax)
!!$    !$omp do
!!$    do i=domain_imin,domain_imax
!!$       do k=domain_kmin,domain_kmax
!!$         
!!$          lower_ilim = max(domain_imin,i-int(floor(  (window_width-1)/2.d0)))
!!$          upper_ilim = min(domain_imax,i+int(ceiling((window_width-1)/2.d0)))
!!$          lower_klim = max(domain_kmin,k-int(floor(  (window_width-1)/2.d0)))
!!$          upper_klim = min(domain_kmax,k+int(ceiling((window_width-1)/2.d0)))
!!$
!!$          data_window = 0.d0
!!$          win_imin = int((lower_ilim-i)+(x0+1.d0))
!!$          win_imax = win_imin + (upper_ilim-lower_ilim)
!!$          win_kmin = int((lower_klim-k)+(y0+1.d0))
!!$          win_kmax = win_kmin + (upper_klim-lower_klim)
!!$
!!$          data_window( win_imin:win_imax, &
!!$                       win_kmin:win_kmax ) = &
!!$               raw_data(lower_ilim:upper_ilim,lower_klim:upper_klim)
!!$
!!$          ocean_mask_window = 0
!!$          ocean_mask_window(win_imin:win_imax, &
!!$                       win_kmin:win_kmax ) = &
!!$               ocean_mask(lower_ilim:upper_ilim,lower_klim:upper_klim)
!!$
!!$          ! try determining the value at (i,k) determined by a local
!!$          ! linear approximation in the data_window
!!$
!!$          if (ocean_mask(i,k) == 1) then
!!$
!!$             matrix(1,1) = sum( xmat * xmat, mask=(ocean_mask_window==1))
!!$             matrix(1,2) = sum( xmat * ymat, mask=(ocean_mask_window==1))
!!$             matrix(1,3) = sum( xmat       , mask=(ocean_mask_window==1))
!!$             matrix(2,2) = sum( ymat * ymat, mask=(ocean_mask_window==1))
!!$             matrix(2,3) = sum( ymat       , mask=(ocean_mask_window==1))
!!$             matrix(3,3) = sum(ocean_mask_window)
!!$             matrix(2,1) = matrix(1,2)
!!$             matrix(3,1) = matrix(1,3)
!!$             matrix(3,2) = matrix(2,3)
!!$
!!$             rhs(1) = sum( xmat * data_window, mask=(ocean_mask_window==1))
!!$             rhs(2) = sum( ymat * data_window, mask=(ocean_mask_window==1))
!!$             rhs(3) = sum( data_window       , mask=(ocean_mask_window==1))
!!$
!!$             ! do Gaussian elim on matrix
!!$             if (matrix(1,1) == 0.d0) then
!!$                print *, i,k
!!$                print *, xmat(:,3)
!!$                print *, xmat(:,2)
!!$                print *, xmat(:,1)
!!$                print *, 'zero entry in 1,1'
!!$                stop 1
!!$             end if
!!$             multiplier = matrix(2,1)/matrix(1,1)
!!$             matrix(2,:) = matrix(2,:) - multiplier*matrix(1,:)
!!$             rhs(2)      = rhs(2) -      multiplier*rhs(1)
!!$
!!$             if (matrix(1,1) == 0.d0) then
!!$                print *, i,k
!!$                print *, matrix(1,:)
!!$                print *, matrix(2,:)
!!$                print *, matrix(3,:)
!!$                print *, 'zero entry in 1,1'
!!$                stop 1
!!$             end if
!!$             multiplier = matrix(3,1)/matrix(1,1)
!!$             matrix(3,:) = matrix(3,:) - multiplier*matrix(1,:)
!!$             rhs(3)      = rhs(3) -      multiplier*rhs(1)
!!$
!!$             if (matrix(2,2) == 0.d0) then
!!$                print *, i,k
!!$                print *, matrix(1,:)
!!$                print *, matrix(2,:)
!!$                print *, matrix(3,:)
!!$                print *, 'zero entry in 2,2'
!!$                stop 1
!!$             end if
!!$             multiplier = matrix(3,2)/matrix(2,2)
!!$             matrix(3,:) = matrix(3,:) - multiplier*matrix(2,:)
!!$             rhs(3)      = rhs(3)      - multiplier*rhs(2)
!!$
!!$             if (matrix(3,3) == 0.d0) then
!!$                print *, i,k
!!$                print *, matrix(1,:)
!!$                print *, matrix(2,:)
!!$                print *, matrix(3,:)
!!$                print *, 'zero entry in 3,3'
!!$                stop 1
!!$             end if
!!$
!!$             long_waves(i,k) = rhs(3) / matrix(3,3)
!!$             bs(i,k) = (rhs(2) - long_waves(i,k)*matrix(2,3)) / matrix(2,2)
!!$             as(i,k) = (rhs(1) - long_waves(i,k)*matrix(1,3) - &
!!$                                    bs(i,k)*matrix(2,1)) / matrix(1,1)
!!$
!!$          else
!!$
!!$             long_waves(i,k) = 0.d0
!!$
!!$          end if
!!$
!!$       end do
!!$    end do
!!$
!!$    !$omp end do
!!$    !$omp end parallel
!!$
!!$    short_waves = 0.d0
!!$    where(ocean_mask == 1)
!!$       short_waves = raw_data - long_waves
!!$    end where
!!$
!!$  end subroutine local_linear

end program nc_regrid
