program last_non_nan_time

  use netcdf

  implicit none

  character(len=1014) :: nc_file_in

  call main()

contains

  subroutine main()

    character(len=512) :: gen_usage = "nc_find_last_nonnan_time <nc_file_in>"
    character(len=512) :: argstr
    integer :: ntime, N

  !local variables  
  integer :: nc_id
  integer :: time_dimid,time_varid,x_dimid,y_dimid,level_dimid
  integer :: x_varid,y_varid,level_varid
  integer :: nx, ny, nt, nlevel
!  integer :: thck_varid,topog_varid,kinbcmask_varid
  integer :: uvelhom_varid,vvelhom_varid,temp_varid
  character(len=NF90_MAX_NAME) :: x1_name,y1_name,t_name,level_name
  integer :: last_good_time = 0

  real(kind=kind(1.d0)),dimension(:,:,:),allocatable :: vels

    if (command_argument_count() .ne. 1) then
     write(*,*) "Incorrect number of arguments.  Usage: ", trim(gen_usage)
     stop 1
   end if

    call get_command_argument(1,argstr)
    read(argstr,'(a)') nc_file_in
    !write(*,*) 'nc_file_in:', trim(nc_file_in)

  call check(nf90_open(trim(nc_file_in), NF90_NOWRITE, nc_id))  

  call check(nf90_inq_dimid(nc_id, 'x1', x_dimid))
  call check(nf90_inq_dimid(nc_id, 'y1', y_dimid))
  call check(nf90_inq_dimid(nc_id, 'time', time_dimid))
  call check(nf90_inq_dimid(nc_id, 'level', level_dimid))

  call check(nf90_inq_varid(nc_id, 'x1', x_varid))
  call check(nf90_inq_varid(nc_id, 'y1', y_varid))
  call check(nf90_inq_varid(nc_id, 'time', time_varid))
  call check(nf90_inq_varid(nc_id, 'level', level_varid))

  call check(nf90_inq_varid(nc_id, 'uvelhom', uvelhom_varid))
  call check(nf90_inq_varid(nc_id, 'vvelhom', vvelhom_varid))

  call check(nf90_inquire_dimension(nc_id, x_dimid, x1_name, nx))
  call check(nf90_inquire_dimension(nc_id, y_dimid, y1_name, ny))
  call check(nf90_inquire_dimension(nc_id, time_dimid, t_name, nt))
  call check(nf90_inquire_dimension(nc_id, level_dimid, level_name, nlevel))

  call check(nf90_inquire_dimension(nc_id, time_dimid, t_name, N))

  allocate(vels(nlevel,ny,nx))

  do n=1,N
     !load the velocities at the current time index
     call check(nf90_get_var(nc_id, uvelhom_varid, vels, &
                             start= (/ 1,1,1,n /), &
                             count= (/ (nx-1),(ny-1),nlevel,1 /) ))
     if (.not.(any(isnan(vels)))) then
      call check(nf90_get_var(nc_id, vvelhom_varid, vels, &
                             start= (/ 1,1,1,n /), &
                             count= (/ (nx-1),(ny-1),nlevel,1 /) ))
      if (.not.(any(isnan(vels)))) then
         last_good_time = last_good_time + 1
      else
         exit
      end if
     else
        exit
     end if


  end do 


  deallocate(vels)

  print *, last_good_time-1

end subroutine main

subroutine check(status_code)

    integer,intent(in) :: status_code

    if (status_code /= 0) then
       write(*,*) 'Fatal netcdf error:',nf90_strerror(status_code)
       stop 1
    end if

end subroutine check

end program last_non_nan_time

