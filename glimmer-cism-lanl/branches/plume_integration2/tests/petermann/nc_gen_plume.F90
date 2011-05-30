program nc_gen_plume

  use netcdf

  implicit none

  real,parameter :: pi = 3.1415926
  real,parameter :: g = 9.81
  !real,parameter :: rhoi_rhow = 1030.d0 / 920.d0

  character(len=512) :: gen_usage = "nc_gen_plume  &
       &<type_code> = [fb] &
       &nc_filename [params] "

  character(len=512) :: fixed_bmlt_params = "<fname> <nx> <ny> <hx> <hy>"

  character (len=2) :: type_code

  !variables to be used in all cases and in netcdf writing
  character(len=1024) :: argstr
  character (len=2056) :: fname
  integer :: nx,ny

  !data arrays
  real,dimension(:),allocatable :: xs,ys
  real,dimension(:,:),allocatable :: bmlt, dummy_zeros
  integer,dimension(:,:),allocatable :: dummy_zeros_int

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

  if (type_code == 'fb') then

     call make_fixed_bmlt_plume()

  else

     write(*,*) "Unrecognized type_code: ", type_code

  end if

contains

  subroutine make_fixed_bmlt_plume()

    ! local variables
    integer :: i,j,j_interp
    integer,parameter :: NCONTROL = 5

    real :: hx,hy
    real :: bmelt_interp

    real,dimension(NCONTROL) :: y_control = (/ 0.d0, 8.d0, 28.d0,47.d0,70.d0 /)*1000.d0
    real,dimension(NCONTROL) :: bmlt_control = (/ 0.d0,24.d0,1.d0,5.d0,0.d0 /)

    character (len=512) :: fixed_bmlt_params = "<fname> <nx> <ny> <hx> <hy>"

    if (command_argument_count() < 6) then
       write(*,*) "Incorrect number of parameters. Fixed bmlt requires:"
       write(*,*) trim(fixed_bmlt_params)
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
    read(argstr,'(f18.12)') hx
    write(*,*) 'hx',hx

    call get_command_argument(6,argstr)
    read(argstr,'(f18.12)') hy
    write(*,*) 'hy',hy

    allocate(xs(nx),ys(ny))
    allocate(bmlt(nx,ny),dummy_zeros(nx,ny))
    allocate(dummy_zeros_int(nx,ny))

    !now populate the dimension variables
    xs = (/ ( (real(i)-3.5d0)*hx,i=1,nx ) /)
    ys = (/ ( (real(j)-3.5d0)*hy,j=1,ny ) /)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! define the bmlt      !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    dummy_zeros = 0.d0
    dummy_zeros_int = 0

    bmlt = 0.d0

    ! Note: j=3.5 corresponds to y = 0, which is the domain boundary

    do j=4,ny-4
       do i=1,(NCONTROL-1)
          if ((y_control(i+1) > ys(j)) .and. (y_control(i) <= ys(j))) then
             j_interp = i
          end if
       end do
       
       bmelt_interp = bmlt_control(j_interp) + &
            (ys(j)-y_control(j_interp))/(y_control(j_interp+1)-y_control(j_interp))* &
            (bmlt_control(j_interp+1) - bmlt_control(j_interp))

       bmlt(1+3:(nx-3),j) = bmelt_interp

    enddo
    

    call write_nc_file()

    deallocate(xs,ys)
    deallocate(bmlt,dummy_zeros,dummy_zeros_int)

  end subroutine make_fixed_bmlt_plume

  subroutine write_nc_file()

    !local variables
    integer :: nc_id
    integer :: time_dimid,x_dimid,y_dimid
    integer :: x_varid,y_varid,time_varid
    integer :: bmlt_varid
    integer :: jcs_varid,jcw_varid,jcd_u_varid,jcd_v_varid,jcd_fl_varid
    integer :: bpos_varid, ipos_varid, pdep_varid
    integer :: u_varid, v_varid, su_varid, sv_varid
    integer :: bmelt_varid, btemp_varid, bsalt_varid
    integer :: rhop_varid, temp_varid, salt_varid, entr_varid

    call check( nf90_create(fname, NF90_CLOBBER, nc_id) )

    ! defining dimensions
    call check( nf90_def_dim(nc_id,'time',1,time_dimid) )
    call check( nf90_def_dim(nc_id,'x',nx,x_dimid) )
    call check( nf90_def_dim(nc_id,'y',ny,y_dimid) )

    ! define variables
    call check( nf90_def_var(nc_id,'time',NF90_DOUBLE,(/time_dimid/),time_varid) )
    call check( nf90_put_att(nc_id, time_varid, 'long_name', 'time') )
    call check( nf90_put_att(nc_id, time_varid, 'units', 'seconds') )

    call check( nf90_def_var(nc_id,'x',NF90_DOUBLE,(/x_dimid/),x_varid) )
    call check( nf90_put_att(nc_id, x_varid, 'long_name', 'Cartisian x-coordinate') )
    call check( nf90_put_att(nc_id, x_varid, 'units', 'meter') )

    call check( nf90_def_var(nc_id,'y',NF90_DOUBLE,(/y_dimid/),y_varid) )
    call check( nf90_put_att(nc_id, y_varid, 'long_name', 'Cartisian y-coordinate') )
    call check( nf90_put_att(nc_id, y_varid, 'units', 'meter') )

    call check( nf90_def_var(nc_id,'bmelt',NF90_DOUBLE,(/x_dimid,y_dimid,time_dimid/),bmlt_varid) )
    call check( nf90_put_att(nc_id, bmlt_varid, 'long_name', 'basal melt rate') )
    call check( nf90_put_att(nc_id, bmlt_varid, 'units', 'meters / year') )


    call check( nf90_def_var(nc_id,'jcs',NF90_INT,(/x_dimid,y_dimid,time_dimid/),jcs_varid) )
    call check( nf90_def_var(nc_id,'jcw',NF90_INT,(/x_dimid,y_dimid,time_dimid/),jcw_varid) )
    call check( nf90_def_var(nc_id,'jcd_u',NF90_INT,(/x_dimid,y_dimid,time_dimid/),jcd_u_varid) )
    call check( nf90_def_var(nc_id,'jcd_v',NF90_INT,(/x_dimid,y_dimid,time_dimid/),jcd_v_varid) )
    call check( nf90_def_var(nc_id,'jcd_fl',NF90_INT,(/x_dimid,y_dimid,time_dimid/),jcd_fl_varid) )
    call check( nf90_def_var(nc_id,'bpos',NF90_DOUBLE,(/x_dimid,y_dimid,time_dimid/),bpos_varid) )
    call check( nf90_def_var(nc_id,'ipos',NF90_DOUBLE,(/x_dimid,y_dimid,time_dimid/),ipos_varid) )
    call check( nf90_def_var(nc_id,'pdep',NF90_DOUBLE,(/x_dimid,y_dimid,time_dimid/),pdep_varid) )
    call check( nf90_def_var(nc_id,'u',NF90_DOUBLE,(/x_dimid,y_dimid,time_dimid/),u_varid) )
    call check( nf90_def_var(nc_id,'v',NF90_DOUBLE,(/x_dimid,y_dimid,time_dimid/),v_varid) )
    call check( nf90_def_var(nc_id,'su',NF90_DOUBLE,(/x_dimid,y_dimid,time_dimid/),su_varid) )
    call check( nf90_def_var(nc_id,'sv',NF90_DOUBLE,(/x_dimid,y_dimid,time_dimid/),sv_varid) )
    call check( nf90_def_var(nc_id,'btemp',NF90_DOUBLE,(/x_dimid,y_dimid,time_dimid/),btemp_varid) )
    call check( nf90_def_var(nc_id,'bsalt',NF90_DOUBLE,(/x_dimid,y_dimid,time_dimid/),bsalt_varid) )
    call check( nf90_def_var(nc_id,'rhop',NF90_DOUBLE,(/x_dimid,y_dimid,time_dimid/),rhop_varid) )
    call check( nf90_def_var(nc_id,'temp',NF90_DOUBLE,(/x_dimid,y_dimid,time_dimid/),temp_varid) )
    call check( nf90_def_var(nc_id,'salt',NF90_DOUBLE,(/x_dimid,y_dimid,time_dimid/),salt_varid) )
    call check( nf90_def_var(nc_id,'entr',NF90_DOUBLE,(/x_dimid,y_dimid,time_dimid/),entr_varid) )


    call check( nf90_enddef(nc_id) )

    ! populate variables

    call check( nf90_put_var(nc_id,x_varid,xs) )
    call check( nf90_put_var(nc_id,y_varid,ys) )
    call check( nf90_put_var(nc_id,time_varid, (/ 1 /) ) )

    ! write the arrays out to the netcdf file and close it
    call check( nf90_put_var(nc_id, bmlt_varid, bmlt) )

    call check( nf90_put_var(nc_id, jcs_varid, dummy_zeros_int) )
    call check( nf90_put_var(nc_id, jcw_varid, dummy_zeros_int) )
    call check( nf90_put_var(nc_id, jcd_u_varid, dummy_zeros_int) )
    call check( nf90_put_var(nc_id, jcd_v_varid, dummy_zeros_int) )
    call check( nf90_put_var(nc_id, jcd_fl_varid, dummy_zeros_int) )

    call check( nf90_put_var(nc_id, bpos_varid, dummy_zeros) )
    call check( nf90_put_var(nc_id, ipos_varid, dummy_zeros) )
    call check( nf90_put_var(nc_id, pdep_varid, dummy_zeros) )
    call check( nf90_put_var(nc_id, u_varid, dummy_zeros) )
    call check( nf90_put_var(nc_id, v_varid, dummy_zeros) )
    call check( nf90_put_var(nc_id, su_varid, dummy_zeros) )
    call check( nf90_put_var(nc_id, sv_varid, dummy_zeros) )
    call check( nf90_put_var(nc_id, btemp_varid, dummy_zeros) )
    call check( nf90_put_var(nc_id, bsalt_varid, dummy_zeros) )
    call check( nf90_put_var(nc_id, rhop_varid, dummy_zeros) )
    call check( nf90_put_var(nc_id, temp_varid, dummy_zeros) )
    call check( nf90_put_var(nc_id, salt_varid, dummy_zeros) )
    call check( nf90_put_var(nc_id, entr_varid, dummy_zeros) )

    call check( nf90_close(nc_id) )

  end subroutine write_nc_file

  subroutine check(status_code)

    integer,intent(in) :: status_code

    if (status_code /= 0) then
       write(*,*) 'fatal netcdf error:',nf90_strerror(status_code)
       stop 1
    end if

  end subroutine check

end program nc_gen_plume
