module parallel
  use netcdf,&
       distributed_get_var=>nf90_get_var,&
       distributed_put_var=>nf90_put_var,&
       parallel_close=>nf90_close,&
       parallel_create=>nf90_create,&
       parallel_def_dim=>nf90_def_dim,&
       parallel_def_var=>nf90_def_var,&
       parallel_enddef=>nf90_enddef,&
       parallel_get_att=>nf90_get_att,&
       parallel_get_var=>nf90_get_var,&
       parallel_inq_attname=>nf90_inq_attname,&
       parallel_inq_dimid=>nf90_inq_dimid,&
       parallel_inq_varid=>nf90_inq_varid,&
       parallel_inquire=>nf90_inquire,&
       parallel_inquire_dimension=>nf90_inquire_dimension,&
       parallel_inquire_variable=>nf90_inquire_variable,&
       parallel_put_att=>nf90_put_att,&
       parallel_put_var=>nf90_put_var,&
       parallel_open=>nf90_open,&
       parallel_redef=>nf90_redef,&
       parallel_sync=>nf90_sync
  implicit none

  integer,parameter :: lhalo = 0
  integer,parameter :: uhalo = 0
  logical,parameter :: main_task = .true.
  integer,parameter :: this_rank = 0

  integer,parameter :: staggered_lhalo = lhalo
  integer,parameter :: staggered_uhalo = 0

  integer,save :: global_ewn,global_nsn

  interface broadcast
     module procedure broadcast_character
     module procedure broadcast_integer
     module procedure broadcast_integer_1d
     module procedure broadcast_logical
     module procedure broadcast_real8_1d
  end interface

  interface parallel_halo
     module procedure parallel_halo_integer_2d
     module procedure parallel_halo_real8_2d
     module procedure parallel_halo_real8_3d
  end interface

  interface parallel_halo_verify
     module procedure parallel_halo_verify_integer_2d
     module procedure parallel_halo_verify_real8_2d
     module procedure parallel_halo_verify_real8_3d
  end interface

  interface parallel_print
     module procedure parallel_print_integer_2d
     module procedure parallel_print_real8_2d
     module procedure parallel_print_real8_3d
  end interface

contains

  subroutine broadcast_character(c)
    implicit none
    character(len=*) :: c
  end subroutine

  subroutine broadcast_integer(i)
    implicit none
    integer :: i
  end subroutine

  subroutine broadcast_integer_1d(a)
    implicit none
    integer,dimension(:) :: a
  end subroutine

  subroutine broadcast_logical(l)
    implicit none
    logical :: l
  end subroutine

  subroutine broadcast_real8_1d(a)
    implicit none
    real(8),dimension(:) :: a
  end subroutine

  subroutine distributed_grid(ewn,nsn)
    implicit none
    integer :: ewn,nsn
    ! begin
    global_ewn = ewn
    global_nsn = nsn
  end subroutine

  function distributed_owner(ew,ewn,ns,nsn)
    implicit none
    logical :: distributed_owner
    integer :: ew,ewn,ns,nsn
    ! begin
    distributed_owner = .true.
  end function

  subroutine global_sum(x,y)
    implicit none
    real(8) :: x,y
  end subroutine

  subroutine not_parallel(file,line)
    implicit none
    integer :: line
    character(len=*) :: file
    ! begin
    write(0,*) "WARNING: not parallel in ",file," at line ",line
  end subroutine

  subroutine parallel_barrier
    implicit none
  end subroutine

  function parallel_boundary(ew,ewn,ns,nsn)
    implicit none
    logical :: parallel_boundary
    integer :: ew,ewn,ns,nsn
    ! begin
    parallel_boundary = (ew==1.or.ew==ewn.or.ns==1.or.ns==nsn)
  end function

  subroutine parallel_finalise
    implicit none
  end subroutine

  subroutine parallel_halo_integer_2d(a)
    implicit none
    integer,dimension(:,:) :: a
  end subroutine

  subroutine parallel_halo_real8_2d(a)
    implicit none
    real(8),dimension(:,:) :: a
  end subroutine

  subroutine parallel_halo_real8_3d(a)
    implicit none
    real(8),dimension(:,:,:) :: a
  end subroutine

  subroutine parallel_halo_verify_integer_2d(a)
    implicit none
    integer,dimension(:,:) :: a
  end subroutine

  subroutine parallel_halo_verify_real8_2d(a)
    implicit none
    real(8),dimension(:,:) :: a
  end subroutine

  subroutine parallel_halo_verify_real8_3d(a)
    implicit none
    real(8),dimension(:,:,:) :: a
  end subroutine

  subroutine parallel_initialise
    implicit none
  end subroutine

  subroutine parallel_print_integer_2d(name,a)
    implicit none
    character(*) :: name
    integer,dimension(:,:) :: a
    
    integer,parameter :: u = 33
    integer :: i,j
    ! begin
    open(unit=u,file=name//".txt",form="formatted",status="replace")
    do j = 1,size(a,2)
       do i = 1,size(a,1)
          write(u,*) j,i,a(i,j)
       end do
       write(u,'()')
    end do
    close(u)
  end subroutine

  subroutine parallel_print_real8_2d(name,a)
    implicit none
    character(*) :: name
    real(8),dimension(:,:) :: a
    
    integer,parameter :: u = 33
    integer :: i,j
    ! begin
    open(unit=u,file=name//".txt",form="formatted",status="replace")
    do j = 1,size(a,2)
       do i = 1,size(a,1)
          write(u,*) j,i,a(i,j)
       end do
       write(u,'()')
    end do
    close(u)
  end subroutine

  subroutine parallel_print_real8_3d(name,a)
    implicit none
    character(*) :: name
    real(8),dimension(:,:,:) :: a
    
    integer,parameter :: u = 33
    integer :: i,j
    ! begin
    open(unit=u,file=name//".txt",form="formatted",status="replace")
    do j = 1,size(a,3)
       do i = 1,size(a,2)
          write(u,'(2i6,100g15.5e3)') j,i,a(:,i,j)
       end do
       write(u,'()')
    end do
    close(u)
  end subroutine

  subroutine parallel_show_minmax(label,values)
    implicit none
    character(*) :: label
    real(8),dimension(:,:,:) :: values
    ! begin
    print *,label,minval(values),maxval(values)
  end subroutine

  subroutine parallel_stop(file,line)
    implicit none
    integer :: line
    character(len=*) :: file
    ! begin
    write(0,*) "STOP in ",file," at line ",line
    ! stop
    write(0,*) "RUNNING in parallel_single mode, so STOP IGNORED."

  end subroutine

  subroutine parallel_temp_halo(a)
    implicit none
    real(8),dimension(:,:,:) :: a
  end subroutine

  subroutine parallel_velo_halo(a)
    implicit none
    real(8),dimension(:,:) :: a
  end subroutine

end module parallel
