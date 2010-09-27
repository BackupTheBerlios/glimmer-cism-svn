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
  integer,parameter :: main_rank = 0

  logical,save :: main_task
  integer,save :: comm,tasks,this_rank

  ! distributed grid
  integer,save :: global_ewn,global_nsn,local_ewn,local_nsn,own_ewn,own_nsn
  integer,save :: ewlb,ewub,nslb,nsub

  integer,parameter :: staggered_lhalo = lhalo
  integer,parameter :: staggered_uhalo = 0

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
    use mpi
    implicit none
    character(len=*) :: c
    integer :: ierror,n
    ! begin
    n = len(c)
    call mpi_bcast(c,n,mpi_character,main_rank,comm,ierror)
  end subroutine

  subroutine broadcast_integer(i)
    use mpi
    implicit none
    integer :: i,ierror
    ! begin
    call mpi_bcast(i,1,mpi_integer,main_rank,comm,ierror)
  end subroutine

  subroutine broadcast_integer_1d(a)
    implicit none
    integer,dimension(:) :: a

    if (tasks >= 2) then
       call not_parallel(__FILE__, __LINE__)
    endif 
  end subroutine

  subroutine broadcast_logical(l)
    use mpi
    implicit none
    logical :: l
    integer :: ierror
    ! begin
    call mpi_bcast(l,1,mpi_logical,main_rank,comm,ierror)
  end subroutine

  subroutine broadcast_real8_1d(a)
    implicit none
    real(8),dimension(:) :: a

    if (tasks >= 2) then
       call not_parallel(__FILE__, __LINE__)
    endif 
  end subroutine

  subroutine distributed_grid(ewn,nsn)
    implicit none
    integer :: ewn,nsn
    ! begin
    global_ewn = ewn
    global_nsn = nsn

    if (tasks >= 2) then
		ewlb = 1
		ewub = global_ewn
		local_ewn = ewub-ewlb+1
		own_ewn = local_ewn-lhalo-uhalo
		ewn = local_ewn
	
		nslb = 1
		nsub = global_nsn
		local_nsn = nsub-nslb+1
		own_nsn = local_nsn-lhalo-uhalo
		nsn = local_nsn
    endif 
  end subroutine

  function distributed_owner(ew,ewn,ns,nsn)
    implicit none
    logical :: distributed_owner
    integer :: ew,ewn,ns,nsn
    ! begin
    distributed_owner = .true.

    if (tasks >= 2) then
       ! Do nothing
    endif 
  end function

  subroutine global_sum(x,y)
    implicit none
    real(8) :: x,y

    if (tasks >= 2) then
        ! Do nothing.  x and y are two values to be summed across all distributed nodes.
        ! In parallel_single model, their local sums should remain unchanged.
    endif 
  end subroutine

  subroutine not_parallel(file,line)
    implicit none
    integer :: line
    character(len=*) :: file
    ! begin
    write(0,*) "WARNING: not parallel in ",file," at line ",line
    STOP
  end subroutine

  subroutine parallel_barrier
    implicit none

    if (tasks >= 2) then
       call not_parallel(__FILE__, __LINE__)
    endif 
  end subroutine

  function parallel_boundary(ew,ewn,ns,nsn)
    implicit none
    logical :: parallel_boundary
    integer :: ew,ewn,ns,nsn
    ! begin
    parallel_boundary = (ew==1.or.ew==ewn.or.ns==1.or.ns==nsn)

    if (tasks >= 2) then
		! Do nothing.
		! The value of parallel_boundary is correct for parallel_single mode.
    endif 
  end function

  subroutine parallel_finalise
    use mpi
    implicit none
    integer :: ierror 
    ! begin 
    call mpi_finalize(ierror)
  end subroutine

  subroutine parallel_halo_integer_2d(a)
    implicit none
    integer,dimension(:,:) :: a

    if (tasks >= 2) then
		! Do nothing
	endif 
  end subroutine

  subroutine parallel_halo_real8_2d(a)
    implicit none
    real(8),dimension(:,:) :: a

    if (tasks >= 2) then
		! Do nothing
	endif 
  end subroutine

  subroutine parallel_halo_real8_3d(a)
    implicit none
    real(8),dimension(:,:,:) :: a

    if (tasks >= 2) then
		! Do nothing
	endif 
  end subroutine

  subroutine parallel_halo_verify_integer_2d(a)
    implicit none
    integer,dimension(:,:) :: a

    if (tasks >= 2) then
		! Do nothing
	endif 
  end subroutine

  subroutine parallel_halo_verify_real8_2d(a)
    implicit none
    real(8),dimension(:,:) :: a

    if (tasks >= 2) then
		! Do nothing
	endif 
  end subroutine

  subroutine parallel_halo_verify_real8_3d(a)
    implicit none
    real(8),dimension(:,:,:) :: a

    if (tasks >= 2) then
		! Do nothing
	endif 
  end subroutine

  subroutine parallel_initialise
    use mpi 
    implicit none
    integer :: ierror 
    ! begin 
    call mpi_init(ierror)
    comm = mpi_comm_world
    call mpi_comm_size(comm,tasks,ierror)
    call mpi_comm_rank(comm,this_rank,ierror)
    main_task = (this_rank==main_rank)
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

    if (tasks >= 2) then
       call not_parallel(__FILE__, __LINE__)
    endif 
  end subroutine

  subroutine parallel_velo_halo(a)
    implicit none
    real(8),dimension(:,:) :: a

    if (tasks >= 2) then
       call not_parallel(__FILE__, __LINE__)
    endif 
  end subroutine

end module parallel
