module parallel
  use netcdf
  implicit none

  integer,parameter :: lhalo = 2
  integer,parameter :: main_rank = 0
  integer,parameter :: uhalo = 2

  logical,save :: main_task
  integer,save :: comm,tasks,this_rank

  ! distributed grid
  integer,save :: global_ewn,global_nsn,local_ewn,local_nsn,own_ewn,own_nsn
  integer,save :: ewlb,ewub,nslb,nsub
  integer,save :: east,north,south,west

  interface broadcast
     module procedure broadcast_character
     module procedure broadcast_integer
     module procedure broadcast_integer_1d
     module procedure broadcast_logical
     module procedure broadcast_real4
     module procedure broadcast_real4_1d
     module procedure broadcast_real8     
     module procedure broadcast_real8_1d
  end interface

  interface distributed_get_var
     module procedure distributed_get_var_integer_2d
     module procedure distributed_get_var_real4_1d
     module procedure distributed_get_var_real4_2d
     module procedure distributed_get_var_real8_2d
     module procedure distributed_get_var_real8_3d
  end interface

  interface distributed_put_var
     module procedure distributed_put_var_integer_2d
     module procedure distributed_put_var_real4_1d
     module procedure distributed_put_var_real4_2d
     module procedure distributed_put_var_real8_2d
     module procedure distributed_put_var_real8_3d
     module procedure parallel_put_var_real8
  end interface

  interface parallel_def_var
     module procedure parallel_def_var_dimids
     module procedure parallel_def_var_nodimids
  end interface

  interface parallel_get_att
     module procedure parallel_get_att_character
     module procedure parallel_get_att_real4
     module procedure parallel_get_att_real4_1d
     module procedure parallel_get_att_real8
     module procedure parallel_get_att_real8_1d
  end interface

  interface parallel_get_var
     module procedure parallel_get_var_integer_1d
     module procedure parallel_get_var_real4_1d
  end interface

  interface parallel_halo
     module procedure parallel_halo_integer_2d
     module procedure parallel_halo_real8_2d
     module procedure parallel_halo_real8_3d
  end interface

  interface parallel_print
     module procedure parallel_print_integer_2d
     module procedure parallel_print_real8_2d
     module procedure parallel_print_real8_3d
  end interface

  interface parallel_put_att
     module procedure parallel_put_att_character
     module procedure parallel_put_att_real4
     module procedure parallel_put_att_real4_1d
     module procedure parallel_put_att_real8
     module procedure parallel_put_att_real8_1d
  end interface

  interface parallel_put_var
     module procedure parallel_put_var_real4
     module procedure parallel_put_var_real8_1d
  end interface

  ! not used and/or wrong?
  interface parallel_velo_halo
     module procedure parallel_velo_halo_real8_2d
     module procedure parallel_velo_halo_real8_3d
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
    use mpi
    implicit none
    integer,dimension(:) :: a
    integer :: ierror
    ! begin
    call mpi_bcast(a,size(a),mpi_integer,main_rank,comm,ierror)
  end subroutine

  subroutine broadcast_logical(l)
    use mpi
    implicit none
    logical :: l
    integer :: ierror
    ! begin
    call mpi_bcast(l,1,mpi_logical,main_rank,comm,ierror)
  end subroutine

  subroutine broadcast_real4(r)
    use mpi
    implicit none
    integer :: ierror
    real(4) :: r
    ! begin
    call mpi_bcast(r,1,mpi_real4,main_rank,comm,ierror)
  end subroutine

  subroutine broadcast_real4_1d(a)
    use mpi
    implicit none
    real(4),dimension(:) :: a
    integer :: ierror
    ! begin
    call mpi_bcast(a,size(a),mpi_real4,main_rank,comm,ierror)
  end subroutine

  subroutine broadcast_real8(r)
    use mpi
    implicit none
    integer :: ierror
    real(8) :: r
    ! begin
    call mpi_bcast(r,1,mpi_real8,main_rank,comm,ierror)
  end subroutine

  subroutine broadcast_real8_1d(a)
    use mpi
    implicit none
    real(8),dimension(:) :: a
    integer :: ierror
    ! begin
    call mpi_bcast(a,size(a),mpi_real8,main_rank,comm,ierror)
  end subroutine

  function distributed_get_var_integer_2d(ncid,varid,values,start)
    use mpi
    implicit none
    integer :: distributed_get_var_integer_2d,ncid,varid
    integer,dimension(:) :: start
    integer,dimension(:,:) :: values

    integer :: ew,i,ierror,ns
    integer,dimension(4) :: mybounds
    integer,dimension(:),allocatable :: displs,sendcounts
    integer,dimension(:,:),allocatable :: bounds
    integer,dimension(:),allocatable :: sendbuf
    integer,dimension(:,:),allocatable :: global_values,recvbuf

    ! begin

    if (size(values,1)==local_ewn) then
       ew = global_ewn
       ns = global_nsn
    else if (size(values,1)==local_ewn-1) then
       ew = global_ewn-1
       ns = global_nsn-1
    else
       call parallel_stop(__FILE__,__LINE__)
    end if
    mybounds(1) = ewlb
    mybounds(2) = ewub
    mybounds(3) = nslb
    mybounds(4) = nsub
    if (main_task) then
       allocate(bounds(4,tasks))
    else
       allocate(bounds(1,1))
    end if
    call mpi_gather(mybounds,4,mpi_integer,bounds,4,mpi_integer,&
         main_rank,comm,ierror)
    if (main_task) then
       allocate(global_values(minval(bounds(1,:)):maxval(bounds(2,:)),&
            minval(bounds(3,:)):maxval(bounds(4,:))))
       global_values(:,:) = 0
       distributed_get_var_integer_2d = nf90_get_var(ncid,varid,&
            global_values(1:ew,1:ns),start)
       allocate(displs(tasks+1))
       allocate(sendcounts(tasks))
       sendcounts(:) = (bounds(2,:)-bounds(1,:)+1)*(bounds(4,:)-bounds(3,:)+1)
       displs(1) = 0
       do i = 1,tasks
          displs(i+1) = displs(i)+sendcounts(i)
       end do
       allocate(sendbuf(displs(tasks+1)))
       do i = 1,tasks
          sendbuf(displs(i)+1:displs(i+1)) = reshape(&
               global_values(bounds(1,i):bounds(2,i),bounds(3,i):bounds(4,i)),&
               (/displs(i+1)-displs(i)/))
       end do
    else
       allocate(displs(1))
       allocate(sendcounts(1))
       allocate(sendbuf(1))
    end if
    call broadcast(distributed_get_var_integer_2d)
    allocate(recvbuf(local_ewn,local_nsn))
    call mpi_scatterv(sendbuf,sendcounts,displs,mpi_integer,&
         recvbuf,size(recvbuf),mpi_integer,main_rank,comm,ierror)
    values(:,:) = recvbuf(:size(values,1),:size(values,2))
    !automatic deallocation
  end function

  function distributed_get_var_real4_1d(ncid,varid,values,start)
    use mpi
    use netcdf
    implicit none
    integer :: distributed_get_var_real4_1d,ncid,varid
    integer,dimension(:) :: start
    real(4),dimension(:) :: values

    integer :: i,ierror,myn,status,x1id,y1id
    integer,dimension(2) :: mybounds
    integer,dimension(:),allocatable :: displs,sendcounts
    integer,dimension(:,:),allocatable :: bounds
    real(4),dimension(:),allocatable :: global_values,sendbuf

    ! begin

    if (main_task) then
       allocate(bounds(2,tasks))
       status = nf90_inq_varid(ncid,"x1",x1id)
       status = nf90_inq_varid(ncid,"y1",y1id)
    else
       allocate(bounds(1,1))
    end if
    call broadcast(x1id)
    call broadcast(y1id)
    if (varid==x1id) then
       mybounds(1) = ewlb
       mybounds(2) = ewub
       myn = global_ewn
    else if (varid==y1id) then
       mybounds(1) = nslb
       mybounds(2) = nsub
       myn = global_nsn
    else
       call parallel_stop(__FILE__,__LINE__)
    end if
    call mpi_gather(mybounds,2,mpi_integer,bounds,2,mpi_integer,&
         main_rank,comm,ierror)
    if (main_task) then
       allocate(global_values(minval(bounds(1,:)):maxval(bounds(2,:))))
       global_values(:) = 0
       distributed_get_var_real4_1d = &
            nf90_get_var(ncid,varid,global_values(1:myn),start)
       allocate(displs(tasks+1))
       allocate(sendcounts(tasks))
       sendcounts(:) = bounds(2,:)-bounds(1,:)+1
       displs(1) = 0
       do i = 1,tasks
          displs(i+1) = displs(i)+sendcounts(i)
       end do
       allocate(sendbuf(displs(tasks+1)))
       do i = 1,tasks
          sendbuf(displs(i)+1:displs(i+1)) = &
               global_values(bounds(1,i):bounds(2,i))
       end do
    else
       allocate(displs(1))
       allocate(sendcounts(1))
       allocate(sendbuf(1))
    end if
    call broadcast(distributed_get_var_real4_1d)
    call mpi_scatterv(sendbuf,sendcounts,displs,mpi_real4,&
         values,size(values),mpi_real4,main_rank,comm,ierror)
    !automatic deallocation
  end function

  function distributed_get_var_real4_2d(ncid,varid,values,start)
    use mpi
    implicit none
    integer :: distributed_get_var_real4_2d,ncid,varid
    integer,dimension(:) :: start
    real(4),dimension(:,:) :: values

    integer :: ew,i,ierror,ns
    integer,dimension(4) :: mybounds
    integer,dimension(:),allocatable :: displs,sendcounts
    integer,dimension(:,:),allocatable :: bounds
    real(4),dimension(:),allocatable :: sendbuf
    real(4),dimension(:,:),allocatable :: global_values,recvbuf

    ! begin

    if (size(values,1)==local_ewn) then
       ew = global_ewn
       ns = global_nsn
    else if (size(values,1)==local_ewn-1) then
       ew = global_ewn-1
       ns = global_nsn-1
    else
       call parallel_stop(__FILE__,__LINE__)
    end if
    mybounds(1) = ewlb
    mybounds(2) = ewub
    mybounds(3) = nslb
    mybounds(4) = nsub
    if (main_task) then
       allocate(bounds(4,tasks))
    else
       allocate(bounds(1,1))
    end if
    call mpi_gather(mybounds,4,mpi_integer,bounds,4,mpi_integer,&
         main_rank,comm,ierror)
    if (main_task) then
       allocate(global_values(minval(bounds(1,:)):maxval(bounds(2,:)),&
            minval(bounds(3,:)):maxval(bounds(4,:))))
       global_values(:,:) = 0
       distributed_get_var_real4_2d = nf90_get_var(ncid,varid,&
            global_values(1:ew,1:ns),start)
       allocate(displs(tasks+1))
       allocate(sendcounts(tasks))
       sendcounts(:) = (bounds(2,:)-bounds(1,:)+1)*(bounds(4,:)-bounds(3,:)+1)
       displs(1) = 0
       do i = 1,tasks
          displs(i+1) = displs(i)+sendcounts(i)
       end do
       allocate(sendbuf(displs(tasks+1)))
       do i = 1,tasks
          sendbuf(displs(i)+1:displs(i+1)) = reshape(&
               global_values(bounds(1,i):bounds(2,i),bounds(3,i):bounds(4,i)),&
               (/displs(i+1)-displs(i)/))
       end do
    else
       allocate(displs(1))
       allocate(sendcounts(1))
       allocate(sendbuf(1))
    end if
    call broadcast(distributed_get_var_real4_2d)
    allocate(recvbuf(local_ewn,local_nsn))
    call mpi_scatterv(sendbuf,sendcounts,displs,mpi_real4,&
         recvbuf,size(recvbuf),mpi_real4,main_rank,comm,ierror)
    values(:,:) = recvbuf(:size(values,1),:size(values,2))
    !automatic deallocation
  end function

  function distributed_get_var_real8_2d(ncid,varid,values,start)
    use mpi
    implicit none
    integer :: distributed_get_var_real8_2d,ncid,varid
    integer,dimension(:) :: start
    real(8),dimension(:,:) :: values

    integer :: ew,i,ierror,ns
    integer,dimension(4) :: mybounds
    integer,dimension(:),allocatable :: displs,sendcounts
    integer,dimension(:,:),allocatable :: bounds
    real(8),dimension(:),allocatable :: sendbuf
    real(8),dimension(:,:),allocatable :: global_values,recvbuf

    ! begin

    if (size(values,1)==local_ewn) then
       ew = global_ewn
       ns = global_nsn
    else if (size(values,1)==local_ewn-1) then
       ew = global_ewn-1
       ns = global_nsn-1
    else
       call parallel_stop(__FILE__,__LINE__)
    end if
    mybounds(1) = ewlb
    mybounds(2) = ewub
    mybounds(3) = nslb
    mybounds(4) = nsub
    if (main_task) then
       allocate(bounds(4,tasks))
    else
       allocate(bounds(1,1))
    end if
    call mpi_gather(mybounds,4,mpi_integer,bounds,4,mpi_integer,&
         main_rank,comm,ierror)
    if (main_task) then
       allocate(global_values(minval(bounds(1,:)):maxval(bounds(2,:)),&
            minval(bounds(3,:)):maxval(bounds(4,:))))
       global_values(:,:) = 0
       distributed_get_var_real8_2d = nf90_get_var(ncid,varid,&
            global_values(1:ew,1:ns),start)
       allocate(displs(tasks+1))
       allocate(sendcounts(tasks))
       sendcounts(:) = (bounds(2,:)-bounds(1,:)+1)*(bounds(4,:)-bounds(3,:)+1)
       displs(1) = 0
       do i = 1,tasks
          displs(i+1) = displs(i)+sendcounts(i)
       end do
       allocate(sendbuf(displs(tasks+1)))
       do i = 1,tasks
          sendbuf(displs(i)+1:displs(i+1)) = reshape(&
               global_values(bounds(1,i):bounds(2,i),bounds(3,i):bounds(4,i)),&
               (/displs(i+1)-displs(i)/))
       end do
    else
       allocate(displs(1))
       allocate(sendcounts(1))
       allocate(sendbuf(1))
    end if
    call broadcast(distributed_get_var_real8_2d)
    allocate(recvbuf(local_ewn,local_nsn))
    call mpi_scatterv(sendbuf,sendcounts,displs,mpi_real8,&
         recvbuf,size(recvbuf),mpi_real8,main_rank,comm,ierror)
    values(:,:) = recvbuf(:size(values,1),:size(values,2))
    !automatic deallocation
  end function

  function distributed_get_var_real8_3d(ncid,varid,values,start)
    use mpi
    implicit none
    integer :: distributed_get_var_real8_3d,ncid,varid
    integer,dimension(:) :: start
    real(8),dimension(:,:,:) :: values

    integer :: ew,i,ierror,ns
    integer,dimension(4) :: mybounds
    integer,dimension(:),allocatable :: displs,sendcounts
    integer,dimension(:,:),allocatable :: bounds
    real(8),dimension(:),allocatable :: sendbuf
    real(8),dimension(:,:,:),allocatable :: global_values,recvbuf

    ! begin

    if (size(values,1)==local_ewn) then
       ew = global_ewn
       ns = global_nsn
    else if (size(values,1)==local_ewn-1) then
       ew = global_ewn-1
       ns = global_nsn-1
    else
       call parallel_stop(__FILE__,__LINE__)
    end if
    mybounds(1) = ewlb
    mybounds(2) = ewub
    mybounds(3) = nslb
    mybounds(4) = nsub
    if (main_task) then
       allocate(bounds(4,tasks))
    else
       allocate(bounds(1,1))
    end if
    call mpi_gather(mybounds,4,mpi_integer,bounds,4,mpi_integer,&
         main_rank,comm,ierror)
    if (main_task) then
       allocate(global_values(minval(bounds(1,:)):maxval(bounds(2,:)),&
            minval(bounds(3,:)):maxval(bounds(4,:)),size(values,3)))
       global_values(:,:,:) = 0
       distributed_get_var_real8_3d = nf90_get_var(ncid,varid,&
            global_values(1:ew,1:ns,:),start)
       allocate(displs(tasks+1))
       allocate(sendcounts(tasks))
       sendcounts(:) = (bounds(2,:)-bounds(1,:)+1)*&
            (bounds(4,:)-bounds(3,:)+1)*size(values,3)
       displs(1) = 0
       do i = 1,tasks
          displs(i+1) = displs(i)+sendcounts(i)
       end do
       allocate(sendbuf(displs(tasks+1)))
       do i = 1,tasks
          sendbuf(displs(i)+1:displs(i+1)) = reshape(global_values(&
               bounds(1,i):bounds(2,i),bounds(3,i):bounds(4,i),:),&
               (/displs(i+1)-displs(i)/))
       end do
    else
       allocate(displs(1))
       allocate(sendcounts(1))
       allocate(sendbuf(1))
    end if
    call broadcast(distributed_get_var_real8_3d)
    allocate(recvbuf(local_ewn,local_nsn,size(values,3)))
    call mpi_scatterv(sendbuf,sendcounts,displs,mpi_real8,&
         recvbuf,size(recvbuf),mpi_real8,main_rank,comm,ierror)
    values(:,:,:) = recvbuf(:size(values,1),:size(values,2),:)
    !automatic deallocation
  end function

  subroutine distributed_grid(ewn,nsn)
    implicit none
    integer :: ewn,nsn

    integer :: best,i,j,metric
    integer :: ewrank,ewtasks,nsrank,nstasks
    real(8) :: rewtasks,rnstasks

    ! begin

    global_ewn = ewn
    global_nsn = nsn

    ewtasks = 0
    nstasks = 0
    best = huge(best)
    do i = 1,min(tasks,global_ewn)
       j = tasks/i
       if (j<=global_nsn.and.i*j==tasks) then ! try to use all tasks
          metric = abs(i*global_nsn-j*global_ewn) ! zero if ewn/nsn == i/j
          if (metric<best) then
             best = metric
             ewtasks = i
             nstasks = j
          end if
       end if
    end do
    if (ewtasks*nstasks/=tasks) call parallel_stop(__FILE__,__LINE__)

    ewrank = mod(this_rank,ewtasks)
    rewtasks = 1/real(ewtasks,8)
    ewlb = nint(ewrank*global_ewn*rewtasks)+1-lhalo
    ewub = nint((ewrank+1)*global_ewn*rewtasks)+uhalo
    local_ewn = ewub-ewlb+1
    own_ewn = local_ewn-lhalo-uhalo
    ewn = local_ewn

    nsrank = this_rank/ewtasks
    rnstasks = 1/real(nstasks,8)
    nslb = nint(nsrank*global_nsn*rnstasks)+1-lhalo
    nsub = nint((nsrank+1)*global_nsn*rnstasks)+uhalo
    local_nsn = nsub-nslb+1
    own_nsn = local_nsn-lhalo-uhalo
    nsn = local_nsn

    east = this_rank-1
    if ((east/ewtasks<this_rank/ewtasks).or.(east<0)) east = east+ewtasks
    west = this_rank+1
    if (west/ewtasks>this_rank/ewtasks) west = west-ewtasks
    north = this_rank-ewtasks
    if (north<0) north = north+tasks
    south = this_rank+ewtasks
    if (south>=tasks) south = south-tasks
  end subroutine

  function distributed_owner(ew,ewn,ns,nsn)
    implicit none
    logical :: distributed_owner
    integer :: ew,ewn,ns,nsn
    ! begin
    distributed_owner = (ew>lhalo.and.ew<=ewn-uhalo.and.&
         ns>lhalo.and.ns<=nsn-uhalo)
  end function

  function distributed_put_var_integer_2d(ncid,varid,values,start)
    use mpi
    implicit none
    integer :: distributed_put_var_integer_2d,ncid,varid
    integer,dimension(:) :: start
    integer,dimension(:,:) :: values

    integer :: ew,i,ierror,ns
    integer,dimension(4) :: mybounds
    integer,dimension(:),allocatable :: displs,recvcounts
    integer,dimension(:,:),allocatable :: bounds
    integer,dimension(:),allocatable :: recvbuf
    integer,dimension(:,:),allocatable :: global_values,sendbuf

    ! begin

    if (size(values,1)==local_ewn) then
       ew = global_ewn
       ns = global_nsn
    else if (size(values,1)==local_ewn-1) then
       ew = global_ewn-1
       ns = global_nsn-1
    else
       call parallel_stop(__FILE__,__LINE__)
    end if
    mybounds(1) = ewlb+lhalo
    mybounds(2) = ewub-uhalo
    mybounds(3) = nslb+lhalo
    mybounds(4) = nsub-uhalo
    if (main_task) then
       allocate(bounds(4,tasks))
    else
       allocate(bounds(1,1))
    end if
    call mpi_gather(mybounds,4,mpi_integer,bounds,4,mpi_integer,&
         main_rank,comm,ierror)
    if (main_task) then
       allocate(global_values(minval(bounds(1,:)):maxval(bounds(2,:)),&
            minval(bounds(3,:)):maxval(bounds(4,:))))
       global_values(:,:) = 0
       allocate(displs(tasks+1))
       allocate(recvcounts(tasks))
       recvcounts(:) = (bounds(2,:)-bounds(1,:)+1)*(bounds(4,:)-bounds(3,:)+1)
       displs(1) = 0
       do i = 1,tasks
          displs(i+1) = displs(i)+recvcounts(i)
       end do
       allocate(recvbuf(displs(tasks+1)))
    else
       allocate(displs(1))
       allocate(recvcounts(1))
       allocate(recvbuf(1))
    end if
    allocate(sendbuf(mybounds(1):mybounds(2),mybounds(3):mybounds(4)))
    sendbuf(:,:) = values(1+lhalo:local_ewn-uhalo,1+lhalo:local_nsn-uhalo)
    call mpi_gatherv(sendbuf,size(sendbuf),mpi_integer,&
         recvbuf,recvcounts,displs,mpi_integer,main_rank,comm,ierror)
    if (main_task) then
       do i = 1,tasks
          global_values(bounds(1,i):bounds(2,i),bounds(3,i):bounds(4,i)) = &
               reshape(recvbuf(displs(i)+1:displs(i+1)), &
               (/bounds(2,i)-bounds(1,i)+1,bounds(4,i)-bounds(3,i)+1/))
       end do
       distributed_put_var_integer_2d = nf90_put_var(ncid,varid,&
            global_values(1:ew,1:ns),start)
    end if
    call broadcast(distributed_put_var_integer_2d)
    !automatic deallocation
  end function

  function distributed_put_var_real4_1d(ncid,varid,values)
    use mpi
    use netcdf
    implicit none
    integer :: distributed_put_var_real4_1d,ncid,varid
    real(4),dimension(:) :: values

    integer :: i,ierror,myn,status,x0id,x1id,y0id,y1id

    integer,dimension(2) :: mybounds
    integer,dimension(:),allocatable :: displs,recvcounts
    integer,dimension(:,:),allocatable :: bounds
    real(4),dimension(:),allocatable :: global_values,recvbuf

    ! begin

    if (main_task) then
       allocate(bounds(2,tasks))
       status = nf90_inq_varid(ncid,"x0",x0id)
       status = nf90_inq_varid(ncid,"x1",x1id)
       status = nf90_inq_varid(ncid,"y0",y0id)
       status = nf90_inq_varid(ncid,"y1",y1id)
    else
       allocate(bounds(1,1))
    end if
    call broadcast(x0id)
    call broadcast(x1id)
    call broadcast(y0id)
    call broadcast(y1id)
    if (varid==x0id) then
       mybounds(1) = ewlb
       mybounds(2) = ewub-1
       myn = global_ewn-1
    else if (varid==x1id) then
       mybounds(1) = ewlb
       mybounds(2) = ewub
       myn = global_ewn
    else if (varid==y0id) then
       mybounds(1) = nslb
       mybounds(2) = nsub-1
       myn = global_nsn-1
    else if (varid==y1id) then
       mybounds(1) = nslb
       mybounds(2) = nsub
       myn = global_nsn
    else
       call parallel_stop(__FILE__,__LINE__)
    end if
    call mpi_gather(mybounds,2,mpi_integer,bounds,2,mpi_integer,&
         main_rank,comm,ierror)
    if (main_task) then
       allocate(global_values(minval(bounds(1,:)):maxval(bounds(2,:))))
       global_values(:) = 0
       allocate(displs(tasks+1))
       allocate(recvcounts(tasks))
       recvcounts(:) = bounds(2,:)-bounds(1,:)+1
       displs(1) = 0
       do i = 1,tasks
          displs(i+1) = displs(i)+recvcounts(i)
       end do
       allocate(recvbuf(displs(tasks+1)))
    else
       allocate(displs(1))
       allocate(recvcounts(1))
       allocate(recvbuf(1))
    end if
    call mpi_gatherv(values,size(values),mpi_real4,recvbuf,recvcounts,&
         displs,mpi_real4,main_rank,comm,ierror)
    if (main_task) then
       do i = 1,tasks
          global_values(bounds(1,i):bounds(2,i)) = &
               recvbuf(displs(i)+1:displs(i+1))
       end do
       distributed_put_var_real4_1d = &
            nf90_put_var(ncid,varid,global_values(1:myn))
    end if
    call broadcast(distributed_put_var_real4_1d)
    !automatic deallocation
  end function

  function distributed_put_var_real4_2d(ncid,varid,values,start)
    use mpi
    implicit none
    integer :: distributed_put_var_real4_2d,ncid,varid
    integer,dimension(:) :: start
    real(4),dimension(:,:) :: values

    integer :: ew,i,ierror,ns
    integer,dimension(4) :: mybounds
    integer,dimension(:),allocatable :: displs,recvcounts
    integer,dimension(:,:),allocatable :: bounds
    real(4),dimension(:),allocatable :: recvbuf
    real(4),dimension(:,:),allocatable :: global_values,sendbuf

    ! begin

    if (size(values,1)==local_ewn) then
       ew = global_ewn
       ns = global_nsn
    else if (size(values,1)==local_ewn-1) then
       ew = global_ewn-1
       ns = global_nsn-1
    else
       call parallel_stop(__FILE__,__LINE__)
    end if
    mybounds(1) = ewlb+lhalo
    mybounds(2) = ewub-uhalo
    mybounds(3) = nslb+lhalo
    mybounds(4) = nsub-uhalo
    if (main_task) then
       allocate(bounds(4,tasks))
    else
       allocate(bounds(1,1))
    end if
    call mpi_gather(mybounds,4,mpi_integer,bounds,4,mpi_integer,&
         main_rank,comm,ierror)
    if (main_task) then
       allocate(global_values(minval(bounds(1,:)):maxval(bounds(2,:)),&
            minval(bounds(3,:)):maxval(bounds(4,:))))
       global_values(:,:) = 0
       allocate(displs(tasks+1))
       allocate(recvcounts(tasks))
       recvcounts(:) = (bounds(2,:)-bounds(1,:)+1)*(bounds(4,:)-bounds(3,:)+1)
       displs(1) = 0
       do i = 1,tasks
          displs(i+1) = displs(i)+recvcounts(i)
       end do
       allocate(recvbuf(displs(tasks+1)))
    else
       allocate(displs(1))
       allocate(recvcounts(1))
       allocate(recvbuf(1))
    end if
    allocate(sendbuf(mybounds(1):mybounds(2),mybounds(3):mybounds(4)))
    sendbuf(:,:) = values(1+lhalo:local_ewn-uhalo,1+lhalo:local_nsn-uhalo)
    call mpi_gatherv(sendbuf,size(sendbuf),mpi_real4,&
         recvbuf,recvcounts,displs,mpi_real4,main_rank,comm,ierror)
    if (main_task) then
       do i = 1,tasks
          global_values(bounds(1,i):bounds(2,i),bounds(3,i):bounds(4,i)) = &
               reshape(recvbuf(displs(i)+1:displs(i+1)), &
               (/bounds(2,i)-bounds(1,i)+1,bounds(4,i)-bounds(3,i)+1/))
       end do
       distributed_put_var_real4_2d = nf90_put_var(ncid,varid,&
            global_values(1:ew,1:ns),start)
    end if
    call broadcast(distributed_put_var_real4_2d)
    !automatic deallocation
  end function

  function distributed_put_var_real8_2d(ncid,varid,values,start)
    use mpi
    implicit none
    integer :: distributed_put_var_real8_2d,ncid,varid
    integer,dimension(:) :: start
    real(8),dimension(:,:) :: values

    integer :: ew,i,ierror,ns
    integer,dimension(4) :: mybounds
    integer,dimension(:),allocatable :: displs,recvcounts
    integer,dimension(:,:),allocatable :: bounds
    real(8),dimension(:),allocatable :: recvbuf
    real(8),dimension(:,:),allocatable :: global_values,sendbuf

    ! begin

    if (size(values,1)==local_ewn) then
       ew = global_ewn
       ns = global_nsn
    else if (size(values,1)==local_ewn-1) then
       ew = global_ewn-1
       ns = global_nsn-1
    else
       call parallel_stop(__FILE__,__LINE__)
    end if
    mybounds(1) = ewlb+lhalo
    mybounds(2) = ewub-uhalo
    mybounds(3) = nslb+lhalo
    mybounds(4) = nsub-uhalo
    if (main_task) then
       allocate(bounds(4,tasks))
    else
       allocate(bounds(1,1))
    end if
    call mpi_gather(mybounds,4,mpi_integer,bounds,4,mpi_integer,&
         main_rank,comm,ierror)
    if (main_task) then
       allocate(global_values(minval(bounds(1,:)):maxval(bounds(2,:)),&
            minval(bounds(3,:)):maxval(bounds(4,:))))
       global_values(:,:) = 0
       allocate(displs(tasks+1))
       allocate(recvcounts(tasks))
       recvcounts(:) = (bounds(2,:)-bounds(1,:)+1)*(bounds(4,:)-bounds(3,:)+1)
       displs(1) = 0
       do i = 1,tasks
          displs(i+1) = displs(i)+recvcounts(i)
       end do
       allocate(recvbuf(displs(tasks+1)))
    else
       allocate(displs(1))
       allocate(recvcounts(1))
       allocate(recvbuf(1))
    end if
    allocate(sendbuf(mybounds(1):mybounds(2),mybounds(3):mybounds(4)))
    sendbuf(:,:) = values(1+lhalo:local_ewn-uhalo,1+lhalo:local_nsn-uhalo)
    call mpi_gatherv(sendbuf,size(sendbuf),mpi_real8,&
         recvbuf,recvcounts,displs,mpi_real8,main_rank,comm,ierror)
    if (main_task) then
       do i = 1,tasks
          global_values(bounds(1,i):bounds(2,i),bounds(3,i):bounds(4,i)) = &
               reshape(recvbuf(displs(i)+1:displs(i+1)), &
               (/bounds(2,i)-bounds(1,i)+1,bounds(4,i)-bounds(3,i)+1/))
       end do
       distributed_put_var_real8_2d = nf90_put_var(ncid,varid,&
            global_values(1:ew,1:ns),start)
    end if
    call broadcast(distributed_put_var_real8_2d)
    !automatic deallocation
  end function

  function distributed_put_var_real8_3d(ncid,varid,values,start)
    use mpi
    implicit none
    integer :: distributed_put_var_real8_3d,ncid,varid
    integer,dimension(:) :: start
    real(8),dimension(:,:,:) :: values

    integer :: ew,i,ierror,ns,nz
    integer,dimension(4) :: mybounds
    integer,dimension(:),allocatable :: displs,recvcounts
    integer,dimension(:,:),allocatable :: bounds
    real(8),dimension(:),allocatable :: recvbuf
    real(8),dimension(:,:,:),allocatable :: global_values,sendbuf

    ! begin

    nz = size(values,3)
    if (size(values,1)==local_ewn) then
       ew = global_ewn
       ns = global_nsn
    else if (size(values,1)==local_ewn-1) then
       ew = global_ewn-1
       ns = global_nsn-1
    else
       call parallel_stop(__FILE__,__LINE__)
    end if
    mybounds(1) = ewlb+lhalo
    mybounds(2) = ewub-uhalo
    mybounds(3) = nslb+lhalo
    mybounds(4) = nsub-uhalo
    if (main_task) then
       allocate(bounds(4,tasks))
    else
       allocate(bounds(1,1))
    end if
    call mpi_gather(mybounds,4,mpi_integer,bounds,4,mpi_integer,&
         main_rank,comm,ierror)
    if (main_task) then
       allocate(global_values(minval(bounds(1,:)):maxval(bounds(2,:)),&
            minval(bounds(3,:)):maxval(bounds(4,:)),nz))
       global_values(:,:,:) = 0
       allocate(displs(tasks+1))
       allocate(recvcounts(tasks))
       recvcounts(:) = (bounds(2,:)-bounds(1,:)+1)*(bounds(4,:)-bounds(3,:)+1)&
            *nz
       displs(1) = 0
       do i = 1,tasks
          displs(i+1) = displs(i)+recvcounts(i)
       end do
       allocate(recvbuf(displs(tasks+1)))
    else
       allocate(displs(1))
       allocate(recvcounts(1))
       allocate(recvbuf(1))
    end if
    allocate(sendbuf(mybounds(1):mybounds(2),mybounds(3):mybounds(4),nz))
    sendbuf(:,:,:) = values(1+lhalo:local_ewn-uhalo,1+lhalo:local_nsn-uhalo,:)
    call mpi_gatherv(sendbuf,size(sendbuf),mpi_real8,&
         recvbuf,recvcounts,displs,mpi_real8,main_rank,comm,ierror)
    if (main_task) then
       do i = 1,tasks
          global_values(bounds(1,i):bounds(2,i),bounds(3,i):bounds(4,i),:) = &
               reshape(recvbuf(displs(i)+1:displs(i+1)), &
               (/bounds(2,i)-bounds(1,i)+1,bounds(4,i)-bounds(3,i)+1,nz/))
       end do
       distributed_put_var_real8_3d = nf90_put_var(ncid,varid,&
            global_values(1:ew,1:ns,:),start)
    end if
    call broadcast(distributed_put_var_real8_3d)
    !automatic deallocation
  end function

  subroutine global_sum(x,y)
    use mpi
    implicit none
    real(8) :: x,y
    
    integer :: ierror
    real(8),dimension(2) :: recvbuf,sendbuf
    ! begin
    sendbuf(1) = x
    sendbuf(2) = y
    call mpi_allreduce(sendbuf,recvbuf,2,mpi_real8,mpi_sum,comm,ierror)
    x = recvbuf(1)
    y = recvbuf(2)
  end subroutine

  subroutine not_parallel(file,line)
    implicit none
    integer :: line
    character(len=*) :: file
    ! begin
    call parallel_stop(file,line)
  end subroutine

  subroutine parallel_barrier
    use mpi
    implicit none
    integer :: ierror
    ! begin
    call mpi_barrier(comm,ierror)
  end subroutine

  function parallel_boundary(ew,ewn,ns,nsn)
    implicit none
    logical :: parallel_boundary
    integer :: ew,ewn,ns,nsn
    ! begin
    parallel_boundary = (ewlb<1.and.ew==1+lhalo).or.&
         (ewub>global_ewn.and.ew==ewn-uhalo).or.&
         (nslb<1.and.ns==1+lhalo).or.&
         (nsub>global_nsn.and.ns==nsn-uhalo)
  end function

  function parallel_close(ncid)
    implicit none
    integer :: ncid,parallel_close
    ! begin
    if (main_task) parallel_close = nf90_close(ncid)
    call broadcast(parallel_close)
  end function

  function parallel_create(path,cmode,ncid)
    implicit none
    integer :: cmode,ncid,parallel_create
    character(len=*) :: path
    ! begin
    if (main_task) parallel_create = nf90_create(path,cmode,ncid)
    call broadcast(parallel_create)
    call broadcast(ncid)
  end function
    
  function parallel_def_dim(ncid,name,len,dimid)
    use netcdf
    implicit none
    integer :: dimid,len,ncid,parallel_def_dim
    character(len=*) :: name
    ! begin
    if (main_task) parallel_def_dim = nf90_def_dim(ncid,name,len,dimid)
    call broadcast(parallel_def_dim)
    call broadcast(dimid)
  end function

  function parallel_def_var_dimids(ncid,name,xtype,dimids,varid)
    implicit none
    integer :: ncid,parallel_def_var_dimids,varid,xtype
    integer,dimension(:) :: dimids
    character(len=*) :: name
    ! begin
    if (main_task) parallel_def_var_dimids = &
         nf90_def_var(ncid,name,xtype,dimids,varid)
    call broadcast(parallel_def_var_dimids)
    call broadcast(varid)
  end function

  function parallel_def_var_nodimids(ncid,name,xtype,varid)
    implicit none
    integer :: ncid,parallel_def_var_nodimids,varid,xtype
    character(len=*) :: name
    ! begin
    if (main_task) parallel_def_var_nodimids = &
         nf90_def_var(ncid,name,xtype,varid)
    call broadcast(parallel_def_var_nodimids)
    call broadcast(varid)
  end function

  function parallel_enddef(ncid)
    implicit none
    integer :: ncid,parallel_enddef
    ! begin
    if (main_task) parallel_enddef = nf90_enddef(ncid)
    call broadcast(parallel_enddef)
  end function

  subroutine parallel_finalise
    use mpi
    implicit none
    integer :: ierror
    ! begin
    call mpi_finalize(ierror)
  end subroutine

  function parallel_get_att_character(ncid,varid,name,values)
    implicit none
    integer :: ncid,parallel_get_att_character,varid
    character(len=*) :: name,values
    ! begin
    if (main_task) parallel_get_att_character = &
         nf90_get_att(ncid,varid,name,values)
    call broadcast(parallel_get_att_character)
    call broadcast(values)
  end function

  function parallel_get_att_real4(ncid,varid,name,values)
    implicit none
    integer :: ncid,parallel_get_att_real4,varid
    character(len=*) :: name
    real(4) :: values
    ! begin
    if (main_task) parallel_get_att_real4 = &
         nf90_get_att(ncid,varid,name,values)
    call broadcast(parallel_get_att_real4)
    call broadcast(values)
  end function

  function parallel_get_att_real4_1d(ncid,varid,name,values)
    implicit none
    integer :: ncid,parallel_get_att_real4_1d,varid
    character(len=*) :: name
    real(4),dimension(:) :: values
    ! begin
    if (main_task) parallel_get_att_real4_1d = &
         nf90_get_att(ncid,varid,name,values)
    call broadcast(parallel_get_att_real4_1d)
    call broadcast(values)
  end function

  function parallel_get_att_real8(ncid,varid,name,values)
    implicit none
    integer :: ncid,parallel_get_att_real8,varid
    character(len=*) :: name
    real(8) :: values
    ! begin
    if (main_task) parallel_get_att_real8 = &
         nf90_get_att(ncid,varid,name,values)
    call broadcast(parallel_get_att_real8)
    call broadcast(values)
  end function

  function parallel_get_att_real8_1d(ncid,varid,name,values)
    implicit none
    integer :: ncid,parallel_get_att_real8_1d,varid
    character(len=*) :: name
    real(8),dimension(:) :: values
    ! begin
    if (main_task) parallel_get_att_real8_1d = &
         nf90_get_att(ncid,varid,name,values)
    call broadcast(parallel_get_att_real8_1d)
    call broadcast(values)
  end function

  function parallel_get_var_integer_1d(ncid,varid,values)
    implicit none
    integer :: ncid,parallel_get_var_integer_1d,varid
    integer,dimension(:) :: values
    ! begin
    if (main_task) parallel_get_var_integer_1d = &
         nf90_get_var(ncid,varid,values)
    call broadcast(parallel_get_var_integer_1d)
    call broadcast(values)
  end function

  function parallel_get_var_real4_1d(ncid,varid,values)
    implicit none
    integer :: ncid,parallel_get_var_real4_1d,varid
    real(4),dimension(:) :: values
    ! begin
    if (main_task) parallel_get_var_real4_1d = &
         nf90_get_var(ncid,varid,values)
    call broadcast(parallel_get_var_real4_1d)
    call broadcast(values)
  end function

  subroutine parallel_halo_integer_2d(a)
    use mpi
    implicit none
    integer,dimension(:,:) :: a
    
    integer :: erequest,ierror,nrequest,srequest,wrequest
    integer,dimension(lhalo,local_nsn-lhalo-uhalo) :: erecv,wsend
    integer,dimension(uhalo,local_nsn-lhalo-uhalo) :: esend,wrecv
    integer,dimension(local_ewn,lhalo) :: nrecv,ssend
    integer,dimension(local_ewn,uhalo) :: nsend,srecv

    ! begin

    ! staggered grid
    if (size(a,1)==local_ewn-1.and.size(a,2)==local_nsn-1) return

    ! unknown grid
    if (size(a,1)/=local_ewn.or.size(a,2)/=local_nsn) &
         call parallel_stop(__FILE__,__LINE__)

    ! unstaggered grid
    call mpi_irecv(wrecv,size(wrecv),mpi_integer,west,west,&
         comm,wrequest,ierror)
    call mpi_irecv(erecv,size(erecv),mpi_integer,east,east,&
         comm,erequest,ierror)
    call mpi_irecv(srecv,size(srecv),mpi_integer,south,south,&
         comm,srequest,ierror)
    call mpi_irecv(nrecv,size(nrecv),mpi_integer,north,north,&
         comm,nrequest,ierror)

    esend(:,:) = a(1+lhalo:1+lhalo+uhalo-1,1+lhalo:local_nsn-uhalo)
    call mpi_send(esend,size(esend),mpi_integer,east,this_rank,comm,ierror)
    wsend(:,:) = &
         a(local_ewn-lhalo-uhalo+1:local_ewn-uhalo,1+lhalo:local_nsn-uhalo)
    call mpi_send(wsend,size(wsend),mpi_integer,west,this_rank,comm,ierror)

    call mpi_wait(wrequest,mpi_status_ignore,ierror)
    a(local_ewn-uhalo+1:,1+lhalo:local_nsn-uhalo) = wrecv(:,:)
    call mpi_wait(erequest,mpi_status_ignore,ierror)
    a(:lhalo,1+lhalo:local_nsn-uhalo) = erecv(:,:)

    nsend(:,:) = a(:,1+lhalo:1+lhalo+uhalo-1)
    call mpi_send(nsend,size(nsend),mpi_integer,north,this_rank,comm,ierror)
    ssend(:,:) = a(:,local_nsn-lhalo-uhalo+1:local_nsn-uhalo)
    call mpi_send(ssend,size(ssend),mpi_integer,south,this_rank,comm,ierror)

    call mpi_wait(srequest,mpi_status_ignore,ierror)
    a(:,local_nsn-uhalo+1:) = srecv(:,:)
    call mpi_wait(nrequest,mpi_status_ignore,ierror)
    a(:,:lhalo) = nrecv(:,:)
  end subroutine

  subroutine parallel_halo_real8_2d(a)
    use mpi
    implicit none
    real(8),dimension(:,:) :: a
    
    integer :: erequest,ierror,nrequest,srequest,wrequest
    real(8),dimension(lhalo,local_nsn-lhalo-uhalo) :: erecv,wsend
    real(8),dimension(uhalo,local_nsn-lhalo-uhalo) :: esend,wrecv
    real(8),dimension(local_ewn,lhalo) :: nrecv,ssend
    real(8),dimension(local_ewn,uhalo) :: nsend,srecv

    ! begin

    ! staggered grid
    if (size(a,1)==local_ewn-1.and.size(a,2)==local_nsn-1) return

    ! unknown grid
    if (size(a,1)/=local_ewn.or.size(a,2)/=local_nsn) &
         call parallel_stop(__FILE__,__LINE__)

    ! unstaggered grid
    call mpi_irecv(wrecv,size(wrecv),mpi_real8,west,west,&
         comm,wrequest,ierror)
    call mpi_irecv(erecv,size(erecv),mpi_real8,east,east,&
         comm,erequest,ierror)
    call mpi_irecv(srecv,size(srecv),mpi_real8,south,south,&
         comm,srequest,ierror)
    call mpi_irecv(nrecv,size(nrecv),mpi_real8,north,north,&
         comm,nrequest,ierror)

    esend(:,:) = a(1+lhalo:1+lhalo+uhalo-1,1+lhalo:local_nsn-uhalo)
    call mpi_send(esend,size(esend),mpi_real8,east,this_rank,comm,ierror)
    wsend(:,:) = &
         a(local_ewn-lhalo-uhalo+1:local_ewn-uhalo,1+lhalo:local_nsn-uhalo)
    call mpi_send(wsend,size(wsend),mpi_real8,west,this_rank,comm,ierror)

    call mpi_wait(wrequest,mpi_status_ignore,ierror)
    a(local_ewn-uhalo+1:,1+lhalo:local_nsn-uhalo) = wrecv(:,:)
    call mpi_wait(erequest,mpi_status_ignore,ierror)
    a(:lhalo,1+lhalo:local_nsn-uhalo) = erecv(:,:)

    nsend(:,:) = a(:,1+lhalo:1+lhalo+uhalo-1)
    call mpi_send(nsend,size(nsend),mpi_real8,north,this_rank,comm,ierror)
    ssend(:,:) = a(:,local_nsn-lhalo-uhalo+1:local_nsn-uhalo)
    call mpi_send(ssend,size(ssend),mpi_real8,south,this_rank,comm,ierror)

    call mpi_wait(srequest,mpi_status_ignore,ierror)
    a(:,local_nsn-uhalo+1:) = srecv(:,:)
    call mpi_wait(nrequest,mpi_status_ignore,ierror)
    a(:,:lhalo) = nrecv(:,:)
  end subroutine

  subroutine parallel_halo_real8_3d(a)
    use mpi
    implicit none
    real(8),dimension(:,:,:) :: a
    
    integer :: erequest,ierror,one,nrequest,srequest,wrequest
    real(8),dimension(size(a,1),lhalo,local_nsn-lhalo-uhalo) :: erecv,wsend
    real(8),dimension(size(a,1),uhalo,local_nsn-lhalo-uhalo) :: esend,wrecv
    real(8),dimension(size(a,1),local_ewn,lhalo) :: nrecv,ssend
    real(8),dimension(size(a,1),local_ewn,uhalo) :: nsend,srecv

    ! begin

    ! staggered grid
    if (size(a,2)==local_ewn-1.and.size(a,3)==local_nsn-1) return

    ! unknown grid
    if (size(a,2)/=local_ewn.or.size(a,3)/=local_nsn) &
         call parallel_stop(__FILE__,__LINE__)

    ! unstaggered grid
    call mpi_irecv(wrecv,size(wrecv),mpi_real8,west,west,&
         comm,wrequest,ierror)
    call mpi_irecv(erecv,size(erecv),mpi_real8,east,east,&
         comm,erequest,ierror)
    call mpi_irecv(srecv,size(srecv),mpi_real8,south,south,&
         comm,srequest,ierror)
    call mpi_irecv(nrecv,size(nrecv),mpi_real8,north,north,&
         comm,nrequest,ierror)

    esend(:,:,:) = a(:,1+lhalo:1+lhalo+uhalo-1,1+lhalo:local_nsn-uhalo)
    call mpi_send(esend,size(esend),mpi_real8,east,this_rank,comm,ierror)
    wsend(:,:,:) = &
         a(:,local_ewn-lhalo-uhalo+1:local_ewn-uhalo,1+lhalo:local_nsn-uhalo)
    call mpi_send(wsend,size(wsend),mpi_real8,west,this_rank,comm,ierror)

    call mpi_wait(wrequest,mpi_status_ignore,ierror)
    a(:,local_ewn-uhalo+1:,1+lhalo:local_nsn-uhalo) = wrecv(:,:,:)
    call mpi_wait(erequest,mpi_status_ignore,ierror)
    a(:,:lhalo,1+lhalo:local_nsn-uhalo) = erecv(:,:,:)

    nsend(:,:,:) = a(:,:,1+lhalo:1+lhalo+uhalo-1)
    call mpi_send(nsend,size(nsend),mpi_real8,north,this_rank,comm,ierror)
    ssend(:,:,:) = a(:,:,local_nsn-lhalo-uhalo+1:local_nsn-uhalo)
    call mpi_send(ssend,size(ssend),mpi_real8,south,this_rank,comm,ierror)

    call mpi_wait(srequest,mpi_status_ignore,ierror)
    a(:,:,local_nsn-uhalo+1:) = srecv(:,:,:)
    call mpi_wait(nrequest,mpi_status_ignore,ierror)
    a(:,:,:lhalo) = nrecv(:,:,:)
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

  function parallel_inq_attname(ncid,varid,attnum,name)
    implicit none
    integer :: attnum,ncid,parallel_inq_attname,varid
    character(len=*) :: name
    ! begin
    if (main_task) parallel_inq_attname = &
         nf90_inq_attname(ncid,varid,attnum,name)
    call broadcast(parallel_inq_attname)
    call broadcast(name)
  end function

  function parallel_inq_dimid(ncid,name,dimid)
    implicit none
    integer :: dimid,ncid,parallel_inq_dimid
    character(len=*) :: name
    ! begin
    if (main_task) parallel_inq_dimid = nf90_inq_dimid(ncid,name,dimid)
    call broadcast(parallel_inq_dimid)
    call broadcast(dimid)
  end function

  function parallel_inq_varid(ncid,name,varid)
    implicit none
    integer :: ncid,parallel_inq_varid,varid
    character(len=*) :: name
    ! begin
    if (main_task) parallel_inq_varid = nf90_inq_varid(ncid,name,varid)
    call broadcast(parallel_inq_varid)
    call broadcast(varid)
  end function

  function parallel_inquire(ncid,nvariables)
    implicit none
    integer :: ncid,parallel_inquire,nvariables
    ! begin
    if (main_task) parallel_inquire = nf90_inquire(ncid,nvariables=nvariables)
    call broadcast(parallel_inquire)
    call broadcast(nvariables)
  end function

  function parallel_inquire_dimension(ncid,dimid,name,len)
    implicit none
    integer :: dimid,ncid,parallel_inquire_dimension
    integer,optional :: len
    character(len=*),optional :: name
    
    integer :: l
    
    ! begin

    if (present(name)) then
       if (main_task) parallel_inquire_dimension = &
            nf90_inquire_dimension(ncid,dimid,name,len=l)
       call broadcast(name)
    else
       if (main_task) parallel_inquire_dimension = &
            nf90_inquire_dimension(ncid,dimid,len=l)
    end if
    call broadcast(parallel_inquire_dimension)
    if (present(len)) then
       call broadcast(l)
       len = l
    end if
  end function

  function parallel_inquire_variable(ncid,varid,name,ndims,dimids,natts)
    implicit none
    integer :: ncid,parallel_inquire_variable,varid
    integer,optional :: ndims,natts
    character(len=*),optional :: name
    integer,dimension(:),optional :: dimids

    integer :: nd,na
    ! begin
    if (present(name)) then
       if (main_task) parallel_inquire_variable = &
            nf90_inquire_variable(ncid,varid,name=name)
       call broadcast(parallel_inquire_variable)
       call broadcast(name)
       if (parallel_inquire_variable/=nf90_noerr) return
    end if
    if (present(dimids)) then
       if (main_task) parallel_inquire_variable = &
            nf90_inquire_variable(ncid,varid,dimids=dimids)
       call broadcast(parallel_inquire_variable)
       call broadcast(dimids)
       if (parallel_inquire_variable/=nf90_noerr) return
    end if
    if (main_task) parallel_inquire_variable = &
         nf90_inquire_variable(ncid,varid,ndims=nd,natts=na)
    call broadcast(parallel_inquire_variable)
    if (present(ndims)) then
       call broadcast(nd)
       ndims = nd
    end if
    if (present(natts)) then
       call broadcast(na)
       natts = na
    end if
  end function

  function parallel_open(path,mode,ncid)
    implicit none
    integer :: mode,ncid,parallel_open
    character(len=*) :: path
    ! begin
    if (main_task) parallel_open = nf90_open(path,mode,ncid)
    call broadcast(parallel_open)
  end function

  subroutine parallel_print_integer_2d(name,values)
    use mpi
    implicit none
    character(*) :: name
    integer,dimension(:,:) :: values
    
    integer,parameter :: u = 33
    character(3) :: ts
    integer :: i,ierror,j,k
    integer,dimension(4) :: mybounds
    integer,dimension(:),allocatable :: displs,recvcounts
    integer,dimension(:,:),allocatable :: bounds
    integer,dimension(:),allocatable :: recvbuf
    integer,dimension(:,:),allocatable :: global_values,sendbuf

    ! begin
    mybounds(1) = ewlb+lhalo
    mybounds(2) = ewub-uhalo
    mybounds(3) = nslb+lhalo
    mybounds(4) = nsub-uhalo
    if (main_task) then
       allocate(bounds(4,tasks))
    else
       allocate(bounds(1,1))
    end if
    call mpi_gather(mybounds,4,mpi_integer,bounds,4,mpi_integer,&
         main_rank,comm,ierror)
    if (main_task) then
       allocate(global_values(minval(bounds(1,:)):maxval(bounds(2,:)),&
            minval(bounds(3,:)):maxval(bounds(4,:))))
       global_values(:,:) = 0
       allocate(displs(tasks+1))
       allocate(recvcounts(tasks))
       recvcounts(:) = (bounds(2,:)-bounds(1,:)+1)*(bounds(4,:)-bounds(3,:)+1)
       displs(1) = 0
       do i = 1,tasks
          displs(i+1) = displs(i)+recvcounts(i)
       end do
       allocate(recvbuf(displs(tasks+1)))
    else
       allocate(displs(1))
       allocate(recvcounts(1))
       allocate(recvbuf(1))
    end if
    allocate(sendbuf(mybounds(1):mybounds(2),mybounds(3):mybounds(4)))
    sendbuf(:,:) = values(1+lhalo:local_ewn-uhalo,1+lhalo:local_nsn-uhalo)
    call mpi_gatherv(sendbuf,size(sendbuf),mpi_integer,&
         recvbuf,recvcounts,displs,mpi_integer,main_rank,comm,ierror)
    if (main_task) then
       do i = 1,tasks
          global_values(bounds(1,i):bounds(2,i),bounds(3,i):bounds(4,i)) = &
               reshape(recvbuf(displs(i)+1:displs(i+1)), &
               (/bounds(2,i)-bounds(1,i)+1,bounds(4,i)-bounds(3,i)+1/))
       end do
       write(ts,'(i3.3)') tasks
       open(unit=u,file=name//ts//".txt",form="formatted",status="replace")
       if (size(values,1)<local_ewn) then
          do j = lbound(global_values,2),ubound(global_values,2)-1
             do i = lbound(global_values,1),ubound(global_values,1)-1
                write(u,*) j,i,global_values(i,j)
             end do
             write(u,'()')
          end do
       else
          do j = lbound(global_values,2),ubound(global_values,2)
             do i = lbound(global_values,1),ubound(global_values,1)
                write(u,*) j,i,global_values(i,j)
             end do
             write(u,'()')
          end do
       end if
       close(u)
    end if
    ! automatic deallocation
  end subroutine

  subroutine parallel_print_real8_2d(name,values)
    use mpi
    implicit none
    character(*) :: name
    real(8),dimension(:,:) :: values
    
    integer,parameter :: u = 33
    character(3) :: ts
    integer :: i,ierror,j,k
    integer,dimension(4) :: mybounds
    integer,dimension(:),allocatable :: displs,recvcounts
    integer,dimension(:,:),allocatable :: bounds
    real(8),dimension(:),allocatable :: recvbuf
    real(8),dimension(:,:),allocatable :: global_values,sendbuf

    ! begin
    mybounds(1) = ewlb+lhalo
    mybounds(2) = ewub-uhalo
    mybounds(3) = nslb+lhalo
    mybounds(4) = nsub-uhalo
    if (main_task) then
       allocate(bounds(4,tasks))
    else
       allocate(bounds(1,1))
    end if
    call mpi_gather(mybounds,4,mpi_integer,bounds,4,mpi_integer,&
         main_rank,comm,ierror)
    if (main_task) then
       allocate(global_values(minval(bounds(1,:)):maxval(bounds(2,:)),&
            minval(bounds(3,:)):maxval(bounds(4,:))))
       global_values(:,:) = 0
       allocate(displs(tasks+1))
       allocate(recvcounts(tasks))
       recvcounts(:) = (bounds(2,:)-bounds(1,:)+1)*(bounds(4,:)-bounds(3,:)+1)
       displs(1) = 0
       do i = 1,tasks
          displs(i+1) = displs(i)+recvcounts(i)
       end do
       allocate(recvbuf(displs(tasks+1)))
    else
       allocate(displs(1))
       allocate(recvcounts(1))
       allocate(recvbuf(1))
    end if
    allocate(sendbuf(mybounds(1):mybounds(2),mybounds(3):mybounds(4)))
    sendbuf(:,:) = values(1+lhalo:local_ewn-uhalo,1+lhalo:local_nsn-uhalo)
    call mpi_gatherv(sendbuf,size(sendbuf),mpi_real8,&
         recvbuf,recvcounts,displs,mpi_real8,main_rank,comm,ierror)
    if (main_task) then
       do i = 1,tasks
          global_values(bounds(1,i):bounds(2,i),bounds(3,i):bounds(4,i)) = &
               reshape(recvbuf(displs(i)+1:displs(i+1)), &
               (/bounds(2,i)-bounds(1,i)+1,bounds(4,i)-bounds(3,i)+1/))
       end do
       write(ts,'(i3.3)') tasks
       open(unit=u,file=name//ts//".txt",form="formatted",status="replace")
       if (size(values,1)<local_ewn) then
          do j = lbound(global_values,2),ubound(global_values,2)-1
             do i = lbound(global_values,1),ubound(global_values,1)-1
                write(u,*) j,i,global_values(i,j)
             end do
             write(u,'()')
          end do
       else
          do j = lbound(global_values,2),ubound(global_values,2)
             do i = lbound(global_values,1),ubound(global_values,1)
                write(u,*) j,i,global_values(i,j)
             end do
             write(u,'()')
          end do
       end if
       close(u)
    end if
    ! automatic deallocation
  end subroutine

  subroutine parallel_print_real8_3d(name,values)
    use mpi
    implicit none
    character(*) :: name
    real(8),dimension(:,:,:) :: values
    
    integer,parameter :: u = 33
    character(3) :: ts
    integer :: i,ierror,j,k
    integer,dimension(4) :: mybounds
    integer,dimension(:),allocatable :: displs,recvcounts
    integer,dimension(:,:),allocatable :: bounds
    real(8),dimension(:),allocatable :: recvbuf
    real(8),dimension(:,:,:),allocatable :: global_values,sendbuf

    ! begin
    mybounds(1) = ewlb+lhalo
    mybounds(2) = ewub-uhalo
    mybounds(3) = nslb+lhalo
    mybounds(4) = nsub-uhalo
    if (main_task) then
       allocate(bounds(4,tasks))
    else
       allocate(bounds(1,1))
    end if
    call mpi_gather(mybounds,4,mpi_integer,bounds,4,mpi_integer,&
         main_rank,comm,ierror)
    if (main_task) then
       allocate(global_values(size(values,1),minval(bounds(1,:)):maxval(bounds(2,:)),&
            minval(bounds(3,:)):maxval(bounds(4,:))))
       global_values(:,:,:) = 0
       allocate(displs(tasks+1))
       allocate(recvcounts(tasks))
       recvcounts(:) = (bounds(2,:)-bounds(1,:)+1)*(bounds(4,:)-bounds(3,:)+1)*size(values,1)
       displs(1) = 0
       do i = 1,tasks
          displs(i+1) = displs(i)+recvcounts(i)
       end do
       allocate(recvbuf(displs(tasks+1)))
    else
       allocate(displs(1))
       allocate(recvcounts(1))
       allocate(recvbuf(1))
    end if
    allocate(sendbuf(size(values,1),mybounds(1):mybounds(2),mybounds(3):mybounds(4)))
    sendbuf(:,:,:) = values(:,1+lhalo:local_ewn-uhalo,1+lhalo:local_nsn-uhalo)
    call mpi_gatherv(sendbuf,size(sendbuf),mpi_real8,&
         recvbuf,recvcounts,displs,mpi_real8,main_rank,comm,ierror)
    if (main_task) then
       do i = 1,tasks
          global_values(:,bounds(1,i):bounds(2,i),bounds(3,i):bounds(4,i)) = &
               reshape(recvbuf(displs(i)+1:displs(i+1)), &
               (/size(values,1),bounds(2,i)-bounds(1,i)+1,bounds(4,i)-bounds(3,i)+1/))
       end do
       write(ts,'(i3.3)') tasks
       open(unit=u,file=name//ts//".txt",form="formatted",status="replace")
       if (size(values,2)<local_ewn) then
          do j = lbound(global_values,3),ubound(global_values,3)-1
             do i = lbound(global_values,2),ubound(global_values,2)-1
                write(u,'(2i6,100g15.5e3)') j,i,global_values(:,i,j)
             end do
             write(u,'()')
          end do
       else
          do j = lbound(global_values,3),ubound(global_values,3)
             do i = lbound(global_values,2),ubound(global_values,2)
                write(u,'(2i6,100g15.5e3)') j,i,global_values(:,i,j)
             end do
             write(u,'()')
          end do
       end if
       close(u)
    end if
    ! automatic deallocation
  end subroutine

  function parallel_put_att_character(ncid,varid,name,values)
    implicit none
    integer :: ncid,parallel_put_att_character,varid
    character(len=*) :: name,values
    ! begin
    if (main_task) parallel_put_att_character = nf90_put_att(ncid,varid,name,values)
    call broadcast(parallel_put_att_character)
  end function

  function parallel_put_att_real4(ncid,varid,name,values)
    implicit none
    integer :: ncid,parallel_put_att_real4,varid
    character(len=*) :: name
    real(4) :: values
    ! begin
    if (main_task) parallel_put_att_real4 = nf90_put_att(ncid,varid,name,values)
    call broadcast(parallel_put_att_real4)
  end function

  function parallel_put_att_real4_1d(ncid,varid,name,values)
    implicit none
    integer :: ncid,parallel_put_att_real4_1d,varid
    character(len=*) :: name
    real(4),dimension(:) :: values
    ! begin
    if (main_task) parallel_put_att_real4_1d = nf90_put_att(ncid,varid,name,values)
    call broadcast(parallel_put_att_real4_1d)
  end function

  function parallel_put_att_real8(ncid,varid,name,values)
    implicit none
    integer :: ncid,parallel_put_att_real8,varid
    character(len=*) :: name
    real(8) :: values
    ! begin
    if (main_task) parallel_put_att_real8 = nf90_put_att(ncid,varid,name,values)
    call broadcast(parallel_put_att_real8)
  end function

  function parallel_put_att_real8_1d(ncid,varid,name,values)
    implicit none
    integer :: ncid,parallel_put_att_real8_1d,varid
    character(len=*) :: name
    real(8),dimension(:) :: values
    ! begin
    if (main_task) parallel_put_att_real8_1d = nf90_put_att(ncid,varid,name,values)
    call broadcast(parallel_put_att_real8_1d)
  end function

  function parallel_put_var_real4(ncid,varid,values,start)
    implicit none
    integer :: ncid,parallel_put_var_real4,varid
    integer,dimension(:) :: start
    real(4) :: values
    ! begin
    if (main_task) parallel_put_var_real4 = &
         nf90_put_var(ncid,varid,values,start)
    call broadcast(parallel_put_var_real4)
  end function

  function parallel_put_var_real8(ncid,varid,values,start)
    implicit none
    integer :: ncid,parallel_put_var_real8,varid
    integer,dimension(:) :: start
    real(8) :: values
    ! begin
    if (main_task) parallel_put_var_real8 = &
         nf90_put_var(ncid,varid,values,start)
    call broadcast(parallel_put_var_real8)
  end function

  function parallel_put_var_real8_1d(ncid,varid,values,start)
    implicit none
    integer :: ncid,parallel_put_var_real8_1d,varid
    integer,dimension(:),optional :: start
    real(8),dimension(:) :: values
    ! begin
    if (main_task) then
       if (present(start)) then
          parallel_put_var_real8_1d = nf90_put_var(ncid,varid,values,start)
       else
          parallel_put_var_real8_1d = nf90_put_var(ncid,varid,values)
       end if
    end if
    call broadcast(parallel_put_var_real8_1d)
  end function

  function parallel_redef(ncid)
    implicit none
    integer :: ncid,parallel_redef
    ! begin
    if (main_task) parallel_redef = nf90_redef(ncid)
    call broadcast(parallel_redef)
  end function

  subroutine parallel_show_minmax(label,values)
    use mpi
    implicit none
    character(*) :: label
    real(8),dimension(:,:,:) :: values
    
    integer :: ierror
    real(8) :: allmin,allmax,mymin,mymax
    ! begin
    mymin = minval(values(:,lhalo:size(values,2)-uhalo,&
         lhalo:size(values,3)-uhalo))
    mymax = maxval(values(:,lhalo:size(values,2)-uhalo,&
         lhalo:size(values,3)-uhalo))
    call mpi_reduce(mymin,allmin,1,mpi_real8,mpi_min,main_rank,comm,ierror)
    call mpi_reduce(mymax,allmax,1,mpi_real8,mpi_max,main_rank,comm,ierror)
    if (main_task) print *,label,allmin,allmax
  end subroutine

  subroutine parallel_stop(file,line)
    use mpi
    implicit none
    integer :: line
    character(len=*) :: file
    integer :: ierror
    ! begin
    if (main_task) write(0,*) "PARALLEL STOP in ",file," at line ",line
    call mpi_finalize(ierror)
    stop "PARALLEL STOP"
  end subroutine

  function parallel_sync(ncid)
    implicit none
    integer :: ncid,parallel_sync
    ! begin
    if (main_task) parallel_sync = nf90_sync(ncid)
    call broadcast(parallel_sync)
  end function

  subroutine parallel_temp_halo(a)
    use mpi
    implicit none
    real(8),dimension(:,:,:) :: a !(:,local_ewn+2,local_nsn+2)
    
    integer :: erequest,ierror,one,nrequest,srequest,wrequest
    real(8),dimension(size(a,1),lhalo+1,local_nsn-lhalo-uhalo) :: erecv,wsend
    real(8),dimension(size(a,1),uhalo+1,local_nsn-lhalo-uhalo) :: esend,wrecv
    real(8),dimension(size(a,1),local_ewn+2,lhalo+1) :: nrecv,ssend
    real(8),dimension(size(a,1),local_ewn+2,uhalo+1) :: nsend,srecv
    ! begin
    call mpi_irecv(wrecv,size(wrecv),mpi_real8,west,west,&
         comm,wrequest,ierror)
    call mpi_irecv(erecv,size(erecv),mpi_real8,east,east,&
         comm,erequest,ierror)
    call mpi_irecv(srecv,size(srecv),mpi_real8,south,south,&
         comm,srequest,ierror)
    call mpi_irecv(nrecv,size(nrecv),mpi_real8,north,north,&
         comm,nrequest,ierror)

    esend(:,:,:) = a(:,2+lhalo:2+lhalo+uhalo,2+lhalo:local_nsn-uhalo+1)
    call mpi_send(esend,size(esend),mpi_real8,east,this_rank,comm,ierror)
    wsend(:,:,:) = a(:,local_ewn-lhalo-uhalo+1:local_ewn-uhalo+2,&
         2+lhalo:local_nsn-uhalo+1)
    call mpi_send(wsend,size(wsend),mpi_real8,west,this_rank,comm,ierror)

    call mpi_wait(wrequest,mpi_status_ignore,ierror)
    a(:,local_ewn-uhalo+2:,2+lhalo:local_nsn-uhalo+1) = wrecv(:,:,:)
    call mpi_wait(erequest,mpi_status_ignore,ierror)
    a(:,:lhalo+1,2+lhalo:local_nsn-uhalo+1) = erecv(:,:,:)

    nsend(:,:,:) = a(:,:,1+lhalo+1:1+lhalo+uhalo+1)
    call mpi_send(nsend,size(nsend),mpi_real8,north,this_rank,comm,ierror)
    ssend(:,:,:) = a(:,:,local_nsn-lhalo-uhalo+1:local_nsn-uhalo+2)
    call mpi_send(ssend,size(ssend),mpi_real8,south,this_rank,comm,ierror)

    call mpi_wait(srequest,mpi_status_ignore,ierror)
    a(:,:,local_nsn-uhalo+2:) = srecv(:,:,:)
    call mpi_wait(nrequest,mpi_status_ignore,ierror)
    a(:,:,:lhalo+1) = nrecv(:,:,:)
  end subroutine

  subroutine parallel_velo_halo_real8_2d(a)
    use mpi
    implicit none
    real(8),dimension(local_ewn-1,local_nsn-1) :: a

    integer :: erequest,ierror,nrequest,srequest,wrequest
    integer,dimension(lhalo,local_nsn-lhalo-uhalo) :: erecv,wsend
    integer,dimension(uhalo-1,local_nsn-lhalo-uhalo) :: esend,wrecv
    integer,dimension(local_ewn-1,lhalo) :: nrecv,ssend
    integer,dimension(local_ewn-1,uhalo-1) :: nsend,srecv
    ! begin
    call mpi_irecv(wrecv,size(wrecv),mpi_real8,west,west,&
         comm,wrequest,ierror)
    call mpi_irecv(erecv,size(erecv),mpi_real8,east,east,&
         comm,erequest,ierror)
    call mpi_irecv(srecv,size(srecv),mpi_real8,south,south,&
         comm,srequest,ierror)
    call mpi_irecv(nrecv,size(nrecv),mpi_real8,north,north,&
         comm,nrequest,ierror)

    esend(:,:) = a(1+lhalo:1+lhalo+uhalo-2,1+lhalo:local_nsn-uhalo)
    call mpi_send(esend,size(esend),mpi_real8,east,this_rank,comm,ierror)
    wsend(:,:) = &
         a(local_ewn-lhalo-uhalo+1:local_ewn-uhalo,1+lhalo:local_nsn-uhalo)
    call mpi_send(wsend,size(wsend),mpi_real8,west,this_rank,comm,ierror)

    call mpi_wait(wrequest,mpi_status_ignore,ierror)
    a(local_ewn-uhalo+1:,1+lhalo:local_nsn-uhalo) = wrecv(:,:)
    call mpi_wait(erequest,mpi_status_ignore,ierror)
    a(:lhalo,1+lhalo:local_nsn-uhalo) = erecv(:,:)

    nsend(:,:) = a(:,1+lhalo:1+lhalo+uhalo-2)
    call mpi_send(nsend,size(nsend),mpi_real8,north,this_rank,comm,ierror)
    ssend(:,:) = a(:,local_nsn-lhalo-uhalo+1:local_nsn-uhalo)
    call mpi_send(ssend,size(ssend),mpi_real8,south,this_rank,comm,ierror)

    call mpi_wait(srequest,mpi_status_ignore,ierror)
    a(:,local_nsn-uhalo+1:) = srecv(:,:)
    call mpi_wait(nrequest,mpi_status_ignore,ierror)
    a(:,:lhalo) = nrecv(:,:)
  end subroutine

  subroutine parallel_velo_halo_real8_3d(a) 
    use mpi
    implicit none
    real(8),dimension(:,:,:) :: a !(:,local_ewn-1,local_nsn-1)

    integer :: erequest,ierror,nrequest,srequest,wrequest
    integer,dimension(size(a,1),lhalo,local_nsn-lhalo-uhalo) :: erecv,wsend
    integer,dimension(size(a,1),uhalo-1,local_nsn-lhalo-uhalo) :: esend,wrecv
    integer,dimension(size(a,1),local_ewn-1,lhalo) :: nrecv,ssend
    integer,dimension(size(a,1),local_ewn-1,uhalo-1) :: nsend,srecv
    ! begin
    call mpi_irecv(wrecv,size(wrecv),mpi_real8,west,west,&
         comm,wrequest,ierror)
    call mpi_irecv(erecv,size(erecv),mpi_real8,east,east,&
         comm,erequest,ierror)
    call mpi_irecv(srecv,size(srecv),mpi_real8,south,south,&
         comm,srequest,ierror)
    call mpi_irecv(nrecv,size(nrecv),mpi_real8,north,north,&
         comm,nrequest,ierror)

    esend(:,:,:) = a(:,1+lhalo:1+lhalo+uhalo-2,1+lhalo:local_nsn-uhalo)
    call mpi_send(esend,size(esend),mpi_real8,east,this_rank,comm,ierror)
    wsend(:,:,:) = &
         a(:,local_ewn-lhalo-uhalo+1:local_ewn-uhalo,1+lhalo:local_nsn-uhalo)
    call mpi_send(wsend,size(wsend),mpi_real8,west,this_rank,comm,ierror)

    call mpi_wait(wrequest,mpi_status_ignore,ierror)
    a(:,local_ewn-uhalo+1:,1+lhalo:local_nsn-uhalo) = wrecv(:,:,:)
    call mpi_wait(erequest,mpi_status_ignore,ierror)
    a(:,:lhalo,1+lhalo:local_nsn-uhalo) = erecv(:,:,:)

    nsend(:,:,:) = a(:,:,1+lhalo:1+lhalo+uhalo-2)
    call mpi_send(nsend,size(nsend),mpi_real8,north,this_rank,comm,ierror)
    ssend(:,:,:) = a(:,:,local_nsn-lhalo-uhalo+1:local_nsn-uhalo)
    call mpi_send(ssend,size(ssend),mpi_real8,south,this_rank,comm,ierror)

    call mpi_wait(srequest,mpi_status_ignore,ierror)
    a(:,:,local_nsn-uhalo+1:) = srecv(:,:,:)
    call mpi_wait(nrequest,mpi_status_ignore,ierror)
    a(:,:,:lhalo) = nrecv(:,:,:)
  end subroutine

end module parallel
