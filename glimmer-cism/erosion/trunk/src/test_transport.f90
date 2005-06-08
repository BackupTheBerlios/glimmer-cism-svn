! test_transport.f90
! Magnus Hagdorn, February 2005
!
! testing the transport code

program testtransport
  use glimmer_global, only:rk
  use glide
  use glimmer_log
  use glimmer_config
  use erosion
  use glide_io
  use paramets, only : len0
  implicit none

  type(erosion_type) :: er           ! erosion
  type(glide_global_type) :: model        ! model instance
  type(ConfigSection), pointer :: config  ! configuration stuff
  character(len=50) :: fname   ! name of paramter file
  real(kind=rk) time
  real, dimension(2) :: start,end
  integer, dimension(2):: s,e

  write(*,*) 'Enter name of GLIDE configuration file to be read'
  read(*,*) fname
  
  ! start logging
  call open_log(unit=50)
  
  ! read configuration
  call ConfigRead(fname,config)  

  call glide_initialise(model,config)
  call er_initialise(er,config,model)
  ! fill dimension variables
  call glide_nc_fillall(model)
  call erosion_nc_fillall(model,er)
  time = model%numerics%tstart


  start(:) = 50000./len0
  end(:) = 150000./len0

  s(1) = int(start(1)/er%dew)
  s(2) = int(start(2)/er%dns)
  e(1) = int(end(1)/er%dew)
  e(2) = int(end(2)/er%dns)

  er%seds1 = 0.
  er%seds1(s(1):e(1),s(2):e(2)) = 5.

  write(*,*) 'total concentration ', sum(er%seds1)
  ! set up sparse matrix
  call set_velos(model%velocity%ubas,model%velocity%vbas,-1.d0)
  call calc_lagrange(er, er%trans, er%dt , er%lag_seds1)
  do while(time.le.model%numerics%tend)
     model%numerics%time = time
     call glide_io_writeall(model,model)
     call erosion_io_writeall(er,model)
     if (model%numerics%tinc .gt. mod(model%numerics%time,model%numerics%tinc*er%ndt)) then
        call transport_scalar(er,er%trans,er%seds1,er%lag_seds1)
     end if
     time = time + model%numerics%tinc
  end do
  write(*,*) 'total concentration ', sum(er%seds1)

  ! finalise GLIDE
  call er_finalise(er)
  call glide_finalise(model)
end program testtransport
