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
  implicit none

  type(erosion_type) :: er           ! erosion
  type(glide_global_type) :: model        ! model instance
  type(ConfigSection), pointer :: config  ! configuration stuff
  character(len=50) :: fname   ! name of paramter file
  real(kind=rk) time

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
  time = model%numerics%tstart
  ! set initial sediment distrib
  er%seds1(13:18,23:28) = 5.
  er%seds1(36:40,15:20) = 5.
  write(*,*) 'total concentration ', sum(er%seds1)
  ! set up sparse matrix
  call calc_lagrange(model, er%trans, real(model%numerics%tinc,kind=dp) , er%lag_seds1)
  call print_sparse(er%lag_seds1,101)
  do while(time.le.model%numerics%tstart+10000.)
     model%numerics%time = time
     call transport_scalar(model,er%trans,er%seds1,er%lag_seds1)
     call glide_io_writeall(model,model)
     call erosion_io_writeall(er,model)
     time = time + model%numerics%tinc
  end do

  ! finalise GLIDE
  call er_finalise(er)
  call glide_finalise(model)
end program testtransport
