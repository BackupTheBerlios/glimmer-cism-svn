! ******************************************************************************
! example.f90
! Magnus Hagdorn
!
! simple example climate driver demonstrating how to use the library
! ******************************************************************************
!
! ChangeLog
! 2004-11-10 Magnus Hagdorn
!  * initial version

program example

  ! load various modules
  use glimmer_global, only:rk ! precision of the model
  use glide                   ! main glide module
  use glimmer_log             ! module for logging messages
  use glimmer_config          ! module for handling configuration files
  use paramets                ! for the scales
  use physcon, only : scyr    ! number of seconds in a year
  implicit none

  ! some variables
  type(glide_global_type) :: model    
  ! this variable holds all the data associated with this particular model
  ! instance.
  type(ConfigSection), pointer :: config
  ! this pointer points to the first element of a linked list which contains
  ! all the configuration variables and their values
  character(len=50) :: fname   
  ! name of paramter file
  real(kind=rk) time
  ! current time
  
  integer ew,ns,ewct,nsct  ! loop variable and grid centre
  real grid, dist          ! node spacing and distance

  ! Ask for configuration file
  write(*,*) 'Enter name of GLIDE configuration file to be read'
  read(*,*) fname
  
  ! start logging
  call open_log(unit=50)
  
  ! read configuration
  call ConfigRead(fname,config)
  
  ! initialise GLIDE
  call glide_initialise(model,config)
  ! fill dimension variables
  call glide_nc_fillall(model)
  ! get current time from start time
  time = model%numerics%tstart
  
  ! setup some variables for BC
  ewct = real(model%general%ewn+1) / 2.0 ! the grid centre (x)
  nsct = real(model%general%nsn+1) / 2.0 ! and (y)
  grid = model%numerics%dew * len0       ! node-spacing

  ! loop over times
  do while(time.le.model%numerics%tend)
     ! setup boundary conditions
     ! this example implements the EISMINT-1 moving margin BC
     
     ! mass balance
     ! loop over grid
     do ns = 1,model%general%nsn
        do ew = 1,model%general%ewn
           ! calculate distance from centre
           dist = grid * sqrt((real(ew) - ewct)**2 + (real(ns) - nsct)**2)
           ! set mass balance
           model%climate%acab(ew,ns) = min(0.5, 1.05e-5 * (450.0e3 - dist))
        end do
     end do
     
     ! surface temperature
     model%climate%artm(:,:) = -3.150  -1.e-2 * model%geometry%thck(:,:)* thk0

     ! scale variables
     model%climate%acab(:,:) = model%climate%acab(:,:) / (acc0 * scyr)

     ! calculate temperature and velocity distribution
     call glide_tstep_p1(model,time)
     ! write to netCDF file, move ice
     call glide_tstep_p2(model)
     ! calculate isostatic adjustment
     call glide_tstep_p3(model)
     ! increment time counter
     time = time + model%numerics%tinc
  end do

  ! finalise GLIDE
  call glide_finalise(model)
end program example
