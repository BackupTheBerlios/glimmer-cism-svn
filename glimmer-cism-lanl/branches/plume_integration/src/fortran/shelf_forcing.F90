#ifdef HAVE_CONFIG_H
#include "config.inc"
#endif

module shelf_forcing

  use glimmer_global, only : sp

  type shelf_climate
     ! holds parameters for the shelf climate

     real(kind=sp) :: artm = 0.0
     real(kind=sp) :: accumulation_rate = 0.0
     real(kind=sp) :: eus = 500.0

  end type shelf_climate

  !MAKE_RESTART
#ifdef RESTARTS
#define RST_SIMPLE_FORCING
#include "glimmer_rst_head.inc"
#undef RST_SIMPLE_FORCING
#endif

contains

#ifdef RESTARTS
#define RST_SIMPLE_FORCING
#include "glimmer_rst_body.inc"
#undef RST_SIMPLE_FORCING
#endif

  subroutine shelf_config_initialise(climate_cfg,config)
    !*FD initialise shelf climate model
    use glimmer_paramets, only: thk0, acc0, scyr
    use glimmer_config
    implicit none
    type(shelf_climate) :: climate_cfg         !*FD structure holding climate info
    type(ConfigSection), pointer :: config  !*FD structure holding sections of configuration file   

  
    call shelf_readconfig(climate_cfg,config)
    call shelf_printconfig(climate_cfg)
           
  end subroutine shelf_config_initialise

  subroutine shelf_readconfig(climate_cfg, config)
    !*FD read configuration
    use glimmer_log
    use glimmer_config
    implicit none
    type(shelf_climate) :: climate_cfg         !*FD structure holding climate info
    type(ConfigSection), pointer :: config  !*FD structure holding sections of configuration file   

    ! local variables
    type(ConfigSection), pointer :: section
    real(kind=sp), dimension(:), pointer :: dummy

    call GetSection(config,section,'Petermann shelf')
    if (associated(section)) then
        call GetValue(section,'air_temperature',climate_cfg%artm)
        call GetValue(section,'accumulation_rate',climate_cfg%accumulation_rate)   
        call GetValue(section,'eustatic_sea_level',climate_cfg%eus)
        return
    else
       !log error
       call write_log('No Petermann section',GM_FATAL)
    end if
 
  end subroutine shelf_readconfig

  subroutine shelf_printconfig(climate_cfg)
    !*FD print simple climate configuration
    use glimmer_log
    implicit none
    type(shelf_climate) :: climate_cfg   !*FD structure holding climate info
    character(len=100) :: message

    call write_log_div
    call write_log('Petermann shelf configuration')
    call write_log('------------------------------------')
    write(message,*) 'air temperature  : ',climate_cfg%artm
    call write_log(message)
    write(message,*) 'accumulation rate: ',climate_cfg%accumulation_rate
    call write_log(message)
    write(message,*) 'eustatic sea level:',climate_cfg%eus
    call write_log(message)

  end subroutine shelf_printconfig
  
end module shelf_forcing
