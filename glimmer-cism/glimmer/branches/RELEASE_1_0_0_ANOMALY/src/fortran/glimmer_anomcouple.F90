! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  glimmer_anomcouple.f90 - part of the GLIMMER ice model   + 
! +                                                           +
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! 
! Copyright (C) 2004 GLIMMER contributors - see COPYRIGHT file 
! for list of contributors.
!
! This program is free software; you can redistribute it and/or 
! modify it under the terms of the GNU General Public License as 
! published by the Free Software Foundation; either version 2 of 
! the License, or (at your option) any later version.
!
! This program is distributed in the hope that it will be useful, 
! but WITHOUT ANY WARRANTY; without even the implied warranty of 
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the 
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License 
! along with this program; if not, write to the Free Software 
! Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 
! 02111-1307 USA
!
! GLIMMER is maintained by:
!
! Ian Rutt
! School of Geographical Sciences
! University of Bristol
! University Road
! Bristol
! BS8 1SS
! UK
!
! email: <i.c.rutt@bristol.ac.uk> or <ian.rutt@physics.org>
!
! GLIMMER is hosted on NeSCForge:
!
! http://forge.nesc.ac.uk/projects/glimmer/
!
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#ifdef HAVE_CONFIG_H
#include <config.inc>
#endif

module glimmer_anomcouple

  !*FD This module provides code to handle anomaly coupling. Although
  !*FD written for use with GLINT, it has general applicability. Temperature coupling
  !*FD is done linearly, precipitation proportionally.

  use glimmer_global

  implicit none

  character(4),parameter :: pvarname='prcp'
  character(4),parameter :: tvarname='artm'
  character(8),dimension(5),parameter :: xvars=(/ &
       'latitude', &
       'lat     ', &
       'x0      ', &
       'x1      ', &
       'x       '/)
  character(9),dimension(5),parameter :: yvars=(/ &
       'longitude', &
       'lon      ', &
       'y0       ', &
       'y1       ', &
       'y        '/)
  character(4),dimension(1),parameter :: tvars=(/ &
       'time'/)

  type anomaly_coupling
     logical :: enabled = .false.
     character(fname_length) :: fname_reference !*FD File containing reference climate
     character(fname_length) :: fname_modelclim !*FD File containing mean model climate
     integer :: nslices                         !*FD Number of time-slices in climatologies
     real(rk),dimension(:,:,:),pointer :: temp_ref => null() !*FD Reference climate (temperature)
     real(rk),dimension(:,:,:),pointer :: temp_mod => null() !*FD Model climate (temperature)
     real(rk),dimension(:,:,:),pointer :: prcp_ref => null() !*FD Reference climate (precip)
     real(rk),dimension(:,:,:),pointer :: prcp_mod => null() !*FD Model climate (precip)
     real(rk),dimension(:)    ,pointer :: time     => null() !*FD Time axis (fraction of year)
     integer :: nx,ny !*FD Grid dimensions (for convenience)
  end type anomaly_coupling

  private
  public :: anomaly_coupling, anomaly_init, anomaly_calc

  !MAKE_RESTART
#ifdef RESTARTS
#define RST_GLIMMER_ANOMCOUPLE
#include "glimmer_rst_head.inc"
#undef RST_GLIMMER_ANOMCOUPLE
#endif

contains

#ifdef RESTARTS
#define RST_GLIMMER_ANOMCOUPLE
#include "glimmer_rst_body.inc"
#undef RST_GLIMMER_ANOMCOUPLE
#endif

  subroutine anomaly_init(params,config)

    use glimmer_config

    type(anomaly_coupling),intent(out) :: params !*FD Parameters to be initialised
    type(ConfigSection),pointer        :: config !*FD Configuation file

    call anomaly_readconfig(params,config)
    if (params%enabled) then
       call anomaly_readdata(params)
       call anomaly_printconfig(params)
    end if

  end subroutine anomaly_init

  !------------------------------------------------------------------------------------------

  subroutine anomaly_calc(params,time,rawtemp,rawprcp,anomtemp,anomprcp)

    type(anomaly_coupling),intent(in) :: params !*FD Parameters to be initialised
    real(rk) :: time
    real(rk),dimension(:,:),intent(in)  :: rawtemp, rawprcp
    real(rk),dimension(:,:),intent(out) :: anomtemp,anomprcp

    real(rk),dimension(size(rawtemp,1),size(rawtemp,2)) :: tempm,prcpm,tempr,prcpr
    integer  :: first
    real(sp) :: frac

    if (params%enabled) then
       call anomaly_index(params%time,time,first,frac)
       tempm=(1.0-frac)*params%temp_mod(:,:,first)+frac*params%temp_mod(:,:,first+1)
       prcpm=(1.0-frac)*params%prcp_mod(:,:,first)+frac*params%prcp_mod(:,:,first+1)
       tempr=(1.0-frac)*params%temp_ref(:,:,first)+frac*params%temp_ref(:,:,first+1)
       prcpr=(1.0-frac)*params%prcp_ref(:,:,first)+frac*params%prcp_ref(:,:,first+1)
       anomtemp=rawtemp-tempm+tempr
       anomprcp=rawprcp*prcpr/prcpm
    else
       anomprcp=rawprcp
       anomtemp=rawtemp
    end if

  end subroutine anomaly_calc

  !------------------------------------------------------------------------------------------
  ! PRIVATE subroutines
  !------------------------------------------------------------------------------------------

  subroutine anomaly_readconfig(params,config)

    use glimmer_config

    type(anomaly_coupling),intent(out) :: params !*FD Parameters to be initialised
    type(ConfigSection),pointer        :: config !*FD Configuation file

    ! local variables
    type(ConfigSection), pointer :: section
    
    call GetSection(config,section,'anomaly coupling')

    if (associated(section)) then
       call GetValue(section,'reference',params%fname_reference)
       call GetValue(section,'model',    params%fname_modelclim)
       params%enabled = .true.
    else
       params%enabled = .false.
    end if

  end subroutine anomaly_readconfig

  !------------------------------------------------------------------------------------------

  subroutine anomaly_printconfig(params)

    use glimmer_log

    type(anomaly_coupling),intent(inout) :: params !*FD Parameters to be initialised
    character(100) :: message

    call write_log_div
    
    call write_log('Anomaly coupling')
    call write_log("----------------")
    write(message,*)"Reference climate:",trim(params%fname_reference)
    call write_log(message)
    write(message,*)"Model climate:    ",trim(params%fname_modelclim)
    call write_log(message)
    write(message,*)"Number of slices: ",params%nslices
    call write_log(message)
    call write_log("")

  end subroutine anomaly_printconfig

  !------------------------------------------------------------------------------------------

  subroutine anomaly_readdata(params)

    use glimmer_log

    type(anomaly_coupling),intent(inout) :: params !*FD Parameters to be initialised

    integer,dimension(4) :: nx,ny,nt
    real(rk),dimension(:),pointer :: timemod => null()
    real(rk),dimension(:),pointer :: timeref => null()

    call anomaly_readnc(params%fname_reference,pvarname,params%prcp_ref,timeref,nx(1),ny(1),nt(1))
    call anomaly_readnc(params%fname_reference,tvarname,params%temp_ref,timeref,nx(2),ny(2),nt(2))
    call anomaly_readnc(params%fname_modelclim,pvarname,params%prcp_mod,timemod,nx(3),ny(3),nt(3))
    call anomaly_readnc(params%fname_modelclim,tvarname,params%temp_mod,timemod,nx(4),ny(4),nt(4))

    if (any(nx(1)/=nx(2:4)).or.any(ny(1)/=ny(2:4)).or.any(nt(1)/=nt(2:4))) &
         call write_log("Anomaly coupling: sizes of arrays in climate files do not agree", &
         GM_FATAL,__FILE__,__LINE__)

    params%nx=nx(1)
    params%ny=ny(1)
    params%nslices=nt(1)

    if (.not.all(timemod==timeref)) &
         call write_log("Anomaly coupling: time axes in climate files do not agree",GM_FATAL,__FILE__,__LINE__)

    if (associated(params%time)) then
       deallocate(params%time)
       params%time => null()
    end if
    allocate(params%time(params%nslices+2))

    params%time=timemod

  end subroutine anomaly_readdata

  !------------------------------------------------------------------------------------------

  subroutine anomaly_readnc(fname,varname,data,timeaxis,nx,ny,nt)

    use netcdf
    use glimmer_log

    character(*),                     intent(in)  :: fname
    character(*),                     intent(in)  :: varname
    real(rk),dimension(:,:,:),pointer             :: data
    real(rk),dimension(:),    pointer             :: timeaxis
    integer,                          intent(out) :: nx,ny,nt

    ! Local variables
    integer :: status, ncid, varid, tvarid, ndims, i
    integer,dimension(3) :: dimids
    integer,dimension(3) :: dimnames
    character(30) :: dntemp,timevar
    real(sp) :: interval
    
    ! Open file
    status=nf90_open(fname,NF90_NOWRITE,ncid)
    call nc_errorhandle(__FILE__,__LINE__,status)

    ! Look for required variable
    status=nf90_inq_varid(ncid,varname,varid)
    call nc_errorhandle(__FILE__,__LINE__,status)

    ! Check we have three dimensions
    status=nf90_inquire_variable(ncid,varid,ndims=ndims)
    call nc_errorhandle(__FILE__,__LINE__,status)
    if (ndims/=3) call write_log("Anomaly coupling: file "//trim(fname)//", variable "// &
         trim(varname)//" should have three dimensions",GM_FATAL,__FILE__,__LINE__)
    status=nf90_inquire_variable(ncid,varid,dimids=dimids)
    call nc_errorhandle(__FILE__,__LINE__,status)

    ! Check we have some sensible dimension names in the right order:
    ! must be x,y,t (this is t,y,x in netcdf-speak...)

    status=nf90_inquire_dimension(ncid,dimids(1),dntemp,nx)
    call nc_errorhandle(__FILE__,__LINE__,status)
    if (.not.any(dntemp==xvars)) &
         call write_log("Anomaly coupling: first dimension in climate file "//trim(fname)// &
         " is not x or longitude",GM_FATAL,__FILE__,__LINE__)

    status=nf90_inquire_dimension(ncid,dimids(2),dntemp,ny)
    call nc_errorhandle(__FILE__,__LINE__,status)
    if (.not.any(dntemp==yvars)) &
         call write_log("Anomaly coupling: second dimension in climate file "//trim(fname)// &
         " is not y or latitude",GM_FATAL,__FILE__,__LINE__)

    status=nf90_inquire_dimension(ncid,dimids(3),dntemp,nt)
    call nc_errorhandle(__FILE__,__LINE__,status)
    if (.not.any(dntemp==tvars)) &
         call write_log("Anomaly coupling: third dimension in climate file "//trim(fname)// &
         " is not time",GM_FATAL,__FILE__,__LINE__)
    timevar=dntemp

    ! If we've got this far, we can retrieve the data itself

    if (associated(data)) then
       deallocate(data)
       data => null()
    end if
    allocate(data(nx,ny,nt+2))

    status=nf90_get_var(ncid,varid,data(:,:,2:nt+1))
    call nc_errorhandle(__FILE__,__LINE__,status)

    data(:,:,1)   =data(:,:,nt+1)
    data(:,:,nt+2)=data(:,:,2)

    ! Now we try and get the time data

    if (associated(timeaxis)) then
       deallocate(timeaxis)
       timeaxis => null()
    end if
    allocate(timeaxis(nt+2))

    status=nf90_inq_varid(ncid,timevar,tvarid)
    select case (status)
    case(NF90_NOERR)
       status=nf90_get_var(ncid,tvarid,timeaxis(2:nt+1))
       call nc_errorhandle(__FILE__,__LINE__,status)
    case(NF90_ENOTVAR)
       ! Time variable not found - construct our own
       interval=1.0/real(nt)
       do i=1,nt
          timeaxis(i+1)=(i-1)*interval+interval/2.0
       end do
       call write_log('Anomaly coupling: Created time-axis')
    case default
       ! Some other error - bail out
       call nc_errorhandle(__FILE__,__LINE__,status)
    end select

    ! Fix up time boundaries
    timeaxis(1)   =1.0-timeaxis(nt+1)
    timeaxis(nt+2)=1.0+timeaxis(2)

    ! Close the file
    status=nf90_close(ncid)
    call nc_errorhandle(__FILE__,__LINE__,status)

  end subroutine anomaly_readnc

  subroutine anomaly_index(timeaxis,time,first,frac)

    use glimmer_log

    real(rk),dimension(:),intent(in)  :: timeaxis
    real(rk),             intent(in)  :: time
    integer,              intent(out) :: first
    real(sp),             intent(out) :: frac 

    first=1

    do
       if (time>=timeaxis(first).and.time<timeaxis(first+1)) then
          frac=(time-timeaxis(first))/(timeaxis(first+1)-timeaxis(first))
          exit
       endif
       first=first+1
       if (first==size(timeaxis)) then
          call write_log("Anomaly coupling: Problem indexing time-slices",GM_FATAL,__FILE__,__LINE__)
       end if
    end do

  end subroutine anomaly_index

end module glimmer_anomcouple
