! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  glimmer_CFproj.f90 - part of the GLIMMER ice model       + 
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
! GLIMMER is hosted on berliOS.de:
!
! https://developer.berlios.de/projects/glimmer-cism/
!
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#ifdef HAVE_CONFIG_H
#include "config.inc"
#endif

module glimmer_map_CFproj

  !*FD Holds derived types and subroutines
  !*FD necessary for handling map projections. Most of the component
  !*FD names of the various derived types are self-explanatory.
  !*FD Note that this doesn't currently interface with the proj4
  !*FD library in anyway, it simply handles NetCDF data and projection
  !*FD parameters in an appropriate format.

  use glimmer_map_types

  implicit none

  private
  public  glimmap_CFGetProj,glimmap_CFPutProj

contains
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! public functions
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function glimmap_CFGetProj(ncid)

    !*FD Read projection from a given netCDF file, returning
    !*FD an instance of type \texttt{glimmap\_proj}.
    !*FDRV Derived type instance containing projection parameters

    use netcdf
    use glimmer_log
    use glimmer_map_init

    implicit none

    type(glimmap_proj) :: glimmap_CFGetProj
    integer, intent(in) :: ncid                !*FD Handle of the file to be read.
    
    !local variables
    integer status
    integer nvars, varid
    integer natts, attid
    logical found_map
    character(len=50) :: attname,mapname

    ! getting variables
    status = nf90_inquire(ncid,nvariables=nvars)
    call nc_errorhandle(__FILE__,__LINE__,status)    
    
    ! looping over variables
    found_map=.false.
    do varid=1,nvars
       status = nf90_inquire_variable(ncid,varid,natts=natts)
       ! and loop over attributes
       do attid=1,natts
          status = nf90_inq_attname(ncid,varid,attid,attname)
          if (trim(attname) == 'grid_mapping_name') then
             found_map = .true.
             status = nf90_get_att(ncid,varid,attname,mapname)
             mapname = adjustl(mapname)
             call nc_errorhandle(__FILE__,__LINE__,status)
             exit
          end if
       end do
       if (found_map) exit
    end do

    if (found_map) then
       glimmap_CFGetProj%found = .true.
       if (index(mapname,'lambert_azimuthal_equal_area').ne.0) then
          glimmap_CFGetProj%laea => CFproj_get_laea(ncid,varid)
          call glimmap_laea_init(glimmap_CFGetProj%laea)
       else if (index(mapname,'albers_conical_equal_area').ne.0) then
          glimmap_CFGetProj%aea => CFproj_get_aea(ncid,varid)
          call glimmap_aea_init(glimmap_CFGetProj%aea)
       else if (index(mapname,'lambert_conformal_conic').ne.0) then
          glimmap_CFGetProj%lcc => CFproj_get_lcc(ncid,varid)
          call glimmap_lcc_init(glimmap_CFGetProj%lcc)
       else if (index(mapname,'polar_stereographic').ne.0) then
          glimmap_CFGetProj%stere => CFproj_get_stere_polar(ncid,varid)
          call glimmap_stere_init(glimmap_CFGetProj%stere)
       else if (index(mapname,'stereographic').ne.0) then
          glimmap_CFGetProj%stere => CFproj_get_stere(ncid,varid)
          call glimmap_stere_init(glimmap_CFGetProj%stere)
       else
          glimmap_CFGetProj%found = .false.
          call write_log('Do not know about this projection: '//(mapname),GM_ERROR)
       end if
    else
        glimmap_CFGetProj%found = .false.
       call write_log('No map projection found',GM_WARNING)
    end if
  end function glimmap_CFGetProj

  !-------------------------------------------------------------------------

  subroutine glimmap_CFPutProj(ncid,mapid,proj)

    !*FD write projection to a netCDF file.

    use netcdf
    use glimmer_log

    implicit none

    type(glimmap_proj) :: proj   !*FD Projection to be written.
    integer, intent(in) :: ncid             !*FD Handle of netCDF file.
    integer, intent(in) :: mapid            !*FD Handle of map projection in netCDF file.

    if (.not.glimmap_allocated(proj)) then
       call write_log('No known projection found!',GM_WARNING)
       return
    end if

    if (associated(proj%laea)) then
       call CFproj_put_laea(ncid,mapid,proj%laea)
       return
    else if (associated(proj%aea)) then
       call CFproj_put_aea(ncid,mapid,proj%aea)
       return
    else if (associated(proj%lcc)) then
       call CFproj_put_lcc(ncid,mapid,proj%lcc)
       return
    else if (associated(proj%stere)) then
       call CFproj_put_stere(ncid,mapid,proj%stere)
       return
    else
       call write_log('No known projection found!',GM_WARNING)
    end if
  end subroutine glimmap_CFPutProj

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! private readers
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function CFproj_get_stere(ncid,mapid)
    use netcdf
    implicit none
    type(proj_stere), pointer :: CFproj_get_stere
    integer, intent(in) :: ncid
    integer, intent(in) :: mapid
    
    integer status

    allocate(CFproj_get_stere)
    status = nf90_get_att(ncid,mapid,'false_easting',CFproj_get_stere%false_easting)
    call nc_errorhandle(__FILE__,__LINE__,status)
    status = nf90_get_att(ncid,mapid,'false_northing',CFproj_get_stere%false_northing)
    call nc_errorhandle(__FILE__,__LINE__,status)
    status = nf90_get_att(ncid,mapid,'longitude_of_projection_origin',CFproj_get_stere%longitude_of_central_meridian)
    call nc_errorhandle(__FILE__,__LINE__,status)
    status = nf90_get_att(ncid,mapid,'latitude_of_projection_origin',CFproj_get_stere%latitude_of_projection_origin)
    call nc_errorhandle(__FILE__,__LINE__,status)
    status = nf90_get_att(ncid,mapid,'scale_factor_at_projection_origin',CFproj_get_stere%scale_factor_at_proj_origin)
    call nc_errorhandle(__FILE__,__LINE__,status)

  end function CFproj_get_stere

  function CFproj_get_stere_polar(ncid,mapid)
    use netcdf
    use glimmer_log
    implicit none
    type(proj_stere), pointer :: CFproj_get_stere_polar
    integer, intent(in) :: ncid
    integer, intent(in) :: mapid
    
    integer status
    real dummy

    allocate(CFproj_get_stere_polar)
    status = nf90_get_att(ncid,mapid,'false_easting',CFproj_get_stere_polar%false_easting)
    call nc_errorhandle(__FILE__,__LINE__,status)
    status = nf90_get_att(ncid,mapid,'false_northing',CFproj_get_stere_polar%false_northing)
    call nc_errorhandle(__FILE__,__LINE__,status)
    status = nf90_get_att(ncid,mapid,'straight_vertical_longitude_from_pole',CFproj_get_stere_polar%longitude_of_central_meridian)
    call nc_errorhandle(__FILE__,__LINE__,status)
    CFproj_get_stere_polar%latitude_of_projection_origin=90.0
    status = nf90_get_att(ncid,mapid,'latitude_of_projection_origin',CFproj_get_stere_polar%latitude_of_projection_origin)
    call nc_errorhandle(__FILE__,__LINE__,status)
    if (abs(abs(CFproj_get_stere_polar%latitude_of_projection_origin)-90.0)>0.001) then
       call write_log('Error (polar stereographic projection) latitude of origin must be +-90.0',&
            GM_FATAL,__FILE__,__LINE__)
    end if
    status = nf90_get_att(ncid,mapid,'scale_factor_at_projection_origin',dummy)
    if (status.eq.NF90_NOERR) then
       CFproj_get_stere_polar%scale_factor_at_proj_origin = dummy
    end if
    status = nf90_get_att(ncid,mapid,'standard_parallel',dummy)
    if (status.eq.NF90_NOERR) then
       CFproj_get_stere_polar%standard_parallel = dummy
    end if
    if (CFproj_get_stere_polar%standard_parallel.ne.0 .and. CFproj_get_stere_polar%scale_factor_at_proj_origin.ne.0.) then
       call write_log('Error (stereographic projection), can only handle either standard_parallel or scale_at_orig',&
            GM_FATAL,__FILE__,__LINE__)
    end if
  end function CFproj_get_stere_polar

  function CFproj_get_laea(ncid,mapid)
    use netcdf
    implicit none
    type(proj_laea), pointer :: CFproj_get_laea
    integer, intent(in) :: ncid
    integer, intent(in) :: mapid
    
    integer status
    allocate(CFproj_get_laea)
    status = nf90_get_att(ncid,mapid,'false_easting',CFproj_get_laea%false_easting)
    call nc_errorhandle(__FILE__,__LINE__,status)
    status = nf90_get_att(ncid,mapid,'false_northing',CFproj_get_laea%false_northing)
    call nc_errorhandle(__FILE__,__LINE__,status)
    status = nf90_get_att(ncid,mapid,'longitude_of_projection_origin',CFproj_get_laea%longitude_of_central_meridian)
    call nc_errorhandle(__FILE__,__LINE__,status)
    status = nf90_get_att(ncid,mapid,'latitude_of_projection_origin',CFproj_get_laea%latitude_of_projection_origin)
    call nc_errorhandle(__FILE__,__LINE__,status)
  end function CFproj_get_laea

  function CFproj_get_aea(ncid,mapid)
    use netcdf
    implicit none
    type(proj_aea), pointer :: CFproj_get_aea
    integer, intent(in) :: ncid
    integer, intent(in) :: mapid
    
    integer status
    allocate(CFproj_get_aea)
    status = nf90_get_att(ncid,mapid,'false_easting',CFproj_get_aea%false_easting)
    call nc_errorhandle(__FILE__,__LINE__,status)
    status = nf90_get_att(ncid,mapid,'false_northing',CFproj_get_aea%false_northing)
    call nc_errorhandle(__FILE__,__LINE__,status)
    status = nf90_get_att(ncid,mapid,'longitude_of_central_meridian',CFproj_get_aea%longitude_of_central_meridian)
    call nc_errorhandle(__FILE__,__LINE__,status)
    status = nf90_get_att(ncid,mapid,'latitude_of_projection_origin',CFproj_get_aea%latitude_of_projection_origin)
    call nc_errorhandle(__FILE__,__LINE__,status)
    status = nf90_get_att(ncid,mapid,'standard_parallel',CFproj_get_aea%standard_parallel)
    call nc_errorhandle(__FILE__,__LINE__,status)
  end function CFproj_get_aea

  function CFproj_get_lcc(ncid,mapid)
    use netcdf
    implicit none
    type(proj_lcc), pointer :: CFproj_get_lcc
    integer, intent(in) :: ncid
    integer, intent(in) :: mapid
    
    integer status
    allocate(CFproj_get_lcc)
    status = nf90_get_att(ncid,mapid,'false_easting',CFproj_get_lcc%false_easting)
    call nc_errorhandle(__FILE__,__LINE__,status)
    status = nf90_get_att(ncid,mapid,'false_northing',CFproj_get_lcc%false_northing)
    call nc_errorhandle(__FILE__,__LINE__,status)
    status = nf90_get_att(ncid,mapid,'longitude_of_central_meridian',CFproj_get_lcc%longitude_of_central_meridian)
    call nc_errorhandle(__FILE__,__LINE__,status)
    status = nf90_get_att(ncid,mapid,'latitude_of_projection_origin',CFproj_get_lcc%latitude_of_projection_origin)
    call nc_errorhandle(__FILE__,__LINE__,status)
    status = nf90_get_att(ncid,mapid,'standard_parallel',CFproj_get_lcc%standard_parallel)
    call nc_errorhandle(__FILE__,__LINE__,status)
  end function CFproj_get_lcc

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! private subroutines to write projection info
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine CFproj_put_stere(ncid,mapid,stere)
    use netcdf
    implicit none
    type(proj_stere), pointer :: stere
    integer, intent(in) :: ncid
    integer, intent(in) :: mapid

    integer status

    if (stere%pole/=0) then
       status = nf90_put_att(ncid,mapid,'grid_mapping_name','polar_stereographic')
    else
       status = nf90_put_att(ncid,mapid,'grid_mapping_name','stereographic')
    end if
    call nc_errorhandle(__FILE__,__LINE__,status)
    status = nf90_put_att(ncid,mapid,'false_easting',stere%false_easting)
    call nc_errorhandle(__FILE__,__LINE__,status)
    status = nf90_put_att(ncid,mapid,'false_northing',stere%false_northing)
    call nc_errorhandle(__FILE__,__LINE__,status)
    if (stere%pole/=0) then
       status = nf90_put_att(ncid,mapid,'straight_vertical_longitude_from_pole',stere%longitude_of_central_meridian)
    else
       status = nf90_put_att(ncid,mapid,'longitude_of_projection_origin',stere%longitude_of_central_meridian)
    end if
    call nc_errorhandle(__FILE__,__LINE__,status)
    status = nf90_put_att(ncid,mapid,'latitude_of_projection_origin',stere%latitude_of_projection_origin)
    call nc_errorhandle(__FILE__,__LINE__,status)
    if (stere%pole/=0) then
       if (stere%standard_parallel.ne.0) then
          status = nf90_put_att(ncid,mapid,'standard_parallel',stere%standard_parallel)
       else
          status = nf90_put_att(ncid,mapid,'scale_factor_at_projection_origin',stere%scale_factor_at_proj_origin)
       end if
    else
       status = nf90_put_att(ncid,mapid,'scale_factor_at_projection_origin',stere%scale_factor_at_proj_origin)
    end if
    call nc_errorhandle(__FILE__,__LINE__,status)
  end subroutine CFproj_put_stere

  subroutine CFproj_put_laea(ncid,mapid,laea)
    use netcdf
    implicit none
    type(proj_laea), pointer :: laea
    integer, intent(in) :: ncid
    integer, intent(in) :: mapid

    integer status

    status = nf90_put_att(ncid,mapid,'grid_mapping_name','lambert_azimuthal_equal_area')
    call nc_errorhandle(__FILE__,__LINE__,status)
    status = nf90_put_att(ncid,mapid,'false_easting',laea%false_easting)
    call nc_errorhandle(__FILE__,__LINE__,status)
    status = nf90_put_att(ncid,mapid,'false_northing',laea%false_northing)
    call nc_errorhandle(__FILE__,__LINE__,status)
    status = nf90_put_att(ncid,mapid,'longitude_of_projection_origin',laea%longitude_of_central_meridian)
    call nc_errorhandle(__FILE__,__LINE__,status)
    status = nf90_put_att(ncid,mapid,'latitude_of_projection_origin',laea%latitude_of_projection_origin)
    call nc_errorhandle(__FILE__,__LINE__,status)
  end subroutine CFproj_put_laea

  subroutine CFproj_put_aea(ncid,mapid,aea)
    use netcdf
    implicit none
    type(proj_aea), pointer :: aea
    integer, intent(in) :: ncid
    integer, intent(in) :: mapid

    integer status

    status = nf90_put_att(ncid,mapid,'grid_mapping_name','albers_conical_equal_area')
    call nc_errorhandle(__FILE__,__LINE__,status)
    status = nf90_put_att(ncid,mapid,'false_easting',aea%false_easting)
    call nc_errorhandle(__FILE__,__LINE__,status)
    status = nf90_put_att(ncid,mapid,'false_northing',aea%false_northing)
    call nc_errorhandle(__FILE__,__LINE__,status)
    status = nf90_put_att(ncid,mapid,'longitude_of_central_meridian',aea%longitude_of_central_meridian)
    call nc_errorhandle(__FILE__,__LINE__,status)
    status = nf90_put_att(ncid,mapid,'latitude_of_projection_origin',aea%latitude_of_projection_origin)
    call nc_errorhandle(__FILE__,__LINE__,status)
    status = nf90_put_att(ncid,mapid,'standard_parallel',aea%standard_parallel)
    call nc_errorhandle(__FILE__,__LINE__,status)
  end subroutine CFproj_put_aea

  subroutine CFproj_put_lcc(ncid,mapid,lcc)
    use netcdf
    implicit none
    type(proj_lcc), pointer :: lcc
    integer, intent(in) :: ncid
    integer, intent(in) :: mapid

    integer status

    status = nf90_put_att(ncid,mapid,'grid_mapping_name','lambert_conformal_conic')
    call nc_errorhandle(__FILE__,__LINE__,status)
    status = nf90_put_att(ncid,mapid,'false_easting',lcc%false_easting)
    call nc_errorhandle(__FILE__,__LINE__,status)
    status = nf90_put_att(ncid,mapid,'false_northing',lcc%false_northing)
    call nc_errorhandle(__FILE__,__LINE__,status)
    status = nf90_put_att(ncid,mapid,'longitude_of_central_meridian',lcc%longitude_of_central_meridian)
    call nc_errorhandle(__FILE__,__LINE__,status)
    status = nf90_put_att(ncid,mapid,'latitude_of_projection_origin',lcc%latitude_of_projection_origin)
    call nc_errorhandle(__FILE__,__LINE__,status)
    status = nf90_put_att(ncid,mapid,'standard_parallel',lcc%standard_parallel)
    call nc_errorhandle(__FILE__,__LINE__,status)
  end subroutine CFproj_put_lcc

end module glimmer_map_CFproj
