
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! This code is taken from the Generic Mapping Tools and was 
! converted into Fortran 90 by Ian Rutt.
!
! Original code (in C) Copyright (c) 1991-2003 by P. Wessel and W. H. F. Smith
!
! Partial translation into Fortran 90 (c) 2004 Ian C. Rutt
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
!
! The Generic Mapping Tools are maintained by Paul Wessel and 
! Walter H. F. Smith. The GMT homepage is:
!
! http://gmt.soest.hawaii.edu/
!
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

module gmt

!*FD This code is taken from the Generic Mapping Tools and was 
!*FD converted into Fortran 90 by Ian Rutt.
!*FD 
!*FD Original code (in C) Copyright (c) 1991-2003 by P. Wessel and W. H. F. Smith.
!*FD 
!*FD Partial translation into Fortran 90 (c) 2004 Ian C. Rutt.
!*FD 
!*FD To use this library:
!*FD \begin{itemize}
!*FD \item Include a use gmt statement at the start of the program;
!*FD \item Declare an instance of gmt\_pinf for each instance of the projection you
!*FD    want to use;
!*FD \item initialise the projection by calling \texttt{gmt\_set(latc,lonc,radea,project\_info)}
!*FD    with the centre of the projection (\texttt{latc},\texttt{lonc}) in degrees, the radius of the earth
!*FD    (\texttt{radea}), in metres (though often it's convenient to set this to 1.0), and the
!*FD    instance of \texttt{gmt\_pinf}.
!*FD \end{itemize}
!*FD 
!*FD NB: Longitude and latitude coordinates are in degrees, relative to the global origin,
!*FD     while x and y are relative to the centre of the projection, and given in metres (or whatever
!*FD     units are used for the radius of the Earth).

  use glimmer_global

  implicit none

  type gmt_pinf
    real(dp) :: cosp             = 0.0     !*FD cos ($\phi_p$)
    real(dp) :: sinp             = 0.0     !*FD sin ($\phi_p$)
    real(dp) :: pole             = 0.0     !*FD +90 or -90, depending on hemisphere
    real(dp) :: central_meridian = 0.0     !*FD Central meridian of projection
    real(dp) :: EQ_RAD           = 0.0     !*FD The radius of earth
    real(dp) :: i_EQ_RAD         = 0.0     !*FD 1/radius of earth.
    real(dp) :: Dx               = 0.0     !*FD A fudge factor for scaling the projection.
    real(dp) :: Dy               = 0.0     !*FD A fudge factor for scaling the projection.
    real(dp) :: iDx              = 0.0     !*FD Inverse of Dx.
    real(dp) :: iDy              = 0.0     !*FD Inverse of Dy.
    real(dp) :: s_c              = 0.0     !*FD Not sure what this does\ldots.
    real(dp) :: s_ic             = 0.0     !*FD Not sure what this does, but it's the inverse of {\tt s\_ic}\ldots.
    logical  :: n_polar          = .false. !*FD not needed - indicates if projection pole is north pole.
    logical  :: s_polar          = .false. !*FD not needed - indicates if projection pole is south pole.
    logical  :: north_pole       = .false. !*FD true if projection in northern hemisphere, false otherwise.
    logical  :: polar            = .false. !*FD Is this a polar projection?
  end type gmt_pinf

  real(dp),parameter :: pi=3.141592654          !*FD The value of $\pi$.
  real(dp),parameter :: M_PI_4=pi/4             !*FD The value of $\pi/4$.
  real(dp),parameter :: M_PI_2=pi/2             !*FD The value of $\pi/2$.
  real(dp),parameter :: D2R=pi/180.0            !*FD Degrees-to-radians conversion factor.
  real(dp),parameter :: R2D=180.0/pi            !*FD Radians-to-degrees conversion factor.
  real(dp),parameter :: DBL_MAX=huge(pi)        !*FD Ceiling value for double-precision real kind.
  real(dp),parameter :: GMT_CONV_LIMIT=1.0e-8   !*FD Convergence limit (a small number).
  logical, parameter :: GMT_convert_latitudes=.false. !*FD An unused constant, left over from C-to-f90 conversion.
  real(dp),parameter :: GMT_map_scale_factor = 1.0    !*FD An unused constant, left over from C-to-f90 conversion.
  real(dp),parameter :: SMALL=1.0e-4                  !*FD A small number.

  private sincos,hypot
  private DBL_MAX,GMT_CONV_LIMIT,GMT_convert_latitudes,M_PI_4

contains

  subroutine gmt_set(p_type,latc,lonc,radea,project_info)

    !*FD Initialise a GMT parameter type

    integer,intent(in) :: p_type  !*FD Type of projection to select
                                  !*FD 
                                  !*FD The projections available are:
                                  !*FD \begin{enumerate}
                                  !*FD \item Lambert Equal Area
                                  !*FD \item Spherical polar
                                  !*FD \item Spherical stereographic (oblique)
                                  !*FD \item Spherical stereographic (equatorial)
                                  !*FD \end{enumerate}
    real(dp),intent(in) :: latc   !*FD Centre of projection (Latitude)
    real(dp),intent(in) :: lonc   !*FD Centre of projection (Longitude)
    real(dp),intent(in) :: radea  !*FD Radius of the Earth (m)
    type(gmt_pinf),intent(inout) :: project_info !*FD GMT parameters to be initialised

    call gmt_type_init(project_info)

    project_info%EQ_RAD=radea                 ! radius of earth
    project_info%i_EQ_RAD=1.0/radea           ! 1/radius of earth

    select case(p_type)
    case(1)   ! Lambert equal area 
      call gmt_lambeq_set(latc,lonc,project_info)
    case(2:4)   ! Spherical/stereographic (polar,oblique equatorial) 
      call gmt_map_init_stereo(latc,lonc,project_info)
    case default
      print*,'* ERROR: ',p_type,' is not a valid projection type.'
      stop 
    end select

  end subroutine gmt_set

  subroutine gmt_type_init(project_info)

    !*FD Initialise the GMT projection parameters to their
    !*FD default values. This subroutine is redundant now we have
    !*FD moved to f95, which allows derived type initialisation.

    type(gmt_pinf),intent(out) :: project_info !*FD GMT parameters to be initialised

    project_info%cosp = 0.0
    project_info%sinp = 0.0
    project_info%pole = 0.0
    project_info%central_meridian = 0.0
    project_info%EQ_RAD = 0.0
    project_info%i_EQ_RAD = 0.0
    project_info%Dx = 0.0
    project_info%Dy = 0.0
    project_info%iDx = 0.0
    project_info%iDy = 0.0
    project_info%s_c = 0.0
    project_info%s_ic = 0.0
    project_info%n_polar = .false.
    project_info%s_polar = .false.
    project_info%north_pole = .false.
    project_info%polar = .false.

  end subroutine gmt_type_init

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! PROJECTION CONVERSION CODES
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Lambert equal area
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine gmt_lambeq_set(latc,lonc,project_info)

    !*FD Initialise a Lambert equal-area projection.

    implicit none

    ! Subroutine arguments ------------------------------------------------------

    real(dp),      intent(in)    :: latc          !*FD Location of projection centre, latitude (degrees).
    real(dp),      intent(in)    :: lonc          !*FD Location of projection centre, longitude (degrees).
    type(gmt_pinf),intent(inout) :: project_info  !*FD GMT parameters to be initialised.

    !----------------------------------------------------------------------------

    project_info%cosp=cos(latc*D2R)           ! cos (phi_p)
    project_info%sinp=sin(latc*D2R)           ! sin (phi_p)

    project_info%central_meridian=lonc        ! central meridian of projection
    project_info%pole=latc                    ! centre of map (latitude)

    ! These next four are not used in this incomplete translation to f90

    project_info%Dx=1.0                       ! Fudge factor for scaling the projection 
    project_info%Dy=1.0                       ! ditto
    project_info%iDx=1.0                      ! inverse of Dx
    project_info%iDy=1.0                      ! inverse of Dy

  end subroutine gmt_lambeq_set

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine gmt_lambeq(lon,lat,x,y,project_info)

    !*FD Converts lat-lon coordinates to x-y coordinates, using a
    !*FD Spherical Lambert Azimuthal Equal-Area projection.

    implicit none

    ! Subroutine arguments ------------------------------------------------------

    real(dp),      intent(in)  :: lon !*FD Longitude (degrees).
    real(dp),      intent(in)  :: lat !*FD Latitude (degrees).
    type(gmt_pinf),intent(in)  :: project_info !*FD GMT parameters to be used.
    real(dp),      intent(out) :: x   !*FD $x$ location (m).
    real(dp),      intent(out) :: y   !*FD $y$ location (m).

    ! Internal variables --------------------------------------------------------

    real(dp) :: k,tmp,sin_lat,cos_lat,sin_lon,cos_lon,c,dlon,dlat

    !----------------------------------------------------------------------------

    dlon = lon-project_info%central_meridian

    do while (dlon.lt.-180.0)
      dlon=dlon+360.0
    enddo

    do while (dlon.gt.180.0)
      dlon=dlon-360.0
    enddo

    dlon = dlon*D2R
    dlat = lat*D2R

    call sincos(dlat,sin_lat,cos_lat);
    call sincos(dlon,sin_lon,cos_lon);
    c = cos_lat * cos_lon

    tmp = 1.0 + project_info%sinp * sin_lat + project_info%cosp * c

    if (tmp > 0.0) then
      k = project_info%EQ_RAD * sqrt (2.0 / tmp)
      x = k * cos_lat * sin_lon
      y = k * (project_info%cosp * sin_lat - project_info%sinp * c)
      if (GMT_convert_latitudes) then  ! Gotta fudge abit 
        x = x*project_info%Dx
        y = y*project_info%Dy
      endif
    else
      x = -DBL_MAX ; y=-DBL_MAX
    endif    

  end subroutine gmt_lambeq

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine gmt_ilambeq(lon,lat,x,y,project_info)

    !*FD Convert x-y coordinates to lat-lon coordinates, using
    !*FD a Lambert Azimuthal Equal-Area projection.

    implicit none

    ! Subroutine arguments ------------------------------------------------------

    real(dp),      intent(in)  :: x    !*FD $x$-location (m).
    real(dp),      intent(in)  :: y    !*FD $y$-location (m).
    type(gmt_pinf),intent(in)  :: project_info !*FD GMT parameters to be used.
    real(dp),      intent(out) :: lon  !*FD Longitude (degrees).
    real(dp),      intent(out) :: lat  !*FD Latitude (degrees).

    ! ---------------------------------------------------------------------------

    real(dp) :: rho,c,sin_c,cos_c,xx,yy
  
    xx=x ; yy=y

    if (GMT_convert_latitudes) then  ! Undo effect of fudge factors
      xx = xx*project_info%iDx
      yy = yy*project_info%iDy
    endif

    rho=hypot (xx,yy)
  
    if (abs(rho) < GMT_CONV_LIMIT) then
      lat = project_info%pole
      lon = project_info%central_meridian
    else
       c = 2.0 * asin(0.5 * rho * project_info%i_EQ_RAD)
      call sincos (c, sin_c, cos_c)
      lat = asin (cos_c * project_info%sinp + (yy * sin_c * project_info%cosp / rho)) * R2D
      if (project_info%n_polar) then
        lon = project_info%central_meridian + R2D * atan2 (xx, -yy)
      else if (project_info%s_polar) then
        lon = project_info%central_meridian + R2D * atan2 (xx, yy)
      else
        lon = project_info%central_meridian + &
          R2D * atan2 (xx * sin_c, (rho * project_info%cosp * cos_c - yy * project_info%sinp * sin_c))
      endif
    endif

  end subroutine gmt_ilambeq

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Spherical polar / Polar Stereographic
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine gmt_map_init_stereo(latg,rlong,project_info)

    !*FD Initialises a polar stereographic projection.

    ! Subroutine arguments ------------------------------------------------------

    real(dp) :: latg      !*FD Centre of map projection, latitude (degrees).
    real(dp) :: rlong     !*FD Centre of map projection, longitude (degrees).
    type(gmt_pinf),intent(inout) :: project_info !*FD GMT parameters to be set.

    ! ---------------------------------------------------------------------------

    call gmt_set_polar (latg,project_info)

    ! Equatorial view has a problem with infinite loops.  Untill I find a cure
    !  we set projection center latitude to 0.001 so equatorial works for now 

    if (abs (latg) < SMALL) latg = 0.001

    call gmt_vstereo(rlong,latg,project_info)

  end subroutine gmt_map_init_stereo

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine gmt_set_polar (plat,project_info)

    !*FD Determines if the projection pole is N or S pole.

    implicit none

    ! Subroutine arguments ------------------------------------------------------

    real(dp),      intent(in)    :: plat         !*FD Latitude of projection centre (degrees). 
    type(gmt_pinf),intent(inout) :: project_info !*FD GMT parameters to be set.

    ! ---------------------------------------------------------------------------

    if (abs (abs (plat) - 90.0) < GMT_CONV_LIMIT) then
      project_info%polar = .true.
      project_info%north_pole  = (plat > 0.0)
       project_info%n_polar = project_info%north_pole
      project_info%s_polar = (.not.project_info%n_polar)
    endif

  end subroutine gmt_set_polar

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine gmt_vstereo (rlong0,plat,project_info)

    !*FD Set up a Stereographic transformation.

    implicit none 

    ! Subroutine arguments ------------------------------------------------------

    real(dp),      intent(in)    :: rlong0       !*FD Centre of the projection, longitude (degrees).
    real(dp),      intent(in)    :: plat         !*FD Centre of the projection, latitude (degrees).
    type(gmt_pinf),intent(inout) :: project_info !*FD GMT parameters to be set.

    ! Internal variables --------------------------------------------------------

    real(dp) :: clat

    ! ---------------------------------------------------------------------------

    clat = plat

    project_info%central_meridian = rlong0
    project_info%pole = plat          ! This is always geodetic
    project_info%sinp = sin (clat)   ! These may be conformal
    project_info%cosp = cos (clat)
    project_info%north_pole = (plat > 0.0)
    project_info%s_c = 2.0 * project_info%EQ_RAD * GMT_map_scale_factor
    project_info%s_ic = 1.0 / project_info%s_c

  end subroutine gmt_vstereo

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine gmt_plrs_sph(lon,lat,x,y,project_info)

  !*FD Convert lon-lat coordinates to x-y coordinates using a 
  !*FD spherical polar projection.

    implicit none

    ! Subroutine arguments ------------------------------------------------------

    real(dp),      intent(in)  :: lon          !*FD Longitude (degrees).
    real(dp),      intent(in)  :: lat          !*FD Latitude (degrees).
    real(dp),      intent(out) :: x            !*FD $x$-location (m).
    real(dp),      intent(out) :: y            !*FD $y$-location (m).
    type(gmt_pinf),intent(in)  :: project_info !*FD GMT parameters to be used.

    ! Internal variables --------------------------------------------------------

    real(dp) :: rho, slon, clon,dlon,dlat
      
    ! ---------------------------------------------------------------------------

    dlon = lon-project_info%central_meridian
    dlat = lat

    do 
    if (.not.(dlon < -180.0)) exit
    dlon = dlon+360.0
    enddo 

    do
    if (.not.(dlon > 180.0)) exit
    dlon = dlon-360.0
    enddo  

    dlon = dlon*D2R
    call sincos (lon,slon,clon)

     if (project_info%north_pole) then
      rho = project_info%s_c * tan (M_PI_4 - 0.5 * D2R * dlat)
      y = -rho * clon
      x =  rho * slon
    else
      rho = project_info%s_c * tan (M_PI_4 + 0.5 * D2R * dlat)
      y = rho * clon
       x = rho * slon
    endif

  end subroutine gmt_plrs_sph

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine gmt_iplrs_sph (lon,lat,x,y,project_info)

    !*FD Convert x-y coordinates to lat-lon coordinates using
    !*FD a spherical polar projection.

    implicit none

    ! Subroutine arguments ------------------------------------------------------

    real(dp),intent(out) :: lon !*FD Longitude (degrees).
    real(dp),intent(out) :: lat !*FD Latitude (degrees).
    real(dp),intent(in)  :: x   !*FD $x$-location (m).
    real(dp),intent(in)  :: y   !*FD $y$-location (m).
    type(gmt_pinf),intent(in) :: project_info !*FD GMT parameters to use.

    ! Internal variables --------------------------------------------------------

    real(dp) :: c
    real(dp) :: xx,yy

    ! ---------------------------------------------------------------------------

    xx=x ; yy=y
    
    if (x == 0.0.and.y == 0.0) then
      lon = project_info%central_meridian
      lat = project_info%pole
      return
    endif

    c = 2.0 * atan (hypot (xx, yy) * project_info%s_ic)

    if (project_info%north_pole) then
      lon = project_info%central_meridian + d_atan2 (xx, -yy) * R2D  
      lat = d_asin (cos (c)) * R2D
    else
      lon = project_info%central_meridian + d_atan2 (xx, yy) * R2D  
      lat = d_asin (-cos (c)) * R2D
    endif

  end subroutine gmt_iplrs_sph

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine gmt_stereo1_sph (lon,lat,x,y,project_info)

    !*FD Convert lon-lat coordinates to x-y coordinates using 
    !*FD a spherical stereographic projection, oblique view.
  
    implicit none

    ! Subroutine arguments ------------------------------------------------------

    real(dp),      intent(in)  :: lon !*FD Longitude (degrees).
    real(dp),      intent(in)  :: lat !*FD Latitude (degrees).
    real(dp),      intent(out) :: x   !*FD $x$-location (m).
    real(dp),      intent(out) :: y   !*FD $y$-location (m).
    type(gmt_pinf),intent(in)  :: project_info !*FD GMT parameters to use.

    ! Internal variables --------------------------------------------------------

    real(dp) :: dlon,sin_dlon,cos_dlon,s,c,cc,A,dlat

    ! ---------------------------------------------------------------------------

    dlon = D2R * (lon - project_info%central_meridian)
    dlat = lat*D2R
    call sincos (dlon,sin_dlon,cos_dlon)
    call sincos (dlat,s,c)
    cc=c*cos_dlon
    A = project_info%s_c / (1.0 + project_info%sinp * s + project_info%cosp * cc)
    x = A * c * sin_dlon
    y = A * (project_info%cosp * s - project_info%sinp * cc)

  end subroutine gmt_stereo1_sph

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine gmt_istereo_sph (lon,lat,x,y,project_info)

    !*FD Convert x-y coordinates to lat-lon coordinates using 
    !*FD a spherical stereographic projection.

    implicit none

    ! Subroutine arguments ------------------------------------------------------
 
    real(dp),      intent(out) :: lon !*FD Longitude (degrees).
    real(dp),      intent(out) :: lat !*FD Latitude (degrees).
    real(dp),      intent(in)  :: x   !*FD $x$-location (m).
    real(dp),      intent(in)  :: y   !*FD $y$-location (m).
    type(gmt_pinf),intent(in)  :: project_info !*FD GMT parameters to use.

    ! Internal variables --------------------------------------------------------

    real(dp) :: rho,c,sin_c,cos_c

    ! ---------------------------------------------------------------------------

    if (x == 0.0 .and. y == 0.0) then
      lon = project_info%central_meridian
      lat = project_info%pole
    else 
      rho = hypot(x,y)
      c = 2.0 * atan (rho * project_info%s_ic)
      call sincos (c, sin_c, cos_c)
      lat = d_asin (cos_c * project_info%sinp + (y * sin_c * project_info%cosp / rho)) * R2D
      lon = R2D * atan (x * sin_c / (rho * project_info%cosp * cos_c - y * project_info%sinp * sin_c)) &
             + project_info%central_meridian
    endif

  end subroutine gmt_istereo_sph

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine gmt_stereo2_sph (lon,lat,x,y,project_info)

    !*FD Convert lon-lat coordinates to x-y coordinates using a stereographic 
    !*FD projection, equatorial view.

    implicit none

    ! Subroutine arguments ------------------------------------------------------

    real(dp),      intent(in)  :: lon !*FD Longitude (degrees).
    real(dp),      intent(in)  :: lat !*FD Latitude (degrees).
    real(dp),      intent(out) :: x   !*FD $x$-location (m).
    real(dp),      intent(out) :: y   !*FD $y$-location (m).
    type(gmt_pinf),intent(in)  :: project_info !*FD GMT parameters to use.

    ! Internal variables --------------------------------------------------------

    real(dp) :: dlon,dlat,s,c,clon,slon,A

    ! ---------------------------------------------------------------------------
    
    dlon = lon - project_info%central_meridian
    if (abs (dlon - 180.0) < GMT_CONV_LIMIT) then
      x = 0.0
      y = 0.0
    else 
      dlon = dlon*D2R
      dlat = lat*D2R
      call sincos (dlat,s,c)
      call sincos (dlon,slon,clon)
      A = project_info%s_c / (1.0 + c * clon)
      x = A * c * slon
      y = A * s
    endif

  end subroutine gmt_stereo2_sph

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Utility routines
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine sincos(a,s,c)

    !*FD Calculates the sin and cos of an angle.

    implicit none

    real(dp),intent(in)  :: a !*FD Input value (radians).
    real(dp),intent(out) :: s !*FD sin(\texttt{a})
    real(dp),intent(out) :: c !*FD cos(\texttt{a})

      s = sin (a)
      c = cos (a)

  end subroutine sincos

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  real(dp) function hypot(x,y)

  !*FD This is an implementation of a standard C library function, 
  !*FD returning the value of $\sqrt{x^2+y^2}$.
  !*RV Returns the value of $\sqrt{x^2+y^2}$.

    implicit none

    real(dp),intent(in) :: x !*FD One input value
    real(dp),intent(in) :: y !*FD Another input value

    hypot=sqrt(x*x+y*y)

  end function hypot

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  real(dp) function d_atan2(x,y)

  !*FD A macro-implemented function from GMT
  !*RV If \texttt{x} and \texttt{0} are both zero, zero is returned, else
  !*FD \texttt{atan2 (y, x)} is returned.

    implicit none

    real(dp),intent(in) :: x,y !*FD Input value

    if (x==0.0 .and. y==0.0) then
      d_atan2=0.0
    else
      d_atan2=atan2 (y, x)
    endif

  end function d_atan2

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  real(dp) function d_asin(x)

  !*FD Macro-implemented function from GMT
  !*RV If $\mathtt{x}\geq 1$, \texttt{sign(M\_PI\_2,x)} is returned, else
  !*RV \texttt{asin(x)} is returned.

    implicit none

    real(dp),intent(in) :: x !*FD Input value

    if (abs(x) >= 1.0) then
      d_asin=sign(M_PI_2,x)
    else
      d_asin=asin(x)
    endif

  end function d_asin

end module gmt
