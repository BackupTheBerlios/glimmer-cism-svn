module plume_context

  ! all subroutine related to defining geographic context of the plume,
  ! such as topography and drag

  use plume_global


contains

  subroutine topog_depth_inflow_set(set_inflow)

    ! set depth and inflow data explicitly for simplified cases

    implicit none

    logical :: set_inflow ! if true, set the inflow region

    ! local variables

    integer :: i,k
    integer :: bnd,kdiag,idiag

    real(kind=kdp) :: toty,runtoty,runtotx,totx
    real(kind=kdp) :: cordis,slopdis,totdis,angle,iangle
    real(kind=kdp) :: ampdep,depy
    real(kind=kdp) :: c1,c2,rand

    ampdep = gldep - ifdep
    if ((bathtype.eq.3).or.(bathtype.eq.5)) then
       iangle = pi/4.d0
       idiag = int(rad*cos(iangle))
       kdiag = int(rad*sin(iangle))
    else
       iangle=0.0
       idiag=0
       kdiag=0
    end if

    ! ----------------
    ! inflow region(s)
    ! ----------------
    if (set_inflow) then
       if ((ninfmin.eq.1).and.(ninfmax.eq.1)) then

          ! set up automated single inflow
          intrin((infloa+1):(infloe-1),(knfloa+1):(knfloe-1)) = 1
          depinf((infloa+1):(infloe-1),(knfloa+1):(knfloe-1)) = depinffix

       else

          ! set up multiple inflows by hand
          do k = knfloa + 1,knfloe - 1
             do i = int(0440.d0*1000.d0/hx)+1,int(0490.d0*1000.d0/hx)-1
                intrin(i,k) = 1
                depinf(i,k) = depinffix
             end do
          end do

          do k = knfloa + 1,knfloe - 1
             do i = int(0510.d0*1000.d0/hx)+1,int(0560.d0*1000.d0/hx)-1
                intrin(i,k) = 2
                depinf(i,k) = depinffix
             end do
          end do

       end if
    end if

    ! ----------------
    ! slope topography
    ! ----------------

    if ((bathtype.ge.1).and.(bathtype.le.3)) then

       ! a) ice shelf has base which rises linearly from inflow
       ! ------------------------------------------------------
       ! find total domain size
       toty = 0.d0
       do k = 1,n_grid
          toty = toty + dy(k)
       end do

       ! (depth decreases from gldep to ifdep, 
       ! position increases from wcdep to wcdep + (gldep-ifdep))
       runtoty = 0.d0
       do k = 1,n_grid
          depy = wcdep + ampdep*(runtoty + dy(k)/2)/toty
          do i = 1,m_grid
             bpos(i,k) = depy
          end do
          runtoty = runtoty + dy(k)
       end do

    else

       ! b) ice shelf has base where isobaths are perpendicular to west wall
       ! -------------------------------------------------------------------

       ! 1)for 90-degree corner (slope of base is preserved, not total 
       !   elevation change)

       if (bathtype.eq.4) then

          ! find total length of wall (to calculate gradient of shelf base)
          cordis = pi*dble(rad)/2.d0

          ! do straight section next to inflow
          do k = 1,kcorn
             depy = wcdep + ampdep*dble(k)/dble(n_grid)
             do i = 1,m_grid
                bpos(i,k) = depy
             end do
          end do

          do k = kcorn + 1,n_grid

             ! do curved section in north-east corner 
             ! (find angle and calculate height)
             do i = infloa - rad + 1,m_grid
                angle = datan(dble(k - kcorn)/dble(i - (infloa - rad)))
                depy = wcdep + ampdep*dble(kcorn)/dble(n_grid) &
                     + ampdep*(2.d0*angle/pi)*cordis/dble(n_grid)
                bpos(i,k) = depy
             end do

             ! do straight section in north-west corner
             do i = 1,infloa - rad
                depy = wcdep + ampdep*(dble(kcorn) + cordis)/dble(n_grid) &
                     + ampdep*dble(infloa - rad - i)/dble(n_grid)
                bpos(i,k) = depy
             end do
          end do

       end if

       ! 2)for 45-degree corner (total elevation change is preserved, not slope)

       if (bathtype.eq.5) then

          ! find total length of wall (to calculate gradient of shelf base)
          ! find length of arc
          cordis = iangle*dble(rad)

          ! find length of diagonal section (minimum because intersects 
          ! one wall or other)
          slopdis = dsqrt(2.d0) &
               *dble(min(n_grid - kcorn + kdiag,infloa - rad + idiag))

          totdis = dble(kcorn) + cordis + slopdis

          ! do straight section next to inflow
          do k = 1,kcorn
             depy = wcdep + ampdep*k/totdis
             do i = 1,m_grid
                bpos(i,k) = depy
             end do
          end do

          do k = kcorn + 1,n_grid

             ! do curved section in north-east triangular section and 
             ! isobaths perpendicular to 45 deg line
             do i = infloa - rad + 1,m_grid
                angle = datan(dble(k - kcorn)/dble(i - (infloa - rad)))
                if (angle.le.iangle) then
                   depy = wcdep + ampdep*kcorn/totdis  &
                        + ampdep*(angle/iangle)*cordis/totdis
                   bpos(i,k) = depy
                else
                   depy = wcdep + ampdep*(kcorn + cordis)/totdis &
                        + (((k - kcorn) - (i - (infloa - rad)))/ &
                        (2.d0*dble(min(n_grid - kcorn - idiag,infloa &
                        - rad + idiag)))) &
                        *ampdep*slopdis/totdis
                   bpos(i,k) = depy
                end if
             end do

             do i = 1,infloa - rad
                depy = wcdep + ampdep*(kcorn + cordis)/totdis &
                     + (((k - kcorn) - (i - (infloa - rad)))/ &
                     (2.d0*dble(min(n_grid - kcorn - idiag,infloa - rad &
                     + idiag)))) &
                     *ampdep*slopdis/totdis
                bpos(i,k) = depy
             end do
          end do

       end if

       ! c) simplified fris bathymetries
       ! -------------------------------

       ! 1) carlson inlet 
       ! (90-degree corner, isobaths perpendicular to west wall and flat 
       ! bottom when round corner)

       if (bathtype.eq.6) then

          ! find total length of slope-adjacent wall (to calculate gradient
          ! of sloping section)
          cordis = pi*dble(rad)/2.d0
          totdis = dble(kcorn) + cordis

          ! do straight section next to inflow
          do k = 1,kcorn
             depy = wcdep + ampdep*k/totdis
             do i = 1,m_grid
                bpos(i,k) = depy
             end do
          end do

          do k = kcorn + 1,n_grid

             ! do curved section in north-east corner (find angle and 
             ! calculate height)
             do i = infloa - rad + 1,m_grid
                angle = datan(dble(k - kcorn)/dble(i - (infloa - rad)))
                depy = wcdep + ampdep*kcorn/totdis  &
                     + ampdep*(2.d0*angle/pi)*cordis/totdis
                bpos(i,k) = depy
             end do

             ! do flat section in north-west corner
             do i = 1,infloa - rad
                depy = gldep + wcdep - ifdep
                bpos(i,k) = depy
             end do
          end do

       end if

       ! 2) rutford ice stream
       ! (90-degree corner, isobaths perpendicular to west wall and 
       ! variable slope)

       if (bathtype.eq.7) then

          ! section 1 - inflow wall and corner

          ! elevation and length
          ampdep = 800.d0
          cordis = pi*dble(rad)/2.d0
          totdis = dble(kcorn) + cordis

          ! straight section next to inflow
          do k = 1,kcorn
             depy = wcdep + ampdep*k/totdis
             do i = 1,m_grid
                bpos(i,k) = depy
             end do
          end do

          ! curved section in north-east corner (find angle and calculate
          ! height)
          do k = kcorn + 1,n_grid
             do i = infloa - rad + 1,m_grid
                angle = datan(dble(k - kcorn)/dble(i - (infloa - rad)))
                depy = wcdep + ampdep*kcorn/totdis  &
                     + ampdep*(2.d0*angle/pi)*cordis/totdis
                bpos(i,k) = depy
             end do
          end do

          ! section 2 - section crossing carlson inlet

          do k = kcorn + 1,n_grid
             do i = infloa - rad - 110,infloa - rad
                depy = wcdep + 800.d0  &
                     + 100.d0*dble(infloa - rad - i)/110.d0
                bpos(i,k) = depy
             end do
          end do

          ! section 3 - southern corner of fowler peninsula

          do k = kcorn + 1,n_grid
             do i = infloa - rad - 140,infloa - rad - 109
                depy = wcdep + 900.d0 &
                     + 300.d0*dble(infloa - rad - 109 - i)/30.d0
                bpos(i,k) = depy
             end do
          end do

          ! section 4 -  flat section in northern corner of fowler peninsula
          do k = kcorn + 1,n_grid
             do i = 1,infloa - rad - 141
                depy = gldep + wcdep - ifdep
                bpos(i,k) = depy
             end do
          end do

       end if

       ! 3) carlson inlet test cases
       ! (height increases northwards and then plateaus, with or without wall)

       if ((bathtype.eq.8).or.(bathtype.eq.9)) then

          ! do sloping section in south
          do k = 1,int(190.d0*1000.d0/hx)
             depy = wcdep + (gldep - ifdep)*dble(k)/(190.d0*1000.d0/hx)
             do i = 1,m_grid
                bpos(i,k) = depy
             end do
          end do

          ! do flat section in north
          do k = int(190.d0*1000.d0/hx) + 1,n_grid
             do i = 1,m_grid
                depy = gldep + wcdep - ifdep
                bpos(i,k) = depy
             end do
          end do

       end if

       ! 4) evans ice stream
       ! (90-degree corner, isobaths perpendicular to west wall,
       ! then flat and dead wall when round corner)

       if (bathtype.eq.10) then

          ! find total length of slope-adjacent wall (to calculate 
          ! gradient of sloping section)
          cordis = pi*dble(rad)/2.d0
          totdis = dble(kcorn) + cordis

          ! do straight section next to inflow
          do k = 1,kcorn
             depy = wcdep + ampdep*k/totdis
             do i = 1,m_grid
                bpos(i,k) = depy
             end do
          end do

          do k = kcorn + 1,n_grid

             ! do curved section in north-east corner (find angle and 
             ! calculate height)
             do i = infloa - rad + 1,m_grid
                angle = datan(dble(k - kcorn)/dble(i - (infloa - rad)))
                depy = wcdep + ampdep*kcorn/totdis  &
                     + ampdep*(2.d0*angle/pi)*cordis/totdis
                bpos(i,k) = depy
             end do

             ! do flat section in north-west corner
             do i = 1,infloa - rad
                depy = gldep + wcdep - ifdep
                bpos(i,k) = depy
             end do
          end do

       end if

       ! 5) whole simplified fris
       ! (simple wedge shelf with walls on both sides)

       if (bathtype.eq.11) then
          ! find total domain size
          toty = 0.d0
          do k = 1,n_grid
             toty = toty + dy(k)
          end do

          ! (depth decreases from gldep to ifdep, 
          ! position increases from wcdep to wcdep + (gldep-ifdep))
          runtoty = 0.d0
          do k = 1,n_grid
             depy = wcdep + ampdep*(runtoty + dy(k)/2)/toty
             do i = 1,m_grid
                bpos(i,k) = depy
             end do
             runtoty = runtoty + dy(k)
          end do

       end if

    end if

    if (bathtype.eq.12) then

       ! d) ice shelf has typical analytical base
       ! ----------------------------------------

       ! find total domain size
       toty = 0.d0
       do k = 1,n_grid
          toty = toty + dy(k)
       end do

       ! (depth decreases from gldep to ifdep, 
       ! position increases from wcdep to wcdep + (gldep-ifdep))
       c2 = toty*ifdep**2/(gldep**2 - ifdep**2)
       c1 = gldep*dsqrt(c2)
       runtoty = 0.d0
       !         
       do k = 1,n_grid
          depy = wcdep + gldep  &
               - c1 / dsqrt(c2 + runtoty + dy(k)/2)
          do i = 1,m_grid
             bpos(i,k) = depy
          end do
          runtoty = runtoty + dy(k)
       end do

    end if

    if (bathtype.eq.13) then


       ! Peterman glacier style longitudinal channels

       ! find total domain size
       toty = 0.d0
       do k = 1,n_grid
          toty = toty + dy(k)
       end do

       totx = 0.d0
       do i = 1,m_grid
          totx = totx + dx(i)
       end do

       ! (depth decreases from gldep to ifdep, 
       ! basal position increases from wcdep to wcdep + (gldep-ifdep)) 
       ! with sinusoidal undulations superimposed

       c0 = wcdep + gldep - ifdep
       c1 = - ampdep / toty

       runtoty = 0.d0

       do k = 1,n_grid

          runtoty = runtoty + dy(k)
          runtotx = 0.d0

          do i = 1,m_grid

             runtotx = runtotx + dx(i)
             call random_number(rand)
             bpos(i,k) = c0 + c1*runtoty + &
                  (0.5)*channel_amplitude *(1-(runtoty/toty))** along_slope_deepening_exp &
                  *(1-cos(2*pi*cross_slope_wavenumber * (runtotx/totx))) &
                  + (0.5d0 - rand)*2.d0 *random_amplitude

          end do
       end do

    end if

    ! -------------------
    ! dry land topography
    ! -------------------

    if (set_inflow) then

       ! in this mode it is assumed that land appears 
       ! to the east and west of inflow areas

       if ((ninfmin.eq.1).and.(ninfmax.eq.1)) then


          ! set automated single inflow cells topography
          do k = knfloa + 1,knfloe - 1
             do i = 1,infloa  
                bpos(i,k) = 0.d0
             end do
             do i = infloe,m_grid
                bpos(i,k) = 0.d0
             end do
          end do

       else

          ! set topography around multiple inflows (by hand)
          do k = knfloa + 1,knfloe - 1
             do i = 1,int(0440.d0*1000.d0/hx)
                bpos(i,k) = 0.d0
             end do
             do i = int(0560.d0*1000.d0/hx),m_grid
                bpos(i,k) = 0.d0
             end do
          end do

          do k = knfloa + 1,knfloe - 1
             do i = int(0490.d0*1000.d0/hx),int(0510.d0*1000.d0/hx)
                bpos(i,k) = 0.d0
             end do
          end do

       end if
    end if

    ! set topography of block west of inflow
    if (((bathtype.ge.2).and.(bathtype.le.7)).or.(bathtype.eq.10))then
       do k = 1,kcorn
          do i = 1,infloa
             bpos(i,k) = 0.d0
          end do
       end do

       ! set topography of rounded 90-degree corner
       do k = kcorn + 1,kcorn + rad
          do i = 1,infloa 
             bnd = (infloa - rad)  &
                  + int(dsqrt(dble(rad*rad - (k - kcorn)**2)))
             if (i.le.bnd) bpos(i,k) = 0.d0
          end do
       end do
    end if

    ! additionally set topography of line for 45-degree corner
    if ((bathtype.eq.3).or.(bathtype.eq.5)) then
       do k = kcorn + kdiag,n_grid
          do i = 1,infloa
             bnd = infloa - rad + idiag - (k - (kcorn + kdiag))
             if (i.le.bnd) bpos(i,k) = 0.d0
          end do
       end do
    end if

    ! set topography of wall for rotating test case
    if (bathtype.eq.9) then
       do k = 1,n_grid
          do i = 1,infloa
             bpos(i,k) = 0.d0
          end do
       end do
    end if

    ! set topography of wall after corner for evans case
    if (bathtype.eq.10) then
       do k = 1,n_grid
          do i = 1,int(1000.d0/hx)
             bpos(i,k) = 0.d0
          end do
       end do
    end if

    ! set topography of both walls for simple fris case
    if (bathtype.eq.11) then
       do k = 1,n_grid
          bpos(1,k) = 0.d0
          bpos(2,k) = 0.d0
          bpos(m_grid-1,k) = 0.d0
          bpos(m_grid,k) = 0.d0
       end do
    end if

    if (bathtype.eq.13) then
       ! set land on east, south and west walls
       bpos(1:2,:) = 0.0
       !bpos(:,1:2) = 0.0
       bpos((m_grid-1):m_grid,:) = 0.0
    end if

  end subroutine topog_depth_inflow_set

  subroutine topog_read_edit()

    ! read depth data from file if using real topography, interpolate, and
    ! hand-edit as required

    implicit none

    ! local variables

    integer i,k,nxin,nyin,ilo,ihi,klo,khi,imin,imax,kmin,kmax

    real(kind=kdp) :: hxin = 0.0, hyin = 0.0
    real(kind=kdp) xfac,yfac
    real(kind=kdp) indraft(lxin,lyin)
    real(kind=kdp) draft(m_grid,n_grid)
    real(kind=kdp) cutdepth

    character(len=11) datapath


    ! ---------------------------------------------
    ! set input data path relative to run directory
    ! ---------------------------------------------

    datapath='../../input'

    ! ---------------------
    ! read slope topography
    ! ---------------------
    ! data must be in the form 0.0 = open ocean, 9999.0 = grounded ice and
    ! other value = ice shelf draft (metres)

    ! ++++++++++++++++++++++++++++++++++++
    ! 1) fris (sandhaeger masked with add)
    ! ++++++++++++++++++++++++++++++++++++

    if (context.eq."fris") then

       nxin = 500
       nyin = 500
       hxin = 2000.d0
       hyin = 2000.d0

       open(22,file=datapath//'/frisdraft_sand_add.dat')

       do k = 1,nyin
          read(22,*) ( indraft(i,k) , i = 1,nxin )
       end do

       close(22)

    end if

    ! +++++++++++++++++++++++++
    ! 2) isomip (70-80s, 0-15e)
    ! +++++++++++++++++++++++++

    if (context.eq."isomip") then

       nxin = 600
       nyin = 600
       hxin = 2000.d0
       hyin = 2000.d0

       open(22,file=datapath//'/isomipdraft.dat')

       do k = 1,nyin
          read(22,*) ( indraft(i,k) , i = 1,nxin )
       end do

       close(22)

    end if

    ! +++++++++++++++++++++++++++++++++++
    ! 2) larsen (airborne data over RAMP)
    ! +++++++++++++++++++++++++++++++++++

    if (context.eq."larsen") then

       nxin = 220 
       nyin = 488
       hxin = 1000.d0
       hyin = 1000.d0

       open(22,file=datapath//'/larsendraft.dat')

       do k = 1,nyin
          read(22,*) ( indraft(i,k) , i = 1,nxin )
       end do

       close(22)

    end if

    ! ----------------------------------
    ! interpolate linearly to model grid
    ! ----------------------------------

    do i = 1,m_grid

       ilo = int(hx*dble(i-1)/hxin) + 1
       ihi = ilo + 1
       xfac = (hx*dble(i-1) - hxin*dble(ilo-1))/hxin

       do k = 1,n_grid

          klo = int(hy*dble(k-1)/hyin) + 1
          khi = klo + 1
          yfac = (hy*dble(k-1) - hyin*dble(klo-1))/hyin

          draft(i,k) = (1.d0-yfac)*( &
               (1.d0-xfac)*indraft(ilo,klo) + xfac*indraft(ihi,klo)) &
               + yfac*((1.d0-xfac)*indraft(ilo,khi) + xfac*indraft(ihi,khi))

       end do
    end do

    ! --------------------------------------
    ! set basal topography where appropriate
    ! --------------------------------------

    do i = 1,m_grid
       do k = 1,n_grid
          if (draft(i,k).le.9998.0) &
               bpos(i,k) = gldep + wcdep - draft(i,k)
       end do
    end do

    ! --------------------------------
    ! hand-edit topography if required
    ! --------------------------------

    if (topedit) then

       ! +++++++
       ! 1) fris
       ! +++++++

       if (context.eq."fris") then

          ! fill support force glacier near grounding line
          cutdepth = 1200.d0
          imin = int(0880.d0*1000.d0/hx)
          imax = int(0920.d0*1000.d0/hx)
          kmin = int(0430.d0*1000.d0/hy)
          kmax = int(0465.d0*1000.d0/hy)

          do i = imin,imax
             do k = kmin,kmax
                if (bpos(i,k).ge.(gldep+wcdep-cutdepth)) then
                   bpos(i,k) = gldep+wcdep-cutdepth 
                end if
             end do
          end do

       end if

       ! +++++++++
       ! 2) isomip
       ! +++++++++

       if (context.eq."isomip") then
       end if

       ! +++++++++
       ! 3) larsen
       ! +++++++++

       if (context.eq."larsen") then
       end if

    end if

  end subroutine topog_read_edit

  subroutine topog_smooth()

    ! smooth topography, both in general and, if required, more in certain areas

    implicit none

    ! local variables

    integer :: i,k,l,imin,imax,kmin,kmax,icen,kcen
    integer,dimension(m_grid,n_grid) :: ii,kk

    real(kind=kdp) tmp(m_grid,n_grid),gmsk(m_grid,n_grid),vmsk(m_grid,n_grid)
    real(kind=kdp) weightsum,valuesum,rotfac


    ! --------------------
    ! set up land/sea mask
    ! --------------------

    where((bpos.gt.0.d0).and.(bpos.lt.gldep+wcdep))
       gmsk = 1.d0
    elsewhere
       gmsk = 0.d0
    end where

    ! --------------------------------------------
    ! smooth all topography (background smoothing)
    ! --------------------------------------------

    if (bsmoothit.gt.0) then

       do l = 1,bsmoothit 

          do i = 1,m_grid
             do k = 1,n_grid
                tmp(i,k) = bpos(i,k)
             end do
          end do

          do i = 1,m_grid
             do k = 1,n_grid
                if (gmsk(i,k).eq.1.d0) then

                   weightsum = cweight*gmsk(i,k) + nweight*( &
                        gmsk(i-1,k-1) + gmsk(i-1,k) + gmsk(i-1,k+1) + &
                        gmsk(i  ,k-1) + gmsk(i,k+1) + &
                        gmsk(i+1,k-1) + gmsk(i+1,k) + gmsk(i+1,k+1) )

                   valuesum = cweight*tmp(i,k)*gmsk(i,k) + nweight*( &
                        tmp(i-1,k-1)*gmsk(i-1,k-1)+tmp(i-1,k)*gmsk(i-1,k) + &
                        tmp(i-1,k+1)*gmsk(i-1,k+1)+tmp(i,k-1)*gmsk(i,k-1) + &
                        tmp(i,k+1)*gmsk(i,k+1)+tmp(i+1,k-1)*gmsk(i+1,k-1) + &
                        tmp(i+1,k)*gmsk(i+1,k)+tmp(i+1,k+1)*gmsk(i+1,k+1))

                   bpos(i,k) = valuesum/weightsum

                end if
             end do
          end do

       end do

    end if

    ! --------------------------------
    ! additionally smooth chosen areas
    ! --------------------------------
    ! unless otherwise indicated the smoothing regions are defined such
    ! that all
    ! cells within a rotated box (defined by limits and then rotated
    ! about centre 
    ! according to rotfac) are additionally smoothed according to the 
    ! use's choice

    ! +++++++
    ! 1) fris
    ! +++++++

    if (context.eq."fris") then

       ! rutford ice stream grounding line region

       if ((smflag(01)).and.(smoothit(01).gt.bsmoothit)) then

          ! set up mask defining area to be additionally smoothed

          rotfac = +000.0d0*pi/180.d0
          imin = int(0225.d0*1000.d0/hx)
          imax = int(0400.d0*1000.d0/hx)
          kmin = int(0000.d0*1000.d0/hy)
          kmax = int(0110.d0*1000.d0/hy)

          icen = imin + nint(dble(imax-imin)/2.d0)
          kcen = kmin + nint(dble(kmax-kmin)/2.d0)

          ii = get_ii(m_grid,n_grid,icen,kcen,rotfac)
          kk = get_kk(m_grid,n_grid,icen,kcen,rotfac)

          where((((ii.ge.imin).and.(ii.le.imax)).and. &
               ((kk.ge.kmin).and.(kk.le.kmax))).and. &
               (gmsk(i,k).eq.1.d0) )
             vmsk = 1.d0
          elsewhere
             vmsk = 0.d0
          end where


          ! perform additional smoothing

          do l = 1,smoothit(01) - bsmoothit

             do i = 1,m_grid
                do k = 1,n_grid
                   tmp(i,k) = bpos(i,k)
                end do
             end do

             do i = 1,m_grid
                do k = 1,n_grid
                   if (vmsk(i,k).eq.1.d0) then

                      weightsum = cweight*vmsk(i,k) + nweight*( &
                           vmsk(i-1,k-1) + vmsk(i-1,k) + vmsk(i-1,k+1) + &
                           vmsk(i  ,k-1) + vmsk(i,k+1) + &
                           vmsk(i+1,k-1) + vmsk(i+1,k) + vmsk(i+1,k+1) )

                      valuesum = cweight*tmp(i,k)*vmsk(i,k) + nweight*( &
                           tmp(i-1,k-1)*vmsk(i-1,k-1)+tmp(i-1,k  )*vmsk(i-1,k  ) + &
                           tmp(i-1,k+1)*vmsk(i-1,k+1)+tmp(i  ,k-1)*vmsk(i  ,k-1) + &
                           tmp(i  ,k+1)*vmsk(i  ,k+1)+tmp(i+1,k-1)*vmsk(i+1,k-1) + &
                           tmp(i+1,k  )*vmsk(i+1,k  )+tmp(i+1,k+1)*vmsk(i+1,k+1) )

                      bpos(i,k) = valuesum/weightsum

                   end if
                end do
             end do

          end do

       end if

       ! support force ice stream grounding line region

       if ((smflag(02)).and.(smoothit(02).gt.bsmoothit)) then

          ! set up mask defining area to be additionally smoothed

          rotfac = +045.0d0*pi/180.d0
          imin = int(0890.d0*1000.d0/hx)
          imax = int(0950.d0*1000.d0/hx)
          kmin = int(0435.d0*1000.d0/hy)
          kmax = int(0505.d0*1000.d0/hy)

          icen = imin + nint(dble(imax-imin)/2.d0)
          kcen = kmin + nint(dble(kmax-kmin)/2.d0)

          ii = get_ii(m_grid,n_grid,icen,kcen,rotfac)
          kk = get_kk(m_grid,n_grid,icen,kcen,rotfac)

          where ((((ii.ge.imin).and.(ii.le.imax)).and. &
               ((kk.ge.kmin).and.(kk.le.kmax))).and. &
               (gmsk.eq. 1.d0)              ) 
             vmsk = 1.d0
          elsewhere
             vmsk = 0.d0
          end where

          ! perform additional smoothing

          do l = 1,smoothit(02) - bsmoothit

             do i = 1,m_grid
                do k = 1,n_grid
                   tmp(i,k) = bpos(i,k)
                end do
             end do

             do i = 1,m_grid
                do k = 1,n_grid
                   if (vmsk(i,k).eq.1.d0) then

                      weightsum = cweight*vmsk(i,k) + nweight*( &
                           vmsk(i-1,k-1) + vmsk(i-1,k) + vmsk(i-1,k+1) + &
                           vmsk(i  ,k-1) + vmsk(i,k+1) + &
                           vmsk(i+1,k-1) + vmsk(i+1,k) + vmsk(i+1,k+1) )

                      valuesum = cweight*tmp(i,k)*vmsk(i,k) + nweight*( &
                           tmp(i-1,k-1)*vmsk(i-1,k-1)+tmp(i-1,k  )*vmsk(i-1,k  ) + &
                           tmp(i-1,k+1)*vmsk(i-1,k+1)+tmp(i  ,k-1)*vmsk(i  ,k-1) + &
                           tmp(i  ,k+1)*vmsk(i  ,k+1)+tmp(i+1,k-1)*vmsk(i+1,k-1) + &
                           tmp(i+1,k  )*vmsk(i+1,k  )+tmp(i+1,k+1)*vmsk(i+1,k+1) )

                      bpos(i,k) = valuesum/weightsum

                   end if
                end do
             end do

          end do

       end if

       ! spare slot

       if (smflag(03)) then

          write(*,*) 'error: smoothing zone not defined'

       end if

       ! end fris

    end if

    ! +++++++++
    ! 2) isomip
    ! +++++++++

    if (context.eq."isomip") then

       if ((smflag(01).or.smflag(02)).or.smflag(03)) then

          write(*,*) 'error: smoothing zone not defined'

       end if

    end if

    ! +++++++++
    ! 3) larsen
    ! +++++++++

    if (context.eq."larsen") then

       ! choyce point region

       if ((smflag(01)).and.(smoothit(01).gt.bsmoothit)) then

          ! set up mask defining area to be additionally smoothed

          rotfac = +000.0d0*pi/180.d0
          imin = int(0010.d0*1000.d0/hx)
          imax = int(0022.d0*1000.d0/hx)
          kmin = int(0110.d0*1000.d0/hy)
          kmax = int(0126.d0*1000.d0/hy)

          icen = imin + nint(dble(imax-imin)/2.d0)
          kcen = kmin + nint(dble(kmax-kmin)/2.d0)

          ii = get_ii(m_grid,n_grid,icen,kcen,rotfac)
          kk = get_kk(m_grid,n_grid,icen,kcen,rotfac)

          where ((((ii.ge.imin).and.(ii.le.imax)).and. &
               ((kk.ge.kmin).and.(kk.le.kmax))).and. &
               (gmsk.eq.1.d0)              ) 
             vmsk = 1.d0
          elsewhere
             vmsk = 0.d0
          end where

          ! perform additional smoothing

          do l = 1,smoothit(01) - bsmoothit

             do i = 1,m_grid
                do k = 1,n_grid
                   tmp(i,k) = bpos(i,k)
                end do
             end do

             do i = 1,m_grid
                do k = 1,n_grid
                   if (vmsk(i,k).eq.1.d0) then

                      weightsum = cweight*vmsk(i,k) + nweight*( &
                           vmsk(i-1,k-1) + vmsk(i-1,k) + vmsk(i-1,k+1) + &
                           vmsk(i  ,k-1) + vmsk(i,k+1) + &
                           vmsk(i+1,k-1) + vmsk(i+1,k) + vmsk(i+1,k+1) ) 

                      valuesum = cweight*tmp(i,k)*vmsk(i,k) + nweight*( &
                           tmp(i-1,k-1)*vmsk(i-1,k-1)+tmp(i-1,k  )*vmsk(i-1,k  ) + &
                           tmp(i-1,k+1)*vmsk(i-1,k+1)+tmp(i  ,k-1)*vmsk(i  ,k-1) + &
                           tmp(i  ,k+1)*vmsk(i  ,k+1)+tmp(i+1,k-1)*vmsk(i+1,k-1) + &
                           tmp(i+1,k  )*vmsk(i+1,k  )+tmp(i+1,k+1)*vmsk(i+1,k+1) )

                      bpos(i,k) = valuesum/weightsum

                   end if
                end do
             end do

          end do

       end if

       ! thuronyi bluff region

       if ((smflag(02)).and.(smoothit(02).gt.bsmoothit)) then

          ! set up mask defining area to be additionally smoothed

          rotfac = +000.0d0*pi/180.d0
          imin = int(0026.d0*1000.d0/hx)
          imax = int(0050.d0*1000.d0/hx)
          kmin = int(0213.d0*1000.d0/hy)
          kmax = int(0236.d0*1000.d0/hy)

          icen = imin + nint(dble(imax-imin)/2.d0)
          kcen = kmin + nint(dble(kmax-kmin)/2.d0)

          ii = get_ii(m_grid,n_grid,icen,kcen,rotfac)
          kk = get_kk(m_grid,n_grid,icen,kcen,rotfac)

          where ((((ii.ge.imin).and.(ii.le.imax)).and. &
               ((kk.ge.kmin).and.(kk.le.kmax))).and. &
               (gmsk.eq.1.d0)              ) 
             vmsk = 1.d0
          elsewhere
             vmsk = 0.d0
          end where

          ! perform additional smoothing

          do l = 1,smoothit(02) - bsmoothit

             tmp = bpos

             do i = 1,m_grid
                do k = 1,n_grid
                   if (vmsk(i,k).eq.1.d0) then

                      weightsum = cweight*vmsk(i,k) + nweight*( &
                           vmsk(i-1,k-1) + vmsk(i-1,k) + vmsk(i-1,k+1) + &
                           vmsk(i  ,k-1) + vmsk(i,k+1) + &
                           vmsk(i+1,k-1) + vmsk(i+1,k) + vmsk(i+1,k+1) )

                      valuesum = cweight*tmp(i,k)*vmsk(i,k) + nweight*( &
                           tmp(i-1,k-1)*vmsk(i-1,k-1)+tmp(i-1,k  )*vmsk(i-1,k  ) + &
                           tmp(i-1,k+1)*vmsk(i-1,k+1)+tmp(i  ,k-1)*vmsk(i  ,k-1) + &
                           tmp(i  ,k+1)*vmsk(i  ,k+1)+tmp(i+1,k-1)*vmsk(i+1,k-1) + &
                           tmp(i+1,k  )*vmsk(i+1,k  )+tmp(i+1,k+1)*vmsk(i+1,k+1) )

                      bpos(i,k) = valuesum/weightsum

                   end if
                end do
             end do

          end do

       end if

       ! spare slot

       if (smflag(03)) then

          write(*,*) 'error: smoothing zone not defined'

       end if

       ! end larsen

    end if

  end subroutine topog_smooth

  subroutine inflow_set_isomip()

    implicit none

    ! set isomip inflow depths by hand if using real topography

    ! local variables

    integer i,k

    ! -----------------------------
    ! set up inflow regions by hand
    ! -----------------------------
    ! unless otherwise indicated the inflows are set up by defining a square
    ! region and then setting all ocean domain within that region to inflow if it
    ! is deeper than a given depth.

    ! the inflows are:
    !  01 - whole southern boundary

    ! whole southern boundary
    ! -----------------------

    if (inflag(01)) then

       do i = m_grid/2-5,m_grid/2+5
          do k = n_grid/2-5,n_grid/2+5
             intrin(i,k) = 1
             depinf(i,k) = depinffix
          end do
       end do

    end if

    ! ----------------------------------------------
    ! mask regions with grounding line and ice front
    ! ----------------------------------------------

    do i = 1,m_grid
       do k = 1,n_grid
          if ((bpos(i,k).le.0.d0).or.(bpos(i,k).ge.(gldep+wcdep))) then
             intrin(i,k) = 0
             depinf(i,k) = 0.d0
          end if
       end do
    end do

  end subroutine inflow_set_isomip

  subroutine inflow_set_fris()

    implicit none

    ! set fris inflow depths by hand if using real topography

    ! local variables

    integer i,k,imin,imax,kmin,kmax

    real(kind=kdp) cutdepth


    ! -----------------------------
    ! set up inflow regions by hand
    ! -----------------------------
    ! unless otherwise indicated the inflows are set up by defining a square
    ! region and then setting all ocean domain within that region to inflow if it
    ! is deeper than a given depth.

    ! the inflows are:
    !  01 - evans ice stream
    !  02 - carlson inlet
    !  03 - rutford ice stream
    !  04 - institute ice stream
    !  05 - mollereisstrom
    !  06 - foundation ice stream
    !  07 - support force glacier
    !  08 - recovery glacier
    !  09 - slessor glacier

    ! evans ice stream
    ! ----------------

    if (inflag(01)) then

       cutdepth = 0900.d0
       imin = int(0050.d0*1000.d0/hx)
       imax = int(0150.d0*1000.d0/hx) 
       kmin = int(0200.d0*1000.d0/hy)
       kmax = int(0300.d0*1000.d0/hy)

       do i = imin,imax
          do k = kmin,kmax
             if (bpos(i,k).le.(gldep+wcdep-cutdepth)) then
                intrin(i,k) = 1
                depinf(i,k) = depinffix
             end if
          end do
       end do

    end if

    ! carlson inlet
    ! -------------

    if (inflag(02)) then

       cutdepth = 1400.d0
       imin = int(0150.d0*1000.d0/hx)
       imax = int(0210.d0*1000.d0/hx)
       kmin = int(0050.d0*1000.d0/hy)
       kmax = int(0150.d0*1000.d0/hy)

       do i = imin,imax
          do k = kmin,kmax
             if (bpos(i,k).le.(gldep+wcdep-cutdepth)) then
                intrin(i,k) = 2
                depinf(i,k) = depinffix
             end if
          end do
       end do

    end if

    ! rutford ice stream
    ! ------------------

    if (inflag(03)) then

       cutdepth = 1500.d0
       imin = int(0200.d0*1000.d0/hx)
       imax = int(0300.d0*1000.d0/hx)
       kmin = int(0000.d0*1000.d0/hy)
       kmax = int(0070.d0*1000.d0/hy)

       do i = imin,imax
          do k = kmin,kmax
             if (bpos(i,k).le.(gldep+wcdep-cutdepth)) then
                intrin(i,k) = 3
                depinf(i,k) = depinffix
             end if
          end do
       end do

    end if

    ! institute ice stream
    ! --------------------

    if (inflag(04)) then

       cutdepth = 1000.d0
       imin = int(0450.d0*1000.d0/hx)
       imax = int(0600.d0*1000.d0/hx)
       kmin = int(0100.d0*1000.d0/hy)
       kmax = int(0200.d0*1000.d0/hy)

       do i = imin,imax
          do k = kmin,kmax
             if (bpos(i,k).le.(gldep+wcdep-cutdepth)) then
                intrin(i,k) = 4
                depinf(i,k) = depinffix
             end if
          end do
       end do

    end if

    ! mollereisstrom
    ! --------------

    if (inflag(05)) then

       cutdepth = 1050.d0
       imin = int(0650.d0*1000.d0/hx)
       imax = int(0750.d0*1000.d0/hx)
       kmin = int(0200.d0*1000.d0/hy)
       kmax = int(0300.d0*1000.d0/hy)

       do i = imin,imax
          do k = kmin,kmax
             if (bpos(i,k).le.(gldep+wcdep-cutdepth)) then
                intrin(i,k) = 5
                depinf(i,k) = depinffix
             end if
          end do
       end do

    end if

    ! foundation ice stream
    ! ---------------------

    if (inflag(06)) then

       cutdepth = 1500.d0
       imin = int(0800.d0*1000.d0/hx)
       imax = int(0900.d0*1000.d0/hx)
       kmin = int(0200.d0*1000.d0/hy)
       kmax = int(0300.d0*1000.d0/hy)

       do i = imin,imax
          do k = kmin,kmax
             if (bpos(i,k).le.(gldep+wcdep-cutdepth)) then
                intrin(i,k) = 6
                depinf(i,k) = depinffix
             end if
          end do
       end do

    end if

    ! support force glacier
    ! ---------------------

    if (inflag(07)) then

       cutdepth = 1300.d0
       imin = int(0850.d0*1000.d0/hx)
       imax = int(0950.d0*1000.d0/hx)
       kmin = int(0450.d0*1000.d0/hy)
       kmax = int(0550.d0*1000.d0/hy)

       do i = imin,imax
          do k = kmin,kmax
             if (bpos(i,k).le.(gldep+wcdep-cutdepth)) then
                intrin(i,k) = 7
                depinf(i,k) = depinffix
             end if
          end do
       end do

    end if

    ! recovery glacier
    ! ----------------

    if (inflag(08)) then

       cutdepth = 1000.d0
       imin = int(0850.d0*1000.d0/hx)
       imax = int(0950.d0*1000.d0/hx)
       kmin = int(0630.d0*1000.d0/hy)
       kmax = int(0750.d0*1000.d0/hy)

       do i = imin,imax
          do k = kmin,kmax
             if (bpos(i,k).le.(gldep+wcdep-cutdepth)) then
                intrin(i,k) = 8
                depinf(i,k) = depinffix
             end if
          end do
       end do

    end if

    ! slessor glacier
    ! ---------------

    if (inflag(09)) then

       cutdepth = 1200.d0
       imin = int(0900.d0*1000.d0/hx)
       imax = int(1000.d0*1000.d0/hx)
       kmin = int(0800.d0*1000.d0/hy)
       kmax = int(0900.d0*1000.d0/hy)

       do i = imin,imax
          do k = kmin,kmax
             if (bpos(i,k).le.(gldep+wcdep-cutdepth)) then
                intrin(i,k) = 9
                depinf(i,k) = depinffix
             end if
          end do
       end do

    end if

    ! ----------------------------------------------
    ! mask regions with grounding line and ice front
    ! ----------------------------------------------

    do i = 1,m_grid
       do k = 1,n_grid
          if ((bpos(i,k).le.0.d0).or.(bpos(i,k).ge.(gldep+wcdep))) then
             intrin(i,k) = 0
             depinf(i,k) = 0.d0
          end if
       end do
    end do

  end subroutine inflow_set_fris

  subroutine inflow_set_larsen()

    implicit none

    ! set larsen inflow depths by hand if using real topography

    ! local variables

    integer :: i,k,imin,imax,kmin,kmax

    real(kind=kdp) :: cutdepth

    ! -----------------------------
    ! set up inflow regions by hand
    ! -----------------------------
    ! unless otherwise indicated the inflows are set up by defining a square
    ! region and then setting all ocean domain within that region to inflow if it
    ! is deeper than a given depth.

    ! the inflows are:
    !  01 - attlee glacier
    !  02 - whole larsen b
    !  03 - whole larsen c

    ! attlee glacier
    ! --------------

    if (inflag(01)) then

       cutdepth = 1000.d0
       imin = int(0065.d0*1000.d0/hx)
       imax = int(0080.d0*1000.d0/hx) 
       kmin = int(0290.d0*1000.d0/hy)
       kmax = int(0310.d0*1000.d0/hy)

       do i = imin,imax
          do k = kmin,kmax
             if (bpos(i,k).le.(gldep+wcdep-cutdepth)) then
                intrin(i,k) = 1
                depinf(i,k) = depinffix
             end if
          end do
       end do

    end if

    ! whole larsen B
    ! --------------

    if (inflag(02)) then

       cutdepth = 0400.d0

       ! south
       imin = int(0130.d0*1000.d0/hx)
       imax = int(0160.d0*1000.d0/hx) 
       kmin = int(0320.d0*1000.d0/hy)
       kmax = int(0425.d0*1000.d0/hy)
       do i = imin,imax
          do k = kmin,kmax
             if (bpos(i,k).le.(gldep+wcdep-cutdepth)) then
                intrin(i,k) = 2
                depinf(i,k) = depinffix
             end if
          end do
       end do

       ! north
       imin = int(0130.d0*1000.d0/hx)
       imax = int(0200.d0*1000.d0/hx) 
       kmin = int(0425.d0*1000.d0/hy)
       kmax = int(0440.d0*1000.d0/hy)
       do i = imin,imax
          do k = kmin,kmax
             if (bpos(i,k).le.(gldep+wcdep-cutdepth)) then
                intrin(i,k) = 2
                depinf(i,k) = depinffix
             end if
          end do
       end do

    end if

    ! whole larsen C (set all>cutdepth to inflow, then remove promontories)
    ! ---------------------------------------------------------------------

    if (inflag(03)) then

       cutdepth = 1000.d0

       ! set all to inflow west of jason peninsula
       imin = int(0000.d0*1000.d0/hx)
       imax = int(0150.d0*1000.d0/hx) 
       kmin = int(0000.d0*1000.d0/hy)
       kmax = int(0310.d0*1000.d0/hy)
       do i = imin,imax
          do k = kmin,kmax
             if (bpos(i,k).le.(gldep+wcdep-cutdepth)) then
                intrin(i,k) = 3
                depinf(i,k) = depinffix
             end if
          end do
       end do

       ! unset cape alexander
       imin = int(0120.d0*1000.d0/hx)
       imax = int(0140.d0*1000.d0/hx) 
       kmin = int(0240.d0*1000.d0/hy)
       kmax = int(0270.d0*1000.d0/hy)
       do i = imin,imax
          do k = kmin,kmax
             intrin(i,k) = 0
             depinf(i,k) = 0.d0
          end do
       end do

       ! unset cole peninsula
       imin = int(0075.d0*1000.d0/hx)
       imax = int(0090.d0*1000.d0/hx) 
       kmin = int(0210.d0*1000.d0/hy)
       kmax = int(0240.d0*1000.d0/hy)
       do i = imin,imax
          do k = kmin,kmax
             intrin(i,k) = 0
             depinf(i,k) = 0.d0
          end do
       end do

       ! unset kenyon peninsula
       imin = int(0070.d0*1000.d0/hx)
       imax = int(0140.d0*1000.d0/hx) 
       kmin = int(0020.d0*1000.d0/hy)
       kmax = int(0060.d0*1000.d0/hy)
       do i = imin,imax
          do k = kmin,kmax
             intrin(i,k) = 0
             depinf(i,k) = 0.d0
          end do
       end do

    end if

    ! ----------------------------------------------
    ! mask regions with grounding line and ice front
    ! ----------------------------------------------

    do i = 1,m_grid
       do k = 1,n_grid
          if ((bpos(i,k).le.0.d0).or.(bpos(i,k).ge.(gldep+wcdep))) then
             intrin(i,k) = 0
             depinf(i,k) = 0.d0
          end if
       end do
    end do

  end subroutine inflow_set_larsen


  subroutine drag_set()

    implicit none

    ! set variable drag-coefficient areas by hand 

    ! local variables
    integer :: imin,imax,kmin,kmax,icen,kcen
    integer, dimension(m_grid,n_grid) :: ii,kk
    real(kind=kdp) :: rotfac

    ! ---------------------------
    ! set up drag regions by hand
    ! ---------------------------
    ! unless otherwise indicated the drag regions are defined such that 
    ! all cells 
    ! within a rotated box (defined by limits and then rotated about centre 
    ! according to rotfac) have a different drag to the background value 

    ! east of henry ice rise
    ! ----------------------

    if (drflag(01)) then

       !rotfac = +055.0d0*pi/180.d0
       !imin = int(0560.d0*1000.d0/hx)
       !imax = int(0620.d0*1000.d0/hx)
       !kmin = int(0370.d0*1000.d0/hy)
       !kmax = int(0500.d0*1000.d0/hy)
       rotfac = +000.0d0*pi/180.d0
       imin = int(0225.d0*1000.d0/hx)
       imax = int(0400.d0*1000.d0/hx)
       kmin = int(0000.d0*1000.d0/hy)
       kmax = int(0110.d0*1000.d0/hy)

       icen = imin + nint(dble(imax-imin)/2.d0)
       kcen = kmin + nint(dble(kmax-kmin)/2.d0)

       ii = get_ii(m_grid,n_grid,icen,kcen,rotfac)
       kk = get_kk(m_grid,n_grid,icen,kcen,rotfac)

       where ( (ii.ge.imin).and.(ii.le.imax) &
            .and.(kk.ge.kmin).and.(kk.le.kmax) )
          drag = cdbvar
       end where

    end if

    ! east of berkner island
    ! ----------------------

    if (drflag(02)) then

       rotfac = +045.0d0*pi/180.d0
       imin = int(0600.d0*1000.d0/hx)
       imax = int(0730.d0*1000.d0/hx)
       kmin = int(0700.d0*1000.d0/hy)
       kmax = int(0830.d0*1000.d0/hy)

       icen = imin + nint(dble(imax-imin)/2.d0)
       kcen = kmin + nint(dble(kmax-kmin)/2.d0)

       ii = get_ii(m_grid,n_grid,icen,kcen,rotfac)
       kk = get_kk(m_grid,n_grid,icen,kcen,rotfac)

       where(((ii.ge.imin).and.(ii.le.imax)) &
            .and.((kk.ge.kmin).and.(kk.le.kmax))) 
          drag = cdbvar
       end where

    end if

    ! spare slot
    ! ----------

    if (drflag(03)) then

       write(*,*) 'error: drag area not set'

    end if

    ! mask regions with grounding line and ice front
    ! ----------------------------------------------

    where ((bpos <= 0.d0).or.(bpos >= (gldep+wcdep)))
       drag = 0.d0
    end where

  end subroutine drag_set

  function get_ii(m_grid,n_grid,icen,kcen,rotfac)

    real(kind=kdp) :: rotfac
    integer :: m_grid,n_grid,icen,kcen
    integer,dimension(m_grid,n_grid) :: get_ii

    get_ii = reshape( &
         (/ ((nint(dble(i-icen)*cos(rotfac)) &
         + nint(dble(k-kcen)*sin(rotfac)) + icen, &
         i=1,m_grid),k=1,n_grid) /), &
         (/ m_grid,n_grid /) )

  end function get_ii

  function get_kk(m_grid,n_grid,icen,kcen,rotfac)

    real(kind=kdp) :: rotfac
    integer :: m_grid,n_grid,icen,kcen
    integer,dimension(m_grid,n_grid) :: get_kk

    get_kk = reshape( &
         (/ ((nint(dble(k-kcen)*cos(rotfac)) &
         - nint(dble(i-icen)*sin(rotfac)) + kcen, &
         i=1,m_grid),k=1,n_grid) /), &
         (/ m_grid,n_grid /) )

  end function get_kk


end module plume_context
