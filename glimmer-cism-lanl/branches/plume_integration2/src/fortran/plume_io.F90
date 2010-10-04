module plume_io

  !all subroutines that perform io for the model

  use plume_global
  use plume_functions
  use netcdf

  implicit none

  private
  public :: plume_logging_initialize,plume_logging_finalize
  public :: plume_io_initialize,plume_io_finalize
  public :: io_append_output,io_output_sys_time
  public :: io_write_amb
  public :: io_write_time, io_write_calculated_time
  public :: io_write_surface_output
  public :: io_write_long_step_output
  public :: plume_netcdf_write_vars

  logical :: suppress_ascii_output,suppress_logging

  integer,parameter :: foutput = 11, &
       ffrazil = 14, &
       fheader = 12, &
       fstates = 13, &
       famb    = 41

  integer :: ndays,nhours,nmin,nsec

  integer,dimension(8) :: idt_values,dt_values

  ! format strings for output
  character(len=*),parameter :: fmt6200 = &
       '(a,i2.2,a,i2.2,a,i4.4,a,i2.2,a,i2.2,a,i2.2)'
  character(len=*),parameter :: fmt6100 = '(a11,i4,a6,i2,a11,i2,a15,i7,a1)'
  character(len=*),parameter :: fmt6400 = '(a,i4,a,i2,a,i2,a,i7,a)'
  character(len=*),parameter :: fmt6300 = '(a,i4,a,i2,a,i2,a,i2,a)'
  character(len=*),parameter :: fmt6350 = '(a,i3,a,d10.3,a)'
  character(len=*),parameter :: fmt6500 = '(a,i4,a,i2,a,i2,a,i2,a)'

  character(len=512)::output_dir
  character(len=64) :: jobid

  real(kind=kdp) :: odep,wtoi

!!! NETCDF stuff

  ! variables available:
  ! parameters; lxmin,lyin,lamb,lice,linf,ldr,lsm,ldec,m_grid,n_grid
  ! params: lat, c0,ci,nu0,pr,sc,fta,ftb,ftc,si,nus,kt,ks,tiuniform
  ! frazil stuff: neglected for now
  ! inflow stuff: neglected
  ! 1D arrays: dx,dxu,dy,dyv
  ! 1D scalar: samb, tamb, rhovf
  ! 1D : ahdx,ahdxu,ahdy,ahdyv
  ! 2D arrays: pdep,ipos,bpos,bmelt,btemp,bsalt,ctempd,tint,
  ! 2D masks: jcs,jcw,jcd_u,jcd_v,jcd_fl,jcd_negdep,jcd_fseed
  ! 2D scalars: rhop,temp,tempa,tins,salt,salta,rhoamb
  !             entr,atemp,asalt,drag
  ! 2D: u,v,ua,va,su,sv,
  ! 3D arrays: c_ice, ca_ice, ctot, ctota, tfreeze, 

  integer :: nc_id,time_counter
  integer :: time_dimid,x_dimid,y_dimid
  integer :: x_varid,y_varid,time_varid
  integer :: u_varid,v_varid,su_varid,sv_varid
  integer :: pdep_varid,bpos_varid,ipos_varid
  integer :: bmelt_varid,btemp_varid,bsalt_varid
  integer :: rhop_varid,temp_varid,salt_varid,entr_varid,artf_entr_frac_varid
  integer :: jcs_varid,jcw_varid,jcd_u_varid,jcd_v_varid,jcd_fl_varid,jcd_negdep_varid

  integer :: debug_varid

contains

  subroutine plume_logging_initialize(output_dir_in,&
       jobid_in, &
       suppress_logging_in)

    character(len=*),intent(in) :: output_dir_in,jobid_in
    logical,intent(in) :: suppress_logging_in

    suppress_logging = suppress_logging_in
    output_dir = trim(output_dir_in)
    jobid = trim(jobid_in)

    if (.not. suppress_logging) then
       ! open output file 
       open(foutput,file=trim(output_dir)//trim(jobid)//'_output',action='write',status='replace')   
    end if

  end subroutine plume_logging_initialize

  subroutine plume_logging_finalize()

    if (.not. suppress_logging) close(foutput)

  end subroutine plume_logging_finalize

  subroutine plume_io_initialize(suppress_ascii_output_in, &
       plume_output_nc_file)

!!!! NB: this should be called AFTER plume_initialise so that certain global
    !        constants are sure to have been set.

    character(len=*),intent(in) :: plume_output_nc_file !should be an absolute path
                                                        !or relative to execution directory
    logical,intent(in) :: suppress_ascii_output_in

    suppress_ascii_output = suppress_ascii_output_in

    ! calculate origin depth
    odep=wcdep+gldep

    ! calculate water m/s to ice m/y conversion factor for output
    wtoi=31557600.d0*(rho0/rhoi)

    call plume_netcdf_init(plume_output_nc_file)

  end subroutine plume_io_initialize

  subroutine io_output_sys_time(prefix)

    character(len=*),intent(in) :: prefix

    integer,dimension(8) :: dt_values
    if (suppress_logging) return

    call date_and_time(VALUES=dt_values)

    write(*,fmt6200) prefix,dt_values(3),'/',&
         dt_values(2),'/',dt_values(1),"; time ", &
         & dt_values(5),':',dt_values(6),':',dt_values(7)

    write(foutput,fmt6200) prefix,dt_values(3),'/',&
         dt_values(2),'/',dt_values(1),"; time ", & 
         & dt_values(5),':',dt_values(6),':',dt_values(7)

  end subroutine io_output_sys_time

  subroutine io_output_plume_vol(totvol)

    real(kind=kdp),intent(in) :: totvol

    if (suppress_logging) return

    write(foutput,*)'total plume volume = ',totvol, ' m^3'
    write(*,*)'total plume volume = ',totvol, ' m^3'

  end subroutine io_output_plume_vol

  subroutine io_write_calculated_time(runtim,nstep)

    real(kind=kdp),intent(in) :: runtim
    integer,intent(in) :: nstep

    if (suppress_logging) return

    ndays=int(runtim/86400.d0)
    nhours=int((runtim-ndays*86400.d0)/3600.d0)
    nmin=int((runtim-ndays*86400.d0-nhours*3600.d0)/60.d0)


    write(foutput,fmt6100) 'calculated ',ndays,' days ',nhours,&
         ' hours and ',nmin,' minutes (step',nstep,')'

    write(*,fmt6100) 'calculated ',ndays,' days ',nhours, &
         ' hours and ',nmin,' minutes (step',nstep,')'

  end subroutine io_write_calculated_time

  subroutine io_write_time(tim,prefix)

    integer,intent(in) :: tim
    character(len=*),intent(in)::prefix

    if (suppress_logging) return

    ndays=int(tim/86400.d0)
    nhours=int((tim-ndays*86400.d0)/3600.d0)
    nmin=int((tim-ndays*86400.d0-nhours*3600.d0)/60.d0)
    nsec=int(tim-ndays*86400.d0-nhours*3600.d0-nmin*60.d0)

    write(foutput,fmt6500) prefix ,ndays,&
         ' days ',nhours,' hours ', &
         nmin,' minutes ',nsec, ' seconds'
    write(foutput,*) ' '


    write(*,fmt6500) prefix,ndays, &
         ' days ',nhours,' hours ', &
         nmin,' minutes ',nsec, ' seconds'
    write(*,*) ' '

  end subroutine io_write_time

  subroutine io_write_long_step_output(icalcan,icalcen,kcalcan,kcalcen,&
       varoutrat, negdep)

    integer :: icalcan,icalcen,kcalcan,kcalcen,varoutrat
    integer :: i,k,l
    real(kind=kdp) :: negdep
    real(kind=kdp) :: totvol

    if (suppress_logging) return

    ! write time and elapsed system time
    call io_append_output(' ')
    call io_output_sys_time('current system time ')
    call io_write_time(get_elapsed_sys_time(), ' elapsed system time ')

    ! write totals of plume volume, negative depths, and negative frazil

    totvol = 0.d0
    do i = 1,m_grid
       do k = 1,n_grid
          totvol = totvol + pdep(i,k)*dx(i)*dy(k)
       end do
    end do

    call io_output_plume_vol(totvol)

    if (negdep.ne.0.d0) then
       write(foutput,*) 'warning: negdep = ',negdep,' (cumulative)'  
       write(*,*) 'warning: negdep = ',negdep,' (cumulative)'

    end if

    do l=1,nice
       if (frzcut(l).gt.0.d0) then
          write(foutput,fmt6350) ' warning: frzcut(',l,') = ',frzcut(l), &
               ' (cumulative)'
          write(*,fmt6350) ' warning: frzcut(',l,') = ',frzcut(l), &
               ' (cumulative)'

       end if
    end do

    ! write ascii rendering of plume thickness
    !    
    call io_plotascii(icalcan,icalcen,kcalcan,kcalcen,varoutrat)


  end subroutine io_write_long_step_output

  subroutine io_write_surface_output(runtim,labtim)

    real(kind=kdp)::runtim,labtim

    ndays=int(runtim/86400.d0)
    nhours=int((runtim-ndays*86400.d0)/3600.d0)
    nmin=int((runtim-ndays*86400.d0-nhours*3600.d0)/60.d0)

    call io_output_files(runtim,labtim)

    if (suppress_logging) return

    write(foutput,fmt6400) ' file output at ',ndays,' days ',nhours,&
         ' hours and ',nmin,&
         ' minutes (step ', int(runtim/dt),')'
    write(foutput,*) ' '
    write(*,fmt6400) ' file output at ',ndays,' days ',nhours, &
         ' hours and ',nmin,&
         ' minutes (step ', int(runtim/dt),')'
    write(*,*) ' '

  end subroutine io_write_surface_output

  subroutine io_append_output(output_data)

    character(len=*),intent(in) :: output_data

    if (suppress_logging) return

    write(foutput,*) output_data
    write(*,*) output_data

  end subroutine io_append_output

  subroutine io_write_amb(namb,depth,tamb,ttt,samb,rhovf,rhopot)

    integer,intent(in) :: namb
    real(kind=kdp),dimension(namb),intent(in) :: depth,tamb,ttt,samb,rhovf,rhopot

    integer :: i
    open(famb,file='ambout')    
    do i=1,namb
       write(famb,*) depth(i),tamb(i),ttt(i),samb(i),rhovf(i),rhopot(i)
    end do
    close(famb)

  end subroutine io_write_amb

  subroutine io_output_files(runtim,labtim)                    

    ! produces surface output (single file version)

    implicit none

    real(kind=kdp),intent(in) :: runtim,labtim

    ! local variables

    character(len=4) :: outftim
    character(len=512) :: path_prefix
    integer :: i,k,l,izo,izu
    real(kind=kdp) :: rhopp(m_grid,n_grid),rhoap(m_grid,n_grid)
    real(kind=kdp) :: zc,difu,difo,tambz,sambz


    ! calculate potential density for ambient and plume (to give differences)

    do i = 1,m_grid
       do k = 1,n_grid
          zc =  wcdep + gldep - ipos(i,k) 
          izo = int(zc/dzincr) + 1
          izu = izo + 1
          difu = dble(izo)*dzincr - zc
          difo = dzincr - difu
          tambz = (difu*tamb(izo) + difo*tamb(izu))/dzincr
          sambz = (difu*samb(izo) + difo*samb(izu))/dzincr

          if (rholinear) then
             rhoap(i,k) = rho_func_linear(tambz,sambz)
             rhopp(i,k) = rho_func_linear(temp(i,k),salt(i,k))
          else
             rhoap(i,k) = rho_func_nonlinear(tambz,sambz,0.d0)
             rhopp(i,k) = rho_func_nonlinear(temp(i,k),salt(i,k),0.d0)
          end if

          if (frazil) then
             rhopp(i,k) = (1.d0 - ctot(i,k))*rhopp(i,k) + ctot(i,k)*rhoi
          end if

       end do
    end do

    if (.not. suppress_ascii_output) then

       ! make file name in units of labtim
       write(outftim,'(i4.4)') int(runtim/labtim)      

       ! write header file
       ! -----------------
       path_prefix = trim(output_dir)//'/'//trim(jobid)//'_'//outftim
       open(fheader,file= trim(path_prefix) //'_header')

       ! write switches
       write(fheader,'(l1)') frazil

       ! write dimensions
       write(fheader,'(3i4)') m_grid,n_grid,nice

       ! write grid
       do i = 1,m_grid
          write(fheader,'(d13.7)') dx(i)
       end do
       do k = 1,n_grid
          write(fheader,'(d13.7)') dy(k)
       end do

       close(fheader)


       ! make file name in units of labtim
       write(outftim,'(i4.4)') int(runtim/labtim)      


       ! note output formats differ from model variables:
       !    draft and plume base position are output relative to sea level
       !    all rates are output in m/y and all melt rates are of ice

       open(fstates,file=trim(path_prefix) // '_states')

       do i = 1,m_grid
          do k = 1,n_grid
             write(fstates,'(21("  ",f16.7))') real(jcw(i,k)),su(i,k), &
                  sv(i,k),temp(i,k),salt(i,k), &
                  pdep(i,k),bpos(i,k)-odep,ipos(i,k)-odep,  &
                  rhoamb(i,k),rhop(i,k),rhoap(i,k),rhopp(i,k), &
                  entr(i,k)*31557600.d0,atemp(i,k),asalt(i,k), &
                  bmelt(i,k)*wtoi,btemp(i,k),bsalt(i,k), &
                  tfreeze(i,k),depinf(i,k),real(jcs(i,k))
          end do
       end do

       close(fstates)


       ! write frazil file (if required)
       ! -------------------------------

       ! note output formats differ from model variables:
       !    all ice interaction rates are output in m/y of ice

       if (frazil) then

          open(ffrazil,file=trim(path_prefix) // '_frazil')

          ! write total fields

          do i = 1,m_grid
             do k = 1,n_grid
                write(ffrazil,'(2(" ",d13.7))') ctot(i,k),ctempd(i,k)
             end do
          end do

          ! write individual size-class fields

          do i = 1,m_grid
             do k = 1,n_grid
                do l = 1,nice
                   write(ffrazil,'(4(" ",d13.7))') c_ice(i,k,l),fmelt(i,k,l)*wtoi, &
                        fppn(i,k,l)*wtoi,fnuc(i,k,l)*wtoi
                end do
             end do
          end do

          close(ffrazil)

       end if

    end if

!!! NB: we write out netcdf output even if suppress_output is .true.

    ! this is redundant but I am switching over to using netcdf output
    call plume_netcdf_write_vars(runtim)

  end subroutine io_output_files


  subroutine io_plotascii(icalcan,icalcen,kcalcan,kcalcen,varoutrat)

    ! plot ascii output of a field

    implicit none

    integer,intent(in) :: icalcan,icalcen,kcalcan,kcalcen,varoutrat

    ! local variables
    character(len=1) :: a(m_grid,n_grid),deca(ldec),dry,land,wet,sep,inf
    integer :: i,k,jci
    integer :: minc,ninc,msize,nsize,izo,izu

    real(kind=kdp) :: tmp,xsig,depth,difu,difo,rhoatmp

    character(len=*),parameter :: fmt50 = '(5x,37(i1," "))'
    character(len=*),parameter :: fmt100 = '(1x,i3,1x,74a1)'

    ! data definitions

    data dry,land,wet,sep,inf/'.','l','w','*','>'/
    data deca/'w','a','b','c','d','e','f','g','h','i','j','k', &
         'l','m','n','o','p','q','r','s','t','u','v','x','y','z'/


    if (suppress_logging) return

    xsig = 5.0d-1

    ! if plume is too large for screen then use coarser resolution
    ! (varible grid step in each direction chosen to fit 50/70 box

    minc = 1
    ninc = 1
    msize = icalcen - icalcan + 1
    nsize = kcalcen - kcalcan + 1
    if (msize.gt.50) minc = msize/50 + 1 
    if (nsize.gt.74) ninc = nsize/74 + 1 

    ! choose output resolution fixed solely on extent of plume in one direction

    if (varoutrat.eq.2) then
       minc = ninc
    end if

    if (varoutrat.eq.3) then
       ninc = minc
    end if

    ! calculate array of plume/land/ambient values

    do i = 1,m_grid
       do k = 1,n_grid
          a(i,k) = land
          if (jcs(i,k).le.0) cycle
          a(i,k) = dry
          if (jcw(i,k).le.0) cycle
          a(i,k) = wet
          tmp = pdep(i,k)
          if (abs(tmp).lt. dcr) cycle
          tmp = tmp + sign(xsig,tmp)
          jci = int(abs(tmp)/10.d0) + 1
          jci = max0(jci,1)
          jci = min0(jci,ldec)
          a(i,k) = deca(jci)
       end do
    end do

    ! over-write plume thickness in inflow region

    do i = 1,m_grid
       do k = 1,n_grid
          if (depinf(i,k).gt.0.d0) a(i,k) = inf
       end do
    end do

    ! over-write plume thickness where it has separated

    do i = 1,m_grid
       do k = 1,n_grid
          depth = wcdep + gldep - ipos(i,k) - 5.0d-1*pdep(i,k)
          izo = int(depth/dzincr) + 1
          izu = izo + 1
          difu = dble(izo)*dzincr - depth
          difo = dzincr - difu
          rhoatmp = (difu*rhovf(izo)+difo*rhovf(izu))/dzincr
          if ((jcw(i,k).eq.1).and.((rhoatmp < rhop(i,k) + septol))) &
               a(i,k) = sep
       end do
    end do

    ! draw horizontal scale

    call io_append_output(' ')

    write(*,fmt50) (k/100,k = kcalcan,kcalcen,2*ninc)
    write(*,fmt50) ((k - (k/100*100))/10,k = kcalcan,kcalcen,2*ninc)
    write(*,fmt50) (k - (k/100*100) - (k - (k/100*100))/10*10 &
         ,k = kcalcan,kcalcen,2*ninc)
    write(*,*)' '


    write(foutput,fmt50) (k/100,k = kcalcan,kcalcen,2*ninc)
    write(foutput,fmt50) ((k - (k/100*100))/10,k = kcalcan,kcalcen,2*ninc)
    write(foutput,fmt50) (k - (k/100*100) - (k - (k/100*100))/10*10 &
         ,k = kcalcan,kcalcen,2*ninc)
    write(foutput,*)' '

    ! draw plume and vertical scale
    do i = icalcan,icalcen,minc
       write(*,fmt100) i,(a(i,k),k = kcalcan,kcalcen,ninc)
       write(foutput,fmt100) i,(a(i,k),k = kcalcan,kcalcen,ninc)

    end do

    write(*,*)' '
    write(foutput,*)' '


  end subroutine io_plotascii

  subroutine plume_io_finalize()

    call plume_netcdf_finalize()

  end subroutine plume_io_finalize

  subroutine plume_netcdf_init(ncdf_filename)

    !open a netcdf file and define all variables

    character(len=*),intent(in) :: ncdf_filename

    ! local variables
    integer :: i,j

    time_counter = 1	

    call check( nf90_create(ncdf_filename, NF90_CLOBBER, nc_id) )

    ! defining dimensions
    call check( nf90_def_dim(nc_id,'time',NF90_UNLIMITED,time_dimid) )
    call check( nf90_def_dim(nc_id,'x',m_grid,x_dimid) )
    call check( nf90_def_dim(nc_id,'y',n_grid,y_dimid) )

    call io_append_output('Creating variable time')
    call check( nf90_def_var(nc_id,'time',NF90_DOUBLE,(/time_dimid/),time_varid) )
    call check( nf90_put_att(nc_id, time_varid, 'long_name', 'Model time') )
    call check( nf90_put_att(nc_id, time_varid, 'standard_name', 'time') )
    call check( nf90_put_att(nc_id, time_varid, 'calendar', 'none')  )
    call check( nf90_put_att(nc_id, time_varid, 'units', 'year since 1-1-1 0:0:0') )

    call io_append_output('Creating variable x')
    call check( nf90_def_var(nc_id,'x',NF90_DOUBLE,(/x_dimid/),x_varid) )
    call check( nf90_put_att(nc_id, x_varid, 'long_name', 'Cartisian x-coordinate') )
    call check( nf90_put_att(nc_id, x_varid, 'units', 'meter') )

    call io_append_output('Creating variable y')
    call check( nf90_def_var(nc_id,'y',NF90_DOUBLE,(/y_dimid/),y_varid) )
    call check( nf90_put_att(nc_id, y_varid, 'long_name', 'Cartisian y-coordinate') )
    call check( nf90_put_att(nc_id, y_varid, 'units', 'meter') )

    call io_append_output('Creating variable u')
    call check( nf90_def_var(nc_id,'u',NF90_DOUBLE,(/x_dimid,y_dimid,time_dimid/),u_varid) )
    call check( nf90_put_att(nc_id, u_varid, 'positive', 'east') )
    call check( nf90_put_att(nc_id, u_varid, 'long_name', 'zonal transport') )
    call check( nf90_put_att(nc_id, u_varid, 'standard_name', 'zonal transport') )
    call check( nf90_put_att(nc_id, u_varid, 'units', 'meters^2/second') )

    call io_append_output('Creating variable v')
    call check( nf90_def_var(nc_id,'v',NF90_DOUBLE,(/x_dimid,y_dimid,time_dimid/),v_varid) )
    call check( nf90_put_att(nc_id, v_varid, 'positive', 'north') )
    call check( nf90_put_att(nc_id, v_varid, 'long_name', 'meridional transport') )
    call check( nf90_put_att(nc_id, v_varid, 'standard_name', 'meridional transport') )
    call check( nf90_put_att(nc_id, v_varid, 'units', 'meters^2/second') )

    call io_append_output('Creating variable su')
    call check( nf90_def_var(nc_id,'su',NF90_DOUBLE,(/x_dimid,y_dimid,time_dimid/),su_varid) )
    call check( nf90_put_att(nc_id, su_varid, 'positive', 'east') )
    call check( nf90_put_att(nc_id, su_varid, 'long_name', 'zonal velocity') )
    call check( nf90_put_att(nc_id, su_varid, 'standard_name', 'zonal velocity') )
    call check( nf90_put_att(nc_id, su_varid, 'units', 'meters/second') )

    call io_append_output('Creating variable sv')
    call check( nf90_def_var(nc_id,'sv',NF90_DOUBLE,(/x_dimid,y_dimid,time_dimid/),sv_varid) )
    call check( nf90_put_att(nc_id, sv_varid, 'positive', 'north') )
    call check( nf90_put_att(nc_id, sv_varid, 'long_name', 'meridional velocity') )
    call check( nf90_put_att(nc_id, sv_varid, 'standard_name', 'meridional velocity') )
    call check( nf90_put_att(nc_id, sv_varid, 'units', 'meters/second') )

    call io_append_output('Creating variable pdep')
    call check( nf90_def_var(nc_id,'pdep',NF90_DOUBLE,(/x_dimid,y_dimid,time_dimid/),pdep_varid) )
    call check( nf90_put_att(nc_id, pdep_varid, 'long_name', 'plume depth') )
    call check( nf90_put_att(nc_id, pdep_varid, 'standard_name', 'plume depth') )
    call check( nf90_put_att(nc_id, pdep_varid, 'units', 'meters') )

    call io_append_output('Creating variable bpos')
    call check( nf90_def_var(nc_id,'bpos',NF90_DOUBLE,(/x_dimid,y_dimid,time_dimid/),bpos_varid) )
    call check( nf90_put_att(nc_id, bpos_varid, 'positive', 'up') )
    call check( nf90_put_att(nc_id, bpos_varid, 'long_name', 'ice shelf lower surface height') )
    call check( nf90_put_att(nc_id, bpos_varid, 'standard_name', 'ice shelf lower surface height') )
    call check( nf90_put_att(nc_id, bpos_varid, 'units', 'meters') )

    call io_append_output('Creating variable ipos')
    call check( nf90_def_var(nc_id,'ipos',NF90_DOUBLE,(/x_dimid,y_dimid,time_dimid/),ipos_varid) )
    call check( nf90_put_att(nc_id, ipos_varid, 'positive', 'up') )
    call check( nf90_put_att(nc_id, ipos_varid, 'long_name', 'lower plume surface') )
    call check( nf90_put_att(nc_id, ipos_varid, 'standard_name', 'lower plume surface') )
    call check( nf90_put_att(nc_id, ipos_varid, 'units', 'meters') )

    call io_append_output('Creating variable bmelt')
    call check( nf90_def_var(nc_id,'bmelt',NF90_DOUBLE,(/x_dimid,y_dimid,time_dimid/),bmelt_varid) )
    call check( nf90_put_att(nc_id, bmelt_varid, 'positive', 'thickening plume') )
    call check( nf90_put_att(nc_id, bmelt_varid, 'long_name', 'basal melt rate') )
    call check( nf90_put_att(nc_id, bmelt_varid, 'standard_name', 'basal melt rate') )
    call check( nf90_put_att(nc_id, bmelt_varid, 'units', 'meters/year') )

    call io_append_output('Creating variable btemp')
    call check( nf90_def_var(nc_id,'btemp',NF90_DOUBLE,(/x_dimid,y_dimid,time_dimid/),btemp_varid) )
    call check( nf90_put_att(nc_id, btemp_varid, 'long_name', 'ice shelf interface temperature') )
    call check( nf90_put_att(nc_id, btemp_varid, 'standard_name', 'ice shelf interface temperature') )
    call check( nf90_put_att(nc_id, btemp_varid, 'units', 'K') )

    call io_append_output('Creating variable bsalt')
    call check( nf90_def_var(nc_id,'bsalt',NF90_DOUBLE,(/x_dimid,y_dimid,time_dimid/),bsalt_varid) )
    call check( nf90_put_att(nc_id, bsalt_varid, 'long_name', 'ice shelf interface salinity') )
    call check( nf90_put_att(nc_id, bsalt_varid, 'standard_name', 'ice shelf interface salinity') )
    call check( nf90_put_att(nc_id, bsalt_varid, 'units', 'psu') )

    call io_append_output('Creating variable rhop')
    call check( nf90_def_var(nc_id,'rhop',NF90_DOUBLE,(/x_dimid,y_dimid,time_dimid/),rhop_varid) )
    call check( nf90_put_att(nc_id, rhop_varid, 'long_name', 'plume density') )
    call check( nf90_put_att(nc_id, rhop_varid, 'standard_name', 'plume density') )
    call check( nf90_put_att(nc_id, rhop_varid, 'units', 'kg/meters^3') )

    call io_append_output('Creating variable temp')
    call check( nf90_def_var(nc_id,'temp',NF90_DOUBLE,(/x_dimid,y_dimid,time_dimid/),temp_varid) )
    call check( nf90_put_att(nc_id, temp_varid, 'long_name', 'plume temperature') )
    call check( nf90_put_att(nc_id, temp_varid, 'standard_name', 'plume temperature') )
    call check( nf90_put_att(nc_id, temp_varid, 'units', 'K') )

    call io_append_output('Creating variable salt')
    call check( nf90_def_var(nc_id,'salt',NF90_DOUBLE,(/x_dimid,y_dimid,time_dimid/),salt_varid) )
    call check( nf90_put_att(nc_id, salt_varid, 'long_name', 'plume salinity') )
    call check( nf90_put_att(nc_id, salt_varid, 'standard_name', 'plume salinity') )
    call check( nf90_put_att(nc_id, salt_varid, 'units', 'psu') )

    call io_append_output('Creating variable entr')
    call check( nf90_def_var(nc_id,'entr',NF90_DOUBLE,(/x_dimid,y_dimid,time_dimid/),entr_varid) )
    call check( nf90_put_att(nc_id, entr_varid, 'positive', 'increasing plume thickness') )
    call check( nf90_put_att(nc_id, entr_varid, 'long_name', 'entrainment velocity') )
    call check( nf90_put_att(nc_id, entr_varid, 'standard_name', 'entrainment velocity') )
    call check( nf90_put_att(nc_id, entr_varid, 'units', 'meters/sec') )

    call io_append_output('Creating variable jcs')
    call check( nf90_def_var(nc_id,'jcs',NF90_DOUBLE,(/x_dimid,y_dimid,time_dimid/),jcs_varid) )
    call check( nf90_put_att(nc_id, jcs_varid, 'long_name', 'ocean mask') )
    call check( nf90_put_att(nc_id, jcs_varid, 'standard_name', 'ocean mask') )

    call io_append_output('Creating variable jcw')
    call check( nf90_def_var(nc_id,'jcw',NF90_DOUBLE,(/x_dimid,y_dimid,time_dimid/),jcw_varid) )
    call check( nf90_put_att(nc_id, jcw_varid, 'long_name', 'wet mask') )
    call check( nf90_put_att(nc_id, jcw_varid, 'standard_name', 'wet mask') )

    call io_append_output('Creating variable jcd_u')
    call check( nf90_def_var(nc_id,'jcd_u',NF90_DOUBLE,(/x_dimid,y_dimid,time_dimid/),jcd_u_varid) )
    call check( nf90_put_att(nc_id, jcd_u_varid, 'long_name', 'u mask') )
    call check( nf90_put_att(nc_id, jcd_u_varid, 'standard_name', 'u mask') )

    call io_append_output('Creating variable jcd_v')
    call check( nf90_def_var(nc_id,'jcd_v',NF90_DOUBLE,(/x_dimid,y_dimid,time_dimid/),jcd_v_varid) )
    call check( nf90_put_att(nc_id, jcd_v_varid, 'long_name', 'v mask') )
    call check( nf90_put_att(nc_id, jcd_v_varid, 'standard_name', 'v mask') )

    call io_append_output('Creating variable artf_entr_frac')
    call check( nf90_def_var(nc_id,'artf_entr_frac',NF90_DOUBLE,(/x_dimid,y_dimid,time_dimid/),artf_entr_frac_varid))
    call check( nf90_put_att(nc_id, artf_entr_frac_varid, 'long_name', 'fraction of artificial entr') )
    call check( nf90_put_att(nc_id, artf_entr_frac_varid, 'standard_name', 'fraction of artificial entr') )

    call io_append_output('Creating variable debug')
    call check( nf90_def_var(nc_id,'debug',         NF90_DOUBLE,(/x_dimid,y_dimid,time_dimid/),debug_varid))
   

    !call io_append_output('Creating variable jcd_negdep')
    !call check( nf90_def_var(nc_id,'jcd_negdep',NF90_DOUBLE,(/x_dimid,y_dimid,time_dimid/),jcd_negdep_varid) )
    !call check( nf90_put_att(nc_id, jcd_negdep_varid, 'long_name', 'negative depth mask') )
    !call check( nf90_put_att(nc_id, jcd_negdep_varid, 'standard_name', 'negative depth mask') )

    !call io_append_output('Creating variable jcd_fl')
    !call check( nf90_def_var(nc_id,'jcd_fl',NF90_DOUBLE,(/x_dimid,y_dimid,time_dimid/),jcd_fl_varid) )
    !call check( nf90_put_att(nc_id, jcd_fl_varid, 'long_name', 'flooded mask') )
    !call check( nf90_put_att(nc_id, jcd_fl_varid, 'standard_name', 'flooded mask') )

    call check( nf90_enddef(nc_id) )

    !now populate the dimension variables

    call check( nf90_put_var(nc_id,x_varid,(/ (i*hx, i=1,m_grid) /)) )
    call check( nf90_put_var(nc_id,y_varid,(/ (j*hy, j=1,n_grid) /)) )
   
    call check( nf90_sync(nc_id) )

  end subroutine plume_netcdf_init

  subroutine plume_netcdf_write_vars(time)

    real(kind=kdp),intent(in) :: time !what time (seconds) in the simulation does the state represent

    integer :: status
    real(kind=kdp),parameter :: sec_per_year = 365.25d0 * 24.d0*3600.d0

    call check(nf90_put_var(nc_id, time_varid, time / sec_per_year,(/ time_counter /)))
    call check(nf90_put_var(nc_id, u_varid, utrans, (/1,1,time_counter/)))
    call check(nf90_put_var(nc_id, v_varid, vtrans, (/1,1,time_counter/)))
    call check(nf90_put_var(nc_id, su_varid, su, (/1,1,time_counter/)))     
    call check(nf90_put_var(nc_id, sv_varid, sv, (/1,1,time_counter/)))     
    call check(nf90_put_var(nc_id, pdep_varid, pdep, (/1,1,time_counter/)))
    call check(nf90_put_var(nc_id, ipos_varid, ipos, (/1,1,time_counter/)))
    call check(nf90_put_var(nc_id, bpos_varid, bpos, (/1,1,time_counter/)))
    call check(nf90_put_var(nc_id, btemp_varid, btemp, (/1,1,time_counter/)))
    call check(nf90_put_var(nc_id, bmelt_varid, bmelt*sec_per_year, (/1,1,time_counter/)))
    call check(nf90_put_var(nc_id, bsalt_varid, bsalt, (/1,1,time_counter/)))
    call check(nf90_put_var(nc_id, rhop_varid, rhop, (/1,1,time_counter/)))
    call check(nf90_put_var(nc_id, temp_varid, temp, (/1,1,time_counter/)))
    call check(nf90_put_var(nc_id, salt_varid, salt, (/1,1,time_counter/)))
    call check(nf90_put_var(nc_id, entr_varid, entr, (/1,1,time_counter/)))
    call check(nf90_put_var(nc_id, jcs_varid, jcs, (/1,1,time_counter/)))
    call check(nf90_put_var(nc_id, jcw_varid, jcw, (/1,1,time_counter/)))
    call check(nf90_put_var(nc_id, jcd_u_varid, jcd_u, (/1,1,time_counter/)))
    call check(nf90_put_var(nc_id, jcd_v_varid, jcd_v, (/1,1,time_counter/)))
    call check(nf90_put_var(nc_id, artf_entr_frac_varid, artf_entr_frac, (/1,1,time_counter/)))
   
    call check(nf90_put_var(nc_id, debug_varid, debug, (/1,1,time_counter/)))

    !call check(nf90_put_var(nc_id, jcd_fl_varid, jcd_fl, (/1,1,time_counter/)))
    !call check(nf90_put_var(nc_id, jcd_negdep_varid, jcd_negdep, (/1,1,time_counter/)))
    
    call check(nf90_sync(nc_id))

    time_counter = time_counter + 1

  end subroutine plume_netcdf_write_vars

  subroutine plume_netcdf_finalize()

    call check(nf90_close(nc_id))

  end subroutine plume_netcdf_finalize

  subroutine check(status_code)

    integer,intent(in) :: status_code

    if (status_code .ne. 0) then
       call io_append_output(NF90_STRERROR(status_code))
       call io_append_output('fatal netcdf error')
       stop
    end if

  end subroutine check

end module plume_io
