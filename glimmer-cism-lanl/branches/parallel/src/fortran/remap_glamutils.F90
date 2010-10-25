
module remap_glamutils      

  ! *sfp** contains various subroutines needed when using LANL incremental remapping code
  ! for thickness evolution in glam/glimmer codes

    use glimmer_paramets, only: sp, dp, len0, thk0, tim0, vel0
    use xls

    type remap_glamutils_workspace
        ! *sfp** arrays needed to pass GLAM variables to/from inc. remapping solver
        real (kind = dp), pointer, dimension(:,:,:) ::   &
            thck_ir,            &
            dew_ir,   dns_ir,   &
            dewt_ir,  dnst_ir,  &
            dewu_ir,  dnsu_ir,  &
            hm_ir,    tarea_ir, &
            ubar_ir,  vbar_ir
        real (kind = dp), pointer, dimension(:,:,:,:) :: trace_ir
        real (kind = dp) :: dt_ir

        ! *sfp** mask to apply for domains where initial ice thickness limits are equivalent
        ! to the domain edge (e.g. the locations where bcs are applied)

        real (kind = dp), pointer, dimension(:,:) :: mask_ir

    end type remap_glamutils_workspace
    contains

!----------------------------------------------------------------------

    subroutine horizontal_remap_init (wk, ewn, nsn)

    ! *sfp** initialize variables for use in inc. remapping code   

      implicit none

      type(remap_glamutils_workspace) :: wk

      integer, intent(in) :: ewn, nsn   ! horizontal dimensions

!whl - to do - Set ntrace to the actual number of tracers we want to transport
!              (e.g., ice temperature, ice age)
!              Hardwired ntrace = 1 for now

      integer, parameter  :: ntrace = 1    ! number of tracers

      wk%dt_ir = 0.0_dp     ! time step

      ! allocate arrays/vars 
      allocate( wk%thck_ir(1:ewn-1,1:nsn-1,1) );         wk%thck_ir = 0.0_dp
      allocate( wk%dew_ir(1:ewn-1,1:nsn-1,1) );          wk%dew_ir = 0.0_dp
      allocate( wk%dns_ir(1:ewn-1,1:nsn-1,1) );          wk%dns_ir = 0.0_dp
      allocate( wk%dewt_ir(1:ewn-1,1:nsn-1,1) );         wk%dewt_ir = 0.0_dp
      allocate( wk%dnst_ir(1:ewn-1,1:nsn-1,1) );         wk%dnst_ir = 0.0_dp
      allocate( wk%dewu_ir(1:ewn-1,1:nsn-1,1) );         wk%dewu_ir = 0.0_dp
      allocate( wk%dnsu_ir(1:ewn-1,1:nsn-1,1) );         wk%dnsu_ir = 0.0_dp
      allocate( wk%hm_ir(1:ewn-1,1:nsn-1,1) );           wk%hm_ir = 0.0_dp
      allocate( wk%tarea_ir(1:ewn-1,1:nsn-1,1) );        wk%tarea_ir = 0.0_dp
      allocate( wk%ubar_ir(1:ewn-1,1:nsn-1,1) );         wk%ubar_ir = 0.0_dp
      allocate( wk%vbar_ir(1:ewn-1,1:nsn-1,1) );         wk%vbar_ir = 0.0_dp
      allocate( wk%trace_ir(1:ewn-1,1:nsn-1,ntrace,1) ); wk%trace_ir = 0.0_dp
      allocate( wk%mask_ir(1:ewn,1:nsn) );               wk%mask_ir = 0.0_dp

    end subroutine horizontal_remap_init

!----------------------------------------------------------------------

    subroutine horizontal_remap_final( wk )

      ! *sfp** deallocate variables for use in inc. remapping code   

      implicit none

      type(remap_glamutils_workspace) :: wk

      ! deallocate arrays 
      deallocate( wk%thck_ir )
      deallocate( wk%dew_ir );  deallocate( wk%dns_ir )
      deallocate( wk%dewt_ir ); deallocate( wk%dnst_ir )
      deallocate( wk%dewu_ir ); deallocate( wk%dnsu_ir )
      deallocate( wk%hm_ir )
      deallocate( wk%tarea_ir )
      deallocate( wk%ubar_ir ); deallocate( wk%vbar_ir )
      deallocate( wk%trace_ir )
      deallocate( wk%mask_ir )

    end subroutine horizontal_remap_final

!----------------------------------------------------------------------

    subroutine horizontal_remap_in( wk, dt,       thck,     &
                                    ntrace,   nghost,       &
                                    dew,      dns,          &
                                    uflx,     vflx,         &
                                    stagthck, thklim )

    ! *sfp** get GLAM variables in order for use in inc. remapping code   
    implicit none

    type(remap_glamutils_workspace) :: wk

!whl - to do - add comments describing the arguments

    integer, intent(out) ::   &
         ntrace           ,&! number of tracer arrays to be remapped
         nghost             ! number of ghost cells

    real (kind = dp), dimension(:,:), intent(in) :: thck, uflx, vflx, stagthck
    real (kind = dp), intent(in) :: dew, dns, dt, thklim
    real (kind = dp) :: dt_cfl

    integer :: i, j         ! indices
    integer :: ewn, nsn     ! grid dimensions

!whl - ewn and nsn used below to fill the mask array
    ewn = size(thck,1)
    nsn = size(thck,2)

!whl - Hardwire ntrace and nghost for now
!      Initially, no tracers are actually remapped, but the remapping routine
!       requires ntrace >= 1

    ntrace = 1

!whl - The number of ghost cells is applicable only when we have a parallel code.
!      For remapping, nghost = 2 is ideal because then no boundary updates are needed.
!      
    nghost = 2

!    where( thck > thklim )
!        wk%thck_ir(:,:,1) = thck(:,:)*thk0
!    elsewhere
!        wk%thck_ir(:,:,1) = 0.0d0
!    end where
    
    wk%thck_ir(:,:,1) = thck(:,:)*thk0
    wk%dew_ir(:,:,1)  = dew*len0; wk%dns_ir(:,:,1) = dns*len0
    wk%dewt_ir(:,:,1) = dew*len0; wk%dnst_ir(:,:,1) = dns*len0
    wk%dewu_ir(:,:,1) = dew*len0; wk%dnsu_ir(:,:,1) = dns*len0
    wk%hm_ir(:,:,1) = 1.0_dp
    wk%tarea_ir = 1.0_dp / ( wk%dew_ir * wk%dns_ir )

    call write_xls("uflx.txt", uflx)
    call write_xls("vflx.txt", vflx)

    where( stagthck > 0.0_dp )
        wk%ubar_ir(:,:,1) = uflx/stagthck*vel0;
        wk%vbar_ir(:,:,1) = vflx/stagthck*vel0;
    elsewhere
        wk%ubar_ir(:,:,1) = 0.0_dp
        wk%vbar_ir(:,:,1) = 0.0_dp
    endwhere

!whl - to do - Fill the tracer array with ice temperature and other tracers
    wk%trace_ir(:,:,:,1) = 1.0_dp
    wk%dt_ir = dt * tim0

!whl - bug fix - Commented out next 3 lines because they lead to out-of-bounds error.
!      Note: thck has dims (ewn-1,nsn-1), but mask has dims (ewn,nsn)

!    where( thck > 0.0_dp )
!        wk%mask_ir  = 1.0_dp
!    end where

!whl - This fix works, but it would be better for all horizontal arrays
!      to have the same dimensions.  Work on this later.
    wk%mask_ir(:,:)  = 0.0_dp
    do j = 1, nsn-1
    do i = 1, ewn-1
       if (thck(i,j) > 0.0_dp) then
          wk%mask_ir(i,j) = 1.0_dp
       endif
    enddo
    enddo

    !*tjb** Checks for IR's CFL-like condition.  These only print a warning for now.
    !Use the conservative, "highly divergent" criterion for now.
    dt_cfl = .5 * max(maxval(wk%dew_ir/wk%ubar_ir), maxval(wk%dns_ir/wk%vbar_ir))
    if (dt_cfl < wk%dt_ir) then
        write(*,*) "WARNING: CFL violation in incremental remapping scheme.  Time step should be <= ", dt_cfl
    end if

    end subroutine horizontal_remap_in

!----------------------------------------------------------------------

    subroutine horizontal_remap_out( wk, thck, acab, dt)
    ! *sfp** take output from inc. remapping and put back in GLAM format

    implicit none

    type(remap_glamutils_workspace) :: wk


    real (kind = dp), intent(in) :: dt
    real (kind = sp), intent(in), dimension(:,:) :: acab
    real (kind = dp), dimension(:,:), intent(inout) :: thck

!    integer :: ewn, nsn, ngew, ngns

!whl - bug fix
    integer :: i, j            ! indices
    integer :: ewn_ir, nsn_ir  ! dimensions of IR work arrays
    
    ewn_ir = size(wk%thck_ir,1)
    nsn_ir = size(wk%thck_ir,2)

!whl - debug
    write (50,*) 'ewn_ir, nsn_ir =', ewn_ir, nsn_ir

!whl - bug fix - Changed this because thck_ir has different dimensions from thck.
!      Later, would like for these two fields to have the same dimensions.

    !Map from IR thickness field back to Glide thickness field
!!    thck = wk%thck_ir(:,:,1) / thk0
    do j = 1, nsn_ir
    do i = 1, ewn_ir
       thck(i,j) = wk%thck_ir(i,j,1) / thk0
    enddo
    enddo

    !Apply accumulation
    thck = thck + acab*dt
 
    !Remove thickness from previously masked out locations
    thck = thck * wk%mask_ir
    
    end subroutine horizontal_remap_out

!----------------------------------------------------------------------

end module remap_glamutils

