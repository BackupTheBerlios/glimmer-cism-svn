subroutine output_res( ewn, nsn, upn, uindx, counter, nu, reswrapped, u_or_v)

! This routine can be used for debugging and evalution of the nonlinear 
! numerical scheme (Picard iteration...and later Newton). It unwraps the 
! residual on the grid, either for the u comp (u_or_v = 1) or the 
! v comp (u_or_v = 2). jfl

use glimmer_paramets, only : dp
use netcdf

  implicit none

  integer, intent(in) :: ewn, nsn, upn, counter, nu, u_or_v
  integer, dimension(ewn-1,nsn-1), intent(in) :: uindx
  real (kind = dp), dimension(nu), intent(in) :: reswrapped
  real (kind = dp), dimension(ewn-1,nsn-1, upn) :: resgrid

  integer, dimension(2) :: loc
  integer :: ew, ns, up, upl, nc_id, ns_id, ew_id, up_id, res_id, status
  integer :: dimids(3)

  character filename*30

      do ns = 1,nsn-1 !had to flip indices in matrix because of netcdf
         do ew = 1,ewn-1
            if (uindx(ew,ns) /= 0) then
               loc(:) = (uindx(ew,ns) - 1) * (upn + 2) + 1 + (/ 1, upn /)
! loc is the same as given by the function getlocrange in glam_str2.F90
               resgrid(ew,ns,:) = abs(reswrapped(loc(1):loc(2)))
            else 
               resgrid(ew,ns,:) = -999d0 ! no ice
            end if
         end do
      end do

      if ( u_or_v .eq. 1) then
         write (filename, '("res_ucomp_k",i2.2,".nc")') counter
      elseif ( u_or_v .eq. 2) then
         write (filename, '("res_vcomp_k",i2.2,".nc")') counter
      endif

      call check(nf90_create(filename, NF90_CLOBBER, nc_id))
      call check(nf90_def_dim(nc_id,'n-s',nsn-1,ns_id))
      call check(nf90_def_dim(nc_id,'e-w',ewn-1,ew_id))
      call check(nf90_def_dim(nc_id,'up',upn,up_id))

      dimids(3) = up_id 
      dimids(1) = ew_id
      dimids(2) = ns_id

      call check(nf90_def_var(nc_id,'residual',NF90_DOUBLE, dimids, res_id))
      call check( nf90_enddef(nc_id) )
      call check(nf90_put_var(nc_id,res_id,resgrid))
      call check(nf90_close(nc_id))

end subroutine output_res

subroutine check(status_code)
      use netcdf
      integer,intent(in) :: status_code

      if (status_code /= 0) then
         write(*,*) 'fatal netcdf error:',nf90_strerror(status_code)
         stop
      end if

end subroutine check

