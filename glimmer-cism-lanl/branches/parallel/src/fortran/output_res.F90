subroutine output_res( ewn, nsn, upn, uindx, counter, nu, reswrapped, u_or_v)

! This unwraps the residual on the grid

use glimmer_paramets, only : dp

  implicit none

  integer, intent(in) :: ewn, nsn, upn, counter, nu, u_or_v
  integer, dimension(ewn-1,nsn-1), intent(in) :: uindx
  real (kind = dp), dimension(nu), intent(in) :: reswrapped
  real (kind = dp), dimension(upn,ewn-1,nsn-1) :: resgrid

  integer, dimension(2) :: loc
  integer :: ew, ns, up, upl

  character filename*30

      do ns = 1,nsn-1
         do ew = 1,ewn-1
            if (uindx(ew,ns) /= 0) then
               loc(:) = (uindx(ew,ns) - 1) * (upn + 2) + 1 + (/ 1, upn /)
! loc is the same as given by the function getlocrange in glam_str2.F90
               resgrid(:,ew,ns) = abs(reswrapped(loc(1):loc(2)))
            else 
               resgrid(:,ew,ns) = -0.05d0 !
            end if
         end do
      end do


      do up = 1, upn

         if ( u_or_v .eq. 1 ) then
            write (filename, '("resu_on_grid_k",i2.2,"_l",i3.3,".dat")') &
            counter, up
         elseif ( u_or_v .eq. 2 ) then
            write (filename, '("resv_on_grid_k",i2.2,"_l",i3.3,".dat")') &
            counter, up
         endif
         open (99, file=filename, status='unknown')
         
         do ns = 1,nsn-1
            write(99,10) (resgrid(up ,ew,ns), ew = 1,ewn-1)
         end do
         
         close(99)
         
      enddo
      
 10   format(1x, 1000(f20.16, 1x))

end subroutine output_res

