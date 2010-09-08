
module glam_bnd_cond

use glide_types, only : glide_bnd_cond_params
use glimmer_paramets, only : dp

implicit none

contains

function is_bnd_point(bnd_cond_params, ew, ns)

type(glide_bnd_cond_params),intent(in) :: bnd_cond_params
integer, intent(in) :: ew,ns
logical :: is_bnd_point

if (ew == 1) then
   is_bnd_point = .true.
else
   is_bnd_point = .false.
end if

end function

subroutine get_normal_vector(bnd_cond_params, ew,ns, norm_vector)

type(glide_bnd_cond_params),intent(in) :: bnd_cond_params
integer, intent(in) :: ew,ns
real(kind=dp),dimension(2),intent(out) :: norm_vector

norm_vector(1) = 1.0
norm_vector(2) = 0.0

end subroutine 

end module glam_bnd_cond