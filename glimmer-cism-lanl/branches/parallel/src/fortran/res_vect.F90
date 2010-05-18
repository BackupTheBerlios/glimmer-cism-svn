!  uvec is either u^l-1 or v^k-1 on input and Av-b or Cu-d on output

subroutine res_vect ( matrix, uvec, bvec, nu, counter, g_flag )


use glimmer_paramets, only : dp
use glimmer_sparse_type
use glide_mask

implicit none

integer :: i, j, counter, nu, nele ! nu: size of uvec and bvec
integer, dimension(nu), intent(in) :: g_flag ! g_flag = 1 for ghost cell

type(sparse_matrix_type),  intent(in) :: matrix

real (kind = dp), dimension(nu), intent(in) :: bvec
real (kind = dp), dimension(nu), intent(inout) :: uvec
real (kind = dp), dimension(nu) :: Au_b

! calculate residual vector of the u OR v component

      Au_b = 0d0
      
      do nele = 1, matrix%nonzeros 
 
         i = matrix%row(nele)
         if (g_flag(i) .eq. 0) then
            j = matrix%col(nele)         
            Au_b(i) = Au_b(i) + matrix%val(nele) * uvec(j)
         endif
         
      enddo
      
      do i = 1, nu
         
         if (g_flag(i) .eq. 0) then
            Au_b(i) = Au_b(i) - bvec(i)
         elseif (g_flag(i) .eq. 1) then
            Au_b(i) = 0d0
         endif
         
      enddo

      uvec = Au_b

      return

end subroutine res_vect
