!  uvec is either u^k-1 or v^k-1 on input and Av-b or Cu-d on output

subroutine res_vect ( matrix, uvec, bvec, nu, counter, g_flag, L2square, whatsparse)


use glimmer_paramets, only : dp
use glimmer_sparse_type
use glimmer_sparse
use glide_mask

implicit none

integer :: i, j, counter, nu, nele, whatsparse ! nu: size of uvec and bvec
integer, dimension(nu), intent(in) :: g_flag ! g_flag = 1 for ghost cell

type(sparse_matrix_type),  intent(in) :: matrix

real (kind = dp), dimension(nu), intent(in) :: bvec
real (kind = dp), dimension(nu), intent(inout) :: uvec
real (kind = dp), dimension(nu) :: Au_b_wig
real (kind = dp), intent(out) :: L2square
! 
real (kind = dp) :: scale_ghosts = 0.0d0

! calculate residual vector of the u OR v component

      Au_b_wig = 0d0 ! regular+ghost cells

      if (whatsparse /= STANDALONE_TRILINOS_SOLVER) then

        do nele = 1, matrix%nonzeros 

           i = matrix%row(nele)
           j = matrix%col(nele)
           Au_b_wig(i) = Au_b_wig(i) + matrix%val(nele) * uvec(j)
         
        enddo

      else 
        call matvecwithtrilinos(uvec, Au_b_wig);
      endif 
      
      do i = 1, nu
         Au_b_wig(i) = Au_b_wig(i) - bvec(i)
      enddo

      uvec = Au_b_wig

! AGS: Residual norm includes scaling to decrease importance of ghost values
! By calling it a redefinition of an inner product, it is kosher.
      L2square = 0.0
      do i = 1, nu
         if (g_flag(i) .eq. 0) then
            L2square = L2square + Au_b_wig(i) * Au_b_wig(i)
         else
            L2square = L2square + scale_ghosts * Au_b_wig(i) * Au_b_wig(i)
         endif
      end do

      return

end subroutine res_vect
