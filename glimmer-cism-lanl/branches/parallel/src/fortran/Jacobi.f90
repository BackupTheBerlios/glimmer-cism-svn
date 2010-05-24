!****************************************************************************
!     performs p iterations on Pdu=rhs !rhs = wk1, du = wk2
!****************************************************************************

! JFL to be removed

subroutine Jacobi (matrix, rhs, dutp, nu)

use glimmer_paramets, only : dp
use glimmer_sparse_type

  implicit none
      
  integer :: i, j, nele, nu, p, Npite

  type(sparse_matrix_type),  intent(in) :: matrix

  real (kind = dp), dimension(nu), intent(in) :: rhs
  real (kind = dp), dimension(nu), intent(inout) :: dutp
  real (kind = dp), dimension(nu) :: Pdu, Diag

  Npite = 10 ! nb of precond iterations
  Pdu = 0d0

  do p = 1, Npite

     do nele = 1, matrix%nonzeros 

        i = matrix%row(nele)
        j = matrix%col(nele)

        Pdu(i) = Pdu(i) + matrix%val(nele) * dutp(j)

        if (i .eq. j) Diag(i) = matrix%val(nele)

!        if (Diag(i) .lt. 1d-12) print *, Diag(i), rhs(i),'watchout dude'

     enddo

     do i = 1, nu

        dutp(i) = (rhs(i)-Pdu(i)) / Diag(i)
        print *, p, i, dutp(i), rhs(i), Pdu(i), Diag(i)

     enddo

  enddo

  return
end subroutine Jacobi






