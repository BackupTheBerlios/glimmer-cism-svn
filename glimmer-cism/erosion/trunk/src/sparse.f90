! sparse.f90
! Magnus Hagdorn, March 2003
!
! sparse matrix operations

module sparse
  use glimmer_global
  ! sparse matrix type
  type sparse_matrix
     integer n
     integer, dimension(:,:), pointer :: pos => NULL()
     real(kind=dp), dimension(:), pointer :: val => NULL()
  end type sparse_matrix

  ! size of sparse matrix 
  integer, parameter, private :: chunksize=1000

contains
  subroutine new_sparse_matrix(n,mat)
    implicit none
    integer, intent(in) :: n
    type(sparse_matrix) :: mat
    
    if (.not.associated(mat%pos)) then
       allocate(mat%pos(2,n))
       allocate(mat%val(n))
    end if
    mat%n = 0
  end subroutine new_sparse_matrix

  subroutine grow_sparse_matrix(matrix)
    implicit none
    type(sparse_matrix) :: matrix

    integer, dimension(:,:), pointer :: newpos
    real(kind=dp), dimension(:), pointer :: newval
    integer oldsize

    oldsize = size(matrix%val)
    
    allocate(newpos(2,chunksize+oldsize))
    allocate(newval(chunksize+oldsize))

    newpos(:,1:oldsize) = matrix%pos(:,:)
    newval(1:oldsize) = matrix%val(:)

    deallocate(matrix%pos)
    deallocate(matrix%val)

    matrix%pos => newpos
    matrix%val => newval

  end subroutine grow_sparse_matrix

  subroutine del_sparse_matrix(matrix)
    implicit none
    type(sparse_matrix) :: matrix

    deallocate(matrix%pos)
    deallocate(matrix%val)
  end subroutine del_sparse_matrix

  subroutine print_sparse(matrix, unit)
    implicit none
    type(sparse_matrix) :: matrix
    integer, intent(in) :: unit

    integer i
    do i = 1, matrix%n
       write(unit,*) matrix%pos(1,i), matrix%pos(2,i), matrix%val(i)
    end do
  end subroutine print_sparse

  subroutine sparse_matrix_vec_prod(matrix, vec, res)
    implicit none
    type(sparse_matrix) :: matrix
    real(kind=dp), intent(in), dimension(:) :: vec
    real(kind=dp), intent(out), dimension(:) :: res

    integer i

    res = 0.
    do i=1,matrix%n
       res(matrix%pos(1,i)) = res(matrix%pos(1,i)) + vec(matrix%pos(2,i))*matrix%val(i)
    end do
  end subroutine sparse_matrix_vec_prod

  subroutine add_to_matrix(matrix, i, j, val)
    implicit none
    type(sparse_matrix) :: matrix
    integer, intent(in) :: i,j
    real(kind=dp), intent(in) :: val

    matrix%n =  matrix%n + 1
    matrix%pos(1,matrix%n) = i
    matrix%pos(2,matrix%n) = j
    matrix%val(matrix%n) = val

    if (matrix%n .eq. size(matrix%val)) then
       call grow_sparse_matrix(matrix)
    end if
  end subroutine add_to_matrix
end module sparse
