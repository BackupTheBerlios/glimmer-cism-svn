module sparse_trilinos

    ! RN_20091006: This module contains three subroutines forming a Fortran
    ! interface to Trilinos solvers.
    use sparse_type
    implicit none

contains

    subroutine trilinos_solver_preprocess(matrix)
	integer :: i
	type(sparse_matrix_type) :: matrix

	write(*,*) '--------------------------------------'
	write(*,*) 'IN PREPROCESS'
	
	! Move from 1-based to 0-based indexing
	matrix%row = matrix%row - 1
	matrix%col = matrix%col - 1
	
	! For debugging =)
!	do 10 i = 1, matrix%nnz
!	    write(*,*) matrix%row(i), matrix%col(i), matrix%val(i)
!10	enddo
	write(*,*) '--------------------------------------'

    end subroutine trilinos_solver_preprocess


!    subroutine trilinos_solve(matrix, rhs, solution, err, niters)
    subroutine trilinos_solve(matrix, rhs, solution)
	integer :: i
!	integer :: niters
!	double precision :: err
	double precision, pointer :: rhs(:), solution(:)
	type(sparse_matrix_type) :: matrix

	write(*,*) '--------------------------------------'
	write(*,*) 'IN SOLVE'

	! For debugging =)
!	do 20 i = 1, matrix%order
!	    write(*,*) rhs(i)
!20	enddo

	! Call to Trilinos solvers
!	call solve(matrix%nnz, matrix%order, matrix%row, matrix%col, &
!                   matrix%val, rhs, solution, err, niters)
	call solve(matrix%nnz, matrix%order, matrix%row, matrix%col, &
                   matrix%val, rhs, solution)

	! For degugging =)
!	write(*,*) 'SOLUTION RETURNED'
!	do 30 i = 1, matrix%order
!	    write(*,*) solution(i)
!30	enddo
	write(*,*) '--------------------------------------'

    end subroutine trilinos_solve


    subroutine trilinos_solver_postprocess(matrix)
	integer :: i
	type(sparse_matrix_type) :: matrix

	write(*,*) '--------------------------------------'
	write(*,*) 'IN POSTPROCESS'

	! Move back from 0-based to 1-based indexing
	matrix%row = matrix%row + 1
	matrix%col = matrix%col + 1

	! For debugging =)
!	do 40 i = 1, matrix%nnz
!	    write(*,*) matrix%row(i), matrix%col(i), matrix%val(i)
!40	enddo
	write(*,*) '--------------------------------------'

    end subroutine trilinos_solver_postprocess

end module sparse_trilinos
