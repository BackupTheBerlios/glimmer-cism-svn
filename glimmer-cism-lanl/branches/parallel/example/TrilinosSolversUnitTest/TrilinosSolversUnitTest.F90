program trilinossolversunittest
    ! RN_20091106: The following is a simple unit test that shows how to
    ! interface with Trilinos solvers.
    use sparse_type
    use sparse_trilinos
    implicit none

    integer :: i, nmax, mtxsize, nnz
    parameter (nmax = 1000, mtxsize = 200)
    double precision, target :: val(nmax)
    integer, target :: col(nmax), row(nmax)
    double precision, target :: firstarray(nmax), secondarray(nmax)

    type(sparse_matrix_type) :: matrix
    double precision, pointer :: rhs(:), solution(:)
    !integer :: niters
    !double precision :: err

    ! Create a sparse matrix of the following pattern, where * denotes a
    ! nonzero.
    !
    ! |* 0 0 0 0 0 0 0 0 0 0|
    ! |0 * * 0 0 0 0 0 0 0 0|
    ! |0 0 * 0 * 0 0 0 0 0 0|
    ! |0 0 0 * 0 0 * 0 0 0 0|
    ! |0 0 0 0 * 0 0 0 * 0 0|
    ! |0 0 0 0 0 * 0 0 0 0 *|
    ! |0 0 0 0 0 0 * 0 0 0 0|
    ! |0 0 0 0 0 0 0 * 0 0 0|
    ! |0 0 0 0 0 0 0 0 * 0 0|
    ! |0 0 0 0 0 0 0 0 0 * 0|
    ! |0 0 0 0 0 0 0 0 0 0 *|

    ! Compute the number of nonzero entries
    if (mod(mtxsize,2) .EQ. 0) then
	nnz = mtxsize + mtxsize / 2 - 1
    else
	nnz = mtxsize + (mtxsize - 1) / 2
    endif

    ! Fill in the arrays for the matrix
    ! These values are of no particular reason...
    do 10 i = 1, nnz
	if (i .LT. mtxsize + 1) then
	    row(i) = i
	    col(i) = i
	    val(i) = 2*i/3 - 1.5/i - 1
	else
	    row(i) = i - mtxsize + 1
	    col(i) = i - mtxsize + row(i)
	    val(i) = 2*i/3 - 1.5/i - 1
	endif
10  enddo

    ! Form the matrix
    matrix%nnz = nnz
    matrix%order = mtxsize
    matrix%row => row
    matrix%col => col
    matrix%val => val

    ! For debugging =)
!    write(*,*) 'Number of nonzeros: ', matrix%nnz
!    write(*,*) 'Size of the matrix: ', matrix%order
!    write(*,*) 'Data of the matrix:'
!    do 20 i = 1, nnz
!	write(*,*) matrix%row(i), matrix%col(i), matrix%val(i)
!20  enddo

    ! Create an array for the rhs
    ! These values are of no particular reason...
    do 30 i = 1, mtxsize
	firstarray(i) = 5 / i - 0.33 * i + i * i
30  enddo

    ! Form the rhs
    rhs => firstarray

    ! For debugging =)
!    write(*,*) 'RHS:'
!    do 40 i = 1, mtxsize
!	write(*,*) rhs(i)
!40  enddo

    ! Create an array for the solution -- not necessary!
    do 50 i = 1, mtxsize
	secondarray(i) = 0
50  enddo

    ! Form the rhs
    solution => secondarray

    ! For debugging =)
!    write(*,*) 'Solution:'
!    do 60 i = 1, mtxsize
!	write(*,*) solution(i)
!60  enddo

    write(*,*) '======================================'
    write(*,*) 'START SOLVING'
    write(*,*) '======================================'


    ! Preprocess the matrix before calling to Trilinos solvers
    call trilinos_solver_preprocess(matrix)

    ! Make a call to Trilinos solvers
!    call trilinos_solve(matrix, rhs, solution, err, niters)
    call trilinos_solve(matrix, rhs, solution)

    ! Postprocess the matrix after solving the system
    call trilinos_solver_postprocess(matrix)

    write(*,*) '======================================'
    write(*,*) 'DONE'
    write(*,*) '======================================'

stop
end
