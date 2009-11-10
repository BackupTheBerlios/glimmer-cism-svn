module sparse_type

    ! RN_20091006: A simple sparse matrix type	
    type sparse_matrix_type
	integer :: nnz
	integer :: order
	integer, dimension(:), pointer :: col => NULL()
	integer, dimension(:), pointer :: row => NULL()
	double precision, dimension(:), pointer :: val => NULL()
    end type sparse_matrix_type

end module sparse_type
