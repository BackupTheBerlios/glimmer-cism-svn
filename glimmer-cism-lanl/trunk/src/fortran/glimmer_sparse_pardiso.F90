module glimmer_sparse_pardiso
    !*FD This module builds on the glimmer_sparse module to provide an 'easy'
    !*FD interface to PARDISO.
    
    use glimmer_sparse_type
    use glimmer_global, only: dp
    implicit none

#ifdef PARDISO_64BIT
        integer, parameter :: parint = 8       
#else
        integer, parameter :: parint = 4
#endif


    type pardiso_solver_workspace
        !*FD Memory addresses used internally by PARDISO
        integer(kind=parint), dimension(64) :: pt
        !*FD Storage for passing parameters to and getting info from PARDISO
        integer(kind=parint), dimension(64) :: iparm
        !*FD Error messages from pardisoinit: 0 no error, -10 no license, -11
        !*FD licence expired (?), -12 wrong user or hostname.
        integer(kind=parint) :: error
        !*FD Work array for multi-recursive solver only
        real(kind=dp) :: dparm(64)
    end type pardiso_solver_workspace

    type pardiso_solver_options
        !*FD Matrix type, we are generally real, non-symmetric (11)
        integer(kind=parint) :: mtype
        !*FD Solver, sparse direct is 0. If we ever have symmetric indefinate
        !*FD matrices, we might try 1; 'multi-recusive iterative'
        integer(kind=parint) :: solver
        real(kind=dp) :: tolerance
    end type pardiso_solver_options

contains
    subroutine pardiso_solver_default_options(opt)
        !*FD Populates a sparse_solver_options (defined above) with default
        !*FD options.  This is necessary because different solvers may define
        !*FD different options beyond the required fields defined above.
        !*FD Filling them in this function allows client code to pick "good"
        !*FD values in a generic way.
        implicit none
        type(pardiso_solver_options), intent(inout) :: opt
        opt%mtype = 11
        opt%solver = 0
    end subroutine pardiso_solver_default_options

    subroutine pardiso_allocate_workspace(matrix, options, workspace, max_nonzeros_arg)
        !*FD Allocate solver workspace.  This needs to be done once
        !*FD (when the maximum number of nonzero entries is first known)
        !*FD Note that the max_nonzeros argument must be optional, and if
        !*FD the current number of nonzeroes must be used.
        implicit none
        type(sparse_matrix_type) :: matrix
        type(pardiso_solver_options) :: options
        type(pardiso_solver_workspace) :: workspace
        integer, optional :: max_nonzeros_arg
        !...Check license of the solver and initialize the solver
        call pardisoinit(workspace%pt, options%mtype, options%solver,&
                         workspace%iparm, workspace%dparm, workspace%error)
    
    end subroutine pardiso_allocate_workspace

    subroutine pardiso_solver_preprocess(matrix, options, workspace)
        !*FD Performs any preprocessing needed to be performed on the sparse
        !*FD matrix.  Workspace must have already been allocated. 
        !*FD This function should be safe to call more than once.
        !*FD 
        !*FD It is an error to call this function on a workspace without already
        !*FD allocated memory.
        !*FD
        !*FD In general sparse_allocate_workspace should perform any actions
        !*FD that depend on the *size* of the sparse matrix, and
        !*FD sprase_solver_preprocess should perform any actions that depend
        !*FD upon the *contents* of the sparse matrix.
        type(sparse_matrix_type) :: matrix
        type(pardiso_solver_options) :: options
        type(pardiso_solver_workspace) :: workspace

        ! PARDISO is looking for a compressed row format for the matrix.
        call to_row_format(matrix)

        !Sort by row within each column in the column format.
        !to_column_format does not necessarily do this
        call sort_row_format(matrix)

    end subroutine pardiso_solver_preprocess

    function pardiso_solve(matrix, rhs, solution, options, workspace,err,niters, verbose)
        !*FD Solves the sparse linear system, and reports status information.
        !*FD This function returns an error code that should be zero if the
        !*FD call succeeded and nonzero if it failed.  No additional error codes
        !*FD are defined.  Although this function reports back the final error
        !*FD and the number of iterations needed to converge, these should *not*
        !*FD be relied upon as not every sparse linear solver may report them.
        implicit none
        type(sparse_matrix_type), intent(inout) :: matrix 
        !*FD Sparse matrix to solve.  This is inout because the sparse solver
        !*FD may have to do some re-arranging of the matrix.
        
        real(kind=dp), dimension(:), intent(in) :: rhs 
        !*FD Right hand side of the solution vector
        
        real(kind=dp), dimension(:), intent(inout) :: solution 
        !*FD Solution vector, containing an initial guess.

        type(pardiso_solver_options), intent(in) :: options
        !*FD Options such as convergence criteria
        
        type(pardiso_solver_workspace), intent(inout) :: workspace
        !*FD Internal solver workspace
        
        real(kind=dp), intent(out) :: err
        !*FD Final solution error
        
        integer, intent(out) :: niters
        !*FD Number of iterations required to reach the solution

        logical, optional, intent(in) :: verbose
        !*FD If present and true, this argument may cause diagnostic information
        !*FD to be printed by the solver (not every solver may implement this).
        
        integer :: sparse_solve
        integer :: pardiso_solve

        integer :: iunit !Unit number to print verbose output to (6=stdout, 0=no output)
        integer :: mtype

        iunit=0
        if (present(verbose)) then
            if(verbose) then
                iunit=1
                write(*,*),"Tolerance=",options%tolerance
            end if
        end if


        !Detect symmetric matrix
        if (matrix%symmetric) then
            mtype = 1
        else
            mtype = 11
        end if
        
        !PARDISO arguments
        write(*,*) "PARDISO Begin!"

        !CALL PARDISO(PT, MAXFCT, MNUM, MTYPE, PHASE, N, A, IA, JA,
        !PERM, NRHS, IPARM, MSGLVL, B, X, ERROR, DPARM)

        call pardiso(workspace%pt, 1, 1, mtype, 13, matrix%order, matrix%val, &
                     matrix%col, matrix%row, 0, matrix%order, workspace%iparm, &
                     iunit, rhs, solution, sparse_solve)
        write(*,*) "PARDISO End!"

        pardiso_solve = sparse_solve

        err = 0
        niters = 1
    end function pardiso_solve

    subroutine pardiso_solver_postprocess(matrix, options, workspace)
        type(sparse_matrix_type) :: matrix
        type(pardiso_solver_options) :: options
        type(pardiso_solver_workspace) :: workspace
    end subroutine

    subroutine pardiso_destroy_workspace(matrix, options, workspace)
        !*FD Deallocates all working memory for the sparse linear solver.
        !*FD This need *not* be safe to call of an unallocated workspace
        !*FD No sparse solver should call this automatically.
        type(sparse_matrix_type) :: matrix
        type(pardiso_solver_options) :: options
        type(pardiso_solver_workspace) :: workspace

        !TODO: Call pardiso with flag to delete

    end subroutine pardiso_destroy_workspace

    subroutine pardiso_interpret_error(error_code, error_string)
        !*FD takes an error code output from sparse_solve and interprets it.
        !*FD error_string must be an optional argument.
        !*FD If it is not provided, the error is printed to standard out
        !*FD instead of being put in the string
        integer :: error_code
        character(*), optional, intent(out) :: error_string
        character(256) :: tmp_error_string
        
        select case (error_code)
            case (0)
                tmp_error_string = "No error"
            case (-1)
                tmp_error_string = "Input inconsistent"
            case (-2)
                tmp_error_string = "Not enough memory"
            case (-3)
                tmp_error_string = "Reordering problem"
            case (-4)
                tmp_error_string = "Zero pivot, numerical fact. or iterative refinement problem"
            case (-5)
                tmp_error_string = "Unclassified (internal) errror"
            case (-6)
                tmp_error_string = "Preordering failed (matrix types 11, 13 only)"
            case (-7)
                tmp_error_string = "Diagonal matrix problem"
            case (-8)
                tmp_error_string = "32-bit integer overflow problem"
        end select

        if (present(error_string)) then
            error_string = tmp_error_string
        else
            write(*,*) tmp_error_string
        endif
    end subroutine pardiso_interpret_error
end module glimmer_sparse_pardiso
