program iswplume

  use plume
  use plume_global
  use plume_io

  implicit none

  integer::iarg_count,istep
  integer::command_argument_count
  character(len=128) :: jobid
  character(len=128) :: output_dir,nl_filename,supp_ascii_output_str,supp_log_str
  logical :: supp_ascii_output = .true.
  logical :: supp_logging = .false.

  ! *******************************************************************
  ! *** process command-line arguments ********************************
  ! *******************************************************************

  iarg_count = command_argument_count()

  if (iarg_count > 0) then
     call get_command_argument(1,jobid)
     if (trim(jobid) == "-h") then
        write(*,*)'Usage: plume <jobid> <namelist_file> [<output_dir>]  [<no logging>] [<no ascii output>] '
        stop
     end if
  else
     !get job id
     write(*,*) 'please enter (3-character) job id:'
     read(*,*) jobid
  end if
  write(*,*) 'jobid: ', trim(jobid)

  if (iarg_count > 1) then
     call get_command_argument(2,nl_filename)
  else
     nl_filename = 'plume.nl'
  end if
  write(*,*) 'using namelist file: ', trim(nl_filename)

  if (iarg_count > 2) then
     call get_command_argument(3,output_dir)
  else
     output_dir = './'
  end if
  write(*,*) 'output_dir:', trim(output_dir)

  if (iarg_count > 3) then
     call get_command_argument(4,supp_log_str)
     read(supp_log_str,'(l)') supp_logging
     if (supp_logging) write(*,*) 'Suppressing all screen/file logging for this run'
  end if

  if (iarg_count > 4) then
     call get_command_argument(5,supp_ascii_output_str)
     read(supp_ascii_output_str,'(l)') supp_ascii_output
     if (supp_ascii_output) write(*,*) 'Suppressing old-style ascii output and screen output for this run'
  end if

  call plume_logging_initialize(trim(output_dir),&
       trim(jobid),&
       supp_logging)
  call plume_initialise(trim(nl_filename),&
       supp_ascii_output,supp_logging)
  call plume_io_initialize(supp_ascii_output, &
       trim(output_dir) // '/' // trim(jobid) // '.nc')

  ! main loop 

  do istep = 1,get_nsteps()

     call plume_runstep()

  end do

  call plume_finalise()
  call plume_io_finalize()
  call plume_logging_finalize()

end program iswplume
