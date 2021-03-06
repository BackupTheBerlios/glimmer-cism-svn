module glimmer_writestats_module
  !> F90 wrapper to gc_writestats
  !!
  !! \author Magnus Hagdorn
  !! \date April 2009
contains
  subroutine glimmer_writestats(resname, cfgname,wallTime,numThreads)
    use glimmer_global, only : dp
    implicit none
    character(len=*), intent(in) :: resname    !< name of the output result file
    character(len=*), intent(in) :: cfgname    !< name of configuration file
    real(kind=dp), intent(in)    :: wallTime   !< elapsed wall clock tine in seconds
    integer, intent(in)          :: numThreads !< the number of threads

    call gf_writestats(resname,cfgname,wallTime,numThreads)
  end subroutine glimmer_writestats

end module glimmer_writestats_module
