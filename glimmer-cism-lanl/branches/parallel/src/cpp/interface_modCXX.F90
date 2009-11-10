module interface_mod

      implicit none
      private
      public :: helloc_test

      interface
         subroutine hellocxx()
           end subroutine
      end interface

contains

      subroutine helloc_test()

          print *, 'Not doing anything here'
          call hellocxx()

      end subroutine helloc_test

end module interface_mod
