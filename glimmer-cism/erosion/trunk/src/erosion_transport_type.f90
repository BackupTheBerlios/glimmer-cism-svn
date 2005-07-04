! erosion_trans_2ndmotype.f90
! Magnus Hagdorn, June 2005
!
! module defining er_transport_type for conservation of 2nd moment advection scheme

module erosion_transport_type
  use advect_2ndmo
  use glimmer_global, only : rk
  type er_transport_type
     type(advect_type) :: mo_seds1
     type(advect_type) :: mo_seds2

     real(kind=rk), dimension(:,:), pointer :: velx, vely
  end type er_transport_type
end module erosion_transport_type
