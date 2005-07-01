! erosion_trans_lagtype.f90
! Magnus Hagdorn, June 2005
!
! module defining er_transport_type for lagrangian type advection schemes

module erosion_transport_type
  
  use geometry
  use glimmer_global, only : dp
  use glimmer_sparse
  
  type er_transport_type
     ! private data
     real(kind=dp), dimension(:), pointer :: lin_stuff,lin_stuff2,lin_con
     real(kind=dp) :: half_xstep, half_ystep
     ! for finite volume
     type(coord_point), dimension(:,:), pointer :: patch_strip 
     type(geom_poly) :: patch, patch1, patch2
     ! for interpolation
     real(kind=dp), dimension(:,:), pointer :: dispx => NULL() ! x-displacement field
     real(kind=dp), dimension(:,:), pointer :: dispy => NULL() ! y-displacement field

     type(sparse_matrix_type) :: lag_seds1                 !*FD sparse matrix holding dirty ice layer
     type(sparse_matrix_type) :: lag_seds2                 !*FD sparse matrix holding deformable sediment layer

  end type er_transport_type
end module erosion_transport_type
