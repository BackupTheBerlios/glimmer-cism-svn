! erosion_transport_2ndmo.f90
! Magnus Hagdorn, May 2005
!
! this module transports some scalar quantity c through a 2D velo field using
! the conservation of 2nd order moments algorithm

module erosion_transport
  
contains
  subroutine init_transport(trans,model,erosion)
    use glide_types
    use erosion_types
    implicit none
    type(er_transport_type) :: trans       ! structure holding transport stuff
    type(glide_global_type) :: model       ! model instance
    type(erosion_type) :: erosion          !*FD structure holding erosion data

    call init_advect(trans%mo_seds1,erosion%seds1,erosion%dew,erosion%dns)
    call init_advect(trans%mo_seds2,erosion%seds2,erosion%dew,erosion%dns)

    ! allocate memory for velocities
    call coordsystem_allocate(erosion%coord, trans%velx)
    call coordsystem_allocate(erosion%coord, trans%vely)

  end subroutine init_transport

  subroutine transport_sediments(erosion,model)
    use erosion_types
    use glide_types
    use glimmer_interpolate2d
    implicit none
    type(erosion_type) :: erosion          !*FD structure holding erosion data
    type(glide_global_type) :: model

    integer ns,ew

    ! set velocities
    if (erosion%grid_magnifier.gt.1) then
       call glimmer_interpolate(erosion%velo_seds,model%velocity%ubas,erosion%trans%velx)
       call glimmer_interpolate(erosion%velo_seds,model%velocity%vbas,erosion%trans%vely)       
    else
       do ns=2,model%general%nsn-1
          do ew=2,model%general%ewn-1
             erosion%trans%velx(ew,ns) = 0.25*sum(model%velocity%ubas(ew-1:ew,ns-1:ns))
             erosion%trans%vely(ew,ns) = 0.25*sum(model%velocity%vbas(ew-1:ew,ns-1:ns))
          end do
       end do
    end if

    ! transport in ice base
    call set_species(erosion%trans%mo_seds1,erosion%seds1)
    call advect(erosion%trans%mo_seds1,erosion%trans%velx,erosion%trans%vely,erosion%dt)
    call get_species(erosion%trans%mo_seds1,erosion%seds1)

    ! transport in deformable sediment layer
    erosion%trans%velx = erosion%transport_fac*erosion%trans%velx
    erosion%trans%vely = erosion%transport_fac*erosion%trans%vely
    call set_species(erosion%trans%mo_seds2,erosion%seds2)
    call advect(erosion%trans%mo_seds2,erosion%trans%velx,erosion%trans%vely,erosion%dt)
    call get_species(erosion%trans%mo_seds2,erosion%seds2)
    
  end subroutine transport_sediments

end module erosion_transport
