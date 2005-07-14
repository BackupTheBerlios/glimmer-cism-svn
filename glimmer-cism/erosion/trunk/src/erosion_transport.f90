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

  end subroutine init_transport

  subroutine transport_sediments(erosion,model)
    use erosion_types
    use glide_types
    use erosion_sediment
    
    implicit none
    type(erosion_type) :: erosion          !*FD structure holding erosion data
    type(glide_global_type) :: model

    integer ns,ew


    ! transport in ice base
    call er_interpolate(erosion,model,model%velocity%ubas,erosion%seds2_vx)
    call er_interpolate(erosion,model,model%velocity%vbas,erosion%seds2_vy)
    call set_species(erosion%trans%mo_seds1,erosion%seds1)
    call advect(erosion%trans%mo_seds1,erosion%seds2_vx,erosion%seds2_vy,erosion%dt)
    call get_species(erosion%trans%mo_seds1,erosion%seds1)

    ! transport in deformable sediment layer
    if (erosion%simple_seds) then
       call er_calc_dthick(erosion,model)
       call er_interpolate(erosion,model,erosion%seds2_max_v,erosion%seds2_max)
       erosion%seds2_vx = erosion%transport_fac*erosion%seds2_vx
       erosion%seds2_vy = erosion%transport_fac*erosion%seds2_vy
    else
       call er_sediment_tstep(erosion%sediment,model)
       call er_interpolate(erosion,model,erosion%sediment%za,erosion%seds2_max)
       call er_interpolate(erosion,model,erosion%sediment%velx,erosion%seds2_vx)
       call er_interpolate(erosion,model,erosion%sediment%vely,erosion%seds2_vy)
    end if
    call set_species(erosion%trans%mo_seds2,erosion%seds2)
    call advect(erosion%trans%mo_seds2,erosion%seds2_vx,erosion%seds2_vy,erosion%dt)
    call get_species(erosion%trans%mo_seds2,erosion%seds2)
    
  end subroutine transport_sediments

  subroutine er_calc_dthick(erosion,model)
    !*FD calculate thickness of deformable sediment bed
    use glide_types
    use glide_velo
    use erosion_types
    implicit none
    type(erosion_type) :: erosion          !*FD structure holding erosion data
    type(glide_global_type) :: model       !*FD model instance

    integer ew,ns

    erosion%seds2_max_v = 0.
    do ns=1,model%general%nsn-1
       do ew=1,model%general%ewn-1
          if (abs(model%velocity%ubas(ew,ns))+abs(model%velocity%vbas(ew,ns)) .gt. 0.) then
             erosion%seds2_max_v(ew,ns) = erosion%soft_b + erosion%soft_a * &
                  sqrt(model%velocity%tau_x(ew,ns)**2 + model%velocity%tau_y(ew,ns)**2)
          end if
       end do
    end do    
  end subroutine er_calc_dthick

  subroutine er_interpolate(erosion,model,in,out)
    !*FD interpolate between glimmer velocity grid and erosion grid
    use glimmer_interpolate2d
    use glide_types
    use erosion_types
    implicit none
    type(erosion_type) :: erosion          !*FD structure holding erosion data
    type(glide_global_type) :: model       !*FD model instance
    real(kind=dp), dimension(:,:), pointer :: in, out 

     integer ew,ns

     if (erosion%grid_magnifier.gt.1) then
        call glimmer_interpolate(erosion%velo_seds,in,out)
     else
        do ns=2,model%general%nsn-1
           do ew=2,model%general%ewn-1
              out(ew,ns) = 0.25*sum(in(ew-1:ew,ns-1:ns))
           end do
        end do
     end if
   end subroutine er_interpolate

end module erosion_transport
