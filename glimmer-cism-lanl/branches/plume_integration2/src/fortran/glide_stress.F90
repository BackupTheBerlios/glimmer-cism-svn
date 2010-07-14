! *sfp* module to hold subroutines for calculation of stress components from converged, higher-order
! stress and effective viscosity fields. To be called at the end of HO vel calculation in 'run_ho_diagnostic'

module stress_hom

    use glimmer_paramets, only : dp
    use glide_types

    private
    public :: glide_stress  

    contains

    subroutine glide_stress( model )

        type(glide_global_type) :: model

!        select case( model%options%which_ho_stresscalc )    !*sfp* still need to alter glide_types and glide_setup
                                                             ! to include these options

!        case( HO_STRESSCALC_PP )

            call calcstrsstr(model%general%ewn,  model%general%nsn,  model%general%upn,     &
                       model%numerics%dew,       model%numerics%dns,                        &
                       model%numerics%sigma,     model%numerics%stagsigma,                  & 
                       model%geometry%thck,                                                 &
                       model%geomderv%dusrfdew,   model%geomderv%dusrfdns,                  &
                       model%geomderv%dthckdew,   model%geomderv%dthckdns,                  &
                       model%velocity_hom%uvel,       model%velocity_hom%vvel,              &
                       model%velocity_hom%efvs,                                             &
                       model%velocity_hom%tau%xx,      model%velocity_hom%tau%yy,           &
                       model%velocity_hom%tau%xy,      model%velocity_hom%tau%scalar,       &
                       model%velocity_hom%tau%xz,      model%velocity_hom%tau%yz )

!        case( HO_STRESSCALC_PBJ )

            !*sfp* need to fill this w/ PB&J HO stress calc. scheme and add appropriate
            ! subroutines below.

!        end select

    end subroutine glide_stress


    subroutine calcstrsstr( ewn,  nsn,  upn,  &
                       dew,        dns,       &
                       sigma,      stagsigma, &  
                       thck,                  &
                       dusrfdew,   dusrfdns,  &
                       dthckdew,   dthckdns,  &
                       uvel,       vvel,      &
                       efvs,                  &
                       tauxx,      tauyy,     &
                       tauxy,      tau,       &
                       tauxz,      tauyz )

        use glimmer_paramets, only : tau0_glam

        implicit none

        integer, intent(in) :: ewn, nsn, upn

        real (kind = dp), intent(in) :: dew, dns 
        real (kind = dp), intent(in), dimension(:)     :: sigma, stagsigma
        real (kind = dp), intent(in), dimension(:,:,:) :: efvs, uvel, vvel
        real (kind = dp), intent(in), dimension(:,:) :: thck, dusrfdew, &
                                                dusrfdns, dthckdew, dthckdns

        real (kind = dp), intent(out), dimension(:,:,:) :: tauxx, tauyy, tauxy, &
                                                   tauxz, tauyz, tau
        !*sfp* local vars
        integer :: ew, ns, up
        real (kind = dp), parameter :: f1 = len0 / thk0
        integer, dimension(2) :: mew, mns
        real (kind = dp) :: dew2, dew4, dns2, dns4
        !real (kind = dp), dimension(upn-1) :: dup, dupm        
        real (kind = dp) :: dup, dupm        

        !*sfp* note that these are already defined and used in glam_strs2. If needed by PB&J 
        ! stress calc routine as well, consider moving the up-scope 
        !dup = (/ ( (sigma(up+1)-sigma(up)), up = 1, upn-1) /)
        dup = sigma(2)-sigma(1)
        dupm = - 0.25_dp / dup
        dew2 = 2.0_dp * dew; dns2 = 2.0_dp * dns        ! *sp* 2x the standard grid spacing
        dew4 = 4.0_dp * dew; dns4 = 4.0_dp * dns        ! *sp* 4x the standard grid spacing

        tauxz = 0.0_dp
        tauyz = 0.0_dp
        tauxx = 0.0_dp
        tauyy = 0.0_dp
        tauxy = 0.0_dp

        do ns = 2,nsn-1
            do ew = 2,ewn-1;

            if (thck(ew,ns) > 0.0_dp) then
                tauxz(1:upn-1,ew,ns) = vertideriv(upn, hsum(uvel(:,ew-1:ew,ns-1:ns)), thck(ew,ns), dupm)
                tauyz(1:upn-1,ew,ns) = vertideriv(upn, hsum(vvel(:,ew-1:ew,ns-1:ns)), thck(ew,ns), dupm)
                tauxx(1:upn-1,ew,ns) = horizderiv(upn,  stagsigma,                &
                              0.5_dp * sum(uvel(:,ew-1:ew,ns-1:ns),3), &         !average over n-s neighbours
                              dew, tauxz(:,ew,ns),           &
                              0.25_dp*sum(dusrfdew(ew-1:ew,ns-1:ns)), &           !average over 4 neighbours
                              0.25_dp*sum(dthckdew(ew-1:ew,ns-1:ns)))             !average over 4 neighbours

                tauyy(1:upn-1,ew,ns) = horizderiv(upn,  stagsigma,                &
                              0.5_dp * sum(vvel(:,ew-1:ew,ns-1:ns),2), &        !average over e-w neighbours
                              dns, tauyz(:,ew,ns),           &
                              0.25_dp * sum(dusrfdns(ew-1:ew,ns-1:ns)), &       !average over 4 neighbours
                              0.25_dp * sum(dthckdns(ew-1:ew,ns-1:ns)))         !average over 4 neighbours

                tauxy(1:upn-1,ew,ns) = horizderiv(upn,  stagsigma,      &     ! du/dy
                              0.5_dp * sum(uvel(:,ew-1:ew,ns-1:ns),2), &       !average over e-w neighbours
                              dns, tauxz(:,ew,ns),           &
                              0.25_dp * sum(dusrfdns(ew-1:ew,ns-1:ns)), &      !average over 4 neighbours
                              0.25_dp * sum(dthckdns(ew-1:ew,ns-1:ns)))  &     !average over 4 neighbours
                            +                                          &
                              horizderiv(upn,  stagsigma,                &     ! dv/dx
                              0.5_dp * sum(vvel(:,ew-1:ew,ns-1:ns),3), &       !average over n-s neighbours
                              dew, tauyz(:,ew,ns),           &
                              0.25_dp * sum(dusrfdew(ew-1:ew,ns-1:ns)), &      !average over 4 neighbours
                              0.25_dp * sum(dthckdew(ew-1:ew,ns-1:ns)))        !average over 4 neighbours

            end if

            end do
        end do

        tauxz = f1 * efvs * tauxz     !* tau0_glam
        tauyz = f1 * efvs * tauyz     !* tau0_glam
        tauxx = 2.0_dp * efvs * tauxx !* tau0_glam
        tauyy = 2.0_dp * efvs * tauyy !* tau0_glam
        tauxy = efvs * tauxy          !* tau0_glam

        !*sfp* expanding this in terms of viscosity and velocity gradients, I've confirmed that 
        ! one gets the same thing as if one took Tau_eff = N_eff * Eps_eff, where Eps_eff is the 
        ! 1st order approx. to the 2nd strain-rate invariant (outlined in model description document).

        tau = sqrt(tauxz**2 + tauyz**2 + tauxx**2 + tauyy**2 + tauxx*tauyy + tauxy**2) !* tau0_glam

        return

    end subroutine calcstrsstr


    function vertideriv( upn, varb, thck, dupm )

        implicit none

        integer, intent(in) :: upn
        real (kind = dp), intent(in), dimension(:) :: varb
        real (kind = dp), intent(in) :: thck
        real (kind = dp), intent(in) :: dupm            

        real (kind = dp), dimension(size(varb)-1) :: vertideriv

        !*sfp* 'dupm' defined as -1/(2*del_sigma), in which case it seems like 
        !there should be a '-' in front of this expression ... or, negative sign
        !may be implicit in the vert indices ( "arb(2:upn) - varb(1:upn-1)" ) and
        !the fact that up=1 at the sfc and up=upn at the bed ??? 

        vertideriv(1:upn-1) = dupm * (varb(2:upn) - varb(1:upn-1)) / thck

        return

   end function vertideriv



    function horizderiv( upn,     stagsigma,   &
                         varb,    dhoriz,        &
                         dvarbdz, dusrfdx, dthckdx)

        implicit none

        integer, intent(in) :: upn
        real (kind = dp), dimension(:), intent(in) :: stagsigma
        real (kind = dp), dimension(:,:), intent(in) :: varb
        real (kind = dp), dimension(:), intent(in) :: dvarbdz
        real (kind = dp), intent(in) :: dusrfdx, dthckdx,dhoriz

        real (kind = dp) :: horizderiv(size(varb,1)-1)

        ! *sfp* where does this factor of 1/4 come from ... averaging? 

        ! *cvg* factor of 1/4 has been moved up into callers of this function, since it was 
        ! indeed used to average the four values of dusrfdx and dthckdx

        ! *cvg* This function uses the formulae: 
        !
        !       dv/dx = (del v/ del x) - (del v/del z)*((del h/del x) - sigma*(del H/del x))
        !       dv/dy = similar
        !
        !       where: v = v(x,sigma) is any 3-dimensional variable
        !              h = h(x,y) is the upper surface
        !              H = H(x,y) is the ice thickness

        horizderiv = (varb(1:upn-1,2) + varb(2:upn,2) - varb(1:upn-1,1) - varb(2:upn,1)) / (2.0_dp*dhoriz) - &
                   dvarbdz * (dusrfdx - stagsigma * dthckdx) 

        return

   end function horizderiv


   function hsum(inp)

      implicit none

      real (kind = dp), dimension(:,:,:), intent(in) :: inp
      real (kind = dp), dimension(size(inp,dim=1)) :: hsum

      hsum = sum(sum(inp(:,:,:),dim=3),dim=2)

      return

   end function hsum

end module stress_hom
    