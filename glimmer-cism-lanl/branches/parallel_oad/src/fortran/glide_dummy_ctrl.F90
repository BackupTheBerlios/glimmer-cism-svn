module glide_dummy_ctrl
    use glimmer_global, only : dp 
    real(dp) inCtrl, outCtrl

    contains
      subroutine set_ctrl(model)
        use glide_types
        type(glide_global_type) :: model        ! model instance
        ! stub for automatic differentiation
        print *,"glide_dummy_ctrl:set_ctr original value  :",  model%temper%temp(5,15,15)
        print *,"glide_dummy_ctrl:set_ctr perturbation abs:",  inCtrl, " rel:", inCtrl/model%temper%temp(5,15,15)
        model%temper%temp(5,15,15)= model%temper%temp(5,15,15)+inCtrl
        print *,"glide_dummy_ctrl:set_ctr perturbed value : ", model%temper%temp(5,15,15)
      end subroutine set_ctrl

end module glide_dummy_ctrl
