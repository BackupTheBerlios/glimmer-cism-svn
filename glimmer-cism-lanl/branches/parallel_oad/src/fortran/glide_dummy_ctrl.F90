module glide_dummy_ctrl
    use glimmer_global, only : dp 
    real(dp) inCtrl, outCtrl

    contains
      subroutine set_ctrl(model)
        use glide_types
        type(glide_global_type) :: model        ! model instance
        ! stub for automatic differentiation
        print *,"glide_dummy_ctrl:set_ctr; ", model%tempwk%inittemp(1,1,1), " perturbation:", inCtrl, &
             " pertutrbed:", model%tempwk%inittemp(1,1,1)+inCtrl
        model%tempwk%inittemp(1,1,1)=model%tempwk%inittemp(1,1,1)+inCtrl
      end subroutine set_ctrl

end module glide_dummy_ctrl
