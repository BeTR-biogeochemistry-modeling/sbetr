module betr_ctrl
!
!logical switches for global control of betr functions
implicit none


  logical, public            :: use_pH_data = .false.
  logical, public            :: is_active_betr_bgc = .false.
  logical, public            :: lglb_diffusion_on = .true.
  logical, public            :: lglb_advection_on = .true.
  logical, public            :: lglb_reaction_on  = .true.
end module betr_ctrl
