subroutine scam_use_iop_srf( cam_in )
!-----------------------------------------------------------------------
    use ppgrid,           only: begchunk, endchunk
    use camsrfexch,       only: cam_in_t    
    use physconst,   only: stebol, latvap
    use scamMod
    use control_mod, only : scm_iop_srf_prop_se

    implicit none
    save

    type(cam_in_t), intent(INOUT) :: cam_in(begchunk:endchunk)
    ! local
    integer :: c    ! Chunk index
    integer :: ncol ! Number of columns
    !
    ! Replace surface fluxes with observed values for IOP forcing if
    ! requested by switch settings in the GUI
    !

!    do c=begchunk,endchunk
!      ncol = cam_in(c)%ncol
!      
!      cam_in(c)%shf(:ncol) = 15.0
!      cam_in(c)%lhf(:ncol) = 115.0
!      cam_in(c)%cflx(:ncol,1) = 115.0/latvap
!      cam_in(c)%ts(:ncol) = 292.5
!      cam_in(c)%lwup(:ncol) = stebol*(292.5)**4
!    end do
    
    if (scm_iop_srf_prop) then
       do c=begchunk,endchunk
          ncol = cam_in(c)%ncol
          if(have_lhflx) then
             cam_in(c)%lhf(:ncol) = lhflxobs(1)
             cam_in(c)%cflx(:ncol,1) = lhflxobs(1)/latvap
          endif
          if(have_shflx) cam_in(c)%shf(:ncol) = shflxobs(1)
          if(have_tg) then
             cam_in(c)%ts(:ncol) = tground(1)
             cam_in(c)%lwup(:ncol) = stebol * tground(1)**4
          endif
       end do
    endif

end subroutine scam_use_iop_srf
