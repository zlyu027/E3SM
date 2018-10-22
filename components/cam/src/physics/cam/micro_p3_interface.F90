module micro_p3_interface

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!
  !! Interface between E3SM and P3 microphysics
  !!
  !! Author: Peter Caldwell
  !!
  !! Last updated: 2018-09-12
  !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  use shr_kind_mod,   only: r8=>shr_kind_r8
  use ppgrid,         only: pcols,pver,pverp

!comment: I think Kai added handle_errmsg. It would be better to 
!use standard E3SM libraries if possible.
  use error_messages, only: handle_errmsg

  use physics_types,  only: physics_state, &
                            physics_ptend, &
                            physics_ptend_init
  use physconst,      only: mwdry, cpair, mwh2o
  use constituents,   only: cnst_add, pcnst, sflxnam, apcnst, bpcnst, pcnst,&
                            cnst_name, cnst_get_ind,cnst_longname
  use physics_buffer, only: physics_buffer_desc, dtype_r8, col_type_subcol, &
                            pbuf_get_field, pbuf_add_field,dyn_time_lvls,dtype_i4, &
                            pbuf_register_subcol
  use ref_pres,       only: top_lev=>trop_cloud_top_lev
  use phys_control,   only: phys_getopts
!  use cam_debug,      only: l_debug, l_summary_debug
  use subcol_utils,   only: subcol_get_scheme
  use cam_abortutils, only: endrun
  use spmd_utils,     only: masterproc
  use cam_logfile,    only: iulog
       
  implicit none

  public :: micro_p3_init, micro_p3_register, micro_p3_tend, &
            micro_p3_init_cnst, micro_p3_implements_cnst

  private

  !Define indices for state%q constituents at module level so
  !defining them in micro_p3_register makes them permanently 
  !available.
  logical :: use_subcol_microp  ! If true, then are using subcolumns in microphysics

  integer, public ::    &
       ixcldliq = -1,   & ! cloud liquid amount index
       ixcldice = -1,      & ! ice index
       ixnumliq = -1,   & ! cloud liquid number index
       ixnumice = -1,   & ! cloud ice number index
       ixrain   = -1,   & ! rain index
       ixnumrain= -1,   & ! rain number index
       ixcldrim = -1,      & ! rime index ??
       ixrimvol  = -1,  & ! rime volume index ??
       ixqirim  = -1      ! ?? index ??

!! pbuf 
   integer :: &
      cldo_idx,           &
      qme_idx,            &
      prain_idx,          &
      nevapr_idx,         &
      rate1_cw2pr_st_idx, &
      dei_idx,            &
      mu_idx,             &
      lambdac_idx,        &
      rei_idx,            &
      rel_idx,            &
      ls_flxprc_idx,      &
      ls_flxsnw_idx,      &
      ls_reffrain_idx,    &
      ls_reffsnow_idx,    &
      cv_reffliq_idx,     &
      cv_reffice_idx,     &
      prer_evap_idx,      &
      cmeliq_idx,         &
      relvar_idx,         &
      accre_enhan_idx,    &
      iciwpst_idx,        &
      iclwpst_idx,        &
      cc_t_idx,           &
      cc_qv_idx,          &
      cc_ql_idx,          &
      cc_qi_idx,          &
      cc_nl_idx,          &
      cc_ni_idx,          &
      cc_qlst_idx
   ! Index fields for precipitation efficiency.
   integer :: &
       acpr_idx = -1, &
       acgcme_idx = -1, &
       acnum_idx = -1

! Physics buffer indices for fields registered by other modules
   integer :: &
      ast_idx = -1,            &
      cld_idx = -1,            &
      concld_idx = -1
! Pbuf fields needed for subcol_SILHS
   integer :: &
      qrain_idx=-1, &
      nrain_idx=-1

   integer :: &
      naai_idx = -1,           &
      naai_hom_idx = -1,       &
      npccn_idx = -1,          &
      rndst_idx = -1,          &
      nacon_idx = -1,          &
      prec_str_idx = -1,       &
      prec_pcw_idx = -1,       &
      prec_sed_idx = -1,       &
      snow_str_idx = -1,       &
      snow_pcw_idx = -1,       &
      snow_sed_idx = -1

! pbuf fields for heterogeneous freezing
   integer :: &
      frzimm_idx = -1, &
      frzcnt_idx = -1, &
     frzdep_idx = -1

   logical :: &
      allow_sed_supersat  ! allow supersaturated conditions after sedimentation loop



   real(r8) :: &
      micro_mg_accre_enhan_fac = huge(1.0_r8), & !Accretion enhancement factor from namelist
      prc_coef1_in             = huge(1.0_r8), &
      prc_exp_in               = huge(1.0_r8), &
      prc_exp1_in              = huge(1.0_r8), &
      cld_sed_in               = huge(1.0_r8), & !scale fac for cloud sedimentation velocity
      nccons                   = huge(1.0_r8), &
      nicons                   = huge(1.0_r8)

   integer :: ncnst

  !Define th (potential temperature) and (water vapor mixing 
  !ratio) qv at module level so "_old" values can easily be 
  !assigned at the beginning of each step from the value at the end 
  !of the step before.

  real(r8) :: th(pcols,pver)
  real(r8) :: qv(pcols,pver)

  character(len=8), parameter :: &      ! Constituent names
     cnst_names(8) = (/'CLDLIQ', 'CLDICE','NUMLIQ','NUMICE', &
                     'RAINQM', 'SNOWQM','NUMRAI','NUMSNO'/)
  contains
  !================================================================================================

  subroutine micro_p3_register

  logical :: prog_modal_aero
  logical :: save_subcol_microp ! If true, then need to store sub-columnized fields in pbuf

  if (masterproc) write(iulog,'(A20)') ' P3_REG Start'
!  if(l_summary_debug) write(6,*) 'micro_p3_register - 001 -'

   call phys_getopts(use_subcol_microp_out = use_subcol_microp, &
                    prog_modal_aero_out   = prog_modal_aero )
   ncnst = 0
    ! Register Microphysics Constituents 
    ! (i.e. members of state%q) and save indices.
    ! TODO make sure the cnst_names match what we think they are here.
    !================
   call cnst_add(cnst_names(1), mwdry, cpair, 0._r8, ixcldliq, &
         longname='Grid box averaged cloud liquid amount', &
         is_convtran1=.true.)
   ncnst = ncnst + 1
   call cnst_add(cnst_names(2), mwdry, cpair, 0._r8, ixcldice, &
         longname='Grid box averaged cloud ice amount', &
         is_convtran1=.true.)
   ncnst = ncnst + 1
   call cnst_add(cnst_names(3), mwh2o, cpair, 0._r8, ixnumliq, &
         longname='Grid box averaged cloud liquid number', &
         is_convtran1=.true.)
   ncnst = ncnst + 1
   call cnst_add(cnst_names(4), mwh2o, cpair, 0._r8, ixnumice, &
         longname='Grid box averaged cloud ice number', &
         is_convtran1=.true.)
   ncnst = ncnst + 1
   call cnst_add(cnst_names(5), mwh2o, cpair, 0._r8, ixrain, &
         longname='Grid box averaged rain amount', &
         is_convtran1=.true.)
   ncnst = ncnst + 1
   call cnst_add(cnst_names(6), mwh2o, cpair, 0._r8, ixcldrim, &
         longname='Grid box averaged riming amount', &
         is_convtran1=.true.)
   ncnst = ncnst + 1
   call cnst_add(cnst_names(7), mwh2o, cpair, 0._r8, ixnumrain, &
         longname='Grid box averaged rain number', &
         is_convtran1=.true.)
   ncnst = ncnst + 1
   call cnst_add(cnst_names(8), mwh2o, cpair, 0._r8, ixrimvol, &
         longname='Grid box averaged riming volume', &
         is_convtran1=.true.)
   ncnst = ncnst + 1
!+++ Aaron
!    call cnst_add(cnst_names(8), mwh2o, cpair, 0._r8, ixqirim, &
!         longname='Grid box averaged riming ???', &
!         is_convtran1=.true.)  ! TODO what is this?

    ! Add Variables to Pbuf
    !================
    !! module radiation_data & module cloud_rad_props
!    call pbuf_add_field('DEI',        'physpkg',dtype_r8,(/pcols,pver/), dei_idx)
!    call pbuf_add_field('MU',         'physpkg',dtype_r8,(/pcols,pver/), mu_idx)
!    call pbuf_add_field('LAMBDAC',    'physpkg',dtype_r8,(/pcols,pver/), lambdac_idx)


!  if(l_summary_debug) write(6,*) 'micro_p3_register - 002 -'

  !! module microp_aero
  call pbuf_add_field('CLDO','global', dtype_r8,(/pcols,pver,dyn_time_lvls/),cldo_idx)

  !! module wetdep 
  call pbuf_add_field('QME',  'physpkg',dtype_r8,(/pcols,pver/), qme_idx)
  call pbuf_add_field('PRAIN','physpkg',dtype_r8,(/pcols,pver/), prain_idx)
  call pbuf_add_field('NEVAPR','physpkg',dtype_r8,(/pcols,pver/), nevapr_idx)

  !! module aero_model
  if (prog_modal_aero) then
     call pbuf_add_field('RATE1_CW2PR_ST','physpkg',dtype_r8,(/pcols,pver/), rate1_cw2pr_st_idx)
  endif

  !! module clubb_intr
  call pbuf_add_field('PRER_EVAP',  'global', dtype_r8,(/pcols,pver/), prer_evap_idx)

  !! module radiation_data & module cloud_rad_props
  call pbuf_add_field('DEI',        'physpkg',dtype_r8,(/pcols,pver/), dei_idx)
  call pbuf_add_field('MU',         'physpkg',dtype_r8,(/pcols,pver/), mu_idx)
  call pbuf_add_field('LAMBDAC',    'physpkg',dtype_r8,(/pcols,pver/), lambdac_idx)

  !! module cospsimulator_intr
  call pbuf_add_field('REL',        'physpkg',dtype_r8,(/pcols,pver/), rel_idx)
  call pbuf_add_field('REI',        'physpkg',dtype_r8,(/pcols,pver/), rei_idx)
  call pbuf_add_field('LS_FLXPRC',  'physpkg',dtype_r8,(/pcols,pverp/), ls_flxprc_idx)
  call pbuf_add_field('LS_FLXSNW',  'physpkg',dtype_r8,(/pcols,pverp/), ls_flxsnw_idx)
  call pbuf_add_field('LS_REFFRAIN','physpkg',dtype_r8,(/pcols,pver/), ls_reffrain_idx)
  call pbuf_add_field('LS_REFFSNOW','physpkg',dtype_r8,(/pcols,pver/), ls_reffsnow_idx)
  call pbuf_add_field('CV_REFFLIQ', 'physpkg',dtype_r8,(/pcols,pver/), cv_reffliq_idx)
  call pbuf_add_field('CV_REFFICE', 'physpkg',dtype_r8,(/pcols,pver/), cv_reffice_idx)

  !! module macrop_driver
  call pbuf_add_field('CC_T',     'global',  dtype_r8,(/pcols,pver,dyn_time_lvls/), cc_t_idx)
  call pbuf_add_field('CC_qv',    'global',  dtype_r8,(/pcols,pver,dyn_time_lvls/), cc_qv_idx)
  call pbuf_add_field('CC_ql',    'global',  dtype_r8,(/pcols,pver,dyn_time_lvls/), cc_ql_idx)
  call pbuf_add_field('CC_qi',    'global',  dtype_r8,(/pcols,pver,dyn_time_lvls/), cc_qi_idx)
  call pbuf_add_field('CC_nl',    'global',  dtype_r8,(/pcols,pver,dyn_time_lvls/), cc_nl_idx)
  call pbuf_add_field('CC_ni',    'global',  dtype_r8,(/pcols,pver,dyn_time_lvls/), cc_ni_idx)
  call pbuf_add_field('CC_qlst',  'global',  dtype_r8,(/pcols,pver,dyn_time_lvls/), cc_qlst_idx)

  !! (internal) Stratiform only in cloud liquid/ice water path for radiation
  call pbuf_add_field('ICLWPST',    'physpkg',dtype_r8,(/pcols,pver/),iclwpst_idx)
  call pbuf_add_field('ICIWPST',    'physpkg',dtype_r8,(/pcols,pver/),iciwpst_idx)


!  if(l_summary_debug) write(6,*) 'micro_p3_register - 003 -'

  ! Register subcolumn pbuf fields
  if (use_subcol_microp) then

    if (masterproc) write(iulog,'(A20)') '  P3_REG subcol'
    call pbuf_register_subcol('CLDO',        'micro_p3_register', cldo_idx)
    call pbuf_register_subcol('QME',         'micro_mg_cam_register', qme_idx)
    call pbuf_register_subcol('PRAIN',       'micro_mg_cam_register', prain_idx)
    call pbuf_register_subcol('NEVAPR',      'micro_mg_cam_register', nevapr_idx)
    call pbuf_register_subcol('PRER_EVAP',   'micro_p3_register', prer_evap_idx)
    call pbuf_register_subcol('DEI',         'micro_p3_register', dei_idx)
    call pbuf_register_subcol('MU',          'micro_p3_register', mu_idx)
    call pbuf_register_subcol('LAMBDAC',     'micro_p3_register', lambdac_idx)
    call pbuf_register_subcol('REL',         'micro_p3_register', rel_idx)
    call pbuf_register_subcol('REI',         'micro_p3_register', rei_idx)
    call pbuf_register_subcol('LS_FLXPRC',   'micro_p3_register', ls_flxprc_idx)
    call pbuf_register_subcol('LS_FLXSNW',   'micro_p3_register', ls_flxsnw_idx)
    call pbuf_register_subcol('LS_REFFRAIN', 'micro_p3_register', ls_reffrain_idx)
    call pbuf_register_subcol('LS_REFFSNOW', 'micro_p3_register', ls_reffsnow_idx)
    call pbuf_register_subcol('CV_REFFLIQ',  'micro_p3_register', cv_reffliq_idx)
    call pbuf_register_subcol('CV_REFFICE',  'micro_p3_register', cv_reffice_idx)
    call pbuf_register_subcol('CC_T',        'micro_p3_register', cc_t_idx)
    call pbuf_register_subcol('CC_qv',       'micro_p3_register', cc_qv_idx)
    call pbuf_register_subcol('CC_ql',       'micro_p3_register', cc_ql_idx)
    call pbuf_register_subcol('CC_qi',       'micro_p3_register', cc_qi_idx)
    call pbuf_register_subcol('CC_nl',       'micro_p3_register', cc_nl_idx)
    call pbuf_register_subcol('CC_ni',       'micro_p3_register', cc_ni_idx)
    call pbuf_register_subcol('CC_qlst',     'micro_p3_register', cc_qlst_idx)
    call pbuf_register_subcol('ICIWPST',     'micro_p3_register', iciwpst_idx)
    call pbuf_register_subcol('ICLWPST',     'micro_p3_register', iclwpst_idx)

    if (prog_modal_aero) then
      call pbuf_register_subcol('RATE1_CW2PR_ST', 'micro_p3_register', rate1_cw2pr_st_idx)
    end if

  end if

  !! (internal) Precipitation efficiency fields across timesteps.
  call pbuf_add_field('ACPRECL',    'global',dtype_r8,(/pcols/), acpr_idx)   ! accumulated precip
  call pbuf_add_field('ACGCME',     'global',dtype_r8,(/pcols/), acgcme_idx) ! accumulated condensation
  call pbuf_add_field('ACNUM',      'global',dtype_i4,(/pcols/), acnum_idx)  ! counter for accumulated # timesteps

  !! module clubb_intr
  call pbuf_add_field('RELVAR',     'global',dtype_r8,(/pcols,pver/), relvar_idx)
  call pbuf_add_field('ACCRE_ENHAN','global',dtype_r8,(/pcols,pver/), accre_enhan_idx)

  ! Diagnostic fields needed for subcol_SILHS, need to be grid-only
  if (subcol_get_scheme() == 'SILHS') then
     call pbuf_add_field('QRAIN',   'global',dtype_r8,(/pcols,pver/), qrain_idx)
     call pbuf_add_field('NRAIN',   'global',dtype_r8,(/pcols,pver/), nrain_idx)
  end if

!  if(l_summary_debug) write(6,*) 'micro_p3_register - 004 -'
!--- Aaron
  if (masterproc) write(iulog,'(A20)') '    P3_REG Finished'
  end subroutine micro_p3_register

  !================================================================================================
  function micro_p3_implements_cnst(name)

    ! Return true if specified constituent is implemented by the
    ! microphysics package

    character(len=*), intent(in) :: name        ! constituent name
    logical :: micro_p3_implements_cnst    ! return value

    micro_p3_implements_cnst = any(name == cnst_names)

  end function micro_p3_implements_cnst


  !================================================================================================

  subroutine micro_p3_init_cnst(name, q, gcid)

    ! Initialize the microphysics constituents, if they are
    ! not read from the initial file.

    character(len=*), intent(in) :: name     ! constituent name
    real(r8), intent(out) :: q(:,:)   ! mass mixing ratio (gcol, plev)
    integer,  intent(in)  :: gcid(:)  ! global column id

    if (micro_p3_implements_cnst(name)) q = 0.0_r8

  end subroutine micro_p3_init_cnst

  !================================================================================================

  subroutine micro_p3_init()
    use micro_p3,       only: p3_init
    use cam_history,    only: addfld, add_default, horiz_only

    character(128) :: p3_lookup_dir, errstring
    integer        :: m, mm

    ! PULL DIRECTORY OF LOOKUP TABLE FROM ATM NAMELIST
    ! =============
    call phys_getopts(p3_lookup_dir_out = p3_lookup_dir)

    ! CALL P3 INIT:
    !==============
    !might want to add all E3SM parameter vals to p3_init call...

    !TODO: add errstring and constants to init function  !DONE
    call p3_init(p3_lookup_dir,errstring)
    call handle_errmsg(errstring, subname="micro_p3_init")

    ! INITIALIZE OUTPUT
    !==============
    !TODO: put addfld and add_default calls here.  !Done and testing
    do m = 1, ncnst
       call cnst_get_ind(cnst_names(m), mm)
       if ( any(mm == (/ ixcldliq, ixcldice, ixrain, ixcldrim /)) ) then
          ! mass mixing ratios
          call addfld(cnst_name(mm), (/ 'lev' /), 'A', 'kg/kg', &
            cnst_longname(mm) )
          call addfld(sflxnam(mm), horiz_only, 'A', 'kg/m2/s', &
            trim(cnst_name(mm))//' surface flux')
       else if ( any(mm == (/ ixnumliq, ixnumice, ixnumrain /)) ) then
          ! number concentrations
          call addfld(cnst_name(mm), (/ 'lev' /), 'A', '1/kg', &
            cnst_longname(mm) )
          call addfld(sflxnam(mm), horiz_only, 'A', '1/m2/s', &
            trim(cnst_name(mm))//' surface flux')
       else if ( mm == ixrimvol ) then
          ! number concentrations
          call addfld(cnst_name(mm), (/ 'lev' /), 'A', 'm3/kg', &
            cnst_longname(mm) )
          call addfld(sflxnam(mm), horiz_only, 'A', 'm3/m2/s', &
            trim(cnst_name(mm))//' surface flux')
       else
          call endrun( "micro_p3_acme_init: &
               &Could not call addfld for constituent with unknown units.")
       endif
    end do
    call addfld(apcnst(ixcldliq), (/ 'lev' /), 'A', 'kg/kg', trim(cnst_name(ixcldliq))//' after physics'  )
    call addfld(apcnst(ixcldice), (/ 'lev' /), 'A', 'kg/kg', trim(cnst_name(ixcldice))//' after physics'  )
    call addfld(bpcnst(ixcldliq), (/ 'lev' /), 'A', 'kg/kg', trim(cnst_name(ixcldliq))//' before physics' )
    call addfld(bpcnst(ixcldice), (/ 'lev' /), 'A', 'kg/kg', trim(cnst_name(ixcldice))//' before physics' )
    call addfld(apcnst(ixrain),   (/ 'lev' /), 'A', 'kg/kg', trim(cnst_name(ixrain))//' after physics'  )
    call addfld(bpcnst(ixrain),   (/ 'lev' /), 'A', 'kg/kg', trim(cnst_name(ixrain))//' before physics' )
    call addfld(apcnst(ixcldrim), (/ 'lev' /), 'A', 'kg/kg', trim(cnst_name(ixcldrim))//' after physics'  )
    call addfld(bpcnst(ixcldrim), (/ 'lev' /), 'A', 'kg/kg', trim(cnst_name(ixcldrim))//' before physics' )


    call addfld ('CME', (/ 'lev' /), 'A', 'kg/kg/s', 'Rate of cond-evap within the cloud'                      )
!    call addfld ('PRODPREC', (/ 'lev' /), 'A', 'kg/kg/s', 'Rate of conversion of condensate to precip'              )
!    call addfld ('EVAPPREC', (/ 'lev' /), 'A', 'kg/kg/s', 'Rate of evaporation of falling precip'                   )
!    call addfld ('EVAPSNOW', (/ 'lev' /), 'A', 'kg/kg/s', 'Rate of evaporation of falling snow'                     )
!    call addfld ('HPROGCLD', (/ 'lev' /), 'A', 'W/kg'    , 'Heating from prognostic clouds'                          )
!    call addfld ('FICE', (/ 'lev' /), 'A', 'fraction', 'Fractional ice content within cloud'                     )
    call addfld ('ICWMRST', (/ 'lev' /), 'A', 'kg/kg', 'Prognostic in-stratus water mixing ratio'                )
    call addfld ('ICIMRST', (/ 'lev' /), 'A', 'kg/kg', 'Prognostic in-stratus ice mixing ratio'                  )

   ! MG microphysics diagnostics
!    call addfld ('QCSEVAP', (/ 'lev' /), 'A', 'kg/kg/s', 'Rate of evaporation of falling cloud water'              )
!    call addfld ('QISEVAP', (/ 'lev' /), 'A', 'kg/kg/s', 'Rate of sublimation of falling cloud ice'                )
!    call addfld ('QVRES', (/ 'lev' /), 'A', 'kg/kg/s', 'Rate of residual condensation term'                      )
    call addfld ('CMEIOUT', (/ 'lev' /), 'A', 'kg/kg/s', 'Rate of deposition/sublimation of cloud ice'             )
!    call addfld ('VTRMC', (/ 'lev' /), 'A', 'm/s', 'Mass-weighted cloud water fallspeed'                     )
!    call addfld ('VTRMI', (/ 'lev' /), 'A', 'm/s', 'Mass-weighted cloud ice fallspeed'                       )
    call addfld ('QCSEDTEN', (/ 'lev' /), 'A', 'kg/kg/s', 'Cloud water mixing ratio tendency from sedimentation'    )
    call addfld ('QISEDTEN', (/ 'lev' /), 'A', 'kg/kg/s', 'Cloud ice mixing ratio tendency from sedimentation'      )
!    call addfld ('PRAO', (/ 'lev' /), 'A', 'kg/kg/s', 'Accretion of cloud water by rain'                        )
!    call addfld ('PRCO', (/ 'lev' /), 'A', 'kg/kg/s', 'Autoconversion of cloud water'                           )
!    call addfld ('MNUCCCO', (/ 'lev' /), 'A', 'kg/kg/s', 'Immersion freezing of cloud water'                       )
!    call addfld ('MNUCCTO', (/ 'lev' /), 'A', 'kg/kg/s', 'Contact freezing of cloud water'                         )
!    call addfld ('MNUCCDO', (/ 'lev' /), 'A', 'kg/kg/s', 'Homogeneous and heterogeneous nucleation from vapor'     )
!    call addfld ('MNUCCDOhet', (/ 'lev' /), 'A','kg/kg/s', 'Heterogeneous nucleation from vapor'                     )
!    call addfld ('MSACWIO', (/ 'lev' /), 'A', 'kg/kg/s', 'Conversion of cloud water from rime-splintering'         )
!    call addfld ('PSACWSO', (/ 'lev' /), 'A', 'kg/kg/s', 'Accretion of cloud water by snow'                        )
!    call addfld ('BERGSO', (/ 'lev' /), 'A', 'kg/kg/s', 'Conversion of cloud water to snow from bergeron'         )
!    call addfld ('BERGO', (/ 'lev' /), 'A', 'kg/kg/s', 'Conversion of cloud water to cloud ice from bergeron'    )
!    call addfld ('MELTO', (/ 'lev' /), 'A', 'kg/kg/s', 'Melting of cloud ice'                                    )
!    call addfld ('HOMOO', (/ 'lev' /), 'A', 'kg/kg/s', 'Homogeneous freezing of cloud water'                     )
!    call addfld ('QCRESO', (/ 'lev' /), 'A', 'kg/kg/s', 'Residual condensation term for cloud water'              )
!    call addfld ('PRCIO', (/ 'lev' /), 'A', 'kg/kg/s', 'Autoconversion of cloud ice'                             )
!    call addfld ('PRAIO', (/ 'lev' /), 'A', 'kg/kg/s', 'Accretion of cloud ice by rain'                          )
!    call addfld ('QIRESO', (/ 'lev' /), 'A', 'kg/kg/s', 'Residual deposition term for cloud ice'                  )
!    call addfld ('MNUCCRO', (/ 'lev' /), 'A', 'kg/kg/s', 'Heterogeneous freezing of rain to snow'                  )
!    call addfld ('PRACSO', (/ 'lev' /), 'A', 'kg/kg/s', 'Accretion of rain by snow'                               )
!    call addfld ('MELTSDT', (/ 'lev' /), 'A', 'W/kg', 'Latent heating rate due to melting of snow'              )
!    call addfld ('FRZRDT', (/ 'lev' /), 'A', 'W/kg', 'Latent heating rate due to homogeneous freezing of rain' )
   call addfld ('QRSEDTEN', (/ 'lev' /), 'A', 'kg/kg/s', 'Rain mixing ratio tendency from sedimentation'           )
!    call addfld ('QSSEDTEN', (/ 'lev' /), 'A', 'kg/kg/s', 'Snow mixing ratio tendency from sedimentation'           )

   ! History variables for CAM5 microphysics
!    call addfld ('MPDT', (/ 'lev' /), 'A', 'W/kg', 'Heating tendency - Morrison microphysics'                )
!    call addfld ('MPDQ', (/ 'lev' /), 'A', 'kg/kg/s', 'Q tendency - Morrison microphysics'                      )
!    call addfld ('MPDLIQ', (/ 'lev' /), 'A', 'kg/kg/s', 'CLDLIQ tendency - Morrison microphysics'                 )
!    call addfld ('MPDICE', (/ 'lev' /), 'A', 'kg/kg/s', 'CLDICE tendency - Morrison microphysics'                 )
!    call addfld ('MPDW2V', (/ 'lev' /), 'A', 'kg/kg/s', 'Water <--> Vapor tendency - Morrison microphysics'       )
!    call addfld ('MPDW2I', (/ 'lev' /), 'A', 'kg/kg/s', 'Water <--> Ice tendency - Morrison microphysics'         )
!    call addfld ('MPDW2P', (/ 'lev' /), 'A', 'kg/kg/s', 'Water <--> Precip tendency - Morrison microphysics'      )
!    call addfld ('MPDI2V', (/ 'lev' /), 'A', 'kg/kg/s', 'Ice <--> Vapor tendency - Morrison microphysics'         )
!    call addfld ('MPDI2W', (/ 'lev' /), 'A', 'kg/kg/s', 'Ice <--> Water tendency - Morrison microphysics'         )
!    call addfld ('MPDI2P', (/ 'lev' /), 'A', 'kg/kg/s', 'Ice <--> Precip tendency - Morrison microphysics'        )
    call addfld ('ICWNC', (/ 'lev' /), 'A', 'm-3', 'Prognostic in-cloud water number conc'                   )
    call addfld ('ICINC', (/ 'lev' /), 'A', 'm-3', 'Prognostic in-cloud ice number conc'                     )
!    call addfld ('EFFLIQ_IND', (/ 'lev' /), 'A','Micron', 'Prognostic droplet effective radius (indirect effect)'   )
    call addfld ('CDNUMC', horiz_only,    'A', '1/m2', 'Vertically-integrated droplet concentration'             )
    call addfld ('MPICLWPI', horiz_only,    'A', 'kg/m2', 'Vertically-integrated &
         &in-cloud Initial Liquid WP (Before Micro)' )
    call addfld ('MPICIWPI', horiz_only,    'A', 'kg/m2', 'Vertically-integrated &
         &in-cloud Initial Ice WP (Before Micro)'    )

   ! This is provided as an example on how to write out subcolumn output
   ! NOTE -- only 'I' should be used for sub-column fields as subc-columns could shift from time-step to time-step
!   if (use_subcol_microp) then
!      call addfld('FICE_SCOL', (/'psubcols','lev     '/), 'I', 'fraction', &
!           'Sub-column fractional ice content within cloud', flag_xyfill=.true., fill_value=1.e30_r8)
!   end if

   ! Averaging for cloud particle number and size
   call addfld ('AWNC', (/ 'lev' /), 'A', 'm-3', 'Average cloud water number conc'                         )
   call addfld ('AWNI', (/ 'lev' /), 'A', 'm-3', 'Average cloud ice number conc'                           )
   call addfld ('AREL', (/ 'lev' /), 'A', 'Micron', 'Average droplet effective radius'                        )
   call addfld ('AREI', (/ 'lev' /), 'A', 'Micron', 'Average ice effective radius'                            )
   ! Frequency arrays for above
   call addfld ('FREQL', (/ 'lev' /), 'A', 'fraction', 'Fractional occurrence of liquid'                          )
   call addfld ('FREQI', (/ 'lev' /), 'A', 'fraction', 'Fractional occurrence of ice'                             )

   ! Average cloud top particle size and number (liq, ice) and frequency
!   call addfld ('ACTREL', horiz_only,    'A', 'Micron', 'Average Cloud Top droplet effective radius'              )
!   call addfld ('ACTREI', horiz_only,    'A', 'Micron', 'Average Cloud Top ice effective radius'                  )
!   call addfld ('ACTNL', horiz_only,    'A', 'Micron', 'Average Cloud Top droplet number'                        )
!   call addfld ('ACTNI', horiz_only,    'A', 'Micron', 'Average Cloud Top ice number'                            )

!   call addfld ('FCTL', horiz_only,    'A', 'fraction', 'Fractional occurrence of cloud top liquid'                )
!   call addfld ('FCTI', horiz_only,    'A', 'fraction', 'Fractional occurrence of cloud top ice'                   )

!   call addfld ('LS_FLXPRC', (/ 'ilev' /), 'A', 'kg/m2/s', 'ls stratiform gbm interface rain+snow flux')
!   call addfld ('LS_FLXSNW', (/ 'ilev' /), 'A', 'kg/m2/s', 'ls stratiform gbm interface snow flux')

   call addfld ('REL', (/ 'lev' /), 'A', 'micron', 'MG REL stratiform cloud effective radius liquid')
   call addfld ('REI', (/ 'lev' /), 'A', 'micron', 'MG REI stratiform cloud effective radius ice')
!   call addfld ('LS_REFFRAIN', (/ 'lev' /), 'A', 'micron', 'ls stratiform rain effective radius')
!   call addfld ('LS_REFFSNOW', (/ 'lev' /), 'A', 'micron', 'ls stratiform snow effective radius')
!   call addfld ('CV_REFFLIQ', (/ 'lev' /), 'A', 'micron', 'convective cloud liq effective radius')
!   call addfld ('CV_REFFICE', (/ 'lev' /), 'A', 'micron', 'convective cloud ice effective radius')

!!== KZ_DCS
   call addfld ('DCST',(/ 'lev' /), 'A','m','dcs')
!!== KZ_DCS
   ! diagnostic precip
!   call addfld ('QRAIN',(/ 'lev' /), 'A','kg/kg','Diagnostic grid-mean rain mixing ratio'         )
!   call addfld ('QSNOW',(/ 'lev' /), 'A','kg/kg','Diagnostic grid-mean snow mixing ratio'         )
!   call addfld ('NRAIN',(/ 'lev' /), 'A','m-3','Diagnostic grid-mean rain number conc'         )
!   call addfld ('NSNOW',(/ 'lev' /), 'A','m-3','Diagnostic grid-mean snow number conc'         )

   ! size of precip
!   call addfld ('RERCLD',(/ 'lev' /), 'A','m','Diagnostic effective radius of Liquid Cloud and Rain' )
!   call addfld ('DSNOW',(/ 'lev' /), 'A','m','Diagnostic grid-mean snow diameter'         )

   ! diagnostic radar reflectivity, cloud-averaged
!   call addfld ('REFL',(/ 'lev' /), 'A','DBz','94 GHz radar reflectivity'       )
!   call addfld ('AREFL',(/ 'lev' /), 'A','DBz','Average 94 GHz radar reflectivity'       )
!   call addfld ('FREFL',(/ 'lev' /), 'A','fraction','Fractional occurrence of radar reflectivity'       )

!   call addfld ('CSRFL',(/ 'lev' /), 'A','DBz','94 GHz radar reflectivity (CloudSat thresholds)'       )
!   call addfld ('ACSRFL',(/ 'lev' /), 'A','DBz','Average 94 GHz radar reflectivity (CloudSat thresholds)'       )
!   call addfld ('FCSRFL',(/ 'lev' /), 'A','fraction','Fractional occurrence of radar reflectivity (CloudSat thresholds)' &
!        )

!   call addfld ('AREFLZ',(/ 'lev' /), 'A','mm^6/m^3','Average 94 GHz radar reflectivity'       )

   ! Aerosol information
   call addfld ('NCAL',(/ 'lev' /), 'A','1/m3','Number Concentation Activated for Liquid')
   call addfld ('NCAI',(/ 'lev' /), 'A','1/m3','Number Concentation Activated for Ice')

   ! Average rain and snow mixing ratio (Q), number (N) and diameter (D), with frequency
!   call addfld ('AQRAIN',(/ 'lev' /), 'A','kg/kg','Average rain mixing ratio'         )
!   call addfld ('AQSNOW',(/ 'lev' /), 'A','kg/kg','Average snow mixing ratio'         )
!   call addfld ('ANRAIN',(/ 'lev' /), 'A','m-3','Average rain number conc'         )
!   call addfld ('ANSNOW',(/ 'lev' /), 'A','m-3','Average snow number conc'         )
   call addfld ('ADRAIN',(/ 'lev' /), 'A','Micron','Average rain effective Diameter'         )
!   call addfld ('ADSNOW',(/ 'lev' /), 'A','Micron','Average snow effective Diameter'         )
   call addfld ('FREQR',(/ 'lev' /), 'A','fraction','Fractional occurrence of rain'       )
!   call addfld ('FREQS',(/ 'lev' /), 'A','fraction','Fractional occurrence of snow'       )

   ! precipitation efficiency & other diagnostic fields
!   call addfld('PE'    ,       horiz_only, 'A', '1', 'Stratiform Precipitation Efficiency  (precip/cmeliq)' )
!   call addfld('APRL'  ,     horiz_only, 'A', 'm/s', 'Average Stratiform Precip Rate over efficiency calculation' )
!   call addfld('PEFRAC',       horiz_only, 'A', '1', 'Fraction of timesteps precip efficiency reported' )
!   call addfld('VPRCO' , horiz_only, 'A', 'kg/kg/s', 'Vertical average of autoconversion rate' )
!   call addfld('VPRAO' , horiz_only, 'A', 'kg/kg/s', 'Vertical average of accretion rate' )
!   call addfld('RACAU' , horiz_only, 'A', 'kg/kg/s', 'Accretion/autoconversion ratio from vertical average' )

   call addfld('UMR', (/ 'lev' /), 'A',   'm/s', 'Mass-weighted rain  fallspeed'              )
!      call addfld('UMS', (/ 'lev' /), 'A',   'm/s', 'Mass-weighted snow fallspeed'               )

   ! qc limiter (only output in versions 1.5 and later)
!   if (.not. (micro_mg_version == 1 .and. micro_mg_sub_version == 0)) then
!      call addfld('QCRAT', (/ 'lev' /), 'A', 'fraction', 'Qc Limiter: Fraction of qc tendency applied')
!   end if


  end subroutine micro_p3_init

  !================================================================================================

  subroutine micro_p3_tend(state, ptend, dtime, pbuf)

    use cam_history,    only: outfld
    use time_manager,   only: is_first_step
    use physics_buffer, only: pbuf_col_type_index

    !INPUT/OUTPUT VARIABLES
    type(physics_state),         intent(in)    :: state
    type(physics_ptend),         intent(out)   :: ptend
    real(r8),                    intent(in)    :: dtime
    type(physics_buffer_desc),   pointer       :: pbuf(:)
    logical :: lq(pcnst)   !list of what constituents to update

    !INTERNAL VARIABLES
    real(r8) :: th_old(pcols,pver)     !potential temperature from last step   K
    real(r8) :: qv_old(pcols,pver)     !water vapor from last step             kg/kg
    real(r8) :: ssat(pcols,pver)       !supersaturated mixing ratio            kg/kg
    real(r8) :: dzq(pcols,pver)        !geometric layer thickness              m
    real(r8) :: cldliq(pcols,pver)     !cloud liquid water mixing ratio        kg/kg
    real(r8) :: numliq(pcols,pver)     !cloud liquid water drop concentraiton  #/kg
    real(r8) :: rain(pcols,pver)       !rain water mixing ratio                kg/kg
    real(r8) :: numrain(pcols,pver)    !rain water number concentration        #/kg
    real(r8) :: qv(pcols,pver)         !water vapor mixing ratio               kg/kg
    real(r8) :: ice(pcols,pver)        !total ice water mixing ratio           kg/kg
    real(r8) :: qirim(pcols,pver)      !rime ice mixing ratio                  kg/kg
    real(r8) :: numice(pcols,pver)     !total ice crystal number concentration #/kg
    real(r8) :: rimvol(pcols,pver)     !rime volume mixing ratio               m3/kg
    real(r8) :: temp(pcols,pver)       !potential temperature                  K
    real(r8) :: rim(pcols,pver)        !rime mixing ratio                      kg/kg
    integer :: it                      !timestep counter                       -
    integer :: kts                     !closest level to TOM                   -
    integer :: kte                     !near surface level                     -
    real(r8), pointer :: rel(:,:)      ! Liquid effective drop radius (microns)
    real(r8), pointer :: rei(:,:)      ! Ice effective drop size (microns)

    logical :: log_predictNc           !prognostic droplet concentration or not?
    integer :: col_type ! Flag to store whether accessing grid or sub-columns in pbuf_get_field
    integer :: icol, ncol, k

    ! Associate Pbuf Variables
    !==============

    call pbuf_col_type_index(use_subcol_microp, col_type=col_type)
    !liq effective radius (m)
    call pbuf_get_field(pbuf, rel_idx,         rel, col_type=col_type)
    !ice effective radius (m)
    call pbuf_get_field(pbuf, rei_idx,         rei, col_type=col_type)

    ! INITIALIZE PTEND
    !==============
    !ptend is an output variable. Since not substepping in micro, don't need 
    !a local copy.

    lq           = .false. !initialize all constituents to false.
    lq(1)        = .true.
    lq(ixcldliq) = .true.
    lq(ixcldice) = .true.
    lq(ixnumliq) = .true.
    lq(ixnumice) = .true.
    lq(ixrain)   = .true.
    lq(ixcldrim)    = .true.
    lq(ixnumrain)= .true.
    lq(ixrimvol)  = .true.

    call physics_ptend_init(ptend, state%psetcols, "micro_p3", ls=.true., lq=lq)

    ! HANDLE AEROSOL ACTIVATION
    !==============
    !saving this for later... use prescribed Nd for now.
    log_predictNc = .false.

    ! GET "OLD" VALUES
    !==============
    ! TODO: answer question: is_first_step =1 for each submission, or just for type=initial first step?
    !          does p3 want _old values from last p3 call, or from some other point?
    ! Aaron- is_first_step=1 only for first step of initial run, per:
    !      logical function is_first_step()
    !              Return true on first step of initial run only.
    if ( is_first_step() ) then
       th_old=state%t*state%exner !this is wrong, just using to get started.
       qv_old=state%q(:,:,1)
    else
       th_old = th !use th from end of last p3 step
       qv_old = qv !use qv from end of last p3 step
    end if

    ! CONVERT T TO POTENTIAL TEMPERATURE
    !==============
    ! Someday we may want to make P3 take in T rather than theta... it uses both.


    ! COMPUTE GEOMETRIC THICKNESS OF GRID
    !==============
    ncol = state%ncol
    do icol = 1,ncol
       do k = 1,pver
          dzq(icol,k) = state%zi(icol,k) - state%zi(icol,k+1)
          th(icol,k)  = state%t(icol,k)*state%exner(icol,k)
       end do
    end do

    ! ASSIGN TOP AND BOTTOM INDICES FOR GRID
    !==============
    !kts is closest level to top of model. Instead of 1 (top-of-model), 
    !we use the previously-defined trop_cloud_top_lev to reduce the number of 
    !levels we need to calculate and to avoid upper-atmos regions where this
    !micro-physics is inappropriate. kte is the near-surface level = pver.

    kts=top_lev
    kte=pver 

    ! DEAL WITH SSAT
    !==============
    !ssat (supersaturated mixing ratio measured in kg/kg) can be prognosed
    !or diagnosed in p3 depending on p3's hardcoded log_predictSsat parameter.
    !ssat is an intent inout variable, but when log_predictSsat is false 
    !(as we will set it initially), ssat is overwritten rather than used by p3.
    !Thus it shouldn't matter what ssat is, so we give it -999.

    ssat(:,:) = -999._r8

    ! HANDLE TIMESTEP COUNTER
    !==============
    !p3 wants to know the timestep number (which it calls "it") because it
    !handles things differently on the first step, where it doesn't have values
    !yet. E3SM has a handy function for deciding if this is the first step, so 
    !we hack "it" with "is_first_step()" for now. Eventually, we should replace
    !"it" with a logical.

    if (is_first_step()) then
       it=1
    else
       it=999 !integer
    end if

    ! MAKE LOCAL COPIES OF VARS MODIFIED BY P3
    !==============
    !local copies are needed because state is passed into this routine as intent=in
    !while P3 seeks to modify state variables in-place. Also, we need a copy of 
    !old values in order to back out ptend values later. Traditionally, a local copy 
    !is created by copying the whole state. It is much cheaper to just copy the 
    !variables we need. 
    
    cldliq  = state%q(:,:,ixcldliq)
    numliq  = state%q(:,:,ixnumliq)
    rain    = state%q(:,:,ixrain)
    numrain = state%q(:,:,ixnumrain)
    qv      = state%q(:,:,1)
    ice     = state%q(:,:,ixcldice)
    qirim   = state%q(:,:,ixcldrim) !Aaron, changed ixqirim to ixcldrim to match Kai's code
    numice  = state%q(:,:,ixnumice)
    rimvol  = state%q(:,:,ixrimvol)

    ! CALL P3
    !==============
#if 0
    call p3_main( &
         cldliq, & ! INOUT  cloud, mass mixing ratio         kg kg-1
         numliq, & ! INOUT  cloud, number mixing ratio       #  kg-1
         rain,   & ! INOUT  rain, mass mixing ratio          kg kg-1
         numrain,& ! INOUT  rain, number mixing ratio        #  kg-1
         th_old,                & ! INOUT  beginning of time step theta     K
         th,                    & ! INOUT  potential temperature            K
         qv_old,                & ! INOUT  beginning of time step qv        kg kg-1
         qv,        & ! INOUT  water vapor mixing ratio         kg kg-1
         dtime,                 & ! IN     model time step                  s
         ice,    & ! INOUT  ice, total mass mixing ratio     kg kg-1
         qirim,  & ! INOUT  ice, rime mass mixing ratio      kg kg-1
         numice, & ! INOUT  ice, total number mixing ratio   #  kg-1
         rimvol, & ! INOUT  ice, rime volume mixing ratio    m3 kg-1
         ssat,                  & ! INOUT  supersaturation (i.e., qv-qvs)   kg kg-1
         state%pmid(:,:),       & ! IN     pressure at cell midpoints       Pa
         dzq,                   & ! IN     vertical grid spacing            m
         it,                    & ! IN     time step counter NOTE: starts at 1 for first time step
         prt_liq,               & ! OUT    surface liquid precip rate       m s-1
         prt_sol,               & ! OUT    surface frozen precip rate       m s-1
         its,                   & ! IN     horizontal index lower bound     -
         ite,                   & ! IN     horizontal index upper bound     -
         kts,                   & ! IN     vertical index lower bound       -
         kte,                   & ! IN     vertical index upper bound       -
         diag_ze,               & ! OUT    equivalent reflectivity          dBZ  UNUSED?
         rel,                   & ! OUT    effective radius, cloud          m
         rei,                   & ! OUT    effective radius, ice            m
         diag_vmi,              & ! OUT    mass-weighted fall speed of ice  m s-1
         diag_di,               & ! OUT    mean diameter of ice             m
         diag_rhoi,             & ! OUT    bulk density of ice              kg m-1
         log_predictNc,         & ! IN     .true.=prognostic Nc, .false.=specified Nc
         )
    !MASSAGE OUTPUT TO FIT E3SM EXPECTATIONS
    !============= 
    precl=prt_liq+prt_sol
#endif

    !TODO: figure out what else other E3SM parameterizations need from micro and make sure 
    !they are assigned here. The comments below are a step in that direction.

    !cloud_rad_props needs ice effective diameter, which Kai calculates as below:
    !   dei = rei*diag_rhopo(i,k,iice)/rhows*2._r8
    !where rhopo is bulk ice density from table lookup (taken from f1pr16, but not 
    !done in my ver yet) and rhows=917.0 is a constant parameter.

    !cloud_rad_props also uses snow radiative properties which aren't available from 
    !P3 (perhaps because ice phase in p3 includes *all* ice already?).

    !BACK OUT TENDENCIES FROM STATE CHANGES
    !=============
    !Aaron, imbed these calls inside do loops to avoid a series of "floating
    ! invalid" error messages at runtime
    do icol = 1,ncol
       do k = 1,pver
          temp(icol,k) = th(icol,k)/state%exner(icol,k) !convert theta to 
          ptend%s                   = cpair*(temp(icol,k) - state%t(icol,k))/dtime 
          ptend%q(icol,k,1)         = (qv(icol,k)      - state%q(icol,k,1) )/dtime
          ptend%q(icol,k,ixcldliq)  = (cldliq(icol,k)  - state%q(icol,k,ixcldliq) )/dtime
          ptend%q(icol,k,ixnumliq)  = (numliq(icol,k)  - state%q(icol,k,ixnumliq) )/dtime
          ptend%q(icol,k,ixrain)    = (rain(icol,k)    - state%q(icol,k,ixrain) )/dtime
          ptend%q(icol,k,ixnumrain) = (numrain(icol,k) - state%q(icol,k,ixnumrain) )/dtime
          ptend%q(icol,k,ixcldice)  = (ice(icol,k)     - state%q(icol,k,ixcldice) )/dtime
          ptend%q(icol,k,ixnumice)  = (numice(icol,k)  - state%q(icol,k,ixnumice) )/dtime
          ptend%q(icol,k,ixcldrim)  = (rim(icol,k)     - state%q(icol,k,ixcldrim) )/dtime
          ptend%q(icol,k,ixrimvol)  = (rimvol(icol,k)  - state%q(icol,k,ixrimvol) )/dtime
       end do
    end do

    !note s=cp*T has units J/kg


    !WRITE OUTPUT
    !=============
    !call outfld('P3_QCAUT',   qcaut_grid,  pcols, lchnk)

    !TODO: add other outfld calls. Probably not worth doing until we get the code compiling...

  end subroutine micro_p3_tend

  !================================================================================================

end module micro_p3_interface
