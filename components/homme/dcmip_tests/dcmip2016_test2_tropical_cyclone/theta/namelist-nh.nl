!
! theta: namelist for dcmip2016 test2: tropical cyclone
!_______________________________________________________________________
&ctl_nl
  nthreads          = 1
  partmethod        = 4                         ! mesh parition method: 4 = space filling curve
  topology          = "cube"                    ! mesh type: cubed sphere
  test_case         = "dcmip2016_test2"         ! test identifier
  ne                = 30                         ! number of elements per cube face
  qsize             = 4                         ! num tracer fields
  ndays             = 10
  statefreq         = 10                        ! number of steps between screen dumps
  restartfreq       = -1                        ! don't write restart files if < 0
  runtype           = 0                         ! 0 => new run
  tstep             = 0.5                       ! largest timestep in seconds
  integration       = 'explicit'                ! explicit time integration
  tstep_type        = 5                         ! 1 => default method
  rsplit            = 1
  qsplit            = 1
  nu                = 1e13
  nu_s              = 1e13
  nu_p              = 1e13
  hypervis_order    = 2                         ! 2 = hyperviscosity
  hypervis_subcycle = 1                         ! 1 = no hyperviz subcycling
  theta_hydrostatic_mode = .false.
/
&vert_nl
  vform         = "ccm"
  vfile_mid     = "../vcoord/camm-30.ascii"
  vfile_int     = "../vcoord/cami-30.ascii"
/
&analysis_nl
  output_dir        = "./movies/"               ! destination dir for netcdf file
  output_timeunits  = 2,                   ! 0=timesteps, 1=days, 2=hours, 3=seconds
  output_frequency  = 1 !4
  output_varnames1  ='T','p','ps','pnh','geo','u','v','w','Th','Q','Q2','Q3','Q4'   ! variables to write to file
  interp_type       = 0                         ! 0=native grid, 1=bilinear
  output_type       ='netcdf'                   ! netcdf or pnetcdf
  num_io_procs      = 16
  interp_nlon       = 360
  interp_nlat       = 181
/
&prof_inparm
  profile_outpe_num   = 100
  profile_single_file	= .true.
/