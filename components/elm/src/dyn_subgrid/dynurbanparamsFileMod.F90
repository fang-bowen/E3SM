module dynurbanparamsFileMod

  !---------------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Handle reading of the dataset that specifies transient urban properties
  ! (furban_timeseries), used when do_transient_urban_properties is .true.
  !
  ! Each year, the urban properties present on the file replace the corresponding
  ! values in urbinp (initially read from fsurdat), and the landunit-level urban
  ! parameters and their derived quantities (view factors, traffic factor,
  ! displacement height, roughness length) are recomputed for all urban landunits.
  ! Values jump to each year's value at the start of the year (no interpolation).
  !
  ! All variables are optional: a property that is not on the file keeps its fsurdat
  ! value for the whole run.
  !
  ! The following urban properties are intentionally NOT read from this file and stay
  ! fixed at their fsurdat values:
  !   - WTLUNIT_ROOF, WTROAD_PERV: set the column weights within the urban landunit;
  !     making them dynamic requires column reweighting with water/heat conservation
  !     (future work).
  !   - THICK_ROOF, THICK_WALL, NLEV_IMPROAD: set the urban column layer structure.
  !   - T_BUILDING_MIN, T_BUILDING_MAX: handled elsewhere (e.g., the urbantv stream).
  !
  ! Energy accounting: dynurbanparams_interp is called from dynSubgrid_driver between
  ! dyn_hwcontent_init and dyn_hwcontent_final, so a change in the heat capacity of the
  ! roof, wall and impervious road layers (CV_*) changes the column heat content and is
  ! booked by the dynbal heat flux (temperatures are kept unchanged).
  ! NOTE: the heat stored in the building interior (air and floor, prognostic building
  ! temperature scheme) is not part of the land energy accounting in ELM, so a change in
  ! HT_ROOF changes the building air volume without a corresponding correction.
  !
  ! Only supported for max_topounits == 1.
  !
  ! !USES:
#include "shr_assert.h"
  use shr_log_mod           , only : errMsg => shr_log_errMsg
  use shr_kind_mod          , only : r8 => shr_kind_r8
  use decompMod             , only : bounds_type, BOUNDS_LEVEL_PROC
  use dynFileMod            , only : dyn_file_type
  use dynVarTimeUninterpMod , only : dyn_var_time_uninterp_type
  use elm_varctl            , only : iulog
  use elm_varcon            , only : grlnd
  use elm_varpar            , only : numrad, nlevurb
  use abortutils            , only : endrun
  use spmdMod               , only : masterproc
  use landunit_varcon       , only : numurbl
  use topounit_varcon       , only : max_topounits
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  implicit none
  private
  save
  public :: dynurbanparams_init    ! initialize information read from furban_timeseries dataset
  public :: dynurbanparams_interp  ! update urban properties for the current year, if needed
  !
  ! ! PRIVATE TYPES
  integer, parameter :: kind_scalar = 1  ! data shape [gridcell, numurbl]
  integer, parameter :: kind_band   = 2  ! data shape [gridcell, numurbl, numrad]
  integer, parameter :: kind_layer  = 3  ! data shape [gridcell, numurbl, nlevurb]

  type :: urbprop_var_type
     character(len=32) :: varname = ' '     ! variable name on file
     integer           :: vkind   = 0       ! kind_scalar, kind_band or kind_layer
     logical           :: on_file = .false. ! whether the variable is on the file
     type(dyn_var_time_uninterp_type) :: dynvar
  end type urbprop_var_type

  integer, parameter :: nvars = 21
  type(urbprop_var_type) :: urbprop_vars(nvars)

  type(dyn_file_type), target :: dynurbanparams_file ! information for the file containing transient urban properties

  integer :: time_index_prev = -1   ! time index of the data used in the last update

  character(len=*), parameter, private :: sourcefile = &
       __FILE__
  !---------------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine dynurbanparams_init(bounds, dynurbanparams_filename)
    !
    ! !DESCRIPTION:
    ! Initialize dataset containing transient urban properties (position it to the
    ! right time samples that bound the initial model date)
    !
    ! !USES:
    use ncdio_pio      , only : check_dim, check_var, var_desc_t
    use dynTimeInfoMod , only : YEAR_POSITION_START_OF_TIMESTEP
    !
    ! !ARGUMENTS:
    type(bounds_type) , intent(in) :: bounds                  ! proc-level bounds
    character(len=*)  , intent(in) :: dynurbanparams_filename ! name of file containing transient urban properties
    !
    ! !LOCAL VARIABLES:
    integer          :: num_points   ! number of spatial points
    integer          :: n            ! variable index
    type(var_desc_t) :: vardesc      ! variable descriptor (unused)
    logical          :: readvar      ! whether the variable is on the file

    character(len=*), parameter :: subname = 'dynurbanparams_init'
    !-----------------------------------------------------------------------

    SHR_ASSERT(bounds%level == BOUNDS_LEVEL_PROC, subname // ': argument must be PROC-level bounds')

    if (max_topounits > 1) then
       write(iulog,*) subname//' ERROR: do_transient_urban_properties is not yet supported with max_topounits > 1'
       call endrun(msg=errMsg(sourcefile, __LINE__))
    end if

    if (masterproc) then
       write(iulog,*) 'Attempting to read transient urban properties data .....'
    end if

    ! Get the year from the START of the timestep, consistent with transient urban area.
    dynurbanparams_file = dyn_file_type(dynurbanparams_filename, YEAR_POSITION_START_OF_TIMESTEP)
    call check_dim(dynurbanparams_file, 'numurbl', numurbl)

    ! Morphology
    call set_var( 1, 'CANYON_HWR'     , kind_scalar)
    call set_var( 2, 'HT_ROOF'        , kind_scalar)
    call set_var( 3, 'WIND_HGT_CANYON', kind_scalar)
    ! Radiative
    call set_var( 4, 'EM_ROOF'        , kind_scalar)
    call set_var( 5, 'EM_IMPROAD'     , kind_scalar)
    call set_var( 6, 'EM_PERROAD'     , kind_scalar)
    call set_var( 7, 'EM_WALL'        , kind_scalar)
    call set_var( 8, 'ALB_ROOF_DIR'   , kind_band)
    call set_var( 9, 'ALB_ROOF_DIF'   , kind_band)
    call set_var(10, 'ALB_IMPROAD_DIR', kind_band)
    call set_var(11, 'ALB_IMPROAD_DIF', kind_band)
    call set_var(12, 'ALB_PERROAD_DIR', kind_band)
    call set_var(13, 'ALB_PERROAD_DIF', kind_band)
    call set_var(14, 'ALB_WALL_DIR'   , kind_band)
    call set_var(15, 'ALB_WALL_DIF'   , kind_band)
    ! Thermal
    call set_var(16, 'TK_ROOF'        , kind_layer)
    call set_var(17, 'TK_WALL'        , kind_layer)
    call set_var(18, 'TK_IMPROAD'     , kind_layer)
    call set_var(19, 'CV_ROOF'        , kind_layer)
    call set_var(20, 'CV_WALL'        , kind_layer)
    call set_var(21, 'CV_IMPROAD'     , kind_layer)

    num_points = (bounds%endg - bounds%begg + 1)

    do n = 1, nvars
       call check_var(dynurbanparams_file, trim(urbprop_vars(n)%varname), vardesc, readvar)
       urbprop_vars(n)%on_file = readvar
       if (.not. readvar) then
          if (masterproc) then
             write(iulog,*) subname//': ', trim(urbprop_vars(n)%varname), &
                  ' not on furban_timeseries; keeping fsurdat value'
          end if
          cycle
       end if

       select case (urbprop_vars(n)%vkind)
       case (kind_scalar)
          urbprop_vars(n)%dynvar = dyn_var_time_uninterp_type( &
               dyn_file = dynurbanparams_file, varname=trim(urbprop_vars(n)%varname), &
               dim1name=grlnd, conversion_factor=1._r8, &
               do_check_sums_equal_1=.false., data_shape=[num_points, numurbl])
       case (kind_band)
          call check_dim(dynurbanparams_file, 'numrad', numrad)
          urbprop_vars(n)%dynvar = dyn_var_time_uninterp_type( &
               dyn_file = dynurbanparams_file, varname=trim(urbprop_vars(n)%varname), &
               dim1name=grlnd, conversion_factor=1._r8, &
               do_check_sums_equal_1=.false., data_shape=[num_points, numurbl, numrad])
       case (kind_layer)
          call check_dim(dynurbanparams_file, 'nlevurb', nlevurb)
          urbprop_vars(n)%dynvar = dyn_var_time_uninterp_type( &
               dyn_file = dynurbanparams_file, varname=trim(urbprop_vars(n)%varname), &
               dim1name=grlnd, conversion_factor=1._r8, &
               do_check_sums_equal_1=.false., data_shape=[num_points, numurbl, nlevurb])
       end select
    end do

    time_index_prev = -1

  contains

    subroutine set_var(n, varname, vkind)
      integer          , intent(in) :: n
      character(len=*) , intent(in) :: varname
      integer          , intent(in) :: vkind
      urbprop_vars(n)%varname = varname
      urbprop_vars(n)%vkind   = vkind
    end subroutine set_var

  end subroutine dynurbanparams_init

  !-----------------------------------------------------------------------
  subroutine dynurbanparams_interp(bounds)
    !
    ! !DESCRIPTION:
    ! Update urban properties for the current model year, if the year's data differ
    ! from the last update.
    !
    ! !USES:
    use UrbanParamsType  , only : urbinp, urbanparams_vars, CheckUrban
    use LandunitType     , only : lun_pp
    use ColumnDataType   , only : col_es
    use column_varcon    , only : icol_roof, icol_sunwall, icol_shadewall
    use column_varcon    , only : icol_road_perv, icol_road_imperv
    use landunit_varcon  , only : isturb_MIN
    use GridcellType     , only : grc_pp
    use ColumnType       , only : col_pp
    use domainMod        , only : ldomain
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds  ! proc-level bounds
    !
    ! !LOCAL VARIABLES:
    integer               :: n, g, l, c, t, dindx   ! indices
    integer               :: time_index_cur         ! time index of the current data
    real(r8), allocatable :: data2d(:,:)            ! current data, scalar variables
    real(r8), allocatable :: data3d(:,:,:)          ! current data, band / layer variables
    real(r8), allocatable :: pcturb(:,:,:)          ! current % urban, for CheckUrban

    character(len=*), parameter :: subname = 'dynurbanparams_interp'
    !-----------------------------------------------------------------------

    SHR_ASSERT(bounds%level == BOUNDS_LEVEL_PROC, subname // ': argument must be PROC-level bounds')

    call dynurbanparams_file%time_info%set_current_year()

    ! Only update when the data year changes (urban properties are annual)
    time_index_cur = dynurbanparams_file%time_info%get_time_index_lower()
    if (time_index_cur == time_index_prev) return
    time_index_prev = time_index_cur

    ! ------------------------------------------------------------------------
    ! Update urbinp with the values on file (topounit index 1; max_topounits == 1)
    ! ------------------------------------------------------------------------
    do n = 1, nvars
       if (.not. urbprop_vars(n)%on_file) cycle

       select case (urbprop_vars(n)%vkind)
       case (kind_scalar)
          allocate(data2d(bounds%begg:bounds%endg, numurbl))
          call urbprop_vars(n)%dynvar%get_current_data(data2d)
          select case (trim(urbprop_vars(n)%varname))
          case ('CANYON_HWR')      ; urbinp%canyon_hwr     (:,1,:) = data2d(:,:)
          case ('HT_ROOF')         ; urbinp%ht_roof        (:,1,:) = data2d(:,:)
          case ('WIND_HGT_CANYON') ; urbinp%wind_hgt_canyon(:,1,:) = data2d(:,:)
          case ('EM_ROOF')         ; urbinp%em_roof        (:,1,:) = data2d(:,:)
          case ('EM_IMPROAD')      ; urbinp%em_improad     (:,1,:) = data2d(:,:)
          case ('EM_PERROAD')      ; urbinp%em_perroad     (:,1,:) = data2d(:,:)
          case ('EM_WALL')         ; urbinp%em_wall        (:,1,:) = data2d(:,:)
          end select
          deallocate(data2d)

       case (kind_band)
          allocate(data3d(bounds%begg:bounds%endg, numurbl, numrad))
          call urbprop_vars(n)%dynvar%get_current_data(data3d)
          select case (trim(urbprop_vars(n)%varname))
          case ('ALB_ROOF_DIR')    ; urbinp%alb_roof_dir   (:,1,:,:) = data3d(:,:,:)
          case ('ALB_ROOF_DIF')    ; urbinp%alb_roof_dif   (:,1,:,:) = data3d(:,:,:)
          case ('ALB_IMPROAD_DIR') ; urbinp%alb_improad_dir(:,1,:,:) = data3d(:,:,:)
          case ('ALB_IMPROAD_DIF') ; urbinp%alb_improad_dif(:,1,:,:) = data3d(:,:,:)
          case ('ALB_PERROAD_DIR') ; urbinp%alb_perroad_dir(:,1,:,:) = data3d(:,:,:)
          case ('ALB_PERROAD_DIF') ; urbinp%alb_perroad_dif(:,1,:,:) = data3d(:,:,:)
          case ('ALB_WALL_DIR')    ; urbinp%alb_wall_dir   (:,1,:,:) = data3d(:,:,:)
          case ('ALB_WALL_DIF')    ; urbinp%alb_wall_dif   (:,1,:,:) = data3d(:,:,:)
          end select
          deallocate(data3d)

       case (kind_layer)
          allocate(data3d(bounds%begg:bounds%endg, numurbl, nlevurb))
          call urbprop_vars(n)%dynvar%get_current_data(data3d)
          select case (trim(urbprop_vars(n)%varname))
          case ('TK_ROOF')    ; urbinp%tk_roof   (:,1,:,:) = data3d(:,:,:)
          case ('TK_WALL')    ; urbinp%tk_wall   (:,1,:,:) = data3d(:,:,:)
          case ('TK_IMPROAD') ; urbinp%tk_improad(:,1,:,:) = data3d(:,:,:)
          case ('CV_ROOF')    ; urbinp%cv_roof   (:,1,:,:) = data3d(:,:,:)
          case ('CV_WALL')    ; urbinp%cv_wall   (:,1,:,:) = data3d(:,:,:)
          case ('CV_IMPROAD') ; urbinp%cv_improad(:,1,:,:) = data3d(:,:,:)
          end select
          deallocate(data3d)
       end select
    end do

    ! ------------------------------------------------------------------------
    ! Check validity of the urban properties wherever an urban landunit has area
    ! ------------------------------------------------------------------------
    allocate(pcturb(bounds%begg:bounds%endg, max_topounits, numurbl))
    pcturb(:,:,:) = 0._r8
    do l = bounds%begl, bounds%endl
       if (lun_pp%urbpoi(l)) then
          g = lun_pp%gridcell(l)
          t = lun_pp%topounit(l) - grc_pp%topi(g) + 1
          dindx = lun_pp%itype(l) - isturb_MIN + 1
          pcturb(g,t,dindx) = lun_pp%wttopounit(l) * 100._r8
       end if
    end do
    call CheckUrban(bounds%begg, bounds%endg, pcturb(bounds%begg:bounds%endg,:,:), &
         subname, ldomain%num_tunits_per_grd)
    deallocate(pcturb)

    ! ------------------------------------------------------------------------
    ! Recompute landunit-level urban parameters and derived quantities
    ! ------------------------------------------------------------------------
    do l = bounds%begl, bounds%endl
       if (lun_pp%urbpoi(l)) then
          call urbanparams_vars%SetLunitParams(l)

          ! Urban ground emissivities are only set at initialization (CanopyTemperature
          ! does not update them for urban columns), so update them here.
          do c = lun_pp%coli(l), lun_pp%colf(l)
             if (col_pp%itype(c) == icol_roof       ) col_es%emg(c) = urbanparams_vars%em_roof(l)
             if (col_pp%itype(c) == icol_sunwall    ) col_es%emg(c) = urbanparams_vars%em_wall(l)
             if (col_pp%itype(c) == icol_shadewall  ) col_es%emg(c) = urbanparams_vars%em_wall(l)
             if (col_pp%itype(c) == icol_road_imperv) col_es%emg(c) = urbanparams_vars%em_improad(l)
             if (col_pp%itype(c) == icol_road_perv  ) col_es%emg(c) = urbanparams_vars%em_perroad(l)
          end do
       end if
    end do

  end subroutine dynurbanparams_interp

end module dynurbanparamsFileMod
