module dynurbanFileMod

  !---------------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Handle reading of the dataset that specifies transient areas of the urban landunits
  ! (PCT_URBAN, one value per urban density class) on the flanduse_timeseries file.
  !
  ! Ported from CTSM (ctsm5.1.dev069 / dev083) and adapted to ELM's topounit structure.
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
  use abortutils            , only : endrun
  use spmdMod               , only : masterproc
  use landunit_varcon       , only : numurbl
  use topounit_varcon       , only : max_topounits
  use UrbanParamsType       , only : CheckUrban
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  implicit none
  private
  save
  public :: dynurban_init    ! initialize information read from landuse.timeseries dataset
  public :: dynurban_interp  ! get urban data for the current time step, if needed
  !
  ! ! PRIVATE TYPES
  type(dyn_file_type), target      :: dynurban_file ! information for the file containing transient urban data
  type(dyn_var_time_uninterp_type) :: wturban       ! weight of the urban landunits

  ! Names of variables on file
  character(len=*), parameter :: urban_varname = 'PCT_URBAN'

  character(len=*), parameter, private :: sourcefile = &
       __FILE__
  !---------------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine dynurban_init(bounds, dynurban_filename)
    !
    ! !DESCRIPTION:
    ! Initialize dataset containing transient urban info (position it to the right time
    ! samples that bound the initial model date)
    !
    ! !USES:
    use ncdio_pio      , only : check_dim
    use dynTimeInfoMod , only : YEAR_POSITION_START_OF_TIMESTEP
    !
    ! !ARGUMENTS:
    type(bounds_type) , intent(in) :: bounds            ! proc-level bounds
    character(len=*)  , intent(in) :: dynurban_filename ! name of file containing transient urban information
    !
    ! !LOCAL VARIABLES:
    integer :: num_points     ! number of spatial points
    integer :: wturb_shape(3) ! shape of the wturb data

    character(len=*), parameter :: subname = 'dynurban_init'
    !-----------------------------------------------------------------------

    SHR_ASSERT(bounds%level == BOUNDS_LEVEL_PROC, subname // ': argument must be PROC-level bounds')

    ! Transient urban is currently only supported for surface datasets without multiple
    ! topounits per gridcell.
    if (max_topounits > 1) then
       write(iulog,*) subname//' ERROR: do_transient_urban is not yet supported with max_topounits > 1'
       call endrun(msg=errMsg(sourcefile, __LINE__))
    end if

    if (masterproc) then
       write(iulog,*) 'Attempting to read urban dynamic landuse data .....'
    end if

    ! Get the year from the START of the timestep; this way, we'll update urban areas
    ! starting after the year boundary. This is consistent with the timing of crop updates.
    dynurban_file = dyn_file_type(dynurban_filename, YEAR_POSITION_START_OF_TIMESTEP)
    call check_dim(dynurban_file, 'numurbl', numurbl)

    ! Note: if you want to change transient urban so that it is interpolated, rather
    ! than jumping to each year's value on Jan 1 of that year, simply change wturban to
    ! be of type dyn_var_time_interp_type (rather than dyn_var_time_uninterp_type), and
    ! change the following constructor accordingly.
    num_points = (bounds%endg - bounds%begg + 1)
    wturb_shape = [num_points, max_topounits, numurbl]
    wturban = dyn_var_time_uninterp_type( &
         dyn_file = dynurban_file, varname=urban_varname, &
         dim1name=grlnd, conversion_factor=100._r8, &
         do_check_sums_equal_1=.false., data_shape=wturb_shape)

  end subroutine dynurban_init

  !-----------------------------------------------------------------------
  subroutine dynurban_interp(bounds)
    !
    ! !DESCRIPTION:
    ! Get urban cover for model time, when needed.
    !
    ! Sets lun_pp%wttopounit for the urban landunits (TBD, HD, MD).
    !
    ! Note that urban cover currently jumps to its new value at the start of the year.
    !
    ! !USES:
    use landunit_varcon   , only : isturb_MIN
    use subgridWeightsMod , only : set_landunit_weight
    use GridcellType      , only : grc_pp
    use domainMod         , only : ldomain
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds  ! proc-level bounds
    !
    ! !LOCAL VARIABLES:
    integer               :: g, t, t2, n          ! indices
    real(r8), allocatable :: wturban_cur(:,:,:)   ! current weight of the urban landunits

    character(len=*), parameter :: subname = 'dynurban_interp'
    !-----------------------------------------------------------------------

    SHR_ASSERT(bounds%level == BOUNDS_LEVEL_PROC, subname // ': argument must be PROC-level bounds')

    call dynurban_file%time_info%set_current_year()

    allocate(wturban_cur(bounds%begg:bounds%endg, max_topounits, numurbl))
    call wturban%get_current_data(wturban_cur)
    do g = bounds%begg, bounds%endg
       do t = grc_pp%topi(g), grc_pp%topf(g)
          t2 = t - grc_pp%topi(g) + 1
          do n = 1, numurbl
             ! density index n = 1,2,3 maps to isturb_tbd, isturb_hd, isturb_md
             call set_landunit_weight(t, isturb_MIN + n - 1, wturban_cur(g,t2,n))
          end do
       end do
    end do

    call CheckUrban(bounds%begg, bounds%endg, wturban_cur(bounds%begg:bounds%endg,:,:), &
         subname, ldomain%num_tunits_per_grd)

    deallocate(wturban_cur)

  end subroutine dynurban_interp

end module dynurbanFileMod
