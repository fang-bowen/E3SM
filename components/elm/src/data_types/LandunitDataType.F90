module LandunitDataType

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Landunit data type allocation and initialization
  ! --------------------------------------------------------
  !
  use shr_kind_mod   , only : r8 => shr_kind_r8
  use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
  use elm_varpar     , only : nlevsno, nlevgrnd, nlevlak, nlevurb
  use elm_varcon     , only : spval, ispval
  use elm_varctl     , only : use_vancouver, use_mexicocity
  use histFileMod    , only : hist_addfld1d
  use ncdio_pio      , only : file_desc_t, ncd_double
  use decompMod      , only : bounds_type
  use restUtilMod
  use LandunitType   , only : lun_pp
  use ColumnType     , only : col_pp   
  use column_varcon  , only : icol_roof, icol_sunwall,icol_shadewall  
  use ColumnDataType , only : column_energy_state

  !
  ! !PUBLIC TYPES:
  implicit none
  save
  private

  !-----------------------------------------------------------------------
  ! Define the data structure that holds energy state information at the landunit level.
  !-----------------------------------------------------------------------
  type, public :: landunit_energy_state
    ! temperature variables
    real(r8), pointer :: t_building       (:)   ! internal building temperature (K)
    real(r8), pointer :: t_roof_inner     (:)   ! roof inside surface temperature (K)
    real(r8), pointer :: t_sunw_inner     (:)   ! sunwall inside surface temperature (K)
    real(r8), pointer :: t_shdw_inner     (:)   ! shadewall inside surface temperature (K)
    real(r8), pointer :: t_floor          (:)   ! floor temperature (K)

    real(r8), pointer :: taf              (:)   ! urban canopy air temperature (K)

  contains
    procedure, public :: Init    => lun_es_init
     procedure, public :: Restart => lun_es_restart
    procedure, public :: Clean   => lun_es_clean
  end type landunit_energy_state

  !-----------------------------------------------------------------------
  ! Define the data structure that holds energy flux information at the landunit level.
  !-----------------------------------------------------------------------
  type, public :: landunit_energy_flux
    ! temperature variables
    real(r8), pointer :: eflx_traffic      (:)   ! traffic sensible heat flux (W/m**2)
    real(r8), pointer :: eflx_wasteheat    (:)   ! sensible heat flux from domestic heating/cooling sources of waste heat (W/m**2)
    real(r8), pointer :: eflx_ventilation  (:)   ! lun sensible heat flux from building ventilation (W/m**2)  
    real(r8), pointer :: eflx_urban_ac     (:)   ! urban air conditioning flux (W/m**2) 
    real(r8), pointer :: eflx_urban_ac_sen (:)   ! sensible heat component of air conditioning flux (W/m**2)
    real(r8), pointer :: eflx_heat_from_ac (:)   ! sensible heat flux to be put back into canyon due to removal by AC (W/m**2)
    real(r8), pointer :: eflx_urban_heat   (:)   ! urban heating flux (W/m**2)
    real(r8), pointer :: eflx_building     (:)   ! building heat flux from change in interior building air temperature (W/m**2)

  contains
    procedure, public :: Init    => lun_ef_init
    procedure, public :: Restart    => lun_ef_restart
    procedure, public :: Clean   => lun_ef_clean
  end type landunit_energy_flux

  !-----------------------------------------------------------------------
  ! Define the data structure that holds water state information at the landunit level.
  !-----------------------------------------------------------------------
  type, public :: landunit_water_state
    ! temperature variables
    real(r8), pointer :: qaf               (:)   ! urban canopy air specific humidity (kg H2O/kg moist air)
    real(r8), pointer :: q_building        (:)   ! lun internal building air specific humidity (kg/kg)
    real(r8), pointer :: rh_building       (:)   ! lun internal building air relative humidity (%) 

  contains
    procedure, public :: Init    => lun_ws_init
     procedure, public :: Restart => lun_ws_restart
    procedure, public :: Clean   => lun_ws_clean
  end type landunit_water_state


    !-----------------------------------------------------------------------
  ! Define the data structure that holds water flux information at the landunit level.
  !-----------------------------------------------------------------------
  type, public :: landunit_water_flux
    ! temperature variables
    real(r8), pointer :: qflx_condensate_from_ac    (:)   ! urban canopy air specific humidity (kg H2O/kg moist air)

  contains
    procedure, public :: Init    => lun_wf_init
    procedure, public :: Restart => lun_wf_restart
    procedure, public :: Clean   => lun_wf_clean
  end type landunit_water_flux

  !-----------------------------------------------------------------------
  ! declare the public instances of landunit-level data types
  !-----------------------------------------------------------------------
  type(landunit_energy_state)  , public, target :: lun_es    ! landunit energy state
  type(landunit_energy_flux)  , public, target :: lun_ef    ! landunit energy flux
  type(landunit_water_state)  , public, target :: lun_ws    ! landunit water state
  type(landunit_water_flux)  , public, target :: lun_wf    ! landunit water flux
  !$acc declare create(lun_es)
  !$acc declare create(lun_ef)
  !$acc declare create(lun_ws)
  !$acc declare create(lun_wf)
  !------------------------------------------------------------------------

contains

  !------------------------------------------------------------------------
  ! Subroutines to initialize and clean landunit energy state data structure
  !------------------------------------------------------------------------
  subroutine lun_es_init(this, begl, endl,col_es, is_simple_buildtemp, is_prog_buildtemp)
    !
    ! !ARGUMENTS:
    class(landunit_energy_state) :: this
    integer, intent(in) :: begl,endl
    logical, intent(in) :: is_simple_buildtemp  ! Simple building temp is being used
    logical, intent(in) :: is_prog_buildtemp    ! Prognostic building temp is being used
    type(column_energy_state), intent(in) :: col_es ! column energy state
    !------------------------------------------------------------------------
    ! !LOCAL VARIABLES:
    integer :: l                        ! indices
    integer :: c
    character(100)    :: lname

    !-----------------------------------------------------------------------
    ! allocate for each member of lun_es
    !-----------------------------------------------------------------------
    allocate(this%t_building            (begl:endl))                      ; this%t_building            (:)   = spval
    allocate(this%t_roof_inner          (begl:endl))                      ; this%t_roof_inner          (:)   = spval
    allocate(this%t_sunw_inner          (begl:endl))                      ; this%t_sunw_inner          (:)   = spval
    allocate(this%t_shdw_inner          (begl:endl))                      ; this%t_shdw_inner          (:)   = spval
    allocate(this%t_floor               (begl:endl))                      ; this%t_floor               (:)   = spval
    allocate(this%taf                   (begl:endl))                      ; this%taf                   (:)   = spval

    !-----------------------------------------------------------------------
    ! initialize history fields for select members of lun_es
    !-----------------------------------------------------------------------
    this%t_building(begl:endl) = spval
    if (      is_simple_buildtemp )then
       lname = 'internal urban building temperature'
    else if ( is_prog_buildtemp   )then
       lname = 'internal urban building air temperature'
    end if
    call hist_addfld1d(fname='TBUILD', units='K',  &
         avgflag='A', long_name=lname, &
         ptr_lunit=this%t_building, set_nourb=spval, l2g_scale_type='unity')

    this%t_roof_inner(begl:endl) = spval
    call hist_addfld1d(fname='TROOF_INNER', units='K',  &
            avgflag='A', long_name='roof inside surface temperature', &
            ptr_lunit=this%t_roof_inner, set_nourb=spval, l2g_scale_type='unity', &
            default='inactive')

    this%t_sunw_inner(begl:endl) = spval
    call hist_addfld1d(fname='TSUNW_INNER', units='K',  &
            avgflag='A', long_name='sunwall inside surface temperature', &
            ptr_lunit=this%t_sunw_inner, set_nourb=spval, l2g_scale_type='unity', &
            default='inactive')

    this%t_shdw_inner(begl:endl) = spval
    call hist_addfld1d(fname='TSHDW_INNER', units='K',  &
            avgflag='A', long_name='shadewall inside surface temperature', &
            ptr_lunit=this%t_shdw_inner, set_nourb=spval, l2g_scale_type='unity', &
            default='inactive')   
         
    if ( is_prog_buildtemp )then

       this%t_floor(begl:endl) = spval
       call hist_addfld1d(fname='TFLOOR', units='K',  &
            avgflag='A', long_name='floor temperature', &
            ptr_lunit=this%t_floor, set_nourb=spval, l2g_scale_type='unity', &
            default='inactive')
    end if       

    !-----------------------------------------------------------------------
    ! cold-start initial conditions for lun_es
    !-----------------------------------------------------------------------
    do l = begl, endl
       if (lun_pp%urbpoi(l)) then
          if (use_vancouver) then
             this%taf(l) = 297.56_r8
          else if (use_mexicocity) then
             this%taf(l) = 289.46_r8
          else
             this%taf(l) = 283._r8
          end if
       end if
    end do

    ! Initialize internal building temperature, inner temperatures of building
    ! surfaces, and floor temperature
    if ( is_prog_buildtemp )then
      do l = begl, endl
        do  c = lun_pp%coli(l),lun_pp%colf(l)
          if (col_pp%itype(c) == icol_roof)  then
            this%t_roof_inner(l) = col_es%t_soisno(c,nlevurb)
            this%t_building(l)   = col_es%t_soisno(c,nlevurb)        ! arbitrarily set to roof temperature
            this%t_floor(l)  = col_es%t_soisno(c,nlevurb)        ! arbitrarily set to roof temperature
          else if (col_pp%itype(c) == icol_sunwall) then
            this%t_sunw_inner(l) = col_es%t_soisno(c,nlevurb)
          else if (col_pp%itype(c) == icol_shadewall) then
            this%t_shdw_inner(l) = col_es%t_soisno(c,nlevurb)
          end if
        end do
      end do
   end if


  end subroutine lun_es_init

  !------------------------------------------------------------------------
   subroutine lun_es_restart(this, bounds, ncid, flag, is_simple_buildtemp, is_prog_buildtemp)
     !
     ! !DESCRIPTION:
     ! Read/Write landunit energy state information to/from restart file.
     !
     ! !USES:
     use spmdMod       , only : masterproc
     use elm_varctl    , only : iulog
     !
     ! !ARGUMENTS:
     class(landunit_energy_state) :: this
     type(bounds_type), intent(in)    :: bounds
     type(file_desc_t), intent(inout) :: ncid
     character(len=*) , intent(in)    :: flag
     logical, intent(in) :: is_simple_buildtemp  ! Simple building temp is being used
     logical, intent(in) :: is_prog_buildtemp    ! Prognostic building temp is being used
     !
     ! !LOCAL VARIABLES:
     logical :: readvar   ! determine if variable is on initial file
     !-----------------------------------------------------------------------

     call restartvar(ncid=ncid, flag=flag, varname='taf', xtype=ncd_double, dim1name='landunit',                       &
          long_name='urban canopy air temperature', units='K',                                                         &
          interpinic_flag='interp', readvar=readvar, data=this%taf)
     

     if ( is_prog_buildtemp )then
      ! landunit type physical state variable - t_building
      call restartvar(ncid=ncid, flag=flag, varname='t_building', xtype=ncd_double,  &
           dim1name='landunit', &
           long_name='internal building air temperature', units='K', &
           interpinic_flag='interp', readvar=readvar, data=this%t_building)
      if (flag=='read' .and. .not. readvar) then
         if (masterproc) write(iulog,*) "can't find t_building in initial file..."
         if (masterproc) write(iulog,*) "Initialize t_building to taf"
         this%t_building(bounds%begl:bounds%endl) = this%taf(bounds%begl:bounds%endl)
      end if

      ! landunit type physical state variable - t_roof_inner
      call restartvar(ncid=ncid, flag=flag, varname='t_roof_inner', xtype=ncd_double,  &
           dim1name='landunit', &
           long_name='roof inside surface temperature', units='K', &
           interpinic_flag='interp', readvar=readvar, data=this%t_roof_inner)
      if (flag=='read' .and. .not. readvar) then
         if (masterproc) write(iulog,*) "can't find t_roof_inner in initial file..."
         if (masterproc) write(iulog,*) "Initialize t_roof_inner to taf"
         this%t_roof_inner(bounds%begl:bounds%endl) = this%taf(bounds%begl:bounds%endl)
      end if

      ! landunit type physical state variable - t_sunw_inner
      call restartvar(ncid=ncid, flag=flag, varname='t_sunw_inner', xtype=ncd_double,  &
           dim1name='landunit', &
           long_name='sunwall inside surface temperature', units='K', &
           interpinic_flag='interp', readvar=readvar, data=this%t_sunw_inner)
      if (flag=='read' .and. .not. readvar) then
         if (masterproc) write(iulog,*) "can't find t_sunw_inner in initial file..."
         if (masterproc) write(iulog,*) "Initialize t_sunw_inner to taf"
         this%t_sunw_inner(bounds%begl:bounds%endl) = this%taf(bounds%begl:bounds%endl)
      end if

      ! landunit type physical state variable - t_shdw_inner
      call restartvar(ncid=ncid, flag=flag, varname='t_shdw_inner', xtype=ncd_double,  &
           dim1name='landunit', &
           long_name='shadewall inside surface temperature', units='K', &
           interpinic_flag='interp', readvar=readvar, data=this%t_shdw_inner)
      if (flag=='read' .and. .not. readvar) then
         if (masterproc) write(iulog,*) "can't find t_shdw_inner in initial file..."
         if (masterproc) write(iulog,*) "Initialize t_shdw_inner to taf"
         this%t_shdw_inner(bounds%begl:bounds%endl) = this%taf(bounds%begl:bounds%endl)
      end if

      ! landunit type physical state variable - t_floor
      call restartvar(ncid=ncid, flag=flag, varname='t_floor', xtype=ncd_double,  &
           dim1name='landunit', &
           long_name='floor temperature', units='K', &
           interpinic_flag='interp', readvar=readvar, data=this%t_floor)
      if (flag=='read' .and. .not. readvar) then
         if (masterproc) write(iulog,*) "can't find t_floor in initial file..."
         if (masterproc) write(iulog,*) "Initialize t_floor to taf"
         this%t_floor(bounds%begl:bounds%endl) = this%taf(bounds%begl:bounds%endl)
      end if
   end if
  end subroutine lun_es_restart

  !------------------------------------------------------------------------
  subroutine lun_es_clean(this)
    !
    ! !ARGUMENTS:
    class(landunit_energy_state) :: this
    !------------------------------------------------------------------------
    deallocate(this%t_building)
    deallocate(this%t_roof_inner)
    deallocate(this%t_sunw_inner)
    deallocate(this%t_shdw_inner)
    deallocate(this%t_floor)
    deallocate(this%taf)

  end subroutine lun_es_clean

  !------------------------------------------------------------------------
  ! Subroutines to initialize, restart and clean landunit energy flux data structure
  !------------------------------------------------------------------------
  subroutine lun_ef_init(this, begl, endl, is_simple_buildtemp, is_prog_buildtemp)
    !
    ! !ARGUMENTS:
    class(landunit_energy_flux) :: this
    integer, intent(in) :: begl,endl
    logical           , intent(in) :: is_simple_buildtemp ! If using simple building temp method
    logical           , intent(in) :: is_prog_buildtemp   ! If using prognostic building temp method
    !------------------------------------------------------------------------
    ! !LOCAL VARIABLES:
    integer :: l                        ! indices

    !-----------------------------------------------------------------------
    ! allocate for each member of lun_ef
    !-----------------------------------------------------------------------
    allocate( this%eflx_heat_from_ac   (begl:endl))             ; this%eflx_heat_from_ac   (:)   = spval
    allocate( this%eflx_traffic        (begl:endl))             ; this%eflx_traffic        (:)   = spval
    allocate( this%eflx_wasteheat      (begl:endl))             ; this%eflx_wasteheat      (:)   = spval
    allocate( this%eflx_ventilation    (begl:endl))             ; this%eflx_ventilation    (:)   = spval
    allocate( this%eflx_urban_ac       (begl:endl))             ; this%eflx_urban_ac       (:)   = spval
    allocate( this%eflx_urban_ac_sen   (begl:endl))             ; this%eflx_urban_ac_sen   (:)   = spval
    allocate( this%eflx_urban_heat     (begl:endl))             ; this%eflx_urban_heat     (:)   = spval
    allocate( this%eflx_building       (begl:endl))             ; this%eflx_building       (:)   = spval

    !-----------------------------------------------------------------------
    ! initialize history fields for select members of lun_ef
    !-----------------------------------------------------------------------
    if (is_prog_buildtemp) then
      this%eflx_urban_ac(begl:endl) = spval
      call hist_addfld1d (fname='EFLXBUILD', units='W/m^2',  &
               avgflag='A', long_name='building heat flux from change in interior building air temperature', &
               ptr_lunit=this%eflx_building, set_nourb=0._r8, l2g_scale_type='unity')

      this%eflx_urban_ac(begl:endl) = spval
      call hist_addfld1d (fname='URBAN_AC', units='W/m^2',  &
               avgflag='A', long_name='urban air conditioning flux', &
               ptr_lunit=this%eflx_urban_ac, set_nourb=0._r8, l2g_scale_type='unity')

      this%eflx_urban_ac_sen(begl:endl) = spval
      call hist_addfld1d (fname='URBAN_AC_SEN', units='W/m^2',  &
            avgflag='A', long_name='sensible heat component of urban air conditioning flux', &
            ptr_lunit=this%eflx_urban_ac_sen, set_nourb=0._r8, l2g_scale_type='unity',default='inactive')

      this%eflx_urban_heat(begl:endl) = spval
      call hist_addfld1d (fname='URBAN_HEAT', units='W/m^2',  &
               avgflag='A', long_name='urban heating flux', &
               ptr_lunit=this%eflx_urban_heat, set_nourb=0._r8, l2g_scale_type='unity')
    end if


    
    !-----------------------------------------------------------------------
    ! cold-start initial conditions for lun_ef
    !-----------------------------------------------------------------------
    do l = begl, endl
       if (.not. lun_pp%urbpoi(l)) then
          this%eflx_traffic(l)   = spval
          this%eflx_wasteheat(l) = spval
          this%eflx_ventilation(l) = spval
          if ( is_prog_buildtemp )then
             this%eflx_building(l)   = 0._r8
             this%eflx_urban_ac(l)   = 0._r8
             this%eflx_urban_ac_sen(l)= 0._r8
             this%eflx_urban_heat(l) = 0._r8
          end if
       else
          if ( is_prog_buildtemp )then
             this%eflx_building(l)   = 0._r8
             this%eflx_urban_ac(l)   = 0._r8
             this%eflx_urban_ac_sen(l)= 0._r8
             this%eflx_urban_heat(l) = 0._r8
             this%eflx_ventilation(l)= 0._r8
          end if

       end if
    end do


  end subroutine lun_ef_init

    !------------------------------------------------------------------------
   subroutine lun_ef_restart(this, bounds, ncid, flag, is_simple_buildtemp, is_prog_buildtemp)
     !
     ! !DESCRIPTION:
     ! Read/Write landunit energy flux information to/from restart file.
     !
     ! !USES:
     use decompMod      , only : get_proc_global
     use ncdio_pio      , only : file_desc_t, ncd_double, ncd_inqvdlen 
     use spmdMod       , only : masterproc
     use elm_varctl    , only : iulog
     !
     ! !ARGUMENTS:
     class(landunit_energy_flux) :: this
     type(bounds_type), intent(in)    :: bounds
     type(file_desc_t), intent(inout) :: ncid
     character(len=*) , intent(in)    :: flag
     logical, intent(in) :: is_simple_buildtemp  ! Simple building temp is being used
     logical, intent(in) :: is_prog_buildtemp    ! Prognostic building temp is being used
     !
     ! !LOCAL VARIABLES:
     logical :: readvar   ! determine if variable is on initial file
     logical :: do_io
     integer :: dimlen
     integer :: numl_global
     integer :: err_code
     !-----------------------------------------------------------------------

    ! Restart for building air temperature method
    call get_proc_global(nl=numl_global)
    if ( is_prog_buildtemp )then
       ! landunit urban energy state variable - eflx_urban_ac
       do_io = .true.
       ! On a read, confirm that this variable has the expected size (landunit-level); if not, 
       ! don't read it (instead give it a default value). This is needed to support older initial
       ! conditions for which this variable had a different size (column-level).
       if (flag == 'read') then
          call ncd_inqvdlen(ncid, 'URBAN_AC_L', 1, dimlen, err_code)
          if (dimlen /= numl_global) then
             do_io = .false.
             readvar = .false.
          end if
       end if
       if (do_io) then
          call restartvar(ncid=ncid, flag=flag, varname='URBAN_AC_L', xtype=ncd_double,  &
               dim1name='landunit',&
               long_name='urban air conditioning flux', units='watt/m^2', &
               interpinic_flag='interp', readvar=readvar, data=this%eflx_urban_ac)
       else
          this%eflx_urban_ac = 0.0_r8
       end if
       ! landunit urban energy state variable - eflx_urban_heat
       do_io = .true.
       ! On a read, confirm that this variable has the expected size (landunit-level); if not, 
       ! don't read it (instead give it a default value). This is needed to support older initial
       ! conditions for which this variable had a different size (column-level).
       if (flag == 'read') then
          call ncd_inqvdlen(ncid, 'URBAN_HEAT_L', 1, dimlen, err_code)
          if (dimlen /= numl_global) then
             do_io = .false.
             readvar = .false.
          end if
       end if
       if (do_io) then
          call restartvar(ncid=ncid, flag=flag, varname='URBAN_HEAT_L', xtype=ncd_double,  &
               dim1name='landunit',&
               long_name='urban heating flux', units='watt/m^2', &
               interpinic_flag='interp', readvar=readvar, data=this%eflx_urban_heat)
       else
          this%eflx_urban_heat = 0.0_r8
       end if

       call restartvar(ncid=ncid, flag=flag, varname='EFLX_VENTILATION', xtype=ncd_double, &
           dim1name='landunit', &
           long_name='sensible heat flux from building ventilation', units='watt/m^2', &
           interpinic_flag='interp', readvar=readvar, data=this%eflx_ventilation)
       if (flag=='read' .and. .not. readvar) then
          if (masterproc) write(iulog,*) "can't find EFLX_VENTILATION in initial file..."
          if (masterproc) write(iulog,*) "Initialize EFLX_VENTILATION to zero"
          this%eflx_ventilation(bounds%begl:bounds%endl) = 0._r8
       end if 
    end if

  end subroutine lun_ef_restart

  !------------------------------------------------------------------------
  subroutine lun_ef_clean(this)
    !
    ! !ARGUMENTS:
    class(landunit_energy_flux) :: this
    !------------------------------------------------------------------------
    deallocate(this%eflx_heat_from_ac)
    deallocate(this%eflx_traffic)
    deallocate(this%eflx_wasteheat)
    deallocate(this%eflx_urban_ac)
    deallocate(this%eflx_urban_heat)
    deallocate(this%eflx_building)

  end subroutine lun_ef_clean

  !------------------------------------------------------------------------
  ! Subroutines to initialize and clean landunit water state data structure
  !------------------------------------------------------------------------
  subroutine lun_ws_init(this, begl, endl, is_prog_buildtemp)
    !
    ! !ARGUMENTS:
    class(landunit_water_state) :: this
    integer, intent(in) :: begl,endl
    logical, intent(in) :: is_prog_buildtemp    ! Prognostic building temp is being used
    !------------------------------------------------------------------------
    ! !LOCAL VARIABLES:
    integer :: l                        ! indices

    !-----------------------------------------------------------------------
    ! allocate for each member of lun_ws
    !-----------------------------------------------------------------------
    allocate(this%qaf          (begl:endl))               ; this%qaf         (:)   = spval
    allocate(this%q_building   (begl:endl))               ; this%q_building  (:)   = spval   ! is it necessary to initialize with spval ?
    allocate(this%rh_building  (begl:endl))               ; this%rh_building (:)   = spval     ! ! is it necessary to initialize with spval ? it was NaN cathy CESM code

    !-----------------------------------------------------------------------
    ! initialize history fields for select members of lun_ws
    !-----------------------------------------------------------------------
    if ( is_prog_buildtemp ) then
       this%qaf(begl:endl) = spval     ! it was already initialized with spval above, so may not be neccesary
       call hist_addfld1d ( &
            fname='QAF_U', &
            units='kg/kg',  &
            avgflag='A', &
            long_name='urban canopy specific humidity', &
            ptr_lunit=this%qaf, l2g_scale_type='unity', set_nourb=spval, default='inactive')        


       this%q_building(begl:endl) = spval     ! it was already initialized with spval above, so may not be neccesary
       call hist_addfld1d ( &
            fname='QBUILD', &
            units='kg/kg',  &
            avgflag='A', &
            long_name='internal urban building air specific humidity', &
            ptr_lunit=this%q_building, l2g_scale_type='unity', set_nourb=spval)   !set_nourb=0._r8  instead of set_nourb=spval, so in agreement with URBAN_AC and URBAN_HEAT
       
       this%rh_building(begl:endl) = spval    ! it was already initialized with spval above, so may not be neccesary
       call hist_addfld1d ( &
            fname='RHBUILD', &                     ! in CESM code this variable does not depend on is_prog_buildtemp
            units='%',  &
            avgflag='A', &
            long_name='Internal urban building air relative humidity', &
            ptr_lunit=this%rh_building, set_nourb=spval, default='inactive')

    end if


    !-----------------------------------------------------------------------
    ! cold-start initial conditions for lun_ws
    !-----------------------------------------------------------------------
    do l = begl, endl
       if (lun_pp%urbpoi(l)) then
          if (use_vancouver) then
             this%qaf(l) = 0.0111_r8
          else if (use_mexicocity) then
             this%qaf(l) = 0.00248_r8
          else
             this%qaf(l) = 1.e-4_r8 ! Arbitrary set since forc_q is not yet available
          end if
       end if
    end do

    ! Initialize internal building specific humidity (following example above and t_building_max in TemperatureType.F90)
    if ( is_prog_buildtemp ) then
       do l = begl, endl
          if (lun_pp%urbpoi(l)) then
             this%q_building(l) = this%qaf(l) ! set to urban canopy specific humidity
          end if
       end do
    end if
    


  end subroutine lun_ws_init

  !------------------------------------------------------------------------
   subroutine lun_ws_restart(this, bounds, ncid, flag, is_prog_buildtemp)
     !
     ! !DESCRIPTION:
     ! Read/Write landunit water state information to/from restart file.
     !
     ! !USES:
     use spmdMod       , only : masterproc
     use elm_varctl    , only : iulog
     !
     ! !ARGUMENTS:
     class(landunit_water_state)      :: this
     type(bounds_type), intent(in)    :: bounds
     type(file_desc_t), intent(inout) :: ncid
     character(len=*) , intent(in)    :: flag
     logical, intent(in) :: is_prog_buildtemp    ! Prognostic building temp is being used
     !
     ! !LOCAL VARIABLES:
     logical :: readvar   ! determine if variable is on initial file
     !-----------------------------------------------------------------------
      
     call restartvar(ncid=ncid, flag=flag, varname='qaf', xtype=ncd_double, dim1name='landunit',                       &
          long_name='urban canopy specific humidity', units='kg/kg',                                                   &
          interpinic_flag='interp', readvar=readvar, data=this%qaf)

     if ( is_prog_buildtemp ) then
       ! landunit type physical state variable - q_building
       call restartvar(ncid=ncid, flag=flag, &
            varname='q_building', &
            xtype=ncd_double, dim1name='landunit', &
            long_name='internal building air specific humidity', &
            units='kg/kg', &
            interpinic_flag='interp', readvar=readvar, data=this%q_building)
       if (flag=='read' .and. .not. readvar) then
          if (masterproc) write(iulog,*) "can't find q_building in initial file..."
          if (masterproc) write(iulog,*) "Initialize q_building to qaf"
          this%q_building(bounds%begl:bounds%endl) = this%qaf(bounds%begl:bounds%endl)
       end if
    end if

     end subroutine lun_ws_restart

  !------------------------------------------------------------------------
  subroutine lun_ws_clean(this)
    !
    ! !ARGUMENTS:
    class(landunit_water_state) :: this
    !------------------------------------------------------------------------
    deallocate(this%qaf)
    deallocate(this%q_building)
    deallocate(this%rh_building)

  end subroutine lun_ws_clean




  !------------------------------------------------------------------------
  ! Subroutines to initialize and clean landunit water flux data structure
  !------------------------------------------------------------------------
  subroutine lun_wf_init(this, begl, endl, is_prog_buildtemp)
    !
    ! !ARGUMENTS:
    class(landunit_water_flux) :: this
    integer, intent(in) :: begl,endl
    logical, intent(in) :: is_prog_buildtemp    ! Prognostic building temp is being used
    !------------------------------------------------------------------------
    ! !LOCAL VARIABLES:
    integer :: l                        ! indices

    !-----------------------------------------------------------------------
    ! allocate for each member of lun_wf
    !-----------------------------------------------------------------------
    allocate(this%qflx_condensate_from_ac (begl:endl))         ; this%qflx_condensate_from_ac(:)   = 0.0_r8 ! REMOVE COMMENT !spval     ! ! is it necessary to initialize with spval ? it was NaN cathy CESM code

    !-----------------------------------------------------------------------
    ! initialize history fields for select members of lun_wf
    !-----------------------------------------------------------------------
    if (is_prog_buildtemp) then     
      this%qflx_condensate_from_ac(begl:endl) = 0.0_r8   
      call hist_addfld1d ( &
            fname='QCOND_FROM_AC_LUN', &
            units='mm/s',  &
            avgflag='A', &
            long_name='Condensed water flux from AC dehumidification (lun var)', &
            ptr_lunit=this%qflx_condensate_from_ac, set_nourb=0.0_r8, l2g_scale_type='unity') !, default='inactive')
    end if


    !-----------------------------------------------------------------------
    ! cold-start initial conditions for lun_wf
    !-----------------------------------------------------------------------    


  end subroutine lun_wf_init

  !------------------------------------------------------------------------
   subroutine lun_wf_restart(this, bounds, ncid, flag, is_prog_buildtemp)
     !
     ! !DESCRIPTION:
     ! Read/Write landunit water flux information to/from restart file.
     !
     ! !USES:
     use spmdMod       , only : masterproc
     use elm_varctl    , only : iulog
     !
     ! !ARGUMENTS:
     class(landunit_water_flux)      :: this
     type(bounds_type), intent(in)    :: bounds
     type(file_desc_t), intent(inout) :: ncid
     character(len=*) , intent(in)    :: flag
     logical, intent(in) :: is_prog_buildtemp    ! Prognostic building temp is being used
     !
     ! !LOCAL VARIABLES:
     logical :: readvar   ! determine if variable is on initial file
     !-----------------------------------------------------------------------
      
  end subroutine lun_wf_restart

  !------------------------------------------------------------------------
  subroutine lun_wf_clean(this)
    !
    ! !ARGUMENTS:
    class(landunit_water_flux) :: this
    !------------------------------------------------------------------------
    deallocate(this%qflx_condensate_from_ac)

  end subroutine lun_wf_clean


end module LandunitDataType
