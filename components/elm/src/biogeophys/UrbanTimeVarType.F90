module UrbanTimeVarType

  !------------------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Urban Time Varying Data
  !
  ! !USES:
  use shr_kind_mod    , only : r8 => shr_kind_r8, CL => shr_kind_CL
  use shr_log_mod     , only : errMsg => shr_log_errMsg
  use abortutils      , only : endrun
  use decompMod       , only : bounds_type
  use elm_varctl      , only : iulog
  use landunit_varcon , only : isturb_MIN, isturb_MAX
  use elm_varcon      , only : spval
  use LandunitType    , only : lun_pp
  use TopounitType    , only : top_pp                  
  use GridcellType    , only : grc_pp
  use mct_mod
  use shr_strdata_mod , only : shr_strdata_type
  !
  implicit none
  save
  private
  !
  !

  ! !PUBLIC TYPE
  type, public :: urbantv_type

     real(r8), public, pointer :: t_building_max(:)    ! lun maximum internal building air temperature (K)
     real(r8), public, pointer :: p_ac(:)              ! lun air-conditioning adoption rate (unitless, between 0 and 1)
     type(shr_strdata_type)    :: sdat_urbantv         ! urban time varying input data stream
   contains

     ! !PUBLIC MEMBER FUNCTIONS:
     procedure, public :: Init              ! Allocate and initialize urbantv
     procedure, public :: urbantv_init      ! Initialize urban time varying stream
     procedure, public :: urbantv_interp    ! Interpolate urban time varying stream
     
  end type urbantv_type

  integer      , private              :: stream_varname_MIN       ! minimum index for stream_varnames
  integer      , private              :: stream_varname_MAX       ! maximum index for stream_varnames
  character(15), private, pointer     :: stream_varnames(:)       ! urban time varying variable names

  character(len=*), parameter, private :: sourcefile = &
       __FILE__

  !-----------------------------------------------------------------------
  ! declare the public instance of urban parameters data types
  !-----------------------------------------------------------------------
  type(urbantv_type)          , public, target :: urbantv_vars    ! urban parameters
  !$acc declare create(urbantv_vars)

contains

  !-----------------------------------------------------------------------
  subroutine Init(this, bounds, NLFilename)
    !
    ! Allocate module variables and data structures
    !
    ! !USES:
    use shr_infnan_mod  , only : nan => shr_infnan_nan, assignment(=)
    use histFileMod     , only : hist_addfld1d
    use UrbanParamsType , only : urban_explicit_ac
    !
    ! !ARGUMENTS:
    class(urbantv_type) :: this
    type(bounds_type) , intent(in) :: bounds  
    character(len=*)  , intent(in) :: NLFilename   ! Namelist filename
    !
    ! !LOCAL VARIABLES:
    integer		:: begl, endl
    !---------------------------------------------------------------------

    begl = bounds%begl; endl = bounds%endl

    ! Determine the minimum and maximum indices for stream_varnames
    stream_varname_MIN = 1
    ! Get value for the maximum index for stream_varnames: if using explicit AC adoption scheme, 
    ! then set maximum index to 6 for reading in tbuildmax and p_ac for three urban density classes;
    ! otherwise, set to 3 to only read in tbuildmax for three urban density classes. 
    if (urban_explicit_ac) then
       stream_varname_MAX = 6
    else
       stream_varname_MAX = 3
    end if

    ! Allocate urbantv data structure

    allocate(this%t_building_max(begl:endl))          ; this%t_building_max      (:)   = nan
    allocate(this%p_ac(begl:endl))                    ; this%p_ac(:) = nan
    allocate(stream_varnames(stream_varname_MIN:stream_varname_MAX))

    call this%urbantv_init(bounds, NLFilename)
    call this%urbantv_interp(bounds)

    ! Add history fields
    call hist_addfld1d (fname='TBUILD_MAX', units='K',      &
          avgflag='A', long_name='prescribed maximum interior building temperature',   &
          ptr_lunit=this%t_building_max, default='inactive', set_nourb=spval, &
          l2g_scale_type='unity')
    if (urban_explicit_ac) then
       call hist_addfld1d (fname='P_AC', units='a fraction between 0 and 1',      &
             avgflag='A', long_name='prescribed air-conditioning ownership rate',   &
             ptr_lunit=this%p_ac, default='inactive', set_nourb=spval, &
             l2g_scale_type='unity')
    end if

  end subroutine Init

  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  subroutine urbantv_init(this, bounds, NLFilename)
   !
   ! !DESCRIPTION:
   ! Initialize data stream information for urban time varying data
   !
   ! !USES:
   use elm_varctl       , only : inst_name
   use elm_time_manager , only : get_calendar
   use ncdio_pio        , only : pio_subsystem
   use shr_pio_mod      , only : shr_pio_getiotype
   use elm_nlUtilsMod   , only : find_nlgroup_name
   use ndepStreamMod    , only : elm_domain_mct
   use spmdMod          , only : masterproc, mpicom, comp_id
   use fileutils        , only : getavu, relavu
   use shr_mpi_mod      , only : shr_mpi_bcast
   use shr_string_mod   , only : shr_string_listAppend
   use shr_strdata_mod  , only : shr_strdata_create, shr_strdata_print
   use decompMod        , only : gsmap_lnd_gdc2glo
   use domainMod        , only : ldomain
   use shr_infnan_mod   , only : nan => shr_infnan_nan, assignment(=)
   use landunit_varcon  , only : isturb_TBD, isturb_HD, isturb_MD
   use UrbanParamsType  , only : urban_explicit_ac
   !
   ! !ARGUMENTS:
   implicit none
   class(urbantv_type)           :: this
   type(bounds_type), intent(in) :: bounds
   character(len=*),  intent(in) :: NLFilename   ! Namelist filename
   !
   ! !LOCAL VARIABLES:
   integer            :: n
   integer            :: begl, endl                           ! landunits
   integer            :: ifield                               ! field index
   integer            :: stream_year_first_urbantv            ! first year in urban tv stream to use
   integer            :: stream_year_last_urbantv             ! last year in urban tv stream to use
   integer            :: model_year_align_urbantv             ! align stream_year_first_urbantv 
                                                              !  with this model year
   integer            :: nu_nml                               ! unit for namelist file
   integer            :: nml_error                            ! namelist i/o error flag
   type(mct_ggrid)    :: dom_elm                              ! domain information 
   character(len=CL)  :: stream_fldFileName_urbantv           ! urban tv streams filename
   character(len=CL)  :: urbantvmapalgo = 'nn'                ! mapping alogrithm for urban ac
   character(len=CL)  :: urbantv_tintalgo = 'linear'          ! time interpolation alogrithm
   character(len=CL)  :: fldList                              ! field string
   character(*), parameter :: subName = "('urbantv_init')"
   character(*), parameter :: F00 = "('(urbantv_init) ',4a)"
   !-----------------------------------------------------------------------
   namelist /urbantv_streams/       &
        stream_year_first_urbantv,  &
        stream_year_last_urbantv,   &
        model_year_align_urbantv,   &
        urbantvmapalgo,             &
        stream_fldFileName_urbantv, &
        urbantv_tintalgo
   !-----------------------------------------------------------------------

   begl = bounds%begl; endl = bounds%endl

   ! Default values for namelist
   stream_year_first_urbantv  = 1       ! first year in stream to use
   stream_year_last_urbantv   = 1       ! last  year in stream to use
   model_year_align_urbantv   = 1       ! align stream_year_first_urbantv with this model year
   stream_fldFileName_urbantv = ' '
   
   ! create the field list for these urbantv fields...use in shr_strdata_create
   stream_varnames(:)          = "NOT_SET"
   stream_varnames(1) = "tbuildmax_TBD"
   stream_varnames(2)  = "tbuildmax_HD"
   stream_varnames(3)  = "tbuildmax_MD"
   if (urban_explicit_ac) then
       stream_varnames(4) = "p_ac_TBD"
       stream_varnames(5) = "p_ac_HD"
       stream_varnames(6) = "p_ac_MD"
    end if
   fldList = ""
   do ifield = stream_varname_MIN,stream_varname_MAX
      call shr_string_listAppend( fldList, stream_varnames(ifield) )
   end do

   ! Read urbantv_streams namelist
   if (masterproc) then
      nu_nml = getavu()
      open( nu_nml, file=trim(NLFilename), status='old', iostat=nml_error )
      call find_nlgroup_name(nu_nml, 'urbantv_streams', status=nml_error)
      if (nml_error == 0) then
         read(nu_nml, nml=urbantv_streams,iostat=nml_error)
         if (nml_error /= 0) then
            call endrun(msg='ERROR reading urbantv_streams namelist'//errMsg(sourcefile, __LINE__))
         end if
      end if
      close(nu_nml)
      call relavu( nu_nml )
   endif

   call shr_mpi_bcast(stream_year_first_urbantv, mpicom)
   call shr_mpi_bcast(stream_year_last_urbantv, mpicom)
   call shr_mpi_bcast(model_year_align_urbantv, mpicom)
   call shr_mpi_bcast(stream_fldFileName_urbantv, mpicom)
   call shr_mpi_bcast(urbantv_tintalgo, mpicom)

   if (masterproc) then
      write(iulog,*) ' '
      write(iulog,*) 'urbantv_streams settings:'
      write(iulog,*) '  stream_year_first_urbantv  = ',stream_year_first_urbantv
      write(iulog,*) '  stream_year_last_urbantv   = ',stream_year_last_urbantv
      write(iulog,*) '  model_year_align_urbantv   = ',model_year_align_urbantv
      write(iulog,*) '  stream_fldFileName_urbantv = ',stream_fldFileName_urbantv
      write(iulog,*) '  urbantv_tintalgo           = ',urbantv_tintalgo
      do n = stream_varname_MIN,stream_varname_MAX
          write(iulog,'(a,a)' ) '  stream_varname         = ',trim(stream_varnames(n))
       end do
      write(iulog,*) ' '
   endif

   call elm_domain_mct (bounds, dom_elm)

   write(iulog,*) 'HERE 1'    ! REMOVE
   call shr_strdata_create(this%sdat_urbantv,name="elmurbantv",     &
        pio_subsystem=pio_subsystem,                   &
        pio_iotype=shr_pio_getiotype(inst_name),       &
        mpicom=mpicom, compid=comp_id,                 &
        gsmap=gsmap_lnd_gdc2glo, ggrid=dom_elm,        &
        nxg=ldomain%ni, nyg=ldomain%nj,                &
        yearFirst=stream_year_first_urbantv,           &
        yearLast=stream_year_last_urbantv,             &
        yearAlign=model_year_align_urbantv,            &
        offset=0,                                      &
        domFilePath='',                                &
        domFileName=trim(stream_fldFileName_urbantv),  &
        domTvarName='time',                            &
        domXvarName='lon' ,                            &
        domYvarName='lat' ,                            &
        domAreaName='area',                            &
        domMaskName='LANDMASK',                        &
        filePath='',                                   &
        filename=(/trim(stream_fldFileName_urbantv)/) , &
        fldListFile=fldList,                           &
        fldListModel=fldList,                          &
        fillalgo='none',                               &
        mapalgo=urbantvmapalgo,                        &
        calendar=get_calendar(),                       &
        tintalgo=urbantv_tintalgo,                     &
        taxmode='extend'                                 )

   write(iulog,*) 'HERE 2'     ! REMOVE

   if (masterproc) then
      call shr_strdata_print(this%sdat_urbantv,'urban time varying data')
   endif

   write(iulog,*) 'HERE 3'      ! REMOVE


  end subroutine urbantv_init

  !-----------------------------------------------------------------------
  subroutine urbantv_interp(this, bounds)
  !
  ! !DESCRIPTION:
  ! Interpolate data stream information for urban time varying data.
  !
  ! !USES:
  use elm_time_manager, only : get_curr_date
  use spmdMod         , only : mpicom
  use shr_strdata_mod , only : shr_strdata_advance
  use elm_varsur      , only : urban_valid
  use UrbanParamsType  , only : urban_explicit_ac
  !
  ! !ARGUMENTS:
  class(urbantv_type)           :: this
  type(bounds_type), intent(in) :: bounds
  !
  ! !LOCAL VARIABLES:
  logical :: found
!   integer :: l, glun,tlun, ig, g,t ip ! REMOVE
  integer :: l, glun, ig, g, t, ip, n, topo_ind
  integer :: year    ! year (0, ...) for nstep+1
  integer :: mon     ! month (1, ..., 12) for nstep+1
  integer :: day     ! day of month (1, ..., 31) for nstep+1
  integer :: sec     ! seconds into current date for nstep+1
  integer :: mcdate  ! Current model date (yyyymmdd)
  integer :: lindx   ! landunit index
  integer :: gindx   ! gridcell index
  integer :: tindx   ! topounit index
  !-----------------------------------------------------------------------

   call get_curr_date(year, mon, day, sec)
   mcdate = year*10000 + mon*100 + day

   call shr_strdata_advance(this%sdat_urbantv, mcdate, sec, mpicom, 'urbantvdyn')

   do l = bounds%begl,bounds%endl
      if (lun_pp%urbpoi(l)) then
         glun  = lun_pp%gridcell(l)

         ig = 0
         do g = bounds%begg,bounds%endg
            ig = ig+1
            if (g == glun) exit
         end do

         ! Since we are within an urban land unit, we know that
         ! lun%itype is within [pisturb_MIN, isturb_MAX]
         do n = stream_varname_MIN,stream_varname_MAX 
            if (trim(stream_varnames(lun_pp%itype(l)-isturb_MIN+1)) == stream_varnames(n)) then
               ! Determine vector index corresponding to glun
               ip = mct_aVect_indexRA(this%sdat_urbantv%avs(1),trim(stream_varnames(n)))
               this%t_building_max(l) = this%sdat_urbantv%avs(1)%rAttr(ip,ig)
            end if
            if (urban_explicit_ac) then
               if (trim(stream_varnames(lun_pp%itype(l)-isturb_MIN+4)) == stream_varnames(n)) then
                  ! Determine vector index corresponding to glun
                  ip = mct_aVect_indexRA(this%sdat_urbantv%avs(1),trim(stream_varnames(n)))
                  this%p_ac(l) = this%sdat_urbantv%avs(1)%rAttr(ip,ig)
               end if
            end if
         end do
      else
         this%t_building_max(l) = spval
         this%p_ac(l) = spval
      end if
   end do

   found = .false.
   do l = bounds%begl,bounds%endl
      if (lun_pp%urbpoi(l)) then
         ! g  = lun_pp%gridcell(l)  ! REMOVE MAYBE
         t  = lun_pp%topounit(l)
         topo_ind = top_pp%topo_grc_ind(t)
         !
         ! Determine vector index corresponding to glun
         !
         ! REMOVE 
         ig = 0
         glun  = lun_pp%gridcell(l)
         do g = bounds%begg,bounds%endg
            ig = ig+1
            if (g == glun) exit
         end do
         ! END REMOVE

         if ( .not. urban_valid(g,topo_ind) .or. (this%t_building_max(l) <= 0._r8)) then
            found = .true.
            gindx = g
            tindx = t
            lindx = l
            exit
         else if (urban_explicit_ac .and. (this%p_ac(l) < 0._r8 .or. this%p_ac(l) > 1._r8)) then
             found = .true.
             gindx = g
             tindx = t
             lindx = l
             exit
         end if
      end if
   end do
   if ( found ) then
      write(iulog,*)'ERROR: no valid urban data for g= ',gindx, ' and t= ',tindx
      write(iulog,*)'landunit type:   ',lun_pp%itype(l)
      write(iulog,*)'urban_valid:     ',urban_valid(gindx,topo_ind)
      write(iulog,*)'t_building_max:  ',this%t_building_max(lindx)
      if (urban_explicit_ac)  write(iulog,*)'p_ac:            ',this%p_ac(lindx)
      call endrun(msg=errmsg(sourcefile, __LINE__))
   end if


  end subroutine urbantv_interp

  !-----------------------------------------------------------------------

end module UrbanTimeVarType
