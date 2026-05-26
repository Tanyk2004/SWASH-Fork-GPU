subroutine SwashReadInput ( comput )
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering and Geosciences              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmers: The SWASH team                               |
!   --|-----------------------------------------------------------|--
!
!
!     SWASH (Simulating WAves till SHore); a non-hydrostatic wave-flow model
!     Copyright (C) 2010-2026  Delft University of Technology
!
!     This program is free software: you can redistribute it and/or modify
!     it under the terms of the GNU General Public License as published by
!     the Free Software Foundation, either version 3 of the License, or
!     (at your option) any later version.
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!     GNU General Public License for more details.
!
!     You should have received a copy of the GNU General Public License
!     along with this program. If not, see <http://www.gnu.org/licenses/>.
!
!
!   Authors
!
!    1.00: Marcel Zijlema
!
!   Updates
!
!    1.00, February 2010: New subroutine
!
!   Purpose
!
!   Reading and processing of the user commands describing the model
!
!   Method
!
!   A new line is read, in which the main keyword determines what the
!   command is. The command is read and processed. Common variables are
!   given proper values. After processing the command the routine returns
!   to label 100, to process the next command. This is repeated until the
!   command STOP is found or the end of input file is reached.
!
!   Depending on the commands, the argument variable 'comput' is given a
!   value, depending on which the program will make a computation in some
!   form or not (see meaning of variable 'comput' below).
!
!   Modules used
!
    use ocpcomm1
    use ocpcomm2
    use ocpcomm4
    use SwashCommdata1
    use SwashCommdata2
    use SwashCommdata3
    use SwashCommdata4
    use SwashTimecomm
    use SwashFlowdata
    use SwashSolvedata
    use SwashRigBoddata
    use outp_data
    use m_genarr
    use m_bndspec
    use m_parall
    use SwanGriddata
    use SwanGridobjects
!
    implicit none
!
!   Argument variables
!
    character(4), intent(inout)        :: comput             ! determines the sort of computation to be performed by Swash
                                                             ! ='COMP'; computation requested
                                                             ! ='NOCO'; no computation but output requested
                                                             ! ='STOP'; make computation, output and stop
!
!   Parameter variables
!
    integer, parameter                 :: mptst  = 50        ! maximum number of test points
    integer, parameter                 :: mxincl = 10        ! maximum of include levels
!
!   Local variables
!
    integer                            :: i                  ! loop counter
    integer                            :: iadv               ! to indicate advection term (horizontal/vertical; see command DISCRET)
    integer, save                      :: ient = 0           ! number of entries in this subroutine
    integer                            :: ieqn               ! to indicate type of equation (see command DISCRET)
    integer                            :: igrd               ! grid number for which a scalar or first component of a vector is read (see command INPGRID)
    integer                            :: igr2               ! grid number for which second component of a vector is read (see command INPGRID)
    integer, save                      :: inclev = 1         ! include level, increases at INCL command, decreases at end-of-file
    integer                            :: indl               ! indicator of layer thickness (see command VERT)
    integer                            :: iostat             ! I/O status in call FOR
    integer                            :: ipp                ! an index or pair of indices for test points
    integer                            :: iside              ! selected side of domain where sponge layer is given (see command SPON)
    integer                            :: itmp1              ! auxiliary integer
    integer                            :: itmp2              ! auxiliary integer
    integer                            :: itmp3              ! auxiliary integer
    integer                            :: ivtype             ! type of output quantity
    integer                            :: ix                 ! loop counter in x-direction
    integer                            :: iy                 ! loop counter in y-direction
    integer                            :: j                  ! loop counter
    integer                            :: kc1                ! index of boundary point of consideration (see command SPON)
    integer                            :: kc2                ! index of subsequent boundary point (see command SPON)
    integer                            :: labl               ! label of rigid body
    integer                            :: mxs                ! user-defined number of cells in x-direction
    integer                            :: mys                ! user-defined number of cells in y-direction
    integer                            :: nfen               ! number of fenders per rigid body (see command BODY FENDER)
    integer                            :: nlen               ! linelength of block output
    integer                            :: nmli               ! number of mooring lines per rigid body (see command BODY MLINE)
    integer                            :: nn                 ! number of stems (see command VEGET)
    integer                            :: numlay             ! actual number of vertical layers (see command VERT)
    integer                            :: nump               ! number of boundary points of selected side (see command SPON)
    integer                            :: nvar               ! actual number of output quantities (see command QUANT)
    integer                            :: vm                 ! boundary marker
    !
    integer, dimension(:), allocatable :: iarr               ! integer array to temporarily store indices of test points
    integer, dimension(mxincl), save   :: incnum = 0         ! unit reference numbers of included files
    !
    real                               :: bb                 ! damping coefficient of the mooring line
    real                               :: bk                 ! spring coefficient of the mooring line/fender
    real                               :: bx                 ! x-component of mooring line attachment point at bottom relative to center of gravity
    real                               :: by                 ! y-component of mooring line attachment point at bottom relative to center of gravity
    real                               :: bz                 ! z-component of mooring line attachment point at bottom relative to center of gravity
    real                               :: cc                 ! drag coefficient (see command VEGET)
    real                               :: cgx                ! x-component of center of gravity of rigid body (see command BODY)
    real                               :: cgy                ! y-component of center of gravity of rigid body (see command BODY)
    real                               :: cgz                ! z-component of center of gravity of rigid body (see command BODY)
    real                               :: cosdir             ! cosine of a direction (see command SPON)
    real                               :: crdm               ! largest crdp resulting in a side to be selected (see command SPON)
    real                               :: crdp               ! cosine of direction of the normal to side (see command SPON)
    real                               :: DEGCNV             ! indicates Cartesian or nautical degrees
    real                               :: dd                 ! stem diameter (see command VEGET)
    real*8                             :: diff               ! difference between end and begin times
    real                               :: dirref             ! user-defined direction of the normal to a side w.r.t. positive x-axis (pointing east) (see command SPON)
    real                               :: dirsi              ! user-defined direction of the normal to a side of domain (see command SPON)
    real                               :: dirsid             ! direction of the normal to a side to be selected (see command SPON)
    real                               :: el                 ! equilibrium length of the mooring line
    real                               :: fx                 ! x-component of mooring line/fender attachment point at floating body relative to center of gravity
    real                               :: fy                 ! y-component of mooring line/fender attachment point at floating body relative to center of gravity
    real                               :: fz                 ! z-component of mooring line/fender attachment point at floating body relative to center of gravity
    real                               :: height             ! wave height of regular wave (see command SOURCE)
    real                               :: hh                 ! layer thickness (see command VERT)
    real                               :: mass               ! mass of rigid body (see command BODY)
    real                               :: moix               ! x-component of moment of inertia of rigid body (see command BODY)
    real                               :: moiy               ! y-component of moment of inertia of rigid body (see command BODY)
    real                               :: moiz               ! z-component of moment of inertia of rigid body (see command BODY)
    real                               :: per                ! wave period of regular wave (see command SOURCE)
    real*8                             :: rtmp               ! auxiliary double precision
    real                               :: sindir             ! sine of a direction (see command SPON)
    real, dimension(5)                 :: spparm             ! parameters used for computation of internal-generated spectrum (see command SOURCE)
                                                             ! =1; wave height (significant or rms)
                                                             ! =2; wave period (peak or mean)
                                                             ! =3; incident or peak wave direction w.r.t. normal on the boundary
                                                             ! =4; directional distribution coefficient
                                                             ! =5; cyclic period of time series
    real                               :: spwid              ! width of sponge layer
    real                               :: sumx               ! sum of mesh sizes in x-direction (see command SPON)
    real                               :: sumy               ! sum of mesh sizes in y-direction (see command SPON)
    real                               :: tsmo               ! period for smoothing internal-generated values during cold start (see command SOURCE)
    real                               :: wdir               ! incident or peak wave direction with respect to problem coordinates (see command SOURCE)
    real                               :: xc                 ! x-coordinate of centroids
    real                               :: yc                 ! y-coordinate of centroids
    !
    real,dimension(12)                 :: qr                 ! user-defined parameters for output quantity (see command QUANT)
    !
    character(18)                      :: DTTIWR             ! write time coding
    character(8)                       :: psname             ! name of point set
    character(6)                       :: qovsnm             ! user-defined short name of output quantity (see command QUANT)
    character(40)                      :: qovlnm             ! user-defined long name of output quantity (see command QUANT)
    !
    logical                            :: chgdry             ! value of epsdry changed
    logical                            :: found              ! keyword found
    logical                            :: KEYWIS             ! indicates whether keyword in user manual is found or not
    logical, save                      :: liwg   = .false.   ! indicates whether linked list of parameters for internal-generated waves is initialized or not (see command SOURCE)
    logical, dimension(6), save        :: logcom = .false.   ! indicates which commands have been given to know if all the
                                                             ! information for a certain command is available. Meaning:
                                                             ! (1) no meaning
                                                             ! (2) command CGRID has been carried out
                                                             ! (3) command READINP BOTTOM has been carried out
                                                             ! (4) command READ COOR has been carried out
                                                             ! (5) command READ UNSTRUC has been carried out
                                                             ! (6) arrays s1, u1 and v1 have been allocated
    logical, save                      :: lrib = .false.     ! indicates whether linked list of rigid body parameters is initialized or not (see command BODY)
    logical, save                      :: runmade = .false.  ! indicates that computation has been made or not
    logical                            :: STPNOW             ! indicates that program must stop
    !
    type(bfsdat), pointer              :: iwgtmp             ! list containing parameters for internal wave generation (see command SOURCE)
    type(opsdat), pointer              :: opstmp             ! list containing parameters for grid
    type(ribdat), pointer              :: bodtmp             ! list containing parameters for rigid bodies (see command BODY)
    type(ribdat), save, pointer        :: crib               ! pointer to rigid body data
    !
    type oqpt                                                ! linked list for output quantity (see command QUANT)
       integer               :: i
       type(oqpt), pointer   :: nexti
    end type oqpt
    type(oqpt), target       :: frstq
    type(oqpt), pointer      :: currq, tmpq
    !
    type vertpt                                              ! linked list for vertical layer schematisation (see command VERT)
       integer               :: i
       real                  :: t
       type(vertpt), pointer :: nextv
    end type vertpt
    type(vertpt), target     :: frstv
    type(vertpt), pointer    :: currv, tmpv
    !
    type vegpt                                               ! linked list for vegetation schematisation (see command VEGET)
       integer               :: n
       real                  :: h, d, c
       type(vegpt), pointer  :: nextvg
    end type vegpt
    type(vegpt), target      :: frstvg
    type(vegpt), pointer     :: currvg, tmpvg
    !
    type sponpt                                              ! linked list for sponge layer (see command SPON)
       integer               :: k
       real                  :: w
       type(sponpt), pointer :: nextsp
    end type sponpt
    type(sponpt), target     :: frstsp
    type(sponpt), pointer    :: currsp, tmpsp
    !
    type fenpt                                               ! linked list for fenders of rigid body (see command BODY FENDER)
       real                  :: k, fx, fy, fz
       type(fenpt), pointer  :: nextfe
    end type fenpt
    type(fenpt), target      :: frstfe
    type(fenpt), pointer     :: currfe, tmpfe
    !
    type mlipt                                               ! linked list for mooring lines of rigid body (see command BODY MLINE)
       real                  :: el, k, b, bx, by, bz,  &
                                fx, fy, fz
       type(mlipt), pointer  :: nextml
    end type mlipt
    type(mlipt), target      :: frstml
    type(mlipt), pointer     :: currml, tmpml
!
!   Structure
!
!     call NWLINE for reading new line of user input
!     call INKEYW to read a new command from user input
!     If the command is equal to one of the SWAN commands, then
!         Read and process the rest of the command
!
!   Source text
!
    call strace (ient,'SwashReadInput')
    !
    ! read command
    !
 100 call NWLINE
    !
    if ( ELTYPE == 'EOF' ) then
       !
       ! end-of-file encountered in (included) input file
       ! return to previous input file
       !
       close (incnum(inclev))
       inclev = inclev - 1
       if ( inclev == 0 ) then
          call msgerr (4, 'unexpected end of command input')
          return
       endif
       !
       ELTYPE = 'USED'
       INPUTF = incnum(inclev)
       !
    endif
    !
    call INKEYW ('REQ',' ')
    !
    ! processing of commands
    !
    ! ==========================================================================
    !
    ! STOP
    !
    ! ==========================================================================
    !
    if ( KEYWIS ('STOP') ) then
       !
       if ( runmade ) then
          comput = 'STOP'
       else
          write (PRINTF,*)' ** No computation requested **'
          comput = 'NOCO'
       endif
       return
       !
    endif
    !
    ! ==========================================================================
    !
    ! INCLude  'FILE'
    !
    ! ==========================================================================
    !
    if ( KEYWIS ('INCL') ) then
       !
       if ( inclev == 1 ) incnum(inclev) = INPUTF
       call INCSTR ('FILE' , FILENM , 'REQ', ' ')
       inclev = inclev + 1
       if ( inclev > mxincl ) then
          call msgerr (4, 'too many INCLUDE levels')
          return
       endif
       iostat = 0
       call FOR (incnum(inclev), FILENM, 'OF', iostat)
       INPUTF = incnum(inclev)
       !
       goto 100
       !
    endif
    !
    ! ==========================================================================
    !
    ! PROJect  'NAME'  'NR'
    !
    !          'title1'
    !
    !          'title2'
    !
    !          'title3'
    !
    ! ==========================================================================
    !
    if ( KEYWIS ('PROJ') ) then
       !
       call INCSTR ('NAME'  , PROJID, 'UNC', BLANK)
       call INCSTR ('NR'    , PROJNR, 'REQ', BLANK)
       call NWLINE
       if ( STPNOW() ) return
       call INCSTR ('TITLE1', PROJT1, 'UNC', ' ')
       call NWLINE
       if ( STPNOW() ) return
       call INCSTR ('TITLE2', PROJT2, 'UNC', ' ')
       call NWLINE
       if ( STPNOW() ) return
       call INCSTR ('TITLE3', PROJT3, 'UNC', ' ')
       !
       goto 100
       !
    endif
    !
    ! ==========================================================================
    !
    ! OUTPut OPTions  'comment'  (TABle [field])  (BLOck  [ndec]  [len])
    !
    ! ==========================================================================
    !
    if ( KEYWIS('OUTP') ) then
       !
       call IGNORE ('OPT')
       call INCSTR ('COMMENT', out_comment, 'UNC', ' ')
       !
       call INKEYW ('STA', ' ')
       if ( KEYWIS('TAB') ) then
          call ININTG ('FIELD', fld_table, 'STA', 15)
          if ( fld_table > 16 ) call msgerr (2, '[field] is too large')
          if ( fld_table < 8  ) call msgerr (2, '[field] is too small')
          write (flt_table, 201) fld_table, fld_table-7
          if ( ITEST >= 30 ) write (PRINTF, 202) flt_table
          call INKEYW ('STA', ' ')
       endif
       !
       if ( KEYWIS('BLO') ) then
          call ININTG ('NDEC', dec_block, 'STA', 6)
          if ( dec_block > 9 ) call msgerr (2, '[ndec] is too large')
          call ININTG ('LEN', nlen, 'STA', 200)
          if ( nlen > 9999 ) call msgerr (2, '[len] is too large')
          write (flt_block, 203) nlen, dec_block+8, dec_block
          if ( ITEST >= 30 ) write (PRINTF, 204) flt_block
       endif
       !
       goto 100
       !
    endif
    !
    ! ==========================================================================
    !
    !       |    STAtionary    |    | -> TWODimensional |
    ! MODE <                    >  <                     >   &
    !       | -> NONSTationary |    |    ONEDimensional |
    !
    !       (SKIPMOMentum) (LINear) (not documented)
    !
    ! ==========================================================================
    !
    if ( KEYWIS ('MODE') ) then
       !
       call INKEYW ('STA',' ')
       if ( KEYWIS('NONST') .or. KEYWIS ('DYN') ) then
          !
          if ( nstatm == 0 ) call msgerr (2, 'mode NONST incorrect here')
          nstatm = 1
          nstatc = 1
          !
       else if ( KEYWIS ('STA') ) then
          !
          if ( nstatm == 1 ) call msgerr (2, 'mode STA incorrect here')
          nstatm = 0
          call INKEYW ('STA',' ')
          !
       endif
       !
       call INKEYW ('STA', ' ')
       if ( KEYWIS ('ONED') ) then
          oned = .true.
       else if ( KEYWIS ('TWOD') ) then
          oned = .false.
       endif
       !
       ! with this optional command the user indicates that the momentum and continuity equations are skipped
       ! the current and water level are then given by the input fields
       ! this command is, for example, useful if only the transport equation must be tested
       !
       call INKEYW ('STA', ' ')
       if (KEYWIS ('SKIPMOM')) then
          momskip = .true.
       endif
       !
       ! with this optional command the linearized shallow water equations are to be dealt with
       !
       call INKEYW ('STA', ' ')
       if (KEYWIS ('LIN')) then
          linswe = .true.
       endif
       !
       goto 100
       !
    endif
    !
    ! ==========================================================================
    !
    !            |  STATionary  [time]                      |
    ! COMPute ( <                                            > )
    !            |                    | -> Sec  |           |
    !            |  ([tbegc] [deltc] <     MIn   > [tendc]) |
    !                                 |    HR   |
    !                                 |    DAy  |
    !
    ! ==========================================================================
    !
    IF ( KEYWIS ('COMP') ) THEN
       !
       comput  = 'COMP'
       runmade = .true.
       !
       call INKEYW ('STA','  ')
       if ( nstatm == 0 .or. KEYWIS('STAT') ) then
          if ( nstatm > 0 ) call INCTIM (ITMOPT,'TIME',tinic,'REQ',0d0)
          if ( tinic < timco ) then
             call msgerr (2, '[time] before current time')
             tinic = timco
          endif
          tfinc  = tinic
          timco  = tinic
          dt     = 1.e10
          rdtim  = 0.
          nstatc = 0
          mtc    = 1
       else
          call IGNORE ('NONST')
          if ( timco < -0.9e10 ) then
             call INCTIM (ITMOPT,'TBEGC',tinic,'REQ',0d0)
          else
             call INCTIM (ITMOPT,'TBEGC',tinic,'STA',timco)
          endif
          if ( tinic < timco-1.01*dt ) then
             call msgerr (2, 'start time [tbegc] before current time')
             tinic = timco
          endif
          call INITVD ('DELTC', dt, 'REQ', 0d0)
          call INCTIM (ITMOPT,'TENDC',tfinc,'REQ',0d0)
          nstatc = 1
          !
          ! tfinc must be greater than tinic
          !
          diff = tfinc - tinic
          if ( .not. diff > 0. ) call msgerr (3, 'start time [tbegc] greater or equal end time [tendc]')
          !
          ! the number of time steps is calculated
          !
          rdtim = 1./dt
          mtc = nint ((tfinc - tinic)/dt)
          if ( mod(tfinc-tinic,dt) > 0.01*dt .and. mod(tfinc-tinic,dt) < 0.99*dt ) &
             call msgerr (1, '[deltc] is not a fraction of the simulation period')
          timco = tinic
          !
       endif
       if ( nstatm > 0 ) chtime = DTTIWR(ITMOPT, timco)
       ncompt = ncompt + 1
       if ( ncompt > 300 ) call msgerr (2, 'no more than 300 COMPUTE commands are allowed')
       rcompt(ncompt,1) = dble(nstatc)
       rcompt(ncompt,2) = dble(mtc)
       rcompt(ncompt,3) = tfinc
       rcompt(ncompt,4) = tinic
       rcompt(ncompt,5) = dt
       !
       return
       !
    endif
    !
    ! ==========================================================================
    !
    ! SET  [level] [nor] [depmin] [maxmes] [maxerr] [seed]       &
    !      [grav] [rhowat] [temp] [salinity] [dynvis] [rhoair]   &
    !      [rhosed] [cdcap] [prmean] [backvisc] [kappa]          &
    !      [latitude]                                            &
    !      CARTesian/NAUTical                                    &
    !      [epshu] (not documented)                              &
    !      CC (not documented)                                   &
    !      [printf]  [prtest] (not documented)                   &
    !      [outlev]
    !
    ! ==========================================================================
    !
    if ( KEYWIS ('SET') ) then
       !
       call INREAL ('LEVEL',    swl,    'UNC', 0.)
       call INREAL ('NOR',      dnorth, 'UNC', 0.)
       call INREAL ('DEPMIN',   epsdry, 'UNC', 0.)
       chgdry = CHGVAL
       call ININTG ('MAXMES',   maxmes, 'UNC', 0 )
       call ININTG ('MAXERR',   MAXERR, 'UNC', 0 )
       call ININTG ('SEED',     myseed, 'UNC', 0 )
       call INREAL ('GRAV',     grav,   'UNC', 0.)
       call INREAL ('RHOWAT',   rhow,   'UNC', 0.)
       call INREAL ('TEMP',     tempw,  'UNC', 0.)
       call INREAL ('SALINITY', salw,   'UNC', 0.)
       call INREAL ('DYNVIS',   dynvis, 'UNC', 0.)
       call INREAL ('RHOAIR',   rhoa,   'UNC', 0.)
       call INREAL ('RHOSED',   rhos,   'UNC', 0.)
       call INREAL ('CDCAP',    cdcap,  'UNC', 0.)
       call INREAL ('PRMEAN',   prmean, 'UNC', 0.)
       call INREAL ('BACKVISC', bvisc,  'UNC', 0.)
       call INREAL ('KAPPA',    vonkar, 'UNC', 0.)
       call INREAL ('LATITUDE', ylatc,  'UNC', 0.)
       !
       call INKEYW ('STA',' ')
       if ( KEYWIS ('NAUT') ) then
          bnaut     = .true.
          outpar(2) = 0.
       endif
       if ( KEYWIS ('CART') ) bnaut = .false.
       !
       call INREAL ('EPSHU', epshu, 'UNC', 0.)
       if ( .not.CHGVAL .and. chgdry ) epshu = epsdry
       !
       call INKEYW ('STA',' ')
       if ( KEYWIS ('CC') ) ccn = .false.
       !
!      the unit numbers PRINTF and PRTEST should be set in swashinit,
!      so not documented
       call ININTG ('PRINTF', PRINTF, 'UNC', 0)
       call ININTG ('PRTEST', PRTEST, 'UNC', 0)
       !
       call ININTG ('OUTLEV', iamout, 'UNC', 0)
       !
       goto 100
       !
    endif
    !
    ! ==========================================================================
    !
    ! COORDinates  /  -> CARTesian
    !              \ SPHErical [rearth]  CCM|QC
    !
    ! ==========================================================================
    !
    if ( KEYWIS ('COORD') ) then
       !
       call INKEYW ('STA',' ')
       if ( KEYWIS ('CART') ) then
          kspher = 0
       else if ( KEYWIS ('SPHE') ) then
          kspher = 1
          call INREAL ('REARTH', rearth, 'UNC', 0.)
          lendeg = rearth * pi / 180.
          !
          ! change properties of output quantities Xp and Yp
          !
          ovunit(1) = 'degr'
          ovllim(1) = -200.
          ovulim(1) =  400.
          ovlexp(1) = -180.
          ovhexp(1) =  360.
          ovexcv(1) = -999.
          !
          ovunit(2) = 'degr'
          ovllim(2) = -100.
          ovulim(2) =  100.
          ovlexp(2) = -90.
          ovhexp(2) =  90.
          ovexcv(2) = -999.
          !
          call INKEYW ('STA','CCM')
          if ( KEYWIS ('QC') ) then
             !
             ! quasi-Cartesian projection method
             !
             mproj = 0
             !
          else if ( KEYWIS ('CCM') ) then
             !
             ! uniform Mercator projection (default for spherical coordinates)
             !
             mproj = 1
             !
          endif
       else
          call WRNKEY
       endif
       !
       goto 100
       !
    endif
    !
    ! ==========================================================================
    !
    ! QUANTity  <...>  'short'  'long'  [lexp]  [hexp]  [excv]   &
    !
    !       [ref]                 (for output quantity TSEC)
    !       [dur] SEC|MIN|HR|DAY  (for output quantity HSIG, SETUP, MVEL, MTKE)
    !       [depth]               (for output quantity HRUN)
    !       [delrp]               (for output quantity RUNUP)
    !       [xcom] [ycom] [zcom]  (for moments acting on floating body)
    !       [alpobj]              (for forces and moments acting on rotated body)
    !       PROBLEM/FRAME         (for directions and vectors)
    !
    ! ==========================================================================
    !
    if ( KEYWIS('QUANT') ) then
       !
       nvar    = 0
       frstq%i = 0
       nullify(frstq%nexti)
       currq => frstq
 110   call svartp (ivtype)
       if ( ivtype >= 1 .and. ivtype <= nmovar ) then
          !
          nvar = nvar + 1
          allocate(tmpq)
          tmpq%i = ivtype
          nullify(tmpq%nexti)
          currq%nexti => tmpq
          currq => tmpq
          call INCSTR ('SHORT' , qovsnm, 'STA', ' ')
          call INCSTR ('LONG'  , qovlnm, 'STA', ' ')
          call INREAL ('LEXP'  , qr(1) , 'STA',  999. )
          call INREAL ('HEXP'  , qr(2) , 'STA', -999. )
          call INREAL ('EXCV'  , qr(3) , 'STA',  999. )
          call INCTIM (ITMOPT  , 'REF' , rtmp, 'STA', -9.99d2)
          qr(4) = real(rtmp)
          call ININTV ('DUR'   , qr(5) , 'STA', -999. )
          call INREAL ('DEPTH' , qr(6) , 'STA', -999. )
          call INREAL ('DELRP' , qr(12), 'STA', -999. )
          call INREAL ('XCOM'  , qr(8) , 'STA',    0. )
          call INREAL ('YCOM'  , qr(9) , 'STA',    0. )
          call INREAL ('ZCOM'  , qr(10), 'STA',    0. )
          call INREAL ('ALPOBJ', qr(11), 'STA',    0. )
          call INKEYW ('STA', ' ')
          if ( KEYWIS('PROBLEM') .or. KEYWIS('USER') ) then
             qr(7) = 0.
          else if ( KEYWIS('FRAME') ) then
             qr(7) = 1.
          else
             qr(7) = -999.
          endif
          goto 110
          !
       else if ( ivtype /= 999 ) then
          !
          call msgerr (2, 'unknown output quantity ')
          !
       endif
       !
       if ( nvar > 0 ) then
          !
          currq => frstq%nexti
          do i = 1, nvar
             ivtype = currq%i
             if ( qovsnm /= ' '  ) ovsnam(ivtype) = qovsnm
             if ( qovlnm /= ' '  ) ovlnam(ivtype) = qovlnm
             if ( qr(1) /=  999. ) ovlexp(ivtype) = qr(1)
             if ( qr(2) /= -999. ) ovhexp(ivtype) = qr(2)
             if ( qr(3) /=  999. ) ovexcv(ivtype) = qr(3)
             if ( ivtype == 40 .or. ivtype == 41 ) then
                if ( nstatm == 0 ) call msgerr (2, 'time output asked in stationary mode')
             endif
             if ( ivtype == 41 ) then
                if ( qr(4) /= -999. ) outpar(1) = qr(4)
             endif
             if ( ivtype == 22 .or. ivtype == 23 .or. ivtype == 24 ) then
                if ( nstatm == 0 ) call msgerr (2, 'wave output asked in stationary mode')
                if ( qr(5) /= -999. ) twavoutp = qr(5)
             endif
             if ( ivtype >= 33 .and. ivtype <= 37 ) then
                if ( nstatm == 0 ) call msgerr (2, 'time-averaged current asked in stationary mode')
                if ( qr(5) /= -999. ) tcuroutp = qr(5)
             endif
             if ( ivtype == 42 .or. ivtype == 43 .or. ivtype == 44 ) then
                if ( nstatm == 0 ) call msgerr (2, 'time-averaged constituent asked in stationary mode')
                if ( qr(5) /= -999. ) ttraoutp = qr(5)
             endif
             if ( ivtype == 57 .or. ivtype == 58 .or. ivtype == 59 ) then
                if ( nstatm == 0 ) call msgerr (2, 'time-averaged turbulence quantities per layer asked in stationary mode')
                if ( qr(5) /= -999. ) tturoutp = qr(5)
             endif
             if ( ivtype >= 84 .and. ivtype <= 88 ) then
                if ( nstatm == 0 ) call msgerr (2, 'time-averaged current per layer asked in stationary mode')
                if ( qr(5) /= -999. ) tcuroutp = qr(5)
             endif
             if ( ivtype == 92 .or. ivtype == 93 .or. ivtype == 94 ) then
                if ( nstatm == 0 ) call msgerr (2, 'time-averaged constituent per layer asked in stationary mode')
                if ( qr(5) /= -999. ) ttraoutp = qr(5)
             endif
             if ( ivtype == 21 ) then
                if ( qr(6) /= -999. ) hrunp = qr(6)
             endif
             if ( ivtype >= 104 .and. ivtype <= 106 ) then
                cogx = qr(8)
                cogy = qr(9)
                cogz = qr(10)
             endif
             if ( ivtype >= 101 .and. ivtype <= 106 ) then
                alpobj = qr(11)
                alpobj = pi2 * (alpobj/360. - nint(alpobj/360.))
             endif
             if ( ivtype == 115 ) then
                if ( qr(12) /= -999. ) delrp = qr(12)
             endif
             if ( ovsvty(ivtype) == 2 .or. ovsvty(ivtype) == 3 ) then
                !
                ! direction or vector
                !
                if ( qr(7) /= -999. ) then
                   outpar(2) = qr(7)
                   if ( bnaut .and. ovsvty(ivtype) == 2 ) call msgerr (1, 'option not allowed with Nautical convention')
                endif
             endif
             currq => currq%nexti
          enddo
          !
          deallocate(tmpq)
          !
       endif
       !
       goto 100
       !
    endif
    !
    ! ==========================================================================
    !
    !            |    LINear [k]                 (not documented)
    !            |
    !            |    CONstant [cf]
    !            |
    !            |    CHEZy [cf]
    !            |
    ! FRICtion  <  -> MANNing [cf]
    !            |
    !            |    COLEbrook [h]
    !            |
    !            |            | -> SMOOTH
    !            |    LOGlaw <
    !                         |    ROUGHness [h]
    !
    ! ==========================================================================
    !
    if ( KEYWIS ('FRIC') ) then
       irough = 3
       call INKEYW ('STA','MANN')
       if ( KEYWIS('CON') ) then
          irough = 1
          call INREAL('CF',pbot(1),'STA',0.002)
       else if ( KEYWIS('CHEZ') ) then
          irough = 2
          call INREAL('CF',pbot(1),'STA',65.)
       else if ( KEYWIS('MANN') ) then
          irough = 3
          call INREAL('CF',pbot(1),'STA',0.019)
       else if ( KEYWIS('LOG') ) then
          irough = 4
          call INKEYW ('STA', 'SMOOTH')
          if ( KEYWIS ('ROUGH') ) then
             call INREAL ('H', pbot(2), 'UNC', 0.)
          else
             call IGNORE ('SMOOTH')
             pbot(2) = 0.
          endif
       else if ( KEYWIS('COLE') ) then
          irough = 5
          call INREAL('H',pbot(1),'UNC', 0.)
       else if ( KEYWIS('LIN') ) then
          irough = 11
          call INREAL('K',pbot(1),'UNC', 0.)
       else
          call WRNKEY
       endif
       goto 100
    endif
    !
    ! ==========================================================================
    !
    !                   | -> CONstant [cd]                       |
    !                   |                                        |
    !                   |    CHARNock [beta] [height]            |
    !                   |                                        |
    !                   |    LINear [a1] [a2] [b] [wlow] [whigh] |
    !                   |                                        |
    !                   |    WU                                  |   | REL  [alpha]
    ! WIND [vel] [dir] <                                          > <
    !                   |    GARRatt                             |   | RELW [crest]
    !                   |                                        |
    !                   |    SMIthbanke                          |
    !                   |                                        |
    !                   |    CHEn   (not documented)             |
    !                   |                                        |
    !                   |    FIT                                 |
    !
    ! ==========================================================================
    !
    if ( KEYWIS ('WIND') ) then
       call INREAL('VEL',u10 ,'UNC',0.)
       call INREAL('DIR',wdip,'UNC',0.)
       iwind = 1
       call INKEYW ('STA','CON')
       if ( KEYWIS('CON') ) then
          iwind = 1
          call INREAL('CD',pwnd(1),'STA',0.002)
       else if ( KEYWIS('CHARN') ) then
          iwind = 2
          call INREAL('BETA'  ,pwnd(2),'STA',0.032)
          call INREAL('HEIGHT',pwnd(3),'STA',10.)
       else if ( KEYWIS('LIN') ) then
          iwind = 3
          call INREAL('A1'   ,pwnd(5),'UNC',0.)
          call INREAL('A2'   ,pwnd(6),'UNC',0.)
          call INREAL('B'    ,pwnd(7),'UNC',0.)
          call INREAL('WLOW' ,pwnd(8),'UNC',0.)
          call INREAL('WHIGH',pwnd(9),'UNC',0.)
       else if ( KEYWIS('WU') ) then
          iwind = 3
          pwnd(5) = 0.8
          pwnd(6) = 0.
          pwnd(7) = 0.065
          pwnd(8) = 7.5
          pwnd(9) = 1000.
       else if ( KEYWIS('GARR') ) then
          iwind = 3
          pwnd(5) = 0.75
          pwnd(6) = 0.
          pwnd(7) = 0.067
          pwnd(8) = -1.
          pwnd(9) = 1000.
       else if ( KEYWIS('SMI') ) then
          iwind = 3
          pwnd(5) = 0.63
          pwnd(6) = 0.
          pwnd(7) = 0.066
          pwnd(8) = 3.
          pwnd(9) = 21.
       else if ( KEYWIS('CHE') ) then
          iwind = 3
          pwnd(5) = 0.2
          pwnd(6) = 18.
          pwnd(7) = 0.065
          pwnd(8) = 7.5
          pwnd(9) = 1000.
       else if ( KEYWIS('FIT') ) then
          iwind = 4
       else
          call WRNKEY
       endif
       !
       ! convert wdip from Nautical degrees to Cartesian degrees, if appropriate
       !
       wdip = DEGCNV (wdip)
       !
       ! wdip is made to be between -pi and pi
       !
       wdip = pi2 * (wdip/360. - nint(wdip/360.))
       !
       call INKEYW ('STA',' ')
       if ( KEYWIS ('RELW') ) then
          relwav = .true.
          call INREAL('CREST',pwnd(10),'STA',0.4)
       else if ( KEYWIS ('REL') ) then
          relwnd = .true.
          call INREAL('ALPHA',pwnd(4),'STA',1.)
       endif
       !
       goto 100
    endif
    !
    ! ==========================================================================
    !
    ! BREaking [alpha] [beta] [nufac]
    !
    ! ==========================================================================
    !
    if ( KEYWIS ('BRE') ) then
       !
       !call ININTG ( 'ISURF', isurf, 'STA', 1 )
       isurf = 1
       !
       if ( isurf == 1 ) then
          !
          call INREAL ( 'ALPHA', psurf(1), 'STA',0.6)
          call INREAL ( 'BETA' , psurf(2), 'STA',-1.)
          call INREAL ( 'NUFAC', psurf(3), 'STA',1.0)
          !
          horwinc = .true.
          if ( ihvisc == 0 ) ihvisc = 4
          !
       else if ( isurf == 2 ) then
          !
          call INREAL ( 'ALPHA', psurf(1), 'STA', 0.4 )
          !
          psurf(2) = 1.e+8
          !
       endif
       !
       goto 100
       !
    endif
    !
    ! ==========================================================================
    !
    !                         |    Cell   |
    ! AMBient  [U] [V] [eta] <             >  [theta]
    !                         | -> Stagg  |
    !
    ! ==========================================================================
    !
    if ( KEYWIS ('AMB') ) then
       iamb = 2
       call INREAL('U'  ,uamb,'UNC',0.)
       call INREAL('V'  ,vamb,'UNC',0.)
       call INREAL('eta',eamb,'UNC',0.)
       call INKEYW ('STA','S')
       if ( KEYWIS ('C') ) then
          iamb = 1
       else
          call IGNORE ('S')
          iamb = 2
       endif
       call INREAL( 'THETA', pnums(71), 'STA', 1. )
       goto 100
    endif
    !
    ! ==========================================================================
    !
    ! POROsity [size] [height] [alpha0] [beta0] [wper]
    !
    ! ==========================================================================
    !
    if ( KEYWIS ('PORO') ) then
       iporos = 1
       call INREAL( 'SIZE'  , ppor(1), 'UNC', 0.    )
       call INREAL( 'HEIGHT', ppor(2), 'STA', 99999.)
       call INREAL( 'ALPHA0', ppor(3), 'STA', 200.  )
       call INREAL( 'BETA0' , ppor(4), 'STA', 1.1   )
       call INREAL( 'WPER'  , ppor(5), 'STA', -1.   )
       goto 100
    endif
    !
    ! ==========================================================================
    !
    ! floating rigid body parameters
    !
    ! ==========================================================================
    !
    if ( KEYWIS ('SHIP') .or. KEYWIS ('VESS') .or. KEYWIS ('WEC') .or. KEYWIS ('FLOAT') .or. KEYWIS ('BODY') .or. KEYWIS ('PONT') .or. KEYWIS ('PFLOW') ) then
       !
       call INKEYW ('REQ',' ')
       if ( KEYWIS('METH') .or. KEYWIS('SOLV') ) then
          !
          ! ==========================================================================
          !
          !              |    NEWmark [beta] [gamma] |
          !              |                           |
          !              | -> CH [rho]               |
          ! BODY SOLVer <                             > COUPler [tol] [maxiter] [relax] &
          !              |    HHT [rho]              |
          !              |                           |
          !              |    WBZ [rho]              |
          !
          !      KBC [theta]
          !
          ! ==========================================================================
          !
          call INKEYW ('STA', 'CH')
          if ( KEYWIS('NEW') ) then
             mfoti = 1
             call INREAL ('BETA' , pship(3), 'STA', 0.25)
             call INREAL ('GAMMA', pship(4), 'STA', 0.50)
          else if ( KEYWIS('HHT') ) then
             mfoti = 3
             call INREAL ('RHO', pship(5), 'STA', 1.0)
          else if ( KEYWIS('WBZ') ) then
             mfoti = 4
             call INREAL ('RHO', pship(5), 'STA', 1.0)
          else
             call IGNORE ('CH')
             mfoti = 2
             call INREAL ('RHO', pship(5), 'STA', 1.0)
          endif
          !
          call INKEYW ('STA','  ')
          if ( KEYWIS('COUP') ) then
             call INREAL('TOL'    , pship(8) , 'STA', 0.0001)
             call ININTG('MAXITER', itmp1    , 'STA', 50)
             call INREAL('RELAX'  , pship(10), 'STA', 0.3)
             pship(9) = real(itmp1)
          endif
          !
          call INKEYW ('STA','  ')
          if ( KEYWIS('KBC') ) then
             call INREAL ('THETA', pship(11), 'STA', 0.5)
          endif
          !
       else if ( KEYWIS('DIM') ) then
          !
          ! ==========================================================================
          !
          ! BODY DIMension [l] [mass] [Ix] [Iy] [Iz] [cogx] [cogy] [cogz]             &
          !
          !      DOF SUrge SWay HEave ROll PItch YAw                                  &
          !
          !      ( MLIne < [K] [B] [apbx] [apby] [apbz] [apfx] [apfy] [apfz] [elen] > &
          !
          !              PRETension )
          !
          !      ( FENder < [K] [apfx] [apfy] [apfz] > )
          !
          ! ==========================================================================
          !
          allocate(bodtmp)
          bodtmp%label     = 0
          bodtmp%parm(1:7) = 0.
          bodtmp%dof(1:6)  = .false.
          !
          nmli = 0
          nfen = 0
          !
          call ININTG ('L'   , labl, 'REQ', 0 )
          call INREAL ('MASS', mass, 'REQ', 0.)
          call INREAL ('IX'  , moix, 'REQ', 0.)
          call INREAL ('IY'  , moiy, 'REQ', 0.)
          call INREAL ('IZ'  , moiz, 'REQ', 0.)
          call INREAL ('COGX', cgx , 'REQ', 0.)
          call INREAL ('COGY', cgy , 'REQ', 0.)
          call INREAL ('COGZ', cgz , 'REQ', 0.)
          bodtmp%label   = labl
          bodtmp%parm(1) = mass
          bodtmp%parm(2) = moix
          bodtmp%parm(3) = moiy
          bodtmp%parm(4) = moiz
          bodtmp%parm(5) = cgx
          bodtmp%parm(6) = cgy
          bodtmp%parm(7) = cgz
          !
          call INKEYW ('REQ', ' ')
          if ( KEYWIS('DOF') ) then
             !
             call INKEYW ('STA', ' ')
             if ( KEYWIS ('SU') ) bodtmp%dof(1) = .true.
             call INKEYW ('STA', ' ')
             if ( KEYWIS ('SW') ) bodtmp%dof(2) = .true.
             call INKEYW ('STA', ' ')
             if ( KEYWIS ('HE') ) bodtmp%dof(3) = .true.
             call INKEYW ('STA', ' ')
             if ( KEYWIS ('RO') ) bodtmp%dof(4) = .true.
             call INKEYW ('STA', ' ')
             if ( KEYWIS ('PI') ) bodtmp%dof(5) = .true.
             call INKEYW ('STA', ' ')
             if ( KEYWIS ('YA') ) bodtmp%dof(6) = .true.
             !
             call INKEYW ('STA', ' ')
             if ( KEYWIS ('ALL') ) bodtmp%dof(1:6) = .true.
             !
          endif
          !
          call INKEYW ('REQ', '  ')
          if ( KEYWIS('MLI') ) then
             !
             frstml%k  = 0.
             frstml%b  = 0.
             frstml%bx = 0.
             frstml%by = 0.
             frstml%bz = 0.
             frstml%fx = 0.
             frstml%fy = 0.
             frstml%fz = 0.
             frstml%el = 0.
             nullify(frstml%nextml)
             currml => frstml
 112         call INREAL ('K', bk, 'REP', -1.)
             if ( bk /= -1. ) then
                call INREAL ('B'   , bb, 'REQ', 0.)
                call INREAL ('APBX', bx, 'REQ', 0.)
                call INREAL ('APBY', by, 'REQ', 0.)
                call INREAL ('APBZ', bz, 'REQ', 0.)
                call INREAL ('APFX', fx, 'REQ', 0.)
                call INREAL ('APFY', fy, 'REQ', 0.)
                call INREAL ('APFZ', fz, 'REQ', 0.)
                call INREAL ('ELEN', el, 'REQ', 0.)
                nmli = nmli + 1
                allocate(tmpml)
                tmpml%k  = bk
                tmpml%b  = bb
                tmpml%bx = bx
                tmpml%by = by
                tmpml%bz = bz
                tmpml%fx = fx
                tmpml%fy = fy
                tmpml%fz = fz
                if ( el /= 0. ) then
                   tmpml%el = el
                else
                   tmpml%el = sqrt((fx-bx)**2+(fy-by)**2+(fz-bz)**2)
                endif
                nullify(tmpml%nextml)
                currml%nextml => tmpml
                currml => tmpml
                goto 112
             endif
             if ( nmli == 0 ) call msgerr (2,'no mooring line parameters found')
             !
             bodtmp%nmli = nmli
             allocate(bodtmp%mli(nmli,9))
             currml => frstml%nextml
             do j = 1, nmli
                bodtmp%mli(j,1) = currml%k
                bodtmp%mli(j,2) = currml%b
                bodtmp%mli(j,3) = currml%bx
                bodtmp%mli(j,4) = currml%by
                bodtmp%mli(j,5) = currml%bz
                bodtmp%mli(j,6) = currml%fx
                bodtmp%mli(j,7) = currml%fy
                bodtmp%mli(j,8) = currml%fz
                bodtmp%mli(j,9) = currml%el
                currml => currml%nextml
             enddo
             deallocate(tmpml)
             !
             call INKEYW ('STA', ' ')
             if ( KEYWIS ('PRET') ) then
                bodtmp%pret = .true.
             else
                bodtmp%pret = .false.
             endif
             !
          endif
          mxmli = max( mxmli, nmli )
          !
          call INKEYW ('REQ', '  ')
          if ( KEYWIS('FEN') ) then
             !
             frstfe%k  = 0.
             frstfe%fx = 0.
             frstfe%fy = 0.
             frstfe%fz = 0.
             nullify(frstfe%nextfe)
             currfe => frstfe
 113         call INREAL ('K', bk, 'REP', -1.)
             if ( bk /= -1. ) then
                call INREAL ('APFX', fx, 'REQ', 0.)
                call INREAL ('APFY', fy, 'REQ', 0.)
                call INREAL ('APFZ', fz, 'REQ', 0.)
                nfen = nfen + 1
                allocate(tmpfe)
                tmpfe%k  = bk
                tmpfe%fx = fx
                tmpfe%fy = fy
                tmpfe%fz = fz
                nullify(tmpfe%nextfe)
                currfe%nextfe => tmpfe
                currfe => tmpfe
                goto 113
             endif
             if ( nfen == 0 ) call msgerr (2,'no fender parameters found')
             !
             bodtmp%nfen = nfen
             allocate(bodtmp%fen(nfen,4))
             currfe => frstfe%nextfe
             do j = 1, nfen
                bodtmp%fen(j,1) = currfe%k
                bodtmp%fen(j,2) = currfe%fx
                bodtmp%fen(j,3) = currfe%fy
                bodtmp%fen(j,4) = currfe%fz
                currfe => currfe%nextfe
             enddo
             deallocate(tmpfe)
             !
          endif
          mxfen = max( mxfen, nfen )
          !
          nullify(bodtmp%nextrib)
          mbod = mbod + 1
          if ( .not.lrib ) then
             frigbod = bodtmp
             crib => frigbod
             lrib = .true.
          else
             crib%nextrib => bodtmp
             crib => bodtmp
          endif
          !
       else
          !
          ! ==========================================================================
          !
          ! PFLOW [alpha] [theta]
          !
          ! ==========================================================================
          !
          call INREAL ( 'ALPHA', pship(1), 'STA', 0. )
          call INREAL ( 'THETA', pship(2), 'STA', 1. )
          !
       endif
       !
       goto 100
       !
    endif
    !
    ! ==========================================================================
    !
    !                             | -> CONstant [visc]
    !                             |
    !            | -> Horizontal <     SMAGorinsky [cs]
    !            |                |
    !            |                |    MIXing [lm]
    !            |
    !            |
    ! VISCosity <     Vertical  KEPS [cfk] [cfe]
    !            |
    !            |
    !            |                | -> LINear
    !            |    FULL  KEPS <
    !            |                |    NONLinear
    !
    ! ==========================================================================
    !
    if ( KEYWIS ('VISC') ) then
       call INKEYW ('STA', 'H')
       if ( KEYWIS('V') ) then
          call INKEYW ('STA','KEPS')
          if ( KEYWIS('KEPS') ) then
             iturb = 1
             ltur  = 2
             call INREAL( 'CFK', pturb(2), 'STA', 0.07 )
             call INREAL( 'CFE', pturb(3), 'STA', 0.16 )
          else
             call WRNKEY
          endif
       else if ( KEYWIS('FULL') ) then
          call INKEYW ('STA','KEPS')
          if ( KEYWIS('KEPS') ) then
             ltur = 2
             call INKEYW ('STA','LIN')
             if ( KEYWIS('LIN') ) then
                iturb = 2
             else if ( KEYWIS('NONL') ) then
                iturb = 3
             else
                call WRNKEY
             endif
          else
             call WRNKEY
          endif
       else
          call IGNORE ('H')
          call INKEYW ('STA','CON')
          if ( KEYWIS('CON') ) then
             ihvisc = 1
             call INREAL('VISC',hvisc,'UNC',0.)
          else if ( KEYWIS('SMAG') ) then
             ihvisc = 2
             call INREAL('CS',csmag,'STA',0.2)
          else if ( KEYWIS('MIX') ) then
             ihvisc = 3
             call INREAL('LM',lmix,'UNC',0.)
          else
             call WRNKEY
          endif
       endif
       goto 100
    endif
    !
    ! ==========================================================================
    !
    ! VEGEtation  < [height]  [diamtr]  [nstems]  [drag] > INERtia [cm]  PORO  V
    !
    ! ==========================================================================
    !
    if ( KEYWIS ('VEGE') ) then
       iveg  = 1
       lmax  = 0
       frstvg%h = 0.
       frstvg%d = 0.
       frstvg%n = 0
       frstvg%c = 0.
       nullify(frstvg%nextvg)
       currvg => frstvg
 115   call INREAL ('HEIGHT', hh, 'REP', -1.)
       if ( hh /= -1. ) then
          if ( hh < 0. ) then
             call msgerr (2,'height is negative')
             hh = 0.
          endif
          call INREAL ('DIAMTR', dd, 'REQ', 0.)
          if ( dd < 0. ) then
             call msgerr (2,'stem diameter is negative')
             dd = 0.
          endif
          call ININTG ('NSTEMS', nn, 'REQ', 1 )
          if ( nn <= 0 ) then
             call msgerr (2,'number of stems is negative or zero')
             nn = 1
          endif
          call INREAL ('DRAG', cc, 'REQ', 0.)
          if ( cc < 0. ) then
             call msgerr (2,'drag coefficient is negative')
             cc = 0.
          endif
          lmax = lmax + 1
          allocate(tmpvg)
          tmpvg%h = hh
          tmpvg%d = dd
          tmpvg%n = nn
          tmpvg%c = cc
          nullify(tmpvg%nextvg)
          currvg%nextvg => tmpvg
          currvg => tmpvg
          goto 115
       endif
       !
       if ( KEYWIS ('INER') .or. KEYWIS ('MASS') ) then
          call INREAL ('CM', cvm, 'REQ', 0.)
       endif
       !
       call INKEYW ('STA', ' ')
       if ( KEYWIS ('PORO') ) iveg = 2
       !
       call INKEYW ('STA', ' ')
       if ( KEYWIS ('V') ) ivegw = 1
       !
       if ( lmax == 0 ) call msgerr (2,'no vegetation parameters found')
       !
       if (.not.allocated(hlayv)) allocate(hlayv(lmax))
       if (.not.allocated(bveg )) allocate(bveg (lmax))
       if (.not.allocated(nveg )) allocate(nveg (lmax))
       if (.not.allocated(cdveg)) allocate(cdveg(lmax))
       !
       currvg => frstvg%nextvg
       do i = 1, lmax
          hlayv(i) = currvg%h
          bveg (i) = currvg%d
          nveg (i) = real(currvg%n)
          cdveg(i) = currvg%c
          currvg => currvg%nextvg
       enddo
       deallocate(tmpvg)
       goto 100
    endif
    !
    ! ==========================================================================
    !
    ! CORIolis [fpar] [epsab2] (not documented)
    !
    ! ==========================================================================
    !
    if ( KEYWIS ('CORI') ) then
       coriolis = .true.
       call INREAL ('FPAR'  , pcor(1), 'UNC', 0. )
       call INREAL ('EPSAB2', pcor(2), 'STA', 0.1)
       goto 100
    endif
    !
    ! ==========================================================================
    !
    !                           | -> Sec |    | -> NONCohesive [size]            |
    ! TRANSPort [diff] [retur] <     MIn  >  <                                    > &
    !                           |    HR  |    | COHesive [tauce] [taucd] [erate] |
    !                           |    DAy |
    !
    !                                        | -> Yes
    !           [fall] [snum] [ak]  DENSity <                                       &
    !                                        | No
    !
    !           [alfa] [crsn] [cp] [ek]   (not documented)                          &
    !
    !           ANTICreep  STAndard | SVK | None
    !
    ! ==========================================================================
    !
    if ( KEYWIS ('TRANSP') .or. KEYWIS('CONST') ) then
       !
       itrans = 1
       !
       call INREAL( 'DIFF' , hdiff  , 'UNC', 0. )
       call ININTV( 'RETUR', tcret  , 'UNC', 0. )
       !
       call INKEYW ('STA','NONC')
       if ( KEYWIS('COH') .or. KEYWIS('MUD') ) then
          !
          call INREAL( 'TAUCE', psed(10), 'REQ', 0.)
          if ( .not. psed(10) > 0. ) then
             call msgerr (2,'critical erosion shear stress is negative or zero')
          endif
          call INREAL( 'TAUCD', psed(11), 'REQ', 0.)
          if ( .not. psed(11) > 0. ) then
             call msgerr (2,'critical deposition shear stress is negative or zero')
          endif
          call INREAL( 'ERATE', psed(12), 'REQ', 0.)
          if ( .not. psed(12) > 0. ) then
             call msgerr (2,'entrainment rate is negative or zero')
          endif
          !
       else
          !
          call IGNORE ('NONC')
          call INREAL( 'SIZE', psed(2), 'UNC', 0.)
          !
       endif
       !
       if ( psed(10) > 0. ) then
          call INREAL( 'FALL', psed(1), 'REQ', 0.)
       else
          call INREAL( 'FALL', psed(1), 'UNC', 0.)
       endif
       !
       call INREAL( 'SNUM', psed(3) , 'STA', 0.7)
       !
       call INREAL( 'AK'  , psed(13), 'UNC', 0. )
       !
       call INKEYW ('STA', ' ')
       if ( KEYWIS ('DENS') .or. KEYWIS('MIXT') ) then
          call INKEYW ('STA', 'Y')
          if ( KEYWIS ('N') ) then
             lmixt = .false.
          else if ( KEYWIS ('Y') ) then
             lmixt = .true.
          else
             call WRNKEY
          endif
       endif
       !
       call INREAL( 'ALFA', psed(4), 'STA', 5.5    )
       call INREAL( 'CRSN', psed(5), 'STA', 0.05   )
       call INREAL( 'CP'  , psed(6), 'STA', 0.00033)
       call INREAL( 'EK'  , psed(7), 'UNC', 0.     )
       !
       call INKEYW ('STA', ' ')
       if ( KEYWIS ('ANTIC') ) then
          call INKEYW ('STA', 'N')
          if ( KEYWIS ('N') ) then
             icreep = 0
          else if ( KEYWIS ('STA') .or. KEYWIS ('Y') ) then
             icreep = 1
          else if ( KEYWIS ('SVK') ) then
             icreep = 2
          else
             call WRNKEY
          endif
       endif
       !
       goto 100
       !
    endif
    !
    ! ==========================================================================
    !
    !                  |    STAndard      |
    !                  |                  |
    ! NONHYDrostatic  <  -> BOX            > [theta]                           &
    !                  |                  |
    !                  |    DEPthaveraged |  (not documented)
    !
    !                 SUBGrid [pmax]  REDuced [qlay]                           &
    !
    !                 SOLVer [rhsaccur] [initaccur] [maxiter] [relax] [precfq] &
    !
    !                 PREConditioner  ILUDS|ILUD|ILU|NONE                      &
    !
    !                 PROJection  ITERative [tol] [maxiter]
    !
    ! ==========================================================================
    !
    if ( KEYWIS ('NONHYD') ) then
       !
       call INKEYW ('STA','BOX')
       if ( KEYWIS('BOX') .or. KEYWIS ('KEL') .or. KEYWIS ('HERM') ) then
          !
          ihydro = 1
          !
       else if ( KEYWIS ('STA') ) then
          !
          ihydro = 2
          !
       else if ( KEYWIS ('DEP') ) then
          !
          ihydro = 3
          lsubg  = .true.
          kpmax  = 1
          !
       else
          call WRNKEY
       endif
       !
       ! implicitness factor for time integration
       !
       call INREAL ('THETA', pnums(5), 'STA', -1.)
       !
       call INKEYW ('STA','  ')
       if ( KEYWIS('SUBG') .or. KEYWIS('PD') ) then
          !
          lsubg = .true.
          call ININTG ('PMAX', kpmax, 'RQI', kmax )
          !
          if ( kpmax == 1 ) ihydro = 3
          !
       endif
       !
       call INKEYW ('STA','  ')
       if ( KEYWIS('RED') ) then
          call ININTG('QLAY', qlay, 'STA', 1 )
       endif
       !
       call INKEYW ('STA','  ')
       if ( KEYWIS('SOLV') ) then
          call INREAL('RHSACCUR' , pnums(22), 'STA', 0.01)
          call INREAL('INITACCUR', pnums(23), 'STA', 0.  )
          call ININTG('MAXITER'  , itmp1    , 'STA', 100 )
          call INREAL('RELAX'    , pnums(27), 'STA', -1. )
          call ININTG('PRECFQ'   , itmp2    , 'STA', 1   )
          pnums(25) = real(itmp1)
          pnums(28) = real(itmp2)
       endif
       !
       call INKEYW ('STA','  ')
       if ( KEYWIS('PREC') ) then
          call INKEYW ('STA', '  ')
          if ( KEYWIS('ILUDS') ) then
             icond = 1
          else if ( KEYWIS('ILUD') ) then
             icond = 2
          else if ( KEYWIS('ILU') ) then
             icond = 3
          else if ( KEYWIS('NONE') ) then
             icond = 3
             pnums(28) = -1.
          endif
       endif
       !
       call INKEYW ('STA', ' ')
       if (KEYWIS ('PROJ')) then
          iproj = 2
          !
          call INKEYW ('STA',' ')
          if ( KEYWIS ('ITER') ) lpproj = .true.
          !
          call INREAL('TOL'    , pnums(58), 'STA', 0.0001)
          call ININTG('MAXITER', itmp1    , 'STA', 50    )
          pnums(59) = real(itmp1)
          !
       endif
       !
       goto 100
       !
    endif
    !
    ! ==========================================================================
    !
    !               | North [width] |
    !               | NW    [width] |
    !               | West  [width] |
    !               | SW    [width] |
    ! SPONgelayer  <  South [width]  > | < [k]  [width] >
    !               | SE    [width] |
    !               | East  [width] |
    !               | NE    [width] |
    !
    ! ==========================================================================
    !
    if ( KEYWIS ('SPON') ) then
       if ( optg /= 5 ) then
          !
          ! structured grid
          !
          call INKEYW ('REQ',' ')
          if ( KEYWIS ('NW') ) then
             dirsi = 45.
          else if ( KEYWIS ('SW') ) then
             dirsi = 135.
          else if ( KEYWIS ('SE') ) then
             dirsi = -135.
          else if ( KEYWIS ('NE') ) then
             dirsi = -45.
          else if ( KEYWIS ('N') .or. KEYWIS('UP') ) then
             dirsi = 0.
          else if ( KEYWIS ('W') .or. KEYWIS('LE') ) then
             dirsi = 90.
          else if ( KEYWIS ('S') .or. KEYWIS('LO') ) then
             dirsi = 180.
          else if ( KEYWIS ('E') .or. KEYWIS('RI') ) then
             dirsi = -90.
          else
             call WRNKEY
          endif
          call INREAL('WIDTH',spwid,'REQ',0.)
          !
          ! select side in the chosen direction
          !
          crdm   = -1.e10
          iside  = 0
          if ( oned ) then
             cosdir = cos(pi*(dnorth+dirsi)/180.)
             sindir = sin(pi*(dnorth+dirsi)/180.)
             do i = 1, 4
                sumx = 0.
                sumy = 0.
                nump = 0
                if ( i == 2 ) then
                   if ( kgrpnt(mxc-1,1) > 1 ) then
                      sumx = xcgrid(mxc-1,1)
                      sumy = ycgrid(mxc-1,1)
                      nump = 1
                   endif
                else if ( i == 4 ) then
                   if ( kgrpnt(1,1) > 1 ) then
                      sumx = xcgrid(1,1)
                      sumy = ycgrid(1,1)
                      nump = 1
                   endif
                endif
                if ( nump > 0 ) then
                   crdp = cosdir*sumx + sindir*sumy
                   ! side with largest crdp is the one selected
                   if ( crdp > crdm ) then
                      crdm = crdp
                      iside = i
                   endif
                endif
             enddo
          else
             do i = 1, 4
                sumx = 0.
                sumy = 0.
                nump = 0
                if ( i == 1 ) then
                   do ix = 1, mxc-1
                      kc2 = kgrpnt(ix,1)
                      if ( ix > 1 ) then
                         if ( kc1 > 1 .and. kc2 > 1 ) then
                            sumx = sumx + xcgrid(ix,1) - xcgrid(ix-1,1)
                            sumy = sumy + ycgrid(ix,1) - ycgrid(ix-1,1)
                            nump = nump + 1
                         endif
                      endif
                      kc1 = kc2
                   enddo
                else if ( i == 2 ) then
                   do iy = 1, myc-1
                      kc2 = kgrpnt(mxc-1,iy)
                      if ( iy > 1 ) then
                         if ( kc1 > 1 .and. kc2 > 1 ) then
                            sumx = sumx + xcgrid(mxc-1,iy) - xcgrid(mxc-1,iy-1)
                            sumy = sumy + ycgrid(mxc-1,iy) - ycgrid(mxc-1,iy-1)
                            nump = nump + 1
                         endif
                      endif
                      kc1 = kc2
                   enddo
                else if ( i == 3 ) then
                   do ix = 1, mxc-1
                      kc2 = kgrpnt(ix,myc-1)
                      if ( ix > 1 ) then
                         if ( kc1 > 1 .and. kc2 > 1 ) then
                            sumx = sumx + xcgrid(ix-1,myc-1) - xcgrid(ix,myc-1)
                            sumy = sumy + ycgrid(ix-1,myc-1) - ycgrid(ix,myc-1)
                            nump = nump + 1
                         endif
                      endif
                      kc1 = kc2
                   enddo
                else if ( i == 4 ) then
                   do iy = 1, myc-1
                      kc2 = kgrpnt(1,iy)
                      if ( iy > 1 ) then
                         if ( kc1 > 1 .and. kc2 > 1 ) then
                            sumx = sumx + xcgrid(1,iy-1) - xcgrid(1,iy)
                            sumy = sumy + ycgrid(1,iy-1) - ycgrid(1,iy)
                            nump = nump + 1
                         endif
                      endif
                      kc1 = kc2
                   enddo
                endif
                if ( nump > 0 ) then
                   dirsid = atan2(sumy,sumx)
                   dirref = pi*(dnorth+dirsi)/180.
                   if ( cvleft ) then
                      crdp = cos(dirsid - 0.5*pi - dirref)
                   else
                      crdp = cos(dirsid + 0.5*pi - dirref)
                   endif
                   ! side with largest crdp is the one selected
                   if ( crdp > crdm ) then
                      crdm = crdp
                      iside = i
                   endif
                endif
             enddo
          endif
          if ( iside == 0 ) call msgerr (2, 'no sponge layer found')
          !
          if ( iside == 1 ) then
             spwidb = spwid
          else if ( iside == 2 ) then
             spwidr = spwid
          else if ( iside == 3 ) then
             spwidt = spwid
          else if ( iside == 4 ) then
             spwidl = spwid
          endif
          !
       else
          !
          ! unstructured mesh
          !
          nspl     = 0
          frstsp%k = 0
          frstsp%w = 0.
          nullify(frstsp%nextsp)
          currsp => frstsp
 125      call ININTG ('K', vm, 'REP', -1)
          if ( vm /= -1 ) then
             if ( vm == 0 ) call msgerr( 2, '[k] is not a valid boundary marker' )
             call INREAL ('WIDTH', spwid, 'REQ', 0.)
             nspl = nspl + 1
             allocate(tmpsp)
             tmpsp%k = vm
             tmpsp%w = spwid
             nullify(tmpsp%nextsp)
             currsp%nextsp => tmpsp
             currsp => tmpsp
             goto 125
          endif
          !
          if (.not.allocated(vmspon)) allocate(vmspon(nspl))
          if (.not.allocated(wdspon)) allocate(wdspon(nspl))
          !
          if ( nspl > 0 ) then
             currsp => frstsp%nextsp
             do i = 1, nspl
                vmspon(i) = currsp%k
                wdspon(i) = currsp%w
                currsp => currsp%nextsp
             enddo
             deallocate(tmpsp)
          endif
       endif
       goto 100
    endif
    !
    ! ==========================================================================
    !
    ! VERTical [kmax]  < [thickness] M|PERC >
    !
    ! ==========================================================================
    !
    if ( KEYWIS ('VERT') ) then
       if ( mcgrd > 1 .or. nverts > 0 ) then
          call msgerr (3, 'command READ COOR or READ UNSTRUC must follow this command')
          goto 100
       endif
       call ININTG ('KMAX',kmax,'RQI',0)
       numlay  = 0
       frstv%t = 0.
       frstv%i = 0
       nullify(frstv%nextv)
       currv => frstv
 120   call INREAL ('THICK',hh,'REP',-1.)
       if ( hh /= -1. ) then
          if ( .not. hh > 0. ) then
             call msgerr (2,'layer thickness is negative or zero')
          endif
          call INKEYW ('STA', 'PERC')
          if ( KEYWIS('M') ) then
             indl = 2
          else
             call IGNORE ('PERC')
             indl = 1
             hh = hh / 100.
          endif
          numlay = numlay + 1
          allocate(tmpv)
          tmpv%t = hh
          tmpv%i = indl
          nullify(tmpv%nextv)
          currv%nextv => tmpv
          currv => tmpv
          goto 120
       endif
       !
       if ( numlay > 0 .and. numlay /= kmax ) then
          call msgerr (1,'[kmax] is set to actual number of layers')
          kmax = numlay
       endif
       !
       kpmax = kmax
       !
       if (.not.allocated(hlay  )) allocate(hlay  (kmax))
       if (.not.allocated(indlay)) allocate(indlay(kmax))
       !
       if ( numlay == 0 ) then
          !
          ! assume equidistant distribution of layers
          !
          indlay = 1
          hlay   = 1./real(kmax)
       else
          currv => frstv%nextv
          do i = 1, kmax
             hlay  (i) = currv%t
             indlay(i) = currv%i
             currv => currv%nextv
          enddo
          deallocate(tmpv)
       endif
       goto 100
    endif
    !
    ! ==========================================================================
    !
    !             | -> EXPL [cfllow] [cflhig] &                               |
    !             |                                                           |
    !             |    RK3 | HANCock | EULer (not documented)                 |
    ! TIMEI METH <                                                             > &
    !             |    IMPL [thetac] [thetas] &                               |
    !             |                                                           |
    !             |         SOLVer [tol] [maxiter] [weight]  NEWTon           |
    !
    !       VERTical [thetau] [thetaw] [thetat]
    !
    ! ==========================================================================
    !
    if ( KEYWIS ('TIMEI') ) then
       mtimei = 1
       inewt  = 0
       call INKEYW ('STA','METH')
       if ( KEYWIS ('METH') ) then
          call INKEYW ('STA', 'EXPL')
          if ( KEYWIS('IMPL') ) then
             mtimei = 2
             call INREAL ('THETAC', pnums(1), 'STA', 0.5)
             call INREAL ('THETAS', pnums(4), 'STA', 0.5)
             !
             call INKEYW ('STA','  ')
             if ( KEYWIS('SOLV') ) then
                call INREAL('TOL'    , pnums(21), 'STA', -1.)
                call ININTG('MAXITER', itmp1    , 'STA', 100)
                call INREAL('WEIGHT' , pnums(26), 'STA', -1.)
                pnums(24) = real(itmp1)
             endif
             !
             call INKEYW ('STA',' ')
             if ( KEYWIS ('NEWT') ) inewt = 1
          else
             call IGNORE ('EXPL')
             mtimei = 1
             call INREAL ('CFLLOW', pnums(2), 'STA', 0.4)
             call INREAL ('CFLHIG', pnums(3), 'STA', 0.8)
             !
             call INKEYW ('STA',' ')
             if ( KEYWIS ('EUL') ) then
                jhk = 3
             else if ( KEYWIS ('HANC') ) then
                jhk = 2
             else if ( KEYWIS ('RK3') ) then
                jhk = 1
             endif
          endif
       else
          call WRNKEY
       endif
       !
       call INKEYW ('STA',' ')
       if ( KEYWIS ('VERT') ) then
          call INREAL ('THETAU', pnums(31), 'STA', 0.5)
          call INREAL ('THETAW', pnums(32), 'STA', 0.5)
          call INREAL ('THETAT', pnums(33), 'STA', 0.5)
       endif
       !
       goto 100
    endif
    !
    ! ==========================================================================
    !
    !                 |         | UMOM  MOMentum|HEAD  / -> Horizontal   |
    !                 |         |                      \    Vertical     |
    !                 | UPWind <                                         |
    !                 |         | WMOM  / -> Horizontal                  |
    !                 |         |       \    Vertical                    |
    !                 |                                                  |
    !                 | CORRdep                                          |
    ! DISCRETization <                                                    >    &
    !                 | TRANSPort  / -> Horizontal                       |
    !                 |            \    Vertical                         |
    !                 |                                                  |
    !                 |       | -> U                                     |
    !                 | ACUR <                                           |
    !                 |       | W                                        |
    !
    !
    !                    | NONe                                 |
    !                    |                                      |
    !                    | FIRstorder                           |
    !                    |                                      |
    !                    | HIGherorder [kappa]                  |
    !                    |                                      |
    !                    |          | -> SWEBy [phi]            |
    !                    |          |                           |
    !                    | LIMiter <  RKAPpa [kappa]            |
    !                    |          |                           |
    !                    |          | PLKAPpa [kappa] [mbound]  |
    !                    | FROmm                                |
    !                    |                                      |
    !                    | -> BDF | LUDs                        |
    !                    |                                      |
    !                   <  QUIck                                 >
    !                    |                                      |
    !                    | CUI                                  |
    !                    |                                      |
    !                    | MINMod                               |
    !                    |                                      |
    !                    | SUPerbee                             |
    !                    |                                      |
    !                    | VANLeer                              |
    !                    |                                      |
    !                    | MUScl                                |
    !                    |                                      |
    !                    | KORen                                |
    !                    |                                      |
    !                    | SMArt                                |
    !
    ! ==========================================================================
    !
    if ( KEYWIS ('DISCRET') ) then
       ieqn   = 0
       iadv   = 0
       propsc = -1
       kappa  = 0.
       mbound = 0.
       phieby = 0.
       call INKEYW ('STA',' ')
       if ( KEYWIS ('UPW') .or. KEYWIS ('ADV') ) then
          call INKEYW ('STA', 'UMOM')
          if ( KEYWIS('WMOM') ) then
             ieqn = 3
             horwinc = .true.
             verwinc = .true.
             call INKEYW ('STA', 'H')
             if ( KEYWIS('V') ) then
                iadv = 2
             else
                call IGNORE ('H')
                iadv = 1
             endif
          else
             call IGNORE ('UMOM')
             ieqn = 1
             call INKEYW ('STA','    ')
             if ( KEYWIS ('MOM') ) then
                strictmom = .true.
             else if ( KEYWIS ('HEAD') ) then
                stricthead = .true.
                strictmom  = .false.
             endif
             call INKEYW ('STA', 'H')
             if ( KEYWIS('V') ) then
                iadv = 2
             else
                call IGNORE ('H')
                iadv = 1
             endif
          endif
       else if ( KEYWIS ('CORR') ) then
          ieqn = 2
          corrdep = .true.
       else if ( KEYWIS ('TRANSP') .or. KEYWIS('CONST') ) then
          ieqn = 4
          call INKEYW ('STA', 'H')
          if ( KEYWIS('V') ) then
             iadv = 2
          else
             call IGNORE ('H')
             iadv = 1
          endif
       else if ( KEYWIS ('MIME') ) then
          mimetic = .true.
          goto 100
       else if ( KEYWIS ('ACUR') .or. KEYWIS('BCUR') .or. KEYWIS('CCUR') ) then
          call INKEYW ('STA', 'U')
          if ( KEYWIS('W') ) then
             ieqn = 6
          else
             call IGNORE ('U')
             ieqn = 5
          endif
       endif
       !
       call INKEYW ('STA',' ')
       if ( KEYWIS('NON') ) then
          propsc = 3
          kappa  = 1.
       else if ( KEYWIS('FIR') ) then
          propsc = 1
       else if ( KEYWIS('FRO') ) then
          propsc = 3
          kappa  = 0.
       else if ( KEYWIS('BDF') .or. KEYWIS ('LUD') .or. KEYWIS ('SEC') ) then
          propsc = 3
          kappa  = -1.
       else if ( KEYWIS('QUI') ) then
          propsc = 3
          kappa  = 0.5
       else if ( KEYWIS('CUI') ) then
          propsc = 3
          kappa  = 1./3.
       else if ( KEYWIS('MINM') ) then
          propsc = 4
          phieby = 1.
       else if ( KEYWIS('SUP') ) then
          propsc = 4
          phieby = 2.
       else if ( KEYWIS('VANL') ) then
          propsc = 5
          kappa  = 0.
       else if ( KEYWIS('MUS') ) then
          propsc = 6
          kappa  = 0.
          mbound = 2.
       else if ( KEYWIS('KOR') ) then
          propsc = 6
          kappa  = 1./3.
          mbound = 2.
       else if ( KEYWIS('SMA') ) then
          propsc = 6
          kappa  = 0.5
          mbound = 4.
       else if ( KEYWIS('HIG') ) then
          propsc = 3
          call INREAL ('KAPPA', kappa, 'REQ', 0.)
       else if ( KEYWIS('LIM') ) then
          call INKEYW ('STA', 'SWEB')
          if ( KEYWIS ('SWEB') ) then
             propsc = 4
             call INREAL ('PHI', phieby, 'REQ', 0.)
          else if ( KEYWIS ('RKAP') ) then
             propsc = 5
             call INREAL ('KAPPA', kappa, 'REQ', 0.)
          else if ( KEYWIS ('PLKAP') ) then
             propsc = 6
             call INREAL ('KAPPA' , kappa , 'REQ', 0.)
             call INREAL ('MBOUND', mbound, 'REQ', 0.)
          endif
       else
          goto 100
       endif
       !
       if ( ieqn == 1 .and. iadv == 1 ) then
          pnums( 6) = real(propsc)
          pnums( 7) = kappa
          pnums( 8) = mbound
          pnums( 9) = phieby
       else if ( ieqn == 1 .and. iadv == 2 ) then
          pnums(36) = real(propsc)
          pnums(37) = kappa
          pnums(38) = mbound
          pnums(39) = phieby
       else if ( ieqn == 2 ) then
          pnums(11) = real(propsc)
          pnums(12) = kappa
          pnums(13) = mbound
          pnums(14) = phieby
       else if ( ieqn == 3 .and. iadv == 1 ) then
          pnums(16) = real(propsc)
          pnums(17) = kappa
          pnums(18) = mbound
          pnums(19) = phieby
       else if ( ieqn == 3 .and. iadv == 2 ) then
          pnums(41) = real(propsc)
          pnums(42) = kappa
          pnums(43) = mbound
          pnums(44) = phieby
       else if ( ieqn == 4 .and. iadv == 1 ) then
          pnums(46) = real(propsc)
          pnums(47) = kappa
          pnums(48) = mbound
          pnums(49) = phieby
       else if ( ieqn == 4 .and. iadv == 2 ) then
          pnums(51) = real(propsc)
          pnums(52) = kappa
          pnums(53) = mbound
          pnums(54) = phieby
       else if ( ieqn == 5 ) then
          pnums(61) = real(propsc)
          pnums(62) = kappa
          pnums(63) = mbound
          pnums(64) = phieby
       else if ( ieqn == 6 ) then
          pnums(66) = real(propsc)
          pnums(67) = kappa
          pnums(68) = mbound
          pnums(69) = phieby
       endif
       !
       goto 100
    endif
    !
    ! ==========================================================================
    !
    !           |    MIN
    !           |
    !           | -> MEAN
    !  DPSopt  <
    !           |    MAX
    !           |
    !           |    SHIFt
    !
    ! ==========================================================================
    !
    if ( KEYWIS ('DPS') .or. KEYWIS ('BOTC') ) then
       dpsopt = 2
       call INKEYW ('STA','MEAN')
       if ( KEYWIS('MIN') ) then
          dpsopt = 1
       else if ( KEYWIS('MEAN') ) then
          dpsopt = 2
       else if ( KEYWIS('MAX') ) then
          dpsopt = 3
       else if ( KEYWIS('SHIF') ) then
          dpsopt = 4
       else
          call WRNKEY
       endif
       goto 100
    endif
    !
    ! ==========================================================================
    !
    !        | REGular [xpc] [ypc] [alpc] [xlenc] [ylenc] [mxc] [myc]  |
    ! CGRID <  CURVilinear [mxc] [myc]  [excval]                        >
    !        | UNSTRUCtured                                            |
    !
    !       REPeating X|Y
    !
    ! ==========================================================================
    !
    if ( KEYWIS('CGRID') .or. KEYWIS('GRID') ) then
       !
       logcom(2) = .true.
       call INKEYW ('STA', 'REG')
       if ( KEYWIS('CURV') ) then
          optg  = 3
          xpc   = 0.
          ypc   = 0.
          alpc  = 0.
          xclen = 0.
          yclen = 0.
       else if ( KEYWIS('UNSTRUC') ) then
          optg  = 5
          xpc   = 0.
          ypc   = 0.
          alpc  = 0.
          xclen = 0.
          yclen = 0.
       else
          call IGNORE ('REG')
          optg = 1
          cvleft = .true.
          if ( kspher == 0 ) then
             call READXY ('XPC', 'YPC', xpc, ypc, 'STA', 0., 0.)
             call INREAL ('ALPC',alpc,'UNC',0.)
          else
             call READXY ('XPC', 'YPC', xpc, ypc, 'REQ', 0., 0.)
             call INREAL ('ALPC',alpc,'STA',0.)
          endif
          call INREAL('XLENC',xclen,'RQI',0.)
          if ( oned ) then
             call INREAL('YLENC',yclen,'STA',0.)
             if ( yclen /= 0 ) then
                call msgerr (1, '1D simulation: [ylenc] set to zero !')
             endif
             yclen = 0.
          else
             call INREAL('YLENC',yclen,'RQI',0.)
          endif
          !
          ! alpc is made to be between -pi and pi
          !
          alpc = pi2 * (alpc/360. - nint(alpc/360.))
       endif
       !
       if ( optg == 5 ) goto 130
       !
       ! mxs is the user-defined number of cells in x-direction
       ! mys is the user-defined number of cells in y-direction
       !
       ! mxc is the number of cells in x-direction (including virtual ones)
       ! myc is the number of cells in y-direction (including virtual ones)
       !
       call ININTG ('MXC',mxs,'RQI',0)
       !
       if ( oned ) then
          call ININTG ('MYC',mys,'STA',0)
          if ( mys /= 0 ) then
             call msgerr (1, '1D simulation: [myc] set to zero !')
          endif
          mys = 0
       else
          call ININTG ('MYC',mys,'RQI',-1)
       endif
       !
       call INKEYW ('STA', ' ')
       if ( KEYWIS ('REP') .or. KEYWIS('PER') ) then
          call INKEYW ('STA', 'X')
          if ( KEYWIS ('Y') ) then
             if ( oned ) then
                call msgerr (1, '1D simulation: grid not repeated in y-direction !')
                lrepty = .false.
             else
                lrepty = .true.
             endif
          else if ( KEYWIS ('X') ) then
             lreptx = .true.
          else
             call WRNKEY
          endif
       endif
       !
       mxc = mxs + 2
       dx = xclen/mxs
       !
       myc = mys + 2
       if ( oned ) then
          dy = dx
       else
          dy = yclen/mys
       endif
       !
       ! for curvilinear grid, read exception values for grid point coordinates
       !
       if ( optg == 3 ) then
          call INKEYW ('STA', ' ')
          if ( KEYWIS ('EXC') ) then
             call INREAL ('EXCVAL', excfld(8), 'REQ', 0.)
             call INREAL ('EXCVAL', excfld(9), 'STA', excfld(8))
          endif
       endif
 130   continue
       !
       if ( ITEST >= 20 ) then
          if ( optg == 1 ) write (PRINTF,205)
          if ( optg == 3 ) write (PRINTF,206)
          if ( optg == 5 ) write (PRINTF,207)
          if ( optg /= 5 ) write (PRINTF,208) mxc, myc
          if ( optg /= 5 ) write (PRINTF,209) dx , dy
       endif
       !
       if(.not.allocated(xcgrid)) allocate(xcgrid(mxc,myc))
       if(.not.allocated(ycgrid)) allocate(ycgrid(mxc,myc))
       !
       if ( optg == 1 ) then
          cospc = cos(alpc)
          sinpc = sin(alpc)
          do j = 1, myc
             do i = 1, mxc
                xcgrid(i,j) = xpc + cospc*(i-1)*dx - sinpc*(j-1)*dy
                ycgrid(i,j) = ypc + sinpc*(i-1)*dx + cospc*(j-1)*dy
             enddo
          enddo
       endif
       !
       ! the computational grid is included in output data
       !
       if ( mxc > 0 .and. myc > 0 ) then
          allocate(opstmp)
          opstmp%psname = 'COMPGRID'
          if ( optg == 1 ) then
             opstmp%pstype = 'F'
             opstmp%opr(1) = xpc
             opstmp%opr(2) = ypc
             opstmp%opr(3) = xclen
             opstmp%opr(4) = yclen
             opstmp%opr(5) = alpc
             opstmp%opr(6) = 0.
          else if ( optg == 3 ) then
             opstmp%pstype = 'H'
             opstmp%opr(1) = float(mxc-2)
             opstmp%opr(2) = float(myc-2)
             opstmp%opr(3) = 0.
             opstmp%opr(4) = 0.
             opstmp%opr(5) = alpc
             opstmp%opr(6) = 0.
          endif
          ! store number of grid points in x- and y-directions
          opstmp%opi(1) = mxc-1
          opstmp%opi(2) = myc-1
          allocate(opstmp%xp(0))
          allocate(opstmp%yp(0))
          nullify(opstmp%nextops)
          if ( .not.lops ) then
             fops = opstmp
             cops => fops
             lops = .true.
          else
             cops%nextops => opstmp
             cops => opstmp
          endif
       endif
       !
       goto 100
    endif
    !
    ! ==========================================================================
    !
    ! INPgrid                                                                &
    !    BOTtom/WLEVel/CURrent/VX/VY/FRic/WInd/WX/WY/PRes/CORI/PORO/PSIZ/    &
    !    HSTRUC/NPLAnts/DRAFt/LABel/SALinity/TEMPerature/SEDiment/MWL/ACUR/  &
    !    AVX/AVY                                                             &
    !    | REG [xpinp] [ypinp] [alpinp]  [mxinp] [myinp]  [dxinp] [dyinp] |
    !   <  CURVilinear [stagrx] [stagry] [mxinp] [myinp]                   > &
    !    | UNSTRUCtured                                                   |
    !    (NONSTATionary [tbeginp] [deltinp] SEC/MIN/HR/DAY [tendinp])
    !
    ! ==========================================================================
    !
    if ( KEYWIS('INP') ) then
       !
       call IGNORE ('GRID')
       call INKEYW ('STA', ' ')
       igr2 = 0
       if ( KEYWIS ('BOT') ) then
          igrd   = 1
          psname = 'BOTTGRID'
       else if ( KEYWIS ('CUR') .or. KEYWIS ('VEL') ) then
          igrd   = 2
          igr2   = 3
          inituf = .true.
          initvf = .true.
          psname = 'VXGRID  '
       else if ( KEYWIS ('VX') ) then
          igrd   = 2
          inituf = .true.
          psname = 'VXGRID  '
       else if ( KEYWIS ('VY') ) then
          igrd   = 3
          initvf = .true.
          psname = 'VYGRID  '
       else if ( KEYWIS ('FR') ) then
          igrd   = 4
          psname = 'FRICGRID'
       else if ( KEYWIS ('WI') ) then
          igrd   = 5
          igr2   = 6
          psname = 'WXGRID  '
       else if ( KEYWIS ('WX') ) then
          igrd   = 5
          psname = 'WXGRID  '
       else if ( KEYWIS ('WY') ) then
          igrd   = 6
          psname = 'WYGRID  '
       else if ( KEYWIS ('WLEV') ) then
          igrd   = 7
          initsf = .true.
          psname = 'WLEVGRID'
       else if ( KEYWIS ('PORO') ) then
          igrd   = 10
          psname = 'POROGRID'
       else if ( KEYWIS ('PSIZ') ) then
          igrd   = 11
          psname = 'PSIZGRID'
       else if ( KEYWIS ('HSTRUC') ) then
          igrd   = 12
          psname = 'HSTRGRID'
       else if ( KEYWIS ('PR') ) then
          igrd   = 13
          psname = 'PRESGRID'
       else if ( KEYWIS ('NPLA') ) then
          igrd   = 14
          psname = 'NPLAGRID'
       else if ( KEYWIS ('SAL') ) then
          igrd   = 15
          psname = 'SALGRID'
       else if ( KEYWIS ('TEMP') ) then
          igrd   = 16
          psname = 'TEMPGRID'
       else if ( KEYWIS ('SED') ) then
          igrd   = 17
          psname = 'SEDGRID'
       else if ( KEYWIS ('SHIP') .or. KEYWIS ('VESS') .or. KEYWIS ('WEC') .or. KEYWIS ('FLOAT') .or. KEYWIS ('DRAF') .or. KEYWIS ('PONT') .or. KEYWIS ('PFLOW') ) then
          igrd   = 18
          psname = 'SHIPGRID'
       else if ( KEYWIS ('LAB') .or. KEYWIS ('BODY') ) then
          igrd   = 19
          psname = 'SHIPGRID'
       else if ( KEYWIS ('ACUR') .or. KEYWIS ('BCUR') .or. KEYWIS ('CCUR') ) then
          igrd   = 20
          igr2   = 21
          psname = 'AVXGRID '
       else if ( KEYWIS ('AVX') ) then
          igrd   = 20
          psname = 'AVXGRID '
       else if ( KEYWIS ('AVY') ) then
          igrd   = 21
          psname = 'AVYGRID '
       else if ( KEYWIS ('MWL') ) then
          igrd   = 22
          psname = 'MWLGRID '
       else if ( KEYWIS ('CORI') ) then
          igrd   = 23
          psname = 'CORIGRID'
       else
          igrd   = 1
          psname = 'BOTTGRID'
       endif
       !
       call SwashInputGrid ( igrd, igr2, psname )
       if (STPNOW()) return
       !
       goto 100
    endif
    !
    ! ==========================================================================
    !
    ! READinp  BOTtom/WLevel/CURrent/FRiction/WInd/PRes/COOR/CORI/PORO/PSIZ/   &
    !          HSTRUC/NPLAnts/DRAFt/LABel/SALinity/TEMPerature/SEDiment/MWL/   &
    !          ACURrent                                                        &
    !      [fac]  / 'fname1'        \
    !             \ SERIES 'fname2' /  [idla] [nhedf] ([nhedt]) (nhedvec])     &
    !      FREE / FORMAT 'form' / [idfm] / UNFORMATTED
    !
    ! or
    !
    !                        | -> ADCirc
    ! READgrid UNSTRUCtured <  TRIAngle \
    !                        | EASYmesh / 'fname'
    !
    ! ==========================================================================
    !
    if ( KEYWIS('READ') ) then
       !
       ! check whether computational grid has been defined
       !
       if ( .not.logcom(2) ) then
          call msgerr (3, 'define computational grid before reading coordinates or input grids')
       endif
       !
       call INKEYW ('REQ',' ')
       if ( KEYWIS('UNSTRUC') ) then
          !
          if ( optg /= 5 ) then
             call msgerr (3, 'define computational grid before reading unstructured grid')
          endif
          if ( oned ) then
             call msgerr (4, '1D simulation cannot be done with unstructured grid')
             return
          endif
          if ( parll ) then
             call msgerr (4, 'unstructured grid is not supported in parallel run')
             return
          endif
          !
          ! read unstructured grid
          !
          logcom(5) = .true.
          call INKEYW ('STA', 'TRIA')
          if ( KEYWIS('ADC') ) then
             grid_generator = meth_adcirc
             !
             ! bottom topography will be taken from fort.14, so
             !
             logcom(3) = .true.
             igtype(1) = 3
             leds(1)   = 2
          else if ( KEYWIS('TRIA') ) then
             grid_generator = meth_triangle
             call INCSTR ('FNAME',FILENM,'REQ', ' ')
          else if ( KEYWIS('EASY') ) then
             grid_generator = meth_easy
             call INCSTR ('FNAME',FILENM,'REQ', ' ')
          endif
          call SwanReadGrid ( FILENM, LENFNM )
          if (STPNOW()) return
          goto 140
       endif
       !
       call SwashInputField ( logcom )
       if (STPNOW()) return
       !
       ! initialize arrays for computational grid as soon as both grid coordinates and bottom are available
       !
 140   if ( optg == 1 ) then
          ! rectilinear grid
          if ( logcom(3) .and. .not.logcom(6) ) then
             call SwashInitCompGrid ( logcom )
             if (STPNOW()) return
          endif
       else if ( optg == 3 ) then
          ! cuvilinear grid
          if ( logcom(3) .and. .not.logcom(4) ) then
             call msgerr (3, 'give CGRID command and read curvilinear coordinates before reading the bottom grid')
          else if ( logcom(2) .and. logcom(3) .and. logcom(4) .and. .not.logcom(6) ) then
             call SwashInitCompGrid ( logcom )
             if (STPNOW()) return
          endif
       else if ( optg == 5 ) then
          ! unstructured grid
          if ( .not.logcom(5) ) then
             call msgerr (3, 'read unstructured grid by means of SMS/ADCIRC, Triangle or Easymesh')
             if ( logcom(4) ) call msgerr (3, 'instead of curvilinear coordinates')
          else if ( logcom(5) .and. logcom(2) .and. .not.logcom(4) .and. .not.logcom(6) ) then
             call SwanCreateEdges
             if (STPNOW()) return
             call SwanGridTopology
             if (STPNOW()) return
             call SwashInitCompUgrid ( logcom )
             if (STPNOW()) return
             !
             ! the computational grid is included in output data
             !
             if ( nverts > 0 ) then
                allocate(opstmp)
                opstmp%psname = 'COMPGRID'
                opstmp%pstype = 'U'
                opstmp%mip    = nverts
                allocate(opstmp%xp(nverts))
                allocate(opstmp%yp(nverts))
                do i = 1, nverts
                   opstmp%xp(i) = xcugrd(i)
                   opstmp%yp(i) = ycugrd(i)
                enddo
                nullify(opstmp%nextops)
                if ( .not.lops ) then
                   fops = opstmp
                   cops => fops
                   lops = .true.
                else
                   cops%nextops => opstmp
                   cops => opstmp
                endif
             endif
             !
             ! the element mesh is included in output data
             !
             if ( ncells > 0 ) then
                allocate(opstmp)
                opstmp%psname = 'ELEMESH'
                opstmp%pstype = 'U'
                opstmp%mip    = ncells
                allocate(opstmp%xp(ncells))
                allocate(opstmp%yp(ncells))
                do i = 1, ncells
                   xc = 0.
                   yc = 0.
                   do j = 1, 3
                      xc = xc + xcugrd(kvertc(j,i))
                      yc = yc + ycugrd(kvertc(j,i))
                   enddo
                   xc = xc / 3.
                   yc = yc / 3.
                   opstmp%xp(i) = xc
                   opstmp%yp(i) = yc
                enddo
                nullify(opstmp%nextops)
                if ( .not.lops ) then
                   fops = opstmp
                   cops => fops
                   lops = .true.
                else
                   cops%nextops => opstmp
                   cops => opstmp
                endif
             endif
          endif
       endif
       goto 100
    endif
    !
    ! ==========================================================================
    !
    !           |    X  |
    !           |       |
    ! SOURce   <  -> Y   >  [centre] [width] [depth] [delta]                   &
    !           |       |
    !           |   [k] |
    !
    !           | REGular  [h] [per] [dir]
    !          <                                                               &
    !           | SPECTrum [h] [per] [dir] [dd] [cycle] SEC|MIN|HR|DAY
    !
    !          SMOOthing [period] SEC|MIN|HR|DAY
    !
    ! ==========================================================================
    !
    if ( KEYWIS ('SOUR') .or. KEYWIS ('LINE') ) then
       !
       iwvgen = 1
       alobnd = .true.
       !
       ! initialize seed
       !
       if ( lseed ) then
          if ( myseed < 1000 ) myseed = 1000*myseed
          call random_seed(size=nseed)
          allocate(seedf(nseed))
          allocate(seedt(nseed))
          do i = 1, nseed
             seedf(i) = abs(myseed) +    i - 1
             seedt(i) = abs(myseed) + 2*(i - 1)
          enddo
          lseed = .false.
       endif
       !
       if ( optg /= 5 ) then
          call INKEYW ('STA', 'Y')
          if ( KEYWIS ('X') ) then
             if ( oned ) then
                call msgerr (1, '1D simulation: source line not parallel to x-direction !')
                lsrcfx = .false.
             else
                lsrcfx = .true.
             endif
          else
             call IGNORE ('Y')
             lsrcfy = .true.
          endif
       else
          call ININTG ('K', vm, 'REQ', 0)
          if ( vm == 0 ) call msgerr( 2, '[k] is not a valid boundary marker' )
          piwg(6) = real(vm)
       endif
       !
       call INREAL( 'CENTRE', piwg(1), 'REQ', 0. )
       call INREAL( 'WIDTH' , piwg(2), 'REQ', 0. )
       call INREAL( 'DEPTH' , piwg(3), 'UNC', 0. )
       call INREAL( 'DELTA' , piwg(4), 'STA', 0.5)
       !
       call INKEYW ('REQ',' ')
       !
       if ( KEYWIS('REG') .or. KEYWIS('MONO') ) then
          !
          call INREAL ('H'  , height, 'REQ', 0.)
          if ( .not. height > 0. ) call msgerr (2, 'wave height is less or equal to zero')
          call INREAL ('PER', per, 'REQ', 0.)
          if ( .not. per > 0. ) call msgerr (2, 'wave period is less or equal to zero')
          call INREAL ('DIR', wdir, 'STA', -999.)
          if ( optg /= 1 .and. .not. wdir /= -999. ) call msgerr (2, 'incident wave direction at boundary is required')
          if ( oned .and. wdir /= -999. ) then
             call msgerr (1, 'incident wave angle normal to boundary in case of 1D computation')
             wdir = -999.
          endif
          !
          ! determine incident/peak wave direction with respect to problem coordinates
          !
          if ( wdir /= -999. ) then
             wdir = DEGCNV(wdir)
             wdir = pi2 * (wdir/360. - nint(wdir/360.))
          endif
          !
          allocate(iwgtmp)
          iwgtmp%nbfs      = -999
          iwgtmp%nfreq     = 1
          iwgtmp%spparm(3) = wdir
          !
          allocate(iwgtmp%ampl (1))
          allocate(iwgtmp%omega(1))
          allocate(iwgtmp%phase(1))
          allocate(iwgtmp%theta(1))
          !
          iwgtmp%ampl (1) = height/2.
          iwgtmp%omega(1) = pi2/per
          iwgtmp%phase(1) = 0.
          iwgtmp%theta(1) = 0.
          !
          nullify(iwgtmp%nextbfs)
          if ( .not.liwg ) then
             fbfs = iwgtmp
             cubfs => fbfs
             liwg = .true.
          else
             cubfs%nextbfs => iwgtmp
             cubfs => iwgtmp
          endif
          !
       else if ( KEYWIS('SPECT') ) then
          !
          call INREAL ('H'  , spparm(1), 'REQ', 0.)
          if ( .not. spparm(1) > 0. ) call msgerr (2, 'wave height is less or equal to zero')
          call INREAL ('PER', spparm(2), 'REQ', 0.)
          if ( .not. spparm(2) > 0. ) call msgerr (2, 'wave period is less or equal to zero')
          call INREAL ('DIR', wdir, 'STA', -999.)
          if ( optg /= 1 .and. .not. wdir /= -999. ) call msgerr (2, 'incident wave direction at boundary is required')
          if ( oned .and. wdir /= -999. ) then
             call msgerr (1, 'incident wave angle normal to boundary in case of 1D computation')
             wdir = -999.
          endif
          call INREAL ('DD' , spparm(4), 'STA', 0.)
          if ( oned .and. spparm(4) /= 0. ) then
             call msgerr (1, 'no multi-directional waves in case of 1D computation')
             spparm(4) = 0.
          endif
          if ( spshape(3) == 1 ) then
             if ( spparm(4) < 0. .or. spparm(4) > 360. ) call msgerr (2, 'directional spreading is less than 0 or larger than 360 degrees')
          else
             if ( spparm(4) < 0. ) call msgerr (2, 'power of cosine is less than zero')
          endif
          call ININTV ('CYCLE', spparm(5), 'STA', 1.e10)
          if ( .not. spparm(5) >  0. ) call msgerr (2, 'cyclic period of time record is less or equal to zero')
          !
          ! determine incident/peak wave direction with respect to problem coordinates
          !
          if ( wdir /= -999. ) then
             wdir = DEGCNV(wdir)
             wdir = pi2 * (wdir/360. - nint(wdir/360.))
          endif
          spparm(3) = wdir
          !
          allocate(iwgtmp)
          iwgtmp%nbfs        = -999
          iwgtmp%nfreq       = -1
          iwgtmp%spparm(1:5) = spparm(1:5)
          !
          nullify(iwgtmp%nextbfs)
          if ( .not.liwg ) then
             fbfs = iwgtmp
             cubfs => fbfs
             liwg = .true.
          else
             cubfs%nextbfs => iwgtmp
             cubfs => iwgtmp
          endif
          !
       endif
       !
       ! apply ramp function (optionally)
       !
       call INKEYW ('STA',' ')
       if ( KEYWIS('SMOO') .or. KEYWIS('RAMP') ) then
          lrampf = .true.
          call ININTV ('PERIOD', tsmo, 'REQ', 0.)
       else
          tsmo = 0.
       endif
       piwg(5) = tsmo
       !
       goto 100
    endif
    !
    ! ==========================================================================
    !
    !                               / -> IJ < [i] [j] > | < [k] > \
    ! TEST [itest] [itrace] POINTS <                               > 'fname'
    !                               \    XY < [x] [y] >           /
    !
    ! ==========================================================================
    !
    if ( KEYWIS ('TEST') ) then
       call ININTG ('ITEST',ITEST,'STA',30)
       if ( ITEST >= 30 ) then
          if ( errpts == 0.and.INODE == MASTER ) then
             errpts = 16
             open (errpts, file='ERRPTS', status='unknown', form='formatted')
          endif
       endif
       call ININTG ('ITRACE',itrace,'UNC',0)
       if ( itrace > 0 ) then
          ltrace = .true.
       else
          ltrace = .false.
       endif
       call INKEYW ('STA', ' ')
       if ( KEYWIS('POI') ) then
          ipp = 2
          if ( optg == 5 ) ipp = 1
          if ( .not.allocated(iarr) ) allocate(iarr(ipp*mptst))
          call SwashReadTestpnts ( ipp*mptst, iarr )
          if (STPNOW()) return
          nptsta = max(1,nptst)
          if ( .not.allocated(xytst) ) allocate(xytst(ipp*nptsta))
          call SWCOPI ( iarr, xytst, ipp*nptsta )
          deallocate(iarr)
       endif
       goto 100
    endif
    !
    ! write solution to file for hot-starting
    !
    if ( KEYWIS('REST') .or. KEYWIS ('BACK') .or. KEYWIS('HOTF') .or. KEYWIS('SAVE') ) then
       if ( mxc <= 1 .and. optg /= 5 ) then
          call msgerr (2, 'command CGRID must precede this command')
       else if ( mcgrd <= 1 .and. nverts <= 0 ) then
          call msgerr(2, 'command READ BOT or READ UNSTRUC must precede this command')
       else
          call SwashBackup
          if (STPNOW()) return
       endif
       goto 100
    endif
    !
    ! process initial conditions
    !
    if ( KEYWIS ('INIT') ) then
       call SwashInitCond
       if (STPNOW()) return
       goto 100
    endif
    !
    ! process boundary conditions
    !
    if ( KEYWIS ('BOU') ) then
       itmp1  = mxc
       itmp2  = myc
       itmp3  = mcgrd
       mxc    = MXCGL
       myc    = MYCGL
       mcgrd  = MCGRDGL
       call SwashBounCond ( XGRDGL, YGRDGL, KGRPGL )
       IF (STPNOW()) return
       mxc    = itmp1
       myc    = itmp2
       mcgrd  = itmp3
       goto 100
    endif
    !
    ! in case of an empty line in the command file proceed to next line
    !
    if ( KEYWRD == '    ' ) goto 100
    !
    found = .false.
    !
    ! no grid - single point - is included in output data
    !
    allocate(opstmp)
    opstmp%psname = 'NOGRID'
    opstmp%pstype = 'P'
    opstmp%mip = 1
    allocate(opstmp%xp(1))
    allocate(opstmp%yp(1))
    opstmp%xp(1) = 0.
    opstmp%yp(1) = 0.
    nullify(opstmp%nextops)
    if ( .not.lops ) then
       fops = opstmp
       cops => fops
       lops = .true.
    else
       cops%nextops => opstmp
       cops => opstmp
    endif
    !
    ! process output locations
    !
    call SwashReqOutL ( found )
    if ( STPNOW() ) return
    if ( found ) goto 100
    !
    ! process output quantities
    !
    call SwashReqOutQ ( found )
    if ( STPNOW() ) return
    !
    if ( .not.found ) then
       LEVERR = max( LEVERR, 3 )
       call WRNKEY
    endif
    !
    goto 100
    !
 201 format ('(E', i2, '.', i1, ')')
 202 format (' Format floating point table: ', a)
 203 format ('(', i4, '(1X,E', i2, '.', i1, '))')
 204 format (' Format floating point block: ', a)
 205 format ('GRID: RECTILINEAR')
 206 format ('GRID: CURVILINEAR')
 207 format ('GRID: UNSTRUCTURED')
 208 format (' MXC: ',i6,' MYC: ',i6)
 209 format (' DX: ',e12.4,' DY: ',e12.4)
    !
end subroutine SwashReadInput
