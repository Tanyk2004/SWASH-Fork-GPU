subroutine SwashBounCond ( xcgrid, ycgrid, kgrpnt )
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
!   Specifies boundary locations and conditions
!
!   Modules used
!
    use ocpcomm2
    use ocpcomm4
    use SwashCommdata3
    use SwashCommdata4
    use m_bndspec
    use m_parall, only: nbggl
    use SwanGriddata
    use SwanGridobjects
    use SwanCompdata
!
    implicit none
!
!   Argument variables
!
    integer, dimension(mxc,myc), intent(in)      :: kgrpnt           ! index table containing the address of each (active) grid point
                                                                     ! =1; not active grid point
                                                                     ! >1; active grid point
    !
    real, dimension(mxc,myc)   , intent(in)      :: xcgrid           ! coordinates of computational grid in x-direction
    real, dimension(mxc,myc)   , intent(in)      :: ycgrid           ! coordinates of computational grid in y-direction
!
!   Local variables
!
    integer                                      :: bcimp            ! indicate how boundary condition is imposed
                                                                     ! = 1; Fourier series
                                                                     ! = 2; regular wave
                                                                     ! = 3; bichromatic wave
                                                                     ! = 4; parametric spectrum
                                                                     ! = 5; time series
                                                                     ! = 6; spectrum from file
                                                                     ! = 7; spectra from SWAN file
    integer                                      :: blayk            ! indicate how boundary condition is specified for each layer
                                                                     ! = -2; logarithmic distribution of velocity in the vertical
                                                                     ! = -1; hyperbolic cosine/Stokes wave distribution of velocity in the vertical
                                                                     ! =  0; uniform distribution of boundary condition in the vertical
                                                                     ! >  0; layer number for which boundary condition is applied
    integer, dimension(:), allocatable           :: bpix             ! indices of boundary points in x-direction
    integer, dimension(:), allocatable           :: bpiy             ! indices of boundary points in y-direction
    integer                                      :: btype            ! boundary type
                                                                     ! = 1; closed (not used here)
                                                                     ! = 2; water level opening
                                                                     ! = 3; velocity opening
                                                                     ! = 5; discharge opening
                                                                     ! = 6; Riemann invariant opening
                                                                     ! = 7; weakly reflective opening
                                                                     ! = 8; Sommerfeld radiation condition
                                                                     ! =10; outflow condition
                                                                     ! < 0; =-btype, where both normal and tangential components are described
    integer                                      :: i                ! loop counter
    integer, dimension(:), allocatable           :: iarr1            ! array containing boundary vertices of unstructured mesh
    integer, dimension(:), allocatable           :: iarr2            ! shifted array containing boundary vertices of unstructured mesh
    integer, save                                :: ient = 0         ! number of entries in this subroutine
    integer, dimension(:), allocatable           :: io               ! array storing the order of sorted array indices (Fourier series)
    integer                                      :: iiopt            ! time coding option
    integer                                      :: ish              ! shift index
    integer                                      :: iside            ! selected side of domain where boundary values are imposed
    integer                                      :: istok            ! indicates the order of Stokes wave theory
                                                                     ! = 0; use hyperbolic cosine distribution for velocity (Airy wave theory)
                                                                     ! = 1; first order Stokes wave (Airy wave theory)
                                                                     ! = 2; second order Stokes wave
                                                                     ! = 3; second order sub- and super-harmonic transfer functions
    integer                                      :: itmp             ! temporary stored integer
    integer                                      :: ix               ! index of point in x-direction
    integer                                      :: ix1              ! index of first endpoint in x-direction
    integer                                      :: ix2              ! index of second endpoint in x-direction
    integer                                      :: ixb1             ! index of first boundary endpoint in unstructured mesh
    integer                                      :: ixb2             ! index of second boundary endpoint in unstructured mesh
    integer                                      :: ixi              ! increment
    integer                                      :: iy               ! index of point in y-direction
    integer                                      :: iy1              ! index of first endpoint in y-direction
    integer                                      :: iy2              ! index of second endpoint in y-direction
    integer                                      :: j                ! index / loop counter
    integer                                      :: jbg              ! type of boundary polygon in unstructured mesh
                                                                     ! =1; sea/mainland boundary
                                                                     ! >1; island boundary
    integer                                      :: k                ! counter
    integer                                      :: kc1              ! index of boundary point of consideration
    integer                                      :: kc2              ! index of subsequent boundary point
    integer                                      :: nbps             ! number of boundary grid points per segment / side
    integer                                      :: nbv1             ! auxiliary integer to store nbv2 from previous variable part of segment / side
    integer                                      :: nbv2             ! user-defined number of boundary values
    integer                                      :: nbvl             ! user-defined number of boundary values
    integer                                      :: ndsd             ! unit reference number of file containing data
    integer                                      :: nfreq            ! number of components in Fourier series
    integer*8                                    :: ntsmo            ! enlarged smoothing period to store as integer
    integer                                      :: nump             ! number of boundary points of selected side
    integer                                      :: vm               ! boundary marker
    !
    real                                         :: ampl             ! amplitude of a Fourier component
    real, dimension(:), allocatable              :: arra             ! auxiliary array to store amplitudes of Fourier series temporarily
    real, dimension(:), allocatable              :: arrf             ! auxiliary array to store frequencies of Fourier series temporarily
    real, dimension(:), allocatable              :: arrp             ! auxiliary array to store phases of Fourier series temporarily
    real                                         :: azero            ! amplitude at zero frequency in Fourier series
    real                                         :: cosdir           ! cosine of a direction
    real                                         :: crdm             ! largest crdp resulting in a side to be selected
    real                                         :: crdp             ! cosine of direction of the normal to side
    real                                         :: DEGCNV           ! indicates Cartesian or nautical degrees
    real                                         :: det              ! determinant
    real                                         :: dirref           ! user-defined direction of the normal to a side w.r.t. positive x-axis (pointing east)
    real                                         :: dirsi            ! user-defined direction of the normal to a side of domain
    real                                         :: dirsid           ! direction of the normal to a side to be selected
    real                                         :: dir1             ! first wave direction of bichromatic wave
    real                                         :: dir2             ! second wave direction of bichromatic wave
    real                                         :: dxa              ! an arbitrary mesh size in x-direction
    real                                         :: dya              ! an arbitrary mesh size in y-direction
    real                                         :: fac              ! an auxiliary factor
    real                                         :: height           ! wave height of regular wave
    real                                         :: h1               ! first wave height of bichromatic wave
    real                                         :: h2               ! second wave height of bichromatic wave
    real                                         :: omega            ! angular frequency of a Fourier component
    real                                         :: per              ! wave period of regular wave
    real                                         :: phase            ! phase of a Fourier component
    real                                         :: rdist            ! distance between two points
    real                                         :: rlen1            ! auxiliary real to store rlen2 from previous variable part of segment / side
    real                                         :: rlen2            ! user-defined distance from a point of segment / side to next point
    real                                         :: rtmp             ! temporary stored real
    real                                         :: sindir           ! sine of a direction
    real, dimension(5)                           :: spparm           ! parameters used for computation of incident spectrum and subsequent Fourier series
                                                                     ! =1; wave height (significant or rms)
                                                                     ! =2; wave period (peak or mean)
                                                                     ! =3; incident or peak wave direction w.r.t. normal on the boundary
                                                                     ! =4; directional distribution coefficient
                                                                     ! =5; cyclic period of time series
    real                                         :: sumx             ! sum of mesh sizes in x-direction
    real                                         :: sumy             ! sum of mesh sizes in y-direction
    real                                         :: tcycl            ! cyclic period of time record
    real                                         :: tsmo             ! period for smoothing boundary values during cold start
    real                                         :: t1               ! first wave period of bichromatic wave
    real                                         :: t2               ! second wave period of bichromatic wave
    real                                         :: wdir             ! incident or peak wave direction with respect to problem coordinates
    real                                         :: x1               ! x-coordinate of boundary point of consideration
    real                                         :: x2               ! x-coordinate of subsequent boundary point
    real                                         :: x3               ! x-coordinate of an arbitrary point
    real                                         :: xc               ! broken x-coordinate of boundary point
    real                                         :: xp               ! x-coordinate of grid point
    real                                         :: y1               ! y-coordinate of boundary point of consideration
    real                                         :: y2               ! y-coordinate of subsequent boundary point
    real                                         :: y3               ! y-coordinate of an arbitrary point
    real                                         :: yc               ! broken y-coordinate of boundary point
    real                                         :: yp               ! y-coordinate of grid point
    !
    logical                                      :: BOUNPT           ! indicates whether a point is a boundary point or not
    logical                                      :: ccw              ! indicates counterclockwise or clockwise direction
    logical                                      :: KEYWIS           ! indicates whether keyword in user manual is found or not
    logical                                      :: lfrst1           ! logical indicates the start of a series of actions to be performed
    logical                                      :: lfrst2           ! another logical indicates the start of a series of actions to be performed
    logical                                      :: lfrst3           ! just another logical indicates the start of a series of actions to be performed
    logical                                      :: lfrst4           ! just another logical indicates the start of a series of actions to be performed
    logical                                      :: locgri           ! indicates whether boundary point is given as pair of indices or location
    logical                                      :: lriem = .false.  ! indicates whether linearized Riemann invariant is imposed or not
    logical                                      :: STPNOW           ! indicates that program must stop
    !
    character(80)                                :: msgstr           ! string to pass message
    !
    type(bfldat), pointer                        :: bfltmp           ! list containing parameters for boundary condition file
    type(bfldat), save, pointer                  :: cubfl            ! current item in list of boundary condition file
    !
    type(bfsdat), pointer                        :: bfstmp           ! list containing parameters for Fourier series
    !
    type(bgpdat), pointer                        :: bgptmp           ! list containing parameters for boundary grid points
    !
    type xypt                                                        ! linked list for boundary grid points
       integer               :: ix, iy
       type(xypt), pointer   :: nextxy
    end type xypt
    type(xypt), target       :: frst
    type(xypt), pointer      :: curr, tmp
    !
    type fspt                                                        ! linked list for Fourier series parameters
       real                  :: ampl, omega, phase
       type(fspt), pointer   :: nextfs
    end type fspt
    type(fspt), target       :: frstf
    type(fspt), pointer      :: currf, tmpf
    !
    type(verttype), dimension(:), pointer        :: vert             ! datastructure for vertices with their attributes
!
!   Structure
!
!   Description of the pseudo code
!
!   Source text
!
    call strace (ient,'SwashBounCond')
    !
    ! point to vertex object
    !
    vert => gridobject%vert_grid
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
    if ( optg == 5 ) then
       !
       ! in case of unstructured grid, make list of boundary points in ascending order
       !
       call SwanBpntlist
       !
    endif
    !
    call INKEYW ('REQ',' ')
    if ( KEYWIS ('SHAP') ) then
       !
       ! specification of the spectral shape
       !
       ! ==========================================================================
       !
       !                     | PM                 |    | -> SIGnificant |
       ! BOUndspec  SHAPe   <  -> JONswap [gamma]  >  <                  >   &
       !                     | TMA                |    | RMS            |
       !
       !             | -> PEAK |           | DEGRees  |
       !            <           >   DSPR  <            >
       !             | MEAN    |           | -> POWer |
       !
       ! ==========================================================================
       !
       call INKEYW ('STA', 'JON')
       if ( KEYWIS ('JON') ) then
          spshape(2) = 2
          call INREAL ('GAMMA', gamma, 'STA', 3.3)
       else if ( KEYWIS ('PM') ) then
          spshape(2) = 1
          gamma      = 1.
       else if ( KEYWIS ('TMA') ) then
          spshape(2) = 3
          gamma      = 3.3
       endif
       ! significant or rms wave height
       call INKEYW ('STA', ' ')
       if ( KEYWIS('RMS') ) then
          spshape(1) = 1
       else
          call IGNORE ('SIG')
          spshape(1) = 2
       endif
       ! peak or mean frequency
       call INKEYW ('STA', ' ')
       if ( KEYWIS('MEAN') ) then
          spshape(2) = -spshape(2)
       else
          call IGNORE ('PEAK')
       endif
       ! directional distribution given by degrees or by power
       call IGNORE ('DSPR')
       call INKEYW ('STA', 'POW')
       if ( KEYWIS('DEGR') ) then
          spshape(3) = 1
       else
          call IGNORE ('POW')
          spshape(3) = 2
       endif
       if ( ITEST >= 50 ) write (PRINTF,201) spshape(1), spshape(2), spshape(3)
       !
    else
       !
       ! boundary condition specification
       !
       ! ==========================================================================
       !
       !                        | North |
       !                        | NW    |
       !                        | West  |
       !                        | SW    |          | -> CCW     |
       !            | -> SIDE  <  South  > | [k]  <              >   |
       !            |           | SE    |          | CLOCKWise  |    |
       !            |           | East  |                            |
       !            |           | NE    |                            |
       ! BOUndcond <                                                  >           &
       !            |           | -> XY  < [x] [y] >           |     |
       !            | SEGment  <                                >    |
       !                        |    IJ  < [i] [j] > | < [k] > |
       !
       !   BTYPe WLEV|VEL|VX|VY|DISCH|RIEMann|LRIEman|WEAKrefl|SOMMerfeld|OUTFlow &
       !
       !   LAYer [k] | HYPerbolic | LOGarithmic                                   &
       !
       !   SMOOthing [period] SEC|MIN|HR|DAY                                      &
       !
       !   ADDBoundwave | ADDIG                                                   &
       !
       !               | FOURier     [azero] < [ampl] [omega] [phase] >
       !               | REGular     [h] [per] [dir]
       !               | BICHromatic [h1] [h2] [per1] [per2] [dir1] [dir2]
       !   | UNIForm  <  SPECTrum    [h] [per] [dir] [dd] [cycle] SEC|MIN|HR|DAY
       !   |           | SERIes      'fname' [itmopt]
       !   |           | SPECFile    'fname' [cycle] SEC|MIN|HR|DAY
       !  <                                                                       &
       !   |           | FOURier  < [len] [azero] < [ampl] [omega] [phase] > >
       !   |           | REGular  < [len] [h] [per] [dir] >
       !   |           | BICHrom  < [len] [h1] [h2] [per1] [per2] [dir1] [dir2] >
       !   | VARiable <  SPECTrum < [len] [h] [per] [dir] [dd] [cycle] S|MI|HR|DA >
       !               | SERIes   < [len] 'fname' [itmopt] >
       !               | SPECFile < [len] 'fname' [cycle] SEC|MIN|HR|DAY >
       !               | SPECSwan 'fname' [cycle] SEC|MIN|HR|DAY
       !
       ! ==========================================================================
       !
       if ( mxc <= 1 .and. optg /= 5 ) then
          call msgerr (3, 'command CGRID must precede this command')
          return
       endif
       if ( mcgrd <= 1 .and. nverts <= 0 ) then
          call msgerr (3, 'command READ BOT or READ UNSTRUC must precede this command')
          return
       endif
       if ( .not.alobnd ) then
          nbfils = 0
          nbval  = 0
          nbgrpt = 0
          nbggl  = 0
          alobnd = .true.
       endif
       !
       ! first define side or segment
       !
       nbps    = 0
       frst%ix = 0
       frst%iy = 0
       nullify(frst%nextxy)
       curr => frst
       call INKEYW ('REQ',' ')
       if ( KEYWIS ('SEG') ) then
          !
          ! boundary condition on a segment
          !
          call INKEYW ('STA','XY')
          if ( KEYWIS('XY') .or. KEYWIS ('LOC') ) then
             locgri = .true.
          else if ( KEYWIS('IJ') .or. KEYWIS ('GRI') ) then
             locgri = .false.
          else
             call WRNKEY
          endif
          !
          ix1    = 1
          iy1    = 1
          lfrst1 = .true.
          !
          ! loop over points describing the segment
          !
          do
             if (locgri) then
                call READXY ('XP','YP', xp, yp, 'REP', -1.e10, -1.e10)
                if ( xp < -.9e10 ) goto 100
                if ( optg /= 5 ) then
                   !
                   ! structured grid
                   !
                   call CVMESH ( xp, yp, xc, yc, kgrpnt, xcgrid, ycgrid )
                   ix2 = nint(xc) + 1
                   iy2 = nint(yc) + 1
                   if ( .not.BOUNPT(ix2, iy2, kgrpnt) ) then
                      call msgerr (2, 'invalid boundary point')
                      write (PRTEST, 202) xp+xoffs, yp+yoffs, xc, yc, ix2, iy2
                   endif
                   !
                else
                   !
                   ! unstructured grid
                   !
                   call SwanFindPoint ( xp, yp, ix2 )
                   if ( ix2 < 0 ) then
                      write (msgstr, '(a,f12.4,a,f12.4,a)') 'boundary point (',xp+xoffs,',',yp+yoffs,') not part of computational grid'
                      call msgerr( 2, trim(msgstr) )
                   endif
                   if ( vert(ix2)%atti(VMARKER) /= 1 ) then
                      write (msgstr, '(a,f12.4,a,f12.4,a)') 'vertex (',xp+xoffs,',',yp+yoffs,') is not a valid boundary point'
                      call msgerr( 2, trim(msgstr) )
                   endif
                   !
                endif
             else
                call ININTG ('I' , ix2, 'REP', -1)
                if ( ix2 < 0 ) goto 100
                if ( optg /= 5 ) then
                   !
                   call ININTG ('J' , iy2, 'REQ', 1)
                   !
                else
                   !
                   if ( ix2 <= 0 .or. ix2 > nverts ) then
                      write (msgstr,'(i4,a)') ix2, ' is not a valid vertex index'
                      call msgerr( 2, trim(msgstr) )
                   else if ( vert(ix2)%atti(VMARKER) /= 1 ) then
                      write (msgstr,'(a,i4,a)') 'vertex with index ',ix2, ' is not a valid boundary point'
                      call msgerr( 2, trim(msgstr) )
                   endif
                   !
                endif
             endif
             if ( ITEST >= 80 .and. optg /= 5 ) write (PRTEST, 202) xcgrid(ix2,iy2)+xoffs, ycgrid(ix2,iy2)+yoffs, xc, yc, ix2, iy2
             !
             ! generate intermediate points on the segment
             !
             if ( optg /= 5 ) then
                !
                ! structured grid
                !
                if ( ix2 > 0 .and. ix2 < mxc .and. iy2 > 0 .and. iy2 < myc ) then
                   if ( lfrst1 ) then
                      k = 1
                      lfrst1 = .false.
                   else
                      k = max (abs(ix2-ix1), abs(iy2-iy1))
                      if ( (ix2 == ix1 .and. lreptx) .or. (iy2 == iy1 .and. lrepty) ) call msgerr (3, 'this boundary segment is part of the repeating grid')
                   endif
                   do i = 1, k
                      fac = real(i) / real(k)
                      !
                      if ( .not.oned ) then
                         ix = ix1 + nint(fac*real(ix2-ix1))
                         iy = iy1 + nint(fac*real(iy2-iy1))
                      else
                         ix = ix1 + nint(fac*real(ix2-ix1))
                         iy = iy1
                      endif
                      !
                      if ( ITEST >= 80 ) write (PRTEST, *) ' boundary point ', fac, ix, iy
                      if ( kgrpnt(ix,iy) > 1 ) then
                         nbps = nbps + 1
                         allocate(tmp)
                         tmp%ix = ix
                         tmp%iy = iy
                         nullify(tmp%nextxy)
                         curr%nextxy => tmp
                         curr => tmp
                      endif
                   enddo
                else
                   call msgerr (2, 'boundary point outside computational grid')
                endif
                !
                iy1 = iy2
                !
             else
                !
                ! unstructured grid
                !
                if ( lfrst1 ) then
                   ixb1 = vert(ix2)%atti(BINDX)
                   jbg  = vert(ix2)%atti(BPOL)
                   det  = 1.
                   lfrst1 = .false.
                else
                   ixb1 = vert(ix1)%atti(BINDX)
                   jbg  = vert(ix1)%atti(BPOL)
                   !
                   ! 1) the data along the given segment can be imposed in counterclockwise or clockwise direction
                   ! 2) content of array blist is ordered in counterclockwise manner for sea/mainland boundary (jbg=1) and
                   !    clockwise for island boundary (jbg>1)
                   ! 3) therefore, determine orientation by means of the determinant of two endpoints of the given segment
                   !    and an arbitrary point inside domain
                   !
                   ! first endpoint of segment
                   x1 = vert(ix1)%attr(VERTX)
                   y1 = vert(ix1)%attr(VERTY)
                   !
                   ! second endpoint of segment
                   x2 = vert(ix2)%attr(VERTX)
                   y2 = vert(ix2)%attr(VERTY)
                   !
                   ! an arbitrary internal point
                   dxa = 0.1*mingsiz
                   dya = 0.1*mingsiz
                   do i = 1, 4
                      x3 = x1 - dxa
                      y3 = y1 - dya
                      call SwanFindPoint ( x3, y3, ix )
                      if ( jbg > 1 .and. ix < 0 ) then
                         x3 = x1 + dxa
                         y3 = y1 + dya
                         exit
                      else if ( jbg==1 .and. ix > 0 ) then
                         exit
                      endif
                      if ( mod(i,2) == 0 ) then
                         dxa = -dxa
                      else
                         dya = -dya
                      endif
                   enddo
                   !
                   det = (y3-y1)*(x2-x1) - (y2-y1)*(x3-x1)
                   if ( det > 0. ) then
                      ! take next boundary point in counterclockwise direction
                      ixb1 = mod(ixb1,nbpt(jbg)) + 1
                   else
                      ! take next boundary point in clockwise direction
                      ixb1 = nbpt(jbg) - mod(nbpt(jbg)+1-ixb1,nbpt(jbg))
                   endif
                endif
                ixb2 = vert(ix2)%atti(BINDX)
                !
                ! determine order of counting
                if ( ixb1 > ixb2 ) then
                   if ( det < 0. ) then
                      ixi  = -1
                   else
                      ixi  = 1
                      ixb2 = ixb2 + nbpt(jbg)
                   endif
                else
                   if ( det > 0. ) then
                      ixi  = 1
                   else
                      ixi  = -1
                      ixb1 = ixb1 + nbpt(jbg)
                   endif
                endif
                !
                do i = ixb1, ixb2, ixi
                   j = mod(i,nbpt(jbg))
                   if ( j == 0 ) j = nbpt(jbg)
                   nbps = nbps + 1
                   ix = blist(j,jbg)
                   vert(ix)%atti(VBC) = 1
                   allocate(tmp)
                   tmp%ix = ix
                   nullify(tmp%nextxy)
                   curr%nextxy => tmp
                   curr => tmp
                enddo
             endif
             !
             ix1 = ix2
             !
          enddo
          !
 100      if ( nbps == 0 ) call msgerr (1, 'no points on the boundaries found')
          if ( lfrst1 ) call msgerr (1, 'at least two points needed for a segment')
          !
       else
          !
          ! boundary condition on one side of the computational grid
          !
          call IGNORE ('SIDE')
          !
          if ( optg == 3 ) then
             call msgerr(2, 'keyword SIDE should not be used for curvilinear grid')
          endif
          !
          ! specification of side for which boundary condition is given
          !
          if ( optg /= 5 ) then
             call INKEYW ('REQ',' ')
             if ( KEYWIS ('NW') ) then
                dirsi = 45.
             else if ( KEYWIS ('SW') ) then
                dirsi = 135.
             else if ( KEYWIS ('SE') ) then
                dirsi = -135.
             else if ( KEYWIS ('NE') ) then
                dirsi = -45.
             else if ( KEYWIS ('N') ) then
                dirsi = 0.
             else if ( KEYWIS ('W') ) then
                dirsi = 90.
             else if ( KEYWIS ('S') ) then
                dirsi = 180.
             else if ( KEYWIS ('E') ) then
                dirsi = -90.
             else
                call WRNKEY
             endif
          else
             call ININTG ('K', vm, 'REQ', 0)
          endif
          !
          ! go along boundary clockwise or counterclockwise (default)
          !
          call INKEYW ('STA', 'CCW')
          if ( KEYWIS('CLOCKW') ) then
             ccw = .false.
          else
             call IGNORE ('CCW')
             ccw = .true.
          endif
          !
          ! select side in the chosen direction
          !
          if ( optg /= 5 ) then
             !
             ! structured grid
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
                   if ( ITEST >= 60 ) write (PRTEST, 203) i, nump, sumx, sumy, dirsid*180/pi, dirref*180/pi, crdp, cvleft
                enddo
             endif
             if ( iside == 0 ) call msgerr (2, 'no open boundary found')
             !
             if ( iside == 1 ) then
                ix1 = 2
                iy1 = 1
                ix2 = mxc-1
                iy2 = 1
             else if ( iside == 2 ) then
                ix1 = mxc-1
                iy1 = 2
                if ( oned ) iy1 = 1
                ix2 = mxc-1
                iy2 = myc-1
             else if ( iside == 3 ) then
                ix1 = mxc-1
                iy1 = myc-1
                ix2 = 2
                iy2 = myc-1
             else if ( iside == 4 ) then
                ix1 = 1
                iy1 = myc-1
                ix2 = 1
                iy2 = 2
                if ( oned ) iy2 = 1
             endif
             !
             if ( .not.ccw .eqv. cvleft ) then
                !
                ! swap end points
                !
                ix  = ix1
                iy  = iy1
                ix1 = ix2
                iy1 = iy2
                ix2 = ix
                iy2 = iy
                !
             endif
             !
             if ( ITEST >= 50 ) write (PRINTF, 204) iside, ix1, iy1, xcgrid(ix1,iy1)+xoffs, ycgrid(ix1,iy1)+yoffs, &
                                                           ix2, iy2, xcgrid(ix2,iy2)+xoffs, ycgrid(ix2,iy2)+yoffs
             if ( ix2 == ix1 .and. lreptx ) then
                call msgerr (3, 'this boundary side is periodic in x-direction')
             else if ( iy2 == iy1 .and. lrepty ) then
                call msgerr (3, 'this boundary side is periodic in y-direction')
             endif
             k = max( abs(ix2-ix1), abs(iy2-iy1) )
             do i = 0, k
                if ( k == 0 ) then
                   fac = 0.
                else
                   fac = real(i) / real(k)
                endif
                ix = ix1 + nint(fac*real(ix2-ix1))
                iy = iy1 + nint(fac*real(iy2-iy1))
                if ( kgrpnt(ix,iy) > 1 ) then
                   nbps = nbps + 1
                   allocate(tmp)
                   tmp%ix = ix
                   tmp%iy = iy
                   nullify(tmp%nextxy)
                   curr%nextxy => tmp
                   curr => tmp
                endif
             enddo
             !
          else
             !
             ! unstructured grid
             !
             do jbg = 1, nbpol
                !
                ! first boundary polyogon is assumed an outer one
                ! (sea/mainland boundary) and hence, content of blist
                ! is ordered in counterclockwise manner
                !
                if ( jbg==1 .eqv. ccw ) then
                   ixb1 = 1
                   ixb2 = nbpt(jbg)
                   ixi  = 1
                else
                   ixb1 = nbpt(jbg)
                   ixb2 = 1
                   ixi  = -1
                endif
                !
                allocate(iarr1(sum(nbpt)))
                k = 0
                do i = ixb1, ixb2, ixi
                   ix = blist(i,jbg)
                   if ( vmark(ix) == vm ) then
                      k = k + 1
                      iarr1(k) = i
                   endif
                enddo
                !
                if ( k /= 0 ) then
                   !
                   allocate(iarr2(k))
                   iarr2(1:k) = iarr1(1:k)
                   ish = 0
                   do i = 2, k
                      if ( iarr2(i) /= iarr2(i-1) + ixi ) then
                         ish = i - 1
                         exit
                      endif
                   enddo
                   iarr2 = cshift(iarr2,ish)
                   !
                   do i = 1, k
                      j = iarr2(i)
                      ix = blist(j,jbg)
                      nbps = nbps + 1
                      vert(ix)%atti(VBC) = 1
                      allocate(tmp)
                      tmp%ix = ix
                      nullify(tmp%nextxy)
                      curr%nextxy => tmp
                      curr => tmp
                   enddo
                   deallocate(iarr2)
                   !
                endif
                deallocate(iarr1)
                !
             enddo
             !
          endif
          !
       endif
       !
       ! boundary type
       !
       call INKEYW ('REQ',' ')
       if ( KEYWIS('BTYP') ) then
          call INKEYW('STA', 'VEL')
          if ( KEYWIS('WLEV') ) then
             btype = 2
          else if ( KEYWIS('VEL') ) then
             btype = 3
          else if ( KEYWIS('VX') ) then
             btype = 3
          else if ( KEYWIS('VY') ) then
             btype = 3
          else if ( KEYWIS('DISCH') ) then
             btype = 5
          else if ( KEYWIS('RIEM') ) then
             btype = 6
          else if ( KEYWIS('LRIE') ) then
             btype = 6
             lriem = .true.
          else if ( KEYWIS('WEAK') ) then
             btype = 7
          else if ( KEYWIS('SOMM') .or. KEYWIS('RADI') ) then
             btype = 8
          else if ( KEYWIS('OUTF') ) then
             btype = 10
          else
             call WRNKEY
          endif
       else
          btype = 3
       endif
       !
       ! layer number / vertical distribution (optionally)
       !
       call INKEYW ('STA','LAY')
       if ( KEYWIS ('LAY') ) then
          if ( btype == 8 .or. btype == 10 ) then
             call ININTG ('K', blayk, 'STA', 0)
          else
             call ININTG ('K', blayk, 'REQ', 0)
          endif
          if ( blayk < 0 .or. blayk > kmax ) then
             call msgerr (2, 'invalid layer number')
             blayk = 0
          endif
       else if ( KEYWIS('HYP') ) then
          blayk = -1
       else if ( KEYWIS('LOG') ) then
          if ( kmax == 1 ) then
             call msgerr (2, 'logarithmic distribution not allowed in depth-averaged mode')
             blayk = 0
          else
             blayk = -2
          endif
       else
          blayk = 0
       endif
       !
       if ( blayk >   0 .and. btype /= 3 .and. btype /= 5                  ) call msgerr (3, 'specification of layer number not allowed for this boundary type')
       if ( blayk == -1 .and. btype /= 3 .and. btype /= 7                  ) call msgerr (3, 'hyperbolic cosine distribution not allowed for this boundary type')
       if ( blayk == -2 .and. btype /= 3 .and. btype /= 5 .and. btype /= 7 ) call msgerr (3, 'logarithmic distribution not allowed for this boundary type')
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
       !
       if ( tsmo < 0.25*86400. ) then
          ntsmo = nint(1e5*tsmo)
       else
          ntsmo = nint(1e3*tsmo)
          ntsmo = ntsmo*100
       endif
       !
       ! add bound waves (optionally)
       !
       call INKEYW ('STA',' ')
       if ( KEYWIS('ADDB') ) then
          istok  = 2
          ibound = 0
       else if ( KEYWIS('ADDIG') ) then
          istok  = 1
          ibound = 1
       else
          istok  = 1
          ibound = 0
       endif
       if ( blayk == -1 ) istok = 0
       !
       if ( istok == 2 .or. ibound == 1 ) then
          if ( btype /= 3 .and. btype /= 7 ) then
             call msgerr (1, 'addition of bound waves not applied for this boundary type')
             istok  = 1
             ibound = 0
          endif
          if ( blayk /= -1 .and. blayk /= 0 ) then
             call msgerr (1, 'addition of bound waves not applied for this vertical distribution')
             istok  = 1
             ibound = 0
          endif
       endif
       !
       ! boundary condition given as Fourier series, spectrum or time series (from file)
       !
       curr => frst%nextxy
       if ( btype == 8 .or. btype == 10 ) then
          call INKEYW ('STA','UNIF')
       else
          call INKEYW ('REQ',' ')
       endif
       if ( KEYWIS('UNIF') .or. KEYWIS('CON') ) then
          !
          call INKEYW('STA', 'FOUR')
          if ( KEYWIS('FOUR') ) then
             !
             if ( btype == 8 .or. btype == 10 ) then
                call INREAL ('AZE', azero, 'STA', 0.)
             else
                call INREAL ('AZE', azero, 'REQ', 0.)
             endif
             ! Fourier series interpreted as short waves
             if ( ( istok == 2 .or. ibound == 1 ) .and. blayk == 0 ) blayk = -1
             if ( blayk == -1 ) azero = -1.e10
             !
             nfreq       = 0
             frstf%ampl  = 0.
             frstf%omega = 0.
             frstf%phase = 0.
             nullify(frstf%nextfs)
             currf => frstf
             do
                call INREAL ('AMP', ampl , 'REP', -1.e10)
                if ( ampl < -.9e10 ) goto 110
                call INREAL ('OME', omega, 'STA', 0.)
                call INREAL ('PHA', phase, 'STA', 0.)
                phase = pi2 * (phase/360. - nint(phase/360.))
                !
                nfreq = nfreq + 1
                allocate(tmpf)
                tmpf%ampl  = ampl
                tmpf%omega = omega
                tmpf%phase = phase
                nullify(tmpf%nextfs)
                currf%nextfs => tmpf
                currf => tmpf
             enddo
             !
 110         nbval = nbval + 1
             nbvl  = nbval
             allocate(bfstmp)
             bfstmp%nbfs  = nbval
             bfstmp%nfreq = nfreq
             bfstmp%azero = azero
             bfstmp%spparm(3) = -999.    ! indicent angle is assumed to be normal on the boundary
             !
             ! also include sub- and super-harmonic transfer functions
             if ( istok == 2 .and. nfreq > 1 ) istok = 3
             !
             if ( nfreq > 0 ) then
                !
                allocate(arra(nfreq))
                allocate(arrf(nfreq))
                allocate(arrp(nfreq))
                currf => frstf%nextfs
                do i = 1, nfreq
                   arra(i) = currf%ampl
                   arrf(i) = currf%omega
                   arrp(i) = currf%phase
                   currf => currf%nextfs
                enddo
                deallocate(tmpf)
                !
                ! initialize the order of array indices
                !
                allocate(io(nfreq))
                do i = 1, nfreq
                   io(i) = i
                enddo
                !
                ! in case of bound waves, sort the frequencies
                !
                if ( istok == 2 .or. ibound == 1 ) then
                   do i = 1, nfreq-1
                      do j = i+1, nfreq
                         if ( arrf(i) > arrf(j) ) then
                            rtmp    = arrf(j)
                            itmp    = io(j)
                            arrf(j) = arrf(i)
                            io(j)   = io(i)
                            arrf(i) = rtmp
                            io(i)   = itmp
                         endif
                      enddo
                   enddo
                endif
                !
                allocate(bfstmp%ampl (nfreq))
                allocate(bfstmp%omega(nfreq))
                allocate(bfstmp%phase(nfreq))
                allocate(bfstmp%theta(nfreq))
                do i = 1, nfreq
                   bfstmp%ampl (i) = arra(io(i))
                   bfstmp%omega(i) = arrf(i)
                   bfstmp%phase(i) = arrp(io(i))
                   bfstmp%theta(i) = 0.
                enddo
                !
                deallocate(arra,arrf,arrp,io)
                !
             endif
             !
             nullify(bfstmp%nextbfs)
             if ( .not.lbfs ) then
                fbfs = bfstmp
                cubfs => fbfs
                lbfs = .true.
             else
                cubfs%nextbfs => bfstmp
                cubfs => bfstmp
             endif
             !
          else if ( KEYWIS('REG') .or. KEYWIS('MONO') ) then
             !
             if ( btype /= 3 .and. btype /= 7 ) call msgerr (3, 'regular wave not allowed for this boundary type')
             if ( .not.oned ) btype = -btype
             !
             if ( blayk /= -1 .and. blayk /= 0 ) then
                call msgerr (1, 'only vertical velocity profile derived from Stokes wave theory is allowed in case of regular wave')
             endif
             blayk = -1
             !
             if ( ibound == 1 ) then
                call msgerr (1, 'addition of bound long wave correction not allowed for regular wave')
                ibound = 0
             endif
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
             nbval = nbval + 1
             nbvl  = nbval
             allocate(bfstmp)
             bfstmp%nbfs      = nbval
             bfstmp%nfreq     = 1
             bfstmp%azero     = -1.e10
             bfstmp%spparm(3) = wdir
             !
             allocate(bfstmp%ampl (1))
             allocate(bfstmp%omega(1))
             allocate(bfstmp%phase(1))
             allocate(bfstmp%theta(1))
             !
             bfstmp%ampl (1) = height/2.
             bfstmp%omega(1) = pi2/per
             bfstmp%phase(1) = 0.
             bfstmp%theta(1) = 0.
             !
             nullify(bfstmp%nextbfs)
             if ( .not.lbfs ) then
                fbfs = bfstmp
                cubfs => fbfs
                lbfs = .true.
             else
                cubfs%nextbfs => bfstmp
                cubfs => bfstmp
             endif
             !
          else if ( KEYWIS('BICH') ) then
             !
             if ( btype /= 3 .and. btype /= 7 ) call msgerr (3, 'bichromatic wave not allowed for this boundary type')
             if ( .not.oned ) btype = -btype
             !
             if ( blayk /= -1 .and. blayk /= 0 ) then
                call msgerr (1, 'only vertical velocity profile derived from Stokes wave theory is allowed in case of bichromatic wave')
             endif
             blayk = -1
             ! also include sub- and super-harmonic transfer functions
             if ( istok == 2 ) istok = 3
             !
             call INREAL ('H1'  , h1, 'REQ', 0.)
             call INREAL ('H2'  , h2, 'REQ', 0.)
             if ( .not. h1 > 0. .or. .not. h2 > 0. ) call msgerr (2, 'wave height is less or equal to zero')
             call INREAL ('PER1', t1, 'REQ', 0.)
             call INREAL ('PER2', t2, 'REQ', 0.)
             if ( .not. t1 > 0. .or. .not. t2 > 0. ) call msgerr (2, 'wave period is less or equal to zero')
             call INREAL ('DIR1', dir1, 'STA', -999.)
             if ( oned .and. dir1 /= -999. ) then
                call msgerr (1, 'first incident wave angle normal to boundary in case of 1D computation')
                dir1 = -999.
             endif
             call INREAL ('DIR2', dir2, 'STA', -999.)
             if ( oned .and. dir2 /= -999. ) then
                call msgerr (1, 'second incident wave angle normal to boundary in case of 1D computation')
                dir2 = -999.
             endif
             !
             ! determine incident wave directions with respect to problem coordinates
             !
             if ( dir1 /= -999. ) then
                dir1 = DEGCNV(dir1)
                dir1 = pi2 * (dir1/360. - nint(dir1/360.))
             endif
             if ( dir2 /= -999. ) then
                dir2 = DEGCNV(dir2)
                dir2 = pi2 * (dir2/360. - nint(dir2/360.))
             endif
             !
             nbval = nbval + 1
             nbvl  = nbval
             allocate(bfstmp)
             bfstmp%nbfs      = nbval
             bfstmp%nfreq     = 2
             bfstmp%azero     = -1.e10
             bfstmp%spparm(3) = -888.   ! to indicate bichromatic wave
             !
             allocate(bfstmp%ampl (2))
             allocate(bfstmp%omega(2))
             allocate(bfstmp%phase(2))
             allocate(bfstmp%theta(2))
             !
             bfstmp%ampl (1) = h1/2.
             bfstmp%omega(1) = pi2/t1
             bfstmp%phase(1) = 0.
             bfstmp%theta(1) = dir1
             !
             bfstmp%ampl (2) = h2/2.
             bfstmp%omega(2) = pi2/t2
             bfstmp%phase(2) = 0.
             bfstmp%theta(2) = dir2
             !
             nullify(bfstmp%nextbfs)
             if ( .not.lbfs ) then
                fbfs = bfstmp
                cubfs => fbfs
                lbfs = .true.
             else
                cubfs%nextbfs => bfstmp
                cubfs => bfstmp
             endif
             !
          else if ( KEYWIS('SPECT') ) then
             !
             if ( btype /= 3 .and. btype /= 7 ) call msgerr (3, 'spectrum not allowed for this boundary type')
             if ( .not.oned ) btype = -btype
             !
             if ( blayk /= -1 .and. blayk /= 0 ) then
                call msgerr (1, 'only vertical velocity profile derived from Stokes wave theory is allowed in case of spectrum')
             endif
             blayk = -1
             ! also include sub- and super-harmonic transfer functions
             if ( istok == 2 ) istok = 3
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
             nbval = nbval + 1
             nbvl  = nbval
             allocate(bfstmp)
             bfstmp%nbfs        = nbval
             bfstmp%nfreq       = -1
             bfstmp%azero       = -1.e10
             bfstmp%spparm(1:5) = spparm(1:5)
             nullify(bfstmp%nextbfs)
             if ( .not.lbfs ) then
                fbfs = bfstmp
                cubfs => fbfs
                lbfs = .true.
             else
                cubfs%nextbfs => bfstmp
                cubfs => bfstmp
             endif
             !
          else if ( KEYWIS('FILE') .or. KEYWIS('SERI') ) then
             !
             if ( blayk == -1 ) call msgerr (3, 'time series not allowed for hyperbolic cosine distribution')
             !
             call INCSTR ('FNAME', FILENM, 'REQ', ' ')
             ! generate new set of file data
             nbfils = nbfils + 1
             nbvl   = nbval
             allocate(bfltmp)
             !
             bfltmp%bfiled = 0
             bfltmp%bctime = -999999999.
             !
             ndsd = 0
             call FOR (ndsd, FILENM, 'OF', 0)
             if (STPNOW()) return
             !
             call INKEYW ('STA', ' ')
             ! read time coding option
             call ININTG ('ITMOPT', iiopt, 'STA', 0)
             !
             if ( nstatm == 0 ) call msgerr (3, 'time information not allowed in stationary mode')
             nstatm = 1
             !
             allocate(bfltmp%bcloc(1))
             bfltmp%bcloc(1) = nbval + 1
             nbval = nbval + 1
             !
             ! store file reading parameters in array bfiled
             !
             bfltmp%bfiled(1)  = 1
             bfltmp%bfiled(4)  = 0
             bfltmp%bfiled(5)  = ndsd
             bfltmp%bfiled(6)  = iiopt
             bfltmp%bfiled(7)  = 0
             bfltmp%bfiled(8)  = 1
             ! ordering of data on file
             bfltmp%bfiled(13) = 0
             bfltmp%bfiled(14) = 0
             bfltmp%bfiled(15) = 0
             bfltmp%bfiled(16) = 0
             !
             nullify(bfltmp%nextbfl)
             if ( .not.lbfils ) then
                fbndfil = bfltmp
                cubfl => fbndfil
                lbfils = .true.
             else
                cubfl%nextbfl => bfltmp
                cubfl => bfltmp
             endif
             nbvl = nbvl + 1
             !
          else if ( KEYWIS('SPECF') ) then
             !
             if ( btype /= 3 .and. btype /= 7 ) call msgerr (3, 'spectrum not allowed for this boundary type')
             if ( .not.oned ) btype = -btype
             !
             if ( blayk /= -1 .and. blayk /= 0 ) then
                call msgerr (1, 'only vertical velocity profile derived from Stokes wave theory is allowed in case of spectrum')
             endif
             blayk = -1
             ! also include sub- and super-harmonic transfer functions
             if ( istok == 2 ) istok = 3
             !
             call INCSTR ('FNAME', FILENM, 'REQ', ' ')
             !
             call INKEYW ('STA', ' ')
             ! read cyclic period of time record
             call ININTV ('CYCLE', tcycl, 'STA', 1.e10)
             if ( .not. tcycl >  0. ) call msgerr (2, 'cyclic period of time record is less or equal to zero')
             !
             if (.not.allocated(bpix)) allocate(bpix(0))
             if (.not.allocated(bpiy)) allocate(bpiy(0))
             call SwashBCspecfile ( FILENM, nbps, bpix, bpiy, xcgrid, ycgrid, kgrpnt, tcycl, .false., btype, blayk, tsmo )
             if (STPNOW()) return
             nbvl = nbval
             !
          endif
          !
          do i = 1, nbps
             ix = curr%ix
             if ( optg /= 5 ) iy = curr%iy
             curr => curr%nextxy
             allocate(bgptmp)
             if ( optg /= 5 ) then
                bgptmp%bgp(1) = kgrpnt(ix,iy)
             else
                bgptmp%bgp(1) = ix
             endif
             bgptmp%bgp(2) = btype
             bgptmp%bgp(3) = 1000
             bgptmp%bgp(4) = nbvl
             bgptmp%bgp(5) = 0
             bgptmp%bgp(6) = 1
             bgptmp%bgp(7) = ntsmo
             if ( blayk == -1 ) then
                if ( istok == 0 ) then
                   blayk = -1
                elseif ( istok == 1 ) then
                   blayk = -11
                elseif ( istok == 2 ) then
                   blayk = -21
                elseif ( istok == 3 ) then
                   blayk = -31
                endif
             endif
             bgptmp%bgp(8) = blayk
             if ( .not.lriem ) then
                bgptmp%bgp(9) = ibound
             else
                bgptmp%bgp(9) = 2
             endif
             nullify(bgptmp%nextbgp)
             if ( .not.lbgp ) then
                fbgp = bgptmp
                cubgp => fbgp
                lbgp = .true.
             else
                cubgp%nextbgp => bgptmp
                cubgp => bgptmp
             endif
          enddo
          !
          nbgrpt = nbgrpt + nbps
          !
       else if ( KEYWIS('VAR') ) then
          !
          call INKEYW('STA', 'FOUR')
          if ( KEYWIS('FOUR') ) then
             bcimp = 1
          else if ( KEYWIS('REG') .or. KEYWIS('MONO') ) then
             bcimp = 2
          else if ( KEYWIS('BICH') ) then
             bcimp = 3
          else if ( KEYWIS('SPECT') ) then
             bcimp = 4
          else if ( KEYWIS('FILE') .or. KEYWIS('SERI') ) then
             bcimp = 5
          else if ( KEYWIS('SPECF') ) then
             bcimp = 6
          else if ( KEYWIS('SPECS') ) then
             bcimp = 7
          endif
          !
          if ( bcimp == 2 .or. bcimp == 3 .or. bcimp == 4 .or. bcimp == 6 .or. bcimp == 7 ) then
             if ( btype /= 3 .and. btype /= 7 ) then
                if ( bcimp == 2 ) then
                   call msgerr (3, 'regular wave not allowed for this boundary type')
                else if ( bcimp == 3 ) then
                   call msgerr (3, 'bichromatic wave not allowed for this boundary type')
                else
                   call msgerr (3, 'spectrum not allowed for this boundary type')
                endif
             endif
             if ( .not.oned ) btype = -btype
             if ( blayk /= -1 .and. blayk /= 0 ) then
                if ( bcimp == 2 ) then
                   call msgerr (1, 'only vertical velocity profile derived from Stokes wave theory is allowed in case of regular wave')
                else if ( bcimp == 3 ) then
                   call msgerr (1, 'only vertical velocity profile derived from Stokes wave theory is allowed in case of bichromatic wave')
                else
                   call msgerr (1, 'only vertical velocity profile derived from Stokes wave theory is allowed in case of spectrum')
                endif
             endif
             blayk = -1
             ! also include sub- and super-harmonic transfer functions for irregular waves
             if ( bcimp /= 2 .and. istok == 2 ) istok = 3
          endif
          if ( bcimp == 5 .and. blayk == -1 ) call msgerr (3, 'time series not allowed for hyperbolic cosine distribution')
          !
          if ( bcimp /= 7 ) then
             !
             rlen1  = -1.e20
             k      = 1
             rdist  = 0.
             nbv1   = 1
             lfrst1 = .true.
             lfrst2 = .true.
             lfrst3 = .true.
             !
             if (.not.allocated(bpix)) allocate(bpix(0))
             if (.not.allocated(bpiy)) allocate(bpiy(0))
             !
             do
                if ( lfrst1 ) then
                   call INREAL('LEN', rlen2, 'REQ', 0.)
                   lfrst1 = .false.
                else
                   call INREAL('LEN', rlen2, 'STA', 1.e20)
                endif
                !
                if ( rlen2 < 0.9e20 ) then
                   !
                   if ( k > nbps ) then
                      call msgerr (1, 'length of segment short, boundary values ignored')
                      write (PRINTF, 205) rdist, rlen2
                   endif
                   !
                   if ( bcimp == 1 ) then
                      !
                      call INREAL ('AZE', azero, 'REQ', 0.)
                      if ( blayk == -1 ) azero = -1.e10
                      !
                      nfreq       = 0
                      frstf%ampl  = 0.
                      frstf%omega = 0.
                      frstf%phase = 0.
                      nullify(frstf%nextfs)
                      currf => frstf
                      do
                         call INREAL ('AMP', ampl , 'REP', -1.e10)
                         if ( ampl < -.9e10 ) goto 120
                         call INREAL ('OME', omega, 'REQ', 0.)
                         call INREAL ('PHA', phase, 'REQ', 0.)
                         phase = pi2 * (phase/360. - nint(phase/360.))
                         !
                         nfreq = nfreq + 1
                         allocate(tmpf)
                         tmpf%ampl  = ampl
                         tmpf%omega = omega
                         tmpf%phase = phase
                         nullify(tmpf%nextfs)
                         currf%nextfs => tmpf
                         currf => tmpf
                      enddo
                      !
 120                  nbval = nbval + 1
                      nbv2  = nbval
                      allocate(bfstmp)
                      bfstmp%nbfs  = nbval
                      bfstmp%nfreq = nfreq
                      bfstmp%azero = azero
                      bfstmp%spparm(3) = -999.    ! indicent angle is assumed to be normal on the boundary
                      !
                      if ( nfreq > 0 ) then
                         !
                         allocate(arra(nfreq))
                         allocate(arrf(nfreq))
                         allocate(arrp(nfreq))
                         currf => frstf%nextfs
                         do i = 1, nfreq
                            arra(i) = currf%ampl
                            arrf(i) = currf%omega
                            arrp(i) = currf%phase
                            currf => currf%nextfs
                         enddo
                         deallocate(tmpf)
                         !
                         ! initialize the order of array indices
                         !
                         allocate(io(nfreq))
                         do i = 1, nfreq
                            io(i) = i
                         enddo
                         !
                         ! in case of bound waves, sort the frequencies
                         !
                         if ( istok == 2 .or. ibound == 1 ) then
                            do i = 1, nfreq-1
                               do j = i+1, nfreq
                                  if ( arrf(i) > arrf(j) ) then
                                     rtmp    = arrf(j)
                                     itmp    = io(j)
                                     arrf(j) = arrf(i)
                                     io(j)   = io(i)
                                     arrf(i) = rtmp
                                     io(i)   = itmp
                                  endif
                               enddo
                            enddo
                         endif
                         !
                         allocate(bfstmp%ampl (nfreq))
                         allocate(bfstmp%omega(nfreq))
                         allocate(bfstmp%phase(nfreq))
                         allocate(bfstmp%theta(nfreq))
                         do i = 1, nfreq
                            bfstmp%ampl (i) = arra(io(i))
                            bfstmp%omega(i) = arrf(i)
                            bfstmp%phase(i) = arrp(io(i))
                            bfstmp%theta(i) = 0.
                         enddo
                         !
                         deallocate(arra,arrf,arrp,io)
                         !
                      endif
                      !
                      nullify(bfstmp%nextbfs)
                      if ( .not.lbfs ) then
                         fbfs = bfstmp
                         cubfs => fbfs
                         lbfs = .true.
                      else
                         cubfs%nextbfs => bfstmp
                         cubfs => bfstmp
                      endif
                      !
                   else if ( bcimp == 2 ) then
                      !
                      if ( ibound == 1 ) then
                         call msgerr (1, 'addition of bound long wave correction not allowed for regular wave')
                         ibound = 0
                      endif
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
                      nbval = nbval + 1
                      nbv2  = nbval
                      allocate(bfstmp)
                      bfstmp%nbfs      = nbval
                      bfstmp%nfreq     = 1
                      bfstmp%azero     = -1.e10
                      bfstmp%spparm(3) = wdir
                      !
                      allocate(bfstmp%ampl (1))
                      allocate(bfstmp%omega(1))
                      allocate(bfstmp%phase(1))
                      allocate(bfstmp%theta(1))
                      !
                      bfstmp%ampl (1) = height/2.
                      bfstmp%omega(1) = pi2/per
                      bfstmp%phase(1) = 0.
                      bfstmp%theta(1) = 0.
                      !
                      nullify(bfstmp%nextbfs)
                      if ( .not.lbfs ) then
                         fbfs = bfstmp
                         cubfs => fbfs
                         lbfs = .true.
                      else
                         cubfs%nextbfs => bfstmp
                         cubfs => bfstmp
                      endif
                      !
                   else if ( bcimp == 3 ) then
                      !
                      call INREAL ('H1'  , h1, 'REQ', 0.)
                      call INREAL ('H2'  , h2, 'REQ', 0.)
                      if ( .not. h1 > 0. .or.  .not. h2 > 0. ) call msgerr (2, 'wave height is less or equal to zero')
                      call INREAL ('PER1', t1, 'REQ', 0.)
                      call INREAL ('PER2', t2, 'REQ', 0.)
                      if ( .not. t1 > 0. .or.  .not. t2 > 0. ) call msgerr (2, 'wave period is less or equal to zero')
                      call INREAL ('DIR1', dir1, 'STA', -999.)
                      if ( oned .and. dir1 /= -999. ) then
                         call msgerr (1, 'first incident wave angle normal to boundary in case of 1D computation')
                         dir1 = -999.
                      endif
                      call INREAL ('DIR2', dir2, 'STA', -999.)
                      if ( oned .and. dir2 /= -999. ) then
                         call msgerr (1, 'second incident wave angle normal to boundary in case of 1D computation')
                         dir2 = -999.
                      endif
                      !
                      ! determine incident wave directions with respect to problem coordinates
                      !
                      if ( dir1 /= -999. ) then
                         dir1 = DEGCNV(dir1)
                         dir1 = pi2 * (dir1/360. - nint(dir1/360.))
                      endif
                      if ( dir2 /= -999. ) then
                         dir2 = DEGCNV(dir2)
                         dir2 = pi2 * (dir2/360. - nint(dir2/360.))
                      endif
                      !
                      nbval = nbval + 1
                      nbv2  = nbval
                      allocate(bfstmp)
                      bfstmp%nbfs      = nbval
                      bfstmp%nfreq     = 2
                      bfstmp%azero     = -1.e10
                      bfstmp%spparm(3) = -888.   ! to indicate bichromatic wave
                      !
                      allocate(bfstmp%ampl (2))
                      allocate(bfstmp%omega(2))
                      allocate(bfstmp%phase(2))
                      allocate(bfstmp%theta(2))
                      !
                      bfstmp%ampl (1) = h1/2.
                      bfstmp%omega(1) = pi2/t1
                      bfstmp%phase(1) = 0.
                      bfstmp%theta(1) = dir1
                      !
                      bfstmp%ampl (2) = h2/2.
                      bfstmp%omega(2) = pi2/t2
                      bfstmp%phase(2) = 0.
                      bfstmp%theta(2) = dir2
                      !
                      nullify(bfstmp%nextbfs)
                      if ( .not.lbfs ) then
                         fbfs = bfstmp
                         cubfs => fbfs
                         lbfs = .true.
                      else
                         cubfs%nextbfs => bfstmp
                         cubfs => bfstmp
                      endif
                      !
                   else if ( bcimp == 4 ) then
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
                      nbval = nbval + 1
                      nbv2  = nbval
                      allocate(bfstmp)
                      bfstmp%nbfs        = nbval
                      bfstmp%nfreq       = -1
                      bfstmp%azero       = -1.e10
                      bfstmp%spparm(1:5) = spparm(1:5)
                      nullify(bfstmp%nextbfs)
                      if ( .not.lbfs ) then
                         fbfs = bfstmp
                         cubfs => fbfs
                         lbfs = .true.
                      else
                         cubfs%nextbfs => bfstmp
                         cubfs => bfstmp
                      endif
                      !
                   else if ( bcimp == 5 ) then
                      !
                      if ( lfrst2 ) then
                         call INCSTR ('FNAME', FILENM, 'REQ', ' ')
                         lfrst2 = .false.
                      else
                         call INCSTR ('FNAME', FILENM, 'STA', ' ')
                      endif
                      !
                      if ( FILENM /= '    ' ) then
                         ! generate new set of file data
                         nbfils = nbfils + 1
                         nbv2   = nbval
                         allocate(bfltmp)
                         !
                         bfltmp%bfiled = 0
                         bfltmp%bctime = -999999999.
                         !
                         ndsd = 0
                         call FOR (ndsd, FILENM, 'OF', 0)
                         if (STPNOW()) return
                         !
                         call INKEYW ('STA', ' ')
                         ! read time coding option
                         call ININTG ('ITMOPT', iiopt, 'STA', 0)
                         !
                         if ( nstatm == 0 ) call msgerr (3, 'time information not allowed in stationary mode')
                         nstatm = 1
                         !
                         allocate(bfltmp%bcloc(1))
                         bfltmp%bcloc(1) = nbval + 1
                         nbval = nbval + 1
                         !
                         ! store file reading parameters in array bfiled
                         !
                         bfltmp%bfiled(1)  = 1
                         bfltmp%bfiled(4)  = 0
                         bfltmp%bfiled(5)  = ndsd
                         bfltmp%bfiled(6)  = iiopt
                         bfltmp%bfiled(7)  = 0
                         bfltmp%bfiled(8)  = 1
                         ! ordering of data on file
                         bfltmp%bfiled(13) = 0
                         bfltmp%bfiled(14) = 0
                         bfltmp%bfiled(15) = 0
                         bfltmp%bfiled(16) = 0
                         !
                         nullify(bfltmp%nextbfl)
                         if ( .not.lbfils ) then
                            fbndfil = bfltmp
                            cubfl => fbndfil
                            lbfils = .true.
                         else
                            cubfl%nextbfl => bfltmp
                            cubfl => bfltmp
                         endif
                      endif
                      nbv2 = nbv2 + 1
                      !
                   else if ( bcimp == 6 ) then
                      !
                      if ( lfrst3 ) then
                         call INCSTR ('FNAME', FILENM, 'REQ', ' ')
                         lfrst3 = .false.
                      else
                         call INCSTR ('FNAME', FILENM, 'STA', ' ')
                      endif
                      !
                      if ( FILENM /= '    ' ) then
                         !
                         call INKEYW ('STA', ' ')
                         ! read cyclic period of time record
                         call ININTV ('CYCLE', tcycl, 'STA', 1.e10)
                         if ( .not. tcycl >  0. ) call msgerr (2, 'cyclic period of time record is less or equal to zero')
                         !
                         call SwashBCspecfile ( FILENM, nbps, bpix, bpiy, xcgrid, ycgrid, kgrpnt, tcycl, .false., btype, blayk, tsmo )
                         if (STPNOW()) return
                         nbv2 = nbval
                         !
                      endif
                      !
                   endif
                   !
                else
                   if ( k > nbps ) goto 140
                endif
                !
                lfrst4 = .true.
                do
                   ix = curr%ix
                   if ( optg /= 5 ) then
                      iy = curr%iy
                      x2 = xcgrid(ix,iy)
                      y2 = ycgrid(ix,iy)
                   else
                      x2 = xcugrd(ix)
                      y2 = ycugrd(ix)
                   endif
                   if ( .not.lfrst4 ) then
                      rdist = rdist + sqrt( (x2-x1)**2 + (y2-y1)**2 )
                   endif
                   lfrst4 = .false.
                   x1 = x2
                   y1 = y2
                   if ( rdist > rlen2 ) goto 130
                   fac = ( rlen2 - rdist )/( rlen2 - rlen1 )
                   !
                   allocate(bgptmp)
                   if ( optg /= 5 ) then
                      bgptmp%bgp(1) = kgrpnt(ix,iy)
                   else
                      bgptmp%bgp(1) = ix
                   endif
                   bgptmp%bgp(2) = btype
                   bgptmp%bgp(3) = nint(1000.*fac)
                   bgptmp%bgp(4) = nbv1
                   bgptmp%bgp(5) = nint(1000.*(1.-fac))
                   bgptmp%bgp(6) = nbv2
                   bgptmp%bgp(7) = ntsmo
                   if ( blayk == -1 ) then
                      if ( istok == 0 ) then
                         blayk = -1
                      elseif ( istok == 1 ) then
                         blayk = -11
                      elseif ( istok == 2 ) then
                         blayk = -21
                      elseif ( istok == 3 ) then
                         blayk = -31
                      endif
                   endif
                   bgptmp%bgp(8) = blayk
                   if ( .not.lriem ) then
                      bgptmp%bgp(9) = ibound
                   else
                      bgptmp%bgp(9) = 2
                   endif
                   nullify(bgptmp%nextbgp)
                   if ( .not.lbgp ) then
                      fbgp = bgptmp
                      cubgp => fbgp
                      lbgp = .true.
                   else
                      cubgp%nextbgp => bgptmp
                      cubgp => bgptmp
                   endif
                   !
                   k = k + 1
                   if ( k > nbps ) goto 130
                   if ( .not.associated(curr%nextxy) ) exit
                   curr => curr%nextxy
                enddo
                !
                ! boundary values have been assigned, read new parameters
                !
 130            if ( rlen2 > 0.9e20 ) goto 140
                !
                rlen1 = rlen2
                nbv1  = nbv2
                !
             enddo
             !
          else
             !
             call INCSTR ('FNAME', FILENM, 'REQ', ' ')
             !
             call INKEYW ('STA', ' ')
             ! read cyclic period of time record
             call ININTV ('CYCLE', tcycl, 'STA', 1.e10)
             if ( .not. tcycl >  0. ) call msgerr (2, 'cyclic period of time record is less or equal to zero')
             !
             ! store indices of boundary points
             !
             if ( allocated(bpix) ) deallocate(bpix)
             allocate(bpix(nbps))
             if ( allocated(bpiy) ) deallocate(bpiy)
             allocate(bpiy(nbps))
             !
             do i = 1, nbps
                bpix(i) = curr%ix
                if ( optg /= 5 ) then
                   bpiy(i) = curr%iy
                else
                   bpiy(i) = 0
                endif
                curr => curr%nextxy
             enddo
             !
             if ( blayk == -1 ) then
                if ( istok == 0 ) then
                   blayk = -1
                elseif ( istok == 1 ) then
                   blayk = -11
                elseif ( istok == 2 ) then
                   blayk = -21
                elseif ( istok == 3 ) then
                   blayk = -31
                endif
             endif
             !
             call SwashBCspecfile ( FILENM, nbps, bpix, bpiy, xcgrid, ycgrid, kgrpnt, tcycl, .true., btype, blayk, tsmo )
             if (STPNOW()) return
             !
          endif
          !
 140      nbgrpt = nbgrpt + nbps
          !
       else
          call WRNKEY
       endif
       !
       if ( .not.oned .and. btype == -2 .or. btype == -6 .or. btype == -8 .or. btype == -10 ) write (PRINTF,206) btype
       if ( oned .and. btype < 0 ) write (PRINTF,207) btype
       !
       if (associated(tmp)) deallocate(tmp)
       if (allocated(bpix)) deallocate(bpix)
       if (allocated(bpiy)) deallocate(bpiy)
       !
    endif
    !
 201 format (' shape of spectrum, height:', i2, ' ; freq:', i2, ' ; dir:', i2)
 202 format (' segment point ', 2f10.2, ' grid ', 2f8.2, 2i4)
 203 format (' side ', 2i4, 2(1x,e11.4), 2(1x,f5.0), 2x, f6.3, 2x, l1)
 204 format (' selected side:', i2, ' from ', 2i4, 2f9.0, ' to ', 2i4, 2f9.0)
 205 format (' segment length=', f9.2, '; [len]=', f9.2)
 206 format (' Internal error: wrong sign for btype = ',i3)
 207 format (' Internal error: invalid btype found in 1D mode. Value of btype = ',i3)
    !
end subroutine SwashBounCond
