subroutine SwashOutput
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
!    1.00, March 2010: New subroutine
!
!   Purpose
!
!   Processes output requests
!
!   Method
!
!   Data structure
!
!   Output requests are encoded in type orqdat. These are set by commands TABLE, BLOCK, etc. (see subroutine SwashReqOutQ).
!   Each output request refers to one set of output locations, and to one or more output quantities.
!
!   Data on output locations are encoded in type opsdat. These are set by commands FRAME, POINTS, CURVE, etc. (see subroutine SwashReqOutL).
!   Each set is characterized by its name (psname in the code) pstype is the type of set (i.e. 'F' for Frame etc.)
!
!   Properties of output quantities are in arrays OVSNAM, OVLNAM, OVUNIT, OVSVTY etc. These are set in subroutine SwashInit.
!   Each output quantity is assigned to a fixed number; i.e. 1=Xp, 2=Yp, 6=water level, 7=velocity, etc.
!   Subroutine SVARTP determines the above number from the name of the quantity as it appears in the user command; this is compared with OVKEYW.
!
!   Modules used
!
    use ocpcomm4
    use SwashCommdata1
    use SwashCommdata3, only: ifloat, kmax
    use SwashCommdata4
    use outp_data
    use m_genarr, only: kgrpnt, xcgrid, ycgrid
    use m_parall
!
    implicit none
!
!   Local variables
!
    integer                              :: i        ! loop counter
    integer, save                        :: ient = 0 ! number of entries in this subroutine
    integer, dimension(:), allocatable   :: ionod    ! indicates in which subdomain output points are located
    integer                              :: j        ! loop counter
    integer                              :: mip      ! actual number of output locations
    integer                              :: mxk      ! number of output points in x-direction
    integer                              :: myk      ! number of output points in y-direction
    integer                              :: n        ! loop counter
    integer                              :: nvar     ! number of output quantities
    integer                              :: nvoqk    ! actual number of layer-dependent output quantities
    integer                              :: nvoqp    ! actual number of output quantities
    integer, dimension(nmovar)           :: voqr     ! indicates where an output quantity is stored in array voq
    !
    real, dimension(:), allocatable      :: voq      ! contains output quantities for all output points
    real, dimension(:), allocatable      :: voqk     ! contains layer-dependent output quantities for all layers and output points
    !
    logical                              :: logact   ! indicates whether output action should be taken or not
    logical, dimension(nmovar)           :: oqproc   ! indicates whether or not an output quantity must be processed
    logical                              :: STPNOW   ! indicates that program must stop
    !
    character(8)                         :: psname   ! name of output point set
    character(1)                         :: pstype   ! type of output point set
    character(4)                         :: rtype    ! type of output request
    !
    type(opsdat), pointer                :: cuops    ! current item in list of point set
    type(orqdat), pointer                :: corq     ! current item in list of request output
!
!   Structure
!
!   After the coordinates of all output locations have been determined,
!   values of all output quantities are calculated, and written into
!   2d array voq (one or two columns for each output quantity, one line
!   for each location). Array voqr shows which quantity is written in
!   each column.
!
!   After array voq is filled, the actual output starts; which subroutine
!   is called depends on rtype, i.e. BLOCK or TABLE.
!
!   Source text
!
    if (ltrace) call strace (ient,'SwashOutput')
    !
    ! indicate number of output requests
    !
    if ( nreoq == 0 ) then
       call msgerr (1, 'no output requested')
       return
    endif
    if ( ITEST >= 10 ) write (PRINTF, 201) nreoq
    !
    ! repeat for all output requests:
    !
    corq => forq
    reqloop: do i = 1, nreoq
       !
       nvoqp = 0
       nvoqk = 0
       !
       ! analyse output quantities
       !
       rtype  = corq%rqtype
       psname = corq%psname
       nvar   = corq%oqi(3)
       call SwashDecOutQ ( corq%oqi, corq%oqr, corq%ivtyp, rtype, psname, nvar, nvoqp, nvoqk, oqproc, voqr, ourqt(i), logact )
       if ( .not.logact ) then
          corq => corq%nextorq
          cycle reqloop
       endif
       corq%rqtype = rtype
       !
       if ( SCREEN /= PRINTF .and. INODE == MASTER ) write (SCREEN, 202) i
       if ( ITEST >= 10 ) write (PRINTF, 203) i
       !
       cuops => FOPS
       do
          if ( cuops%psname == psname ) exit
          if ( .not.associated(cuops%nextops) ) then
             call msgerr (3, 'output requested for non-existing points')
             write (PRINTF, 204) psname
             corq => corq%nextorq
             cycle reqloop
          endif
          cuops => cuops%nextops
       enddo
       !
       if ( ITEST >= 80 .or. IOUTES >= 10 ) write (PRTEST, 205) i, nvoqp, nvoqk, rtype, psname
       !
       ! analyse output point set
       !
       pstype = cuops%pstype
       mip    = cuops%mip
       call SwashDecOutL ( cuops%opi, cuops%opr, psname, pstype, mip, mxk, myk )
       !
       ! indicator for interpolating output quantities
       ! =0: interpolation; =1: no interpolation
       !
       if ( pstype=='F' .or. pstype=='H' ) then
          outpar(3) = cuops%opr(6)
       else
          outpar(3) = 0.
       endif
       !
       ! assign memory to array voq (contains output quantities for all output points)
       allocate (voq(mip*nvoqp))
       !
       ! assign memory to array voqk (contains layer-dependent output quantities for all layers and output points)
       allocate (voqk(mip*(kmax+1)*nvoqk))
       !
       ! assign memory to array ionod (indicates in which subdomain output points are located)
       allocate (ionod(mip))
       !
       if ( .not.PARLL ) then
          ionod = MASTER
       else
          ionod = -999
       endif
       !
       ! calculate grid coordinates of the output points
       !
       call SwashCoorOutp ( pstype, cuops%opi, cuops%opr, cuops%xp, cuops%yp, mip, voq(1), voq(1+mip), voq(1+2*mip), voq(1+3*mip), kgrpnt, xcgrid, ycgrid )
       !
       ! calculate requested output quantities ...
       !
       if ( psname(1:7) /= 'ELEMESH' ) then
          !
          ! ... in the output points
          !
          call SwashQuanOutp ( rtype, oqproc, mip, nvoqp, nvoqk, voq(1+2*mip), voq(1+3*mip), voqr, voq, voqk, ionod )
          !
       else
          !
          ! ... in the mesh elements
          !
          call SwashElemOutp ( oqproc, mip, nvoqp, nvoqk, voqr, voq, voqk )
          !
       endif
       !
       ! calculate requested space-independent output quantities as time series
       !
       if ( psname(1:6) == 'NOGRID' ) then
          !
          ! the hydrodynamic loads and/or body motions
          !
          if ( ifloat /= 0 ) call SwashFlobjOutp ( oqproc, nvoqp, voqr, voq )
          !
          ! the runup height
          !
          if ( oqproc(115) ) call SwashRunupHeight ( oqproc, nvoqp, voqr, voq )
          !
       endif
       !
       if ( ITEST >= 100 ) then
          write (PRTEST, 206) (voqr(n), n=1, nmovar)
          do j = 1, min(mip,20)
             write (PRTEST, 207) (voq(j+(n-1)*mip), n=1, nvoqp)
          enddo
       endif
       !
       if ( rtype(1:3) == 'BLK' ) then
          !
          ! block output
          !
          if (rtype == 'BLKV' ) then
             call SWBLKV ( corq%oqi, corq%oqr, corq%ivtyp, mxk, myk, voqr, voq, voqk, pstype, psname, ionod )
          else if ( PARLL .AND. IMRGE.EQ.1 ) then
             call SWBLKP ( corq%oqi, corq%ivtyp, mxk, myk, voqr, voq, voqk, ionod )
          else
             call SWBLOK ( rtype, corq%oqi, corq%ivtyp, corq%fac, psname, mxk, myk, i, voqr, voq, voqk )
          endif
          if (STPNOW()) return
          !
       else if ( rtype(1:3) == 'TAB' ) then
          !
          ! table output
          !
          call SWTABP ( rtype, corq%oqi, corq%ivtyp, psname, mip, voqr, voq, voqk, ionod )
          if (STPNOW()) return
          !
       else
          !
          write (PRINTF, 208) i, nvoqp, rtype, psname, mip
          !
       endif
       !
       deallocate ( voq, voqk, ionod )
       corq => corq%nextorq
       !
    enddo reqloop
    !
 201 format (1x, i3, ' output requests')
 202 format ('+SWASH is processing output request ', i4)
 203 format (' SWASH is processing output request ', i4)
 204 format (' point set: ', a)
 205 format (' test SwashOutput ', 3i6, 2x, a4, 2x, a16)
 206 format (' arrays voqr and voq:', 30i3)
 207 format (12(1x,e10.4))
 208 format (' error in output request ', 2i6, 2x, a4, 2x, a16, i6)
    !
end subroutine SwashOutput
