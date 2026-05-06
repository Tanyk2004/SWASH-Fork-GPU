subroutine SwashUpdPress
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
!    1.00, May 2013: New subroutine
!
!   Purpose
!
!   Initializes / updates atmospheric pressure based on space varying input field
!   Also correct water level on open boundaries, if appropriate
!
!   Modules used
!
    use ocpcomm4
    use SwashCommdata2
    use SwashCommdata3
    use m_genarr
    use m_parall
    use SwashFlowdata
!
    implicit none
!
!   Local variables
!
    integer, save :: ient = 0 ! number of entries in this subroutine
    integer       :: m        ! loop counter
    integer       :: md       ! index of point m-1
    integer       :: n        ! loop counter
    integer       :: nd       ! index of point n-1
    integer       :: ndm      ! pointer to m,n-1
    integer       :: ndmd     ! pointer to m-1,n-1
    integer       :: nfm      ! pointer to m,nf
    integer       :: nfum     ! pointer to m,nfu
    integer       :: nlm      ! pointer to m,nl
    integer       :: nlum     ! pointer to m,nlu
    integer       :: nm       ! pointer to m,n
    integer       :: nmd      ! pointer to m-1,n
    integer       :: nmf      ! pointer to mf,n
    integer       :: nmfu     ! pointer to mfu,n
    integer       :: nml      ! pointer to ml,n
    integer       :: nmlu     ! pointer to mlu,n
    !
    logical       :: STPNOW   ! indicates that program must stop
!
!   Structure
!
!   Description of the pseudo code
!
!   Source text
!
    if (ltrace) call strace (ient,'SwashUpdPress')
    !
    if ( ifldyn(13) == 1 ) then
       !
       ! determine pressure in wl-points
       !
       if ( oned ) then
          !
          ! loop over wl-points in 1D computational grid
          !
          do m = mfu, ml
             !
             md = m - 1
             !
             nm  = kgrpnt(m ,1)
             nmd = kgrpnt(md,1)
             !
             patm(nm) = 0.5 * ( presf(nm,2) + presf(nmd,2) )
             !
          enddo
          !
          ! exchange pressure with neighbouring subdomains
          !
          call SWEXCHG ( patm, kgrpnt, 1, 1 )
          if (STPNOW()) return
          !
          ! at real boundaries, copy them from internal points
          !
          nmf  = kgrpnt(mf ,1)
          nmfu = kgrpnt(mfu,1)
          nml  = kgrpnt(ml ,1)
          nmlu = kgrpnt(mlu,1)
          !
          if ( LMXF ) patm(nmf ) = patm(nmfu)
          if ( LMXL ) patm(nmlu) = patm(nml )
          !
       else
          !
          ! loop over wl-points in 2D computational grid
          !
          do n = nfu, nl
             !
             do m = mfu, ml
                !
                md = m - 1
                nd = n - 1
                !
                nm   = kgrpnt(m ,n )
                nmd  = kgrpnt(md,n )
                ndm  = kgrpnt(m ,nd)
                ndmd = kgrpnt(md,nd)
                !
                ! for permanently dry neighbours, corresponding values will be mirrored
                !
                if ( nmd  == 1 ) nmd  = nm
                if ( ndm  == 1 ) ndm  = nm
                if ( ndmd == 1 ) ndmd = ndm
                !
                patm(nm) = 0.25 * ( presf(nm,2) + presf(nmd,2) + presf(ndm,2) + presf(ndmd,2) )
                !
             enddo
             !
          enddo
          !
          ! set to zero for permanently dry points
          !
          patm(1) = 0.
          !
          ! exchange pressure with neighbouring subdomains
          !
          call SWEXCHG ( patm, kgrpnt, 1, 1 )
          if (STPNOW()) return
          !
          ! at real boundaries, copy them from internal points or apply periodicity
          !
          if ( .not.lreptx ) then
             !
             do n = nfu, nl
                !
                nmf  = kgrpnt(mf ,n)
                nmfu = kgrpnt(mfu,n)
                nml  = kgrpnt(ml ,n)
                nmlu = kgrpnt(mlu,n)
                !
                if ( LMXF ) then
                   !
                   patm(nmf) = patm(nmfu)
                   !
                endif
                !
                if ( LMXL ) then
                   !
                   patm(nmlu) = patm(nml)
                   !
                endif
                !
             enddo
             !
          else
             !
             call periodic ( patm, kgrpnt, 1, 1 )
             !
          endif
          !
          if ( .not.lrepty ) then
             !
             do m = mfu, ml
                !
                nfm  = kgrpnt(m,nf )
                nfum = kgrpnt(m,nfu)
                nlm  = kgrpnt(m,nl )
                nlum = kgrpnt(m,nlu)
                !
                if ( LMYF ) then
                   !
                   patm(nfm) = patm(nfum)
                   !
                endif
                !
                if ( LMYL ) then
                   !
                   patm(nlum) = patm(nlm)
                   !
                endif
                !
             enddo
             !
          else
             !
             call periodic ( patm, kgrpnt, 1, 1 )
             !
          endif
          !
       endif
       !
    endif
    !
    ! correct water level on boundary for local atmospheric pressure, if appropriate
    !
    if ( prmean > 0. ) then
       !
       if ( oned ) then
          !
          if ( ibl(1) == 2 ) s1(nmf ) = s1(nmf ) + ( prmean - patm(nmf ) ) / (rhow*grav)
          if ( ibr(1) == 2 ) s1(nmlu) = s1(nmlu) + ( prmean - patm(nmlu) ) / (rhow*grav)
          !
       else
          !
          if ( .not.lreptx ) then
             !
             do n = nfu, nl
                !
                nmf  = kgrpnt(mf ,n)
                nmlu = kgrpnt(mlu,n)
                !
                if ( ibl(n) == 2 ) s1(nmf ) = s1(nmf ) + ( prmean - patm(nmf ) ) / (rhow*grav)
                if ( ibr(n) == 2 ) s1(nmlu) = s1(nmlu) + ( prmean - patm(nmlu) ) / (rhow*grav)
                !
             enddo
             !
          endif
          !
          if ( .not.lrepty ) then
             !
             do m = mfu, ml
                !
                nfm  = kgrpnt(m,nf )
                nlum = kgrpnt(m,nlu)
                !
                if ( ibb(m) == 2 ) s1(nfm ) = s1(nfm ) + ( prmean - patm(nfm ) ) / (rhow*grav)
                if ( ibt(m) == 2 ) s1(nlum) = s1(nlum) + ( prmean - patm(nlum) ) / (rhow*grav)
                !
             enddo
             !
          endif
          !
       endif
       !
    endif
    !
end subroutine SwashUpdPress
