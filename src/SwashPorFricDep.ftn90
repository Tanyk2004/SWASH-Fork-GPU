subroutine SwashPorFricDep
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
!    1.00, November 2016: New subroutine
!
!   Purpose
!
!   Adapts frictional forces and added mass inside porous structures in depth-averaged mode
!
!   Modules used
!
    use ocpcomm4
    use SwashCommdata2
    use SwashCommdata3
    use m_genarr
    use SwashFlowdata
!
    implicit none
!
!   Local variables
!
    integer, save :: ient = 0 ! number of entries in this subroutine
    integer       :: m        ! loop counter
    integer       :: n        ! loop counter
    integer       :: md       ! index of point m-1
    integer       :: nd       ! index of point n-1
    integer       :: ndm      ! pointer to m,n-1
    integer       :: nm       ! pointer to m,n
    integer       :: nmd      ! pointer to m-1,n
    !
    real          :: dsiz     ! average size of porous material
    real          :: faca     ! multiplication factor to get proper laminar friction factor
    real          :: facb     ! multiplication factor to get proper turbulent friction factor
    real          :: facc     ! multiplication factor to get proper added mass coefficient
    real          :: kcn      ! Keulegan-Carpenter number
    real          :: npom     ! mean porosity
    real          :: wper     ! wave period
!
!   Structure
!
!   Description of the pseudo code
!
!   Source text
!
    if (ltrace) call strace (ient,'SwashPorFricDep')
    !
    dsiz = ppor(1)
    wper = ppor(5)
    !
    do nm = 1, mcgrd
       !
       apomu(nm,1) = aporu(nm)
       bpomu(nm,1) = bporu(nm)
       cpomu(nm,1) = cporu(nm)
       !
       if ( nporu(nm) > 0.1 .and. nporu(nm) < 0.99 ) then
          !
          if ( .not. hporu(nm) > hum(nm) ) then
             !
             ! breakwater is submerged, so
             !
             ! first compute mean porosity
             !
             npom = ( (nporu(nm)-1.)*hporu(nm) + hum(nm) ) / hum(nm)
             !
             ! and then adapt friction factors and added mass coefficient
             !
             faca = ( 1. - npom )**2. * nporu(nm) * nporu(nm) / ( ( 1. - nporu(nm) )**2. * npom * npom )
             facb = ( 1. - npom )     * nporu(nm) * nporu(nm) / ( ( 1. - nporu(nm) )     * npom * npom )
             facc = ( 1. - npom )     * nporu(nm)             / ( ( 1. - nporu(nm) )     * npom        )
             !
             apomu(nm,1) = faca * aporu(nm)
             bpomu(nm,1) = facb * bporu(nm)
             cpomu(nm,1) = facc * cporu(nm)
             !
          endif
          !
       endif
       !
    enddo
    !
    if ( .not.oned ) then
       !
       do nm = 1, mcgrd
          !
          apomv(nm,1) = aporv(nm)
          bpomv(nm,1) = bporv(nm)
          cpomv(nm,1) = cporv(nm)
          !
          if ( nporv(nm) > 0.1 .and. nporv(nm) < 0.99 ) then
             !
             if ( .not. hporv(nm) > hvm(nm) ) then
                !
                ! breakwater is submerged, so
                !
                ! first compute mean porosity
                !
                npom = ( (nporv(nm)-1.)*hporv(nm) + hvm(nm) ) / hvm(nm)
                !
                ! and then adapt friction factors and added mass coefficient
                !
                faca = ( 1. - npom )**2. * nporv(nm) * nporv(nm) / ( ( 1. - nporv(nm) )**2. * npom * npom )
                facb = ( 1. - npom )     * nporv(nm) * nporv(nm) / ( ( 1. - nporv(nm) )     * npom * npom )
                facc = ( 1. - npom )     * nporv(nm)             / ( ( 1. - nporv(nm) )     * npom        )
                !
                apomv(nm,1) = faca * aporv(nm)
                bpomv(nm,1) = facb * bporv(nm)
                cpomv(nm,1) = facc * cporv(nm)
                !
             endif
             !
          endif
          !
       enddo
       !
    endif
    !
    ! adapt turbulent frictional force using the Keulegan-Carpenter number, if appropriate
    ! note: flow velocity here is the pore velocity
    !
    if ( wper > 0. ) then
       !
       if ( oned ) then
          !
          if ( vargs ) then
             !
             do m = mf, ml
                !
                nm = kgrpnt(m,1)
                !
                if ( nporu(nm) > 0.1 .and. nporu(nm) < 0.99 ) then
                   !
                   dsiz = gsizf(nm)
                   !
                   if ( dsiz /= 0. ) then
                      !
                      kcn         = max(abs(u0(nm,1)),0.01) * max(1.,wper) / dsiz
                      bpomu(nm,1) = bpomu(nm,1) * ( 1. + 7.5/kcn )
                      !
                   endif
                   !
                endif
                !
             enddo
             !
          else
             !
             do m = mf, ml
                !
                nm = kgrpnt(m,1)
                !
                if ( nporu(nm) > 0.1 .and. nporu(nm) < 0.99 ) then
                   !
                   if ( dsiz /= 0. ) then
                      !
                      kcn         = max(abs(u0(nm,1)),0.01) * max(1.,wper) / dsiz
                      bpomu(nm,1) = bpomu(nm,1) * ( 1. + 7.5/kcn )
                      !
                   endif
                   !
                endif
                !
             enddo
             !
          endif
          !
       else
          !
          if ( vargs ) then
             !
             ! loop over u-points
             !
             do n = nfu, nl
                do m = mf, ml
                   !
                   nd = n - 1
                   !
                   nm  = kgrpnt(m,n )
                   ndm = kgrpnt(m,nd)
                   !
                   ! for permanently dry neighbours, corresponding values will be mirrored
                   !
                   if ( ndm == 1 ) ndm = nm
                   !
                   if ( nporu(nm) > 0.1 .and. nporu(nm) < 0.99 ) then
                      !
                      dsiz = 0.5 * ( gsizf(nm) + gsizf(ndm) )
                      !
                      if ( dsiz /= 0. ) then
                         !
                         kcn         = max(abs(u0(nm,1)),0.01) * max(1.,wper) / dsiz
                         bpomu(nm,1) = bpomu(nm,1) * ( 1. + 7.5/kcn )
                         !
                      endif
                      !
                   endif
                   !
                enddo
             enddo
             !
             ! loop over v-points
             !
             do m = mfu, ml
                do n = nf, nl
                   !
                   md = m - 1
                   !
                   nm  = kgrpnt(m ,n)
                   nmd = kgrpnt(md,n)
                   !
                   ! for permanently dry neighbours, corresponding values will be mirrored
                   !
                   if ( nmd == 1 ) nmd = nm
                   !
                   if ( nporv(nm) > 0.1 .and. nporv(nm) < 0.99 ) then
                      !
                      dsiz = 0.5 * ( gsizf(nm) + gsizf(nmd) )
                      !
                      if ( dsiz /= 0. ) then
                         !
                         kcn         = max(abs(v0(nm,1)),0.01) * max(1.,wper) / dsiz
                         bpomv(nm,1) = bpomv(nm,1) * ( 1. + 7.5/kcn )
                         !
                      endif
                      !
                   endif
                   !
                enddo
             enddo
             !
          else
             !
             ! loop over u-points
             !
             do n = nfu, nl
                do m = mf, ml
                   !
                   nm = kgrpnt(m,n)
                   !
                   if ( nporu(nm) > 0.1 .and. nporu(nm) < 0.99 ) then
                      !
                      if ( dsiz /= 0. ) then
                         !
                         kcn         = max(abs(u0(nm,1)),0.01) * max(1.,wper) / dsiz
                         bpomu(nm,1) = bpomu(nm,1) * ( 1. + 7.5/kcn )
                         !
                      endif
                      !
                   endif
                   !
                enddo
             enddo
             !
             ! loop over v-points
             !
             do m = mfu, ml
                do n = nf, nl
                   !
                   nm = kgrpnt(m,n)
                   !
                   if ( nporv(nm) > 0.1 .and. nporv(nm) < 0.99 ) then
                      !
                      if ( dsiz /= 0. ) then
                         !
                         kcn         = max(abs(v0(nm,1)),0.01) * max(1.,wper) / dsiz
                         bpomv(nm,1) = bpomv(nm,1) * ( 1. + 7.5/kcn )
                         !
                      endif
                      !
                   endif
                   !
                enddo
             enddo
             !
          endif
          !
       endif
       !
    endif
    !
end subroutine SwashPorFricDep
