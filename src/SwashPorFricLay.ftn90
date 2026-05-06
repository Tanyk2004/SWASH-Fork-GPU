subroutine SwashPorFricLay
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
!   Adapts frictional forces and added mass inside porous structures in multi-layered mode
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
    integer       :: k        ! loop counter over layers
    integer       :: m        ! loop counter
    integer       :: n        ! loop counter
    integer       :: md       ! index of point m-1
    integer       :: nd       ! index of point n-1
    integer       :: ndm      ! pointer to m,n-1
    integer       :: nm       ! pointer to m,n
    integer       :: nmd      ! pointer to m-1,n
    integer       :: ndmd     ! pointer to m-1,n-1
    !
    real          :: dsiz     ! average size of porous material
    real          :: faca     ! multiplication factor to get proper laminar friction factor
    real          :: facb     ! multiplication factor to get proper turbulent friction factor
    real          :: facc     ! multiplication factor to get proper added mass coefficient
    real          :: kcn      ! Keulegan-Carpenter number
    real          :: npok     ! mean porosity within layer
    real          :: wper     ! wave period
    real          :: zkh      ! structure height level
!
!   Structure
!
!   Description of the pseudo code
!
!   Source text
!
    if (ltrace) call strace (ient,'SwashPorFricLay')
    !
    dsiz = ppor(1)
    wper = ppor(5)
    !
    do nm = 1, mcgrd
       !
       apomu(nm,:) = aporu(nm)
       bpomu(nm,:) = bporu(nm)
       cpomu(nm,:) = cporu(nm)
       !
       if ( nporu(nm) > 0.1 .and. nporu(nm) < 0.99 ) then
          !
          ! compute structure height level relative to reference plane (positive upwards)
          !
          zkh = zkum(nm,kmax) + hporu(nm)
          !
          if ( zkh > zkum(nm,0) ) then
             !
             ! breakwater is emerged
             !
             apomu(nm,:) = aporu(nm)
             bpomu(nm,:) = bporu(nm)
             cpomu(nm,:) = cporu(nm)
             !
          else
             !
             ! breakwater is submerged, so
             !
             do k = 1, kmax
                !
                if ( zkh < zkum(nm,k) ) then
                   !
                   ! layer is fluid
                   !
                   apomu(nm,k) = 0.
                   bpomu(nm,k) = 0.
                   cpomu(nm,k) = 0.
                   !
                else if ( zkh > zkum(nm,k-1) ) then
                   !
                   ! layer is completely porous
                   !
                   apomu(nm,k) = aporu(nm)
                   bpomu(nm,k) = bporu(nm)
                   cpomu(nm,k) = cporu(nm)
                   !
                else
                   !
                   ! layer is partly porous, so
                   !
                   ! first compute mean porosity within layer
                   !
                   npok = ( nporu(nm)*( zkh - zkum(nm,k) ) + ( zkum(nm,k-1) - zkh ) ) / ( zkum(nm,k-1) - zkum(nm,k) )
                   !
                   ! and then adapt friction factors and added mass coefficient
                   !
                   faca = ( 1. - npok )**2. * nporu(nm) * nporu(nm) / ( ( 1. - nporu(nm) )**2. * npok * npok )
                   facb = ( 1. - npok )     * nporu(nm) * nporu(nm) / ( ( 1. - nporu(nm) )     * npok * npok )
                   facc = ( 1. - npok )     * nporu(nm)             / ( ( 1. - nporu(nm) )     * npok        )
                   !
                   apomu(nm,k) = faca * aporu(nm)
                   bpomu(nm,k) = facb * bporu(nm)
                   cpomu(nm,k) = facc * cporu(nm)
                   !
                endif
                !
             enddo
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
          apomv(nm,:) = aporv(nm)
          bpomv(nm,:) = bporv(nm)
          cpomv(nm,:) = cporv(nm)
          !
          if ( nporv(nm) > 0.1 .and. nporv(nm) < 0.99 ) then
             !
             ! compute structure height level relative to reference plane (positive upwards)
             !
             zkh = zkvm(nm,kmax) + hporv(nm)
             !
             if ( zkh > zkvm(nm,0) ) then
                !
                ! breakwater is emerged
                !
                apomv(nm,:) = aporv(nm)
                bpomv(nm,:) = bporv(nm)
                cpomv(nm,:) = cporv(nm)
                !
             else
                !
                ! breakwater is submerged, so
                !
                do k = 1, kmax
                   !
                   if ( zkh < zkvm(nm,k) ) then
                      !
                      ! layer is fluid
                      !
                      apomv(nm,k) = 0.
                      bpomv(nm,k) = 0.
                      cpomv(nm,k) = 0.
                      !
                   else if ( zkh > zkvm(nm,k-1) ) then
                      !
                      ! layer is completely porous
                      !
                      apomv(nm,k) = aporv(nm)
                      bpomv(nm,k) = bporv(nm)
                      cpomv(nm,k) = cporv(nm)
                      !
                   else
                      !
                      ! layer is partly porous, so
                      !
                      ! first compute mean porosity within layer
                      !
                      npok = ( nporv(nm)*( zkh - zkvm(nm,k) ) + ( zkvm(nm,k-1) - zkh ) ) / ( zkvm(nm,k-1) - zkvm(nm,k) )
                      !
                      ! and then adapt friction factors and added mass coefficient
                      !
                      faca = ( 1. - npok )**2. * nporv(nm) * nporv(nm) / ( ( 1. - nporv(nm) )**2. * npok * npok )
                      facb = ( 1. - npok )     * nporv(nm) * nporv(nm) / ( ( 1. - nporv(nm) )     * npok * npok )
                      facc = ( 1. - npok )     * nporv(nm)             / ( ( 1. - nporv(nm) )     * npok        )
                      !
                      apomv(nm,k) = faca * aporv(nm)
                      bpomv(nm,k) = facb * bporv(nm)
                      cpomv(nm,k) = facc * cporv(nm)
                      !
                   endif
                   !
                enddo
                !
             endif
             !
          endif
          !
       enddo
       !
    endif
    !
    if ( ihydro == 1 .or. ihydro == 2 ) then
       !
       do nm = 1, mcgrd
          !
          apoks(nm,:) = apors(nm)
          bpoks(nm,:) = bpors(nm)
          cpoks(nm,:) = cpors(nm)
          npoks(nm,:) = npors(nm)
          !
          if ( npors(nm) > 0.1 .and. npors(nm) < 0.99 ) then
             !
             ! compute structure height level relative to reference plane (positive upwards)
             !
             zkh = zks(nm,kmax) + hpors(nm)
             !
             do k = 0, kmax-1
                !
                if ( zkh < zks(nm,k) ) then
                   !
                   ! layer interface is fluid
                   !
                   apoks(nm,k) = 0.
                   bpoks(nm,k) = 0.
                   cpoks(nm,k) = 0.
                   npoks(nm,k) = 1.
                   !
                endif
                !
             enddo
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
                      do k = 1, kmax
                         kcn         = max(abs(u0(nm,k)),0.01) * max(1.,wper) / dsiz
                         bpomu(nm,k) = bpomu(nm,k) * ( 1. + 7.5/kcn )
                      enddo
                      !
                   endif
                   !
                endif
                !
             enddo
             !
             if ( ihydro == 1 .or. ihydro == 2 ) then
                !
                ! loop over wl-points
                !
                do m = mfu, ml
                   !
                   md = m - 1
                   !
                   nm  = kgrpnt(m ,1)
                   nmd = kgrpnt(md,1)
                   !
                   if ( npors(nm) > 0.1 .and. npors(nm) < 0.99 ) then
                      !
                      dsiz = 0.5 * ( gsizf(nm) + gsizf(nmd) )
                      !
                      if ( dsiz /= 0. ) then
                         !
                         do k = 0, kmax-1
                            kcn         = max(abs(w0(nm,k)),0.01) * max(1.,wper) / dsiz
                            bpoks(nm,k) = bpoks(nm,k) * ( 1. + 7.5/kcn )
                         enddo
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
             do m = mf, ml
                !
                nm = kgrpnt(m,1)
                !
                if ( nporu(nm) > 0.1 .and. nporu(nm) < 0.99 ) then
                   !
                   if ( dsiz /= 0. ) then
                      !
                      do k = 1, kmax
                         kcn         = max(abs(u0(nm,k)),0.01) * max(1.,wper) / dsiz
                         bpomu(nm,k) = bpomu(nm,k) * ( 1. + 7.5/kcn )
                      enddo
                      !
                   endif
                   !
                endif
                !
             enddo
             !
             if ( ihydro == 1 .or. ihydro == 2 ) then
                !
                ! loop over wl-points
                !
                do m = mfu, ml
                   !
                   nm = kgrpnt(m,1)
                   !
                   if ( npors(nm) > 0.1 .and. npors(nm) < 0.99 ) then
                      !
                      if ( dsiz /= 0. ) then
                         !
                         do k = 0, kmax-1
                            kcn         = max(abs(w0(nm,k)),0.01) * max(1.,wper) / dsiz
                            bpoks(nm,k) = bpoks(nm,k) * ( 1. + 7.5/kcn )
                         enddo
                         !
                      endif
                      !
                   endif
                   !
                enddo
                !
             endif
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
                         do k = 1, kmax
                            kcn         = max(abs(u0(nm,k)),0.01) * max(1.,wper) / dsiz
                            bpomu(nm,k) = bpomu(nm,k) * ( 1. + 7.5/kcn )
                         enddo
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
                         do k = 1, kmax
                            kcn         = max(abs(v0(nm,k)),0.01) * max(1.,wper) / dsiz
                            bpomv(nm,k) = bpomv(nm,k) * ( 1. + 7.5/kcn )
                         enddo
                         !
                      endif
                      !
                   endif
                   !
                enddo
             enddo
             !
             if ( ihydro == 1 .or. ihydro == 2 ) then
                !
                ! loop over wl-points
                !
                do n = nfu, nl
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
                      if ( npors(nm) > 0.1 .and. npors(nm) < 0.99 ) then
                         !
                         dsiz = 0.25 * ( gsizf(nm) + gsizf(nmd) + gsizf(ndm) + gsizf(ndmd) )
                         !
                         if ( dsiz /= 0. ) then
                            !
                            do k = 0, kmax-1
                               kcn         = max(abs(w0(nm,k)),0.01) * max(1.,wper) / dsiz
                               bpoks(nm,k) = bpoks(nm,k) * ( 1. + 7.5/kcn )
                            enddo
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
                         do k = 1, kmax
                            kcn         = max(abs(u0(nm,k)),0.01) * max(1.,wper) / dsiz
                            bpomu(nm,k) = bpomu(nm,k) * ( 1. + 7.5/kcn )
                         enddo
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
                         do k = 1, kmax
                            kcn         = max(abs(v0(nm,k)),0.01) * max(1.,wper) / dsiz
                            bpomv(nm,k) = bpomv(nm,k) * ( 1. + 7.5/kcn )
                         enddo
                         !
                      endif
                      !
                   endif
                   !
                enddo
             enddo
             !
             if ( ihydro == 1 .or. ihydro == 2 ) then
                !
                ! loop over wl-points
                !
                do n = nfu, nl
                   do m = mfu, ml
                      !
                      nm = kgrpnt(m,n)
                      !
                      if ( npors(nm) > 0.1 .and. npors(nm) < 0.99 ) then
                         !
                         if ( dsiz /= 0. ) then
                            !
                            do k = 0, kmax-1
                               kcn         = max(abs(w0(nm,k)),0.01) * max(1.,wper) / dsiz
                               bpoks(nm,k) = bpoks(nm,k) * ( 1. + 7.5/kcn )
                            enddo
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
    endif
    !
end subroutine SwashPorFricLay
