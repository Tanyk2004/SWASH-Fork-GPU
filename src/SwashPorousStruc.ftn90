subroutine SwashPorousStruc
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
!    1.00, March 2012: New subroutine
!
!   Purpose
!
!   Determines parameters for porous structures
!
!   Method
!
!   Both laminar and turbulent frictional forces inside porous medium
!   are calculated based on the Darcy and Forchheimer formulation
!   The empirical formula's of Van Gent (1995) are applied
!   Added mass effect for accelerating fluid in porous medium is included
!
!   Note: the intrinsic volume-averaged variables are employed
!         (not the Darcy's averaged variables!)
!         the frictional forces are then multiplied by the volumetric porosity
!
!   M.R.A. van Gent (1995)
!   Wave interaction with permeable coastal structures
!   PhD thesis, Delft University of Technology, Delft, the Netherlands
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
!   Parameter variables
!
    real,    parameter :: gamma = 0.34 ! empirical coefficient for added mass
!
!   Local variables
!
    integer, save :: ient = 0 ! number of entries in this subroutine
    integer       :: m        ! loop counter
    integer       :: n        ! loop counter
    integer       :: md       ! index of point m-1
    integer       :: nd       ! index of point n-1
    integer       :: ndm      ! pointer to m,n-1
    integer       :: ndmd     ! pointer to m-1,n-1
    integer       :: nm       ! pointer to m,n
    integer       :: nmd      ! pointer to m-1,n
!
    real          :: alpha0   ! dimensionless constant for laminar friction factor
    real          :: beta0    ! dimensionless constant for turbulent friction factor
    real          :: dsiz     ! grain size
    !
    logical       :: STPNOW   ! indicates that program must stop
!
!   Structure
!
!   Description of the pseudo code
!
!   Source text
!
    if (ltrace) call strace (ient,'SwashPorousStruc')
    !
    dsiz   = ppor(1)
    alpha0 = ppor(3)
    beta0  = ppor(4)
    !
    ! determine porosity in water level and velocity points
    !
    if ( oned ) then
       !
       do m = mfu, ml
          !
          md = m - 1
          !
          nm  = kgrpnt(m ,1)
          nmd = kgrpnt(md,1)
          !
          npors(nm) = 0.5 * ( nporf(nm) + nporf(nmd) )
          !
       enddo
       !
       do m = mf, ml
          !
          nm = kgrpnt(m,1)
          !
          nporu(nm) = nporf(nm)
          !
       enddo
       !
    else
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
             npors(nm) = 0.25 * ( nporf(nm) + nporf(nmd) + nporf(ndm) + nporf(ndmd) )
             !
          enddo
       enddo
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
             nporu(nm) = 0.5 * ( nporf(nm) + nporf(ndm) )
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
             nporv(nm) = 0.5 * ( nporf(nm) + nporf(nmd) )
             !
          enddo
       enddo
       !
       ! set to one for permanently dry points
       !
       npors(1) = 1.
       nporu(1) = 1.
       nporv(1) = 1.
       !
    endif
    !
    ! exchange porosity with neighbouring subdomains
    !
    call SWEXCHG ( npors, kgrpnt, 1, 1 )
    call SWEXCHG ( nporu, kgrpnt, 1, 1 )
    if ( .not.oned ) call SWEXCHG ( nporv, kgrpnt, 1, 1 )
    if (STPNOW()) return
    !
    ! synchronize porosity at appropriate boundaries in case of repeating grid
    !
    call periodic ( npors, kgrpnt, 1, 1 )
    call periodic ( nporu, kgrpnt, 1, 1 )
    if ( .not.oned ) call periodic ( nporv, kgrpnt, 1, 1 )
    !
    ! determine structure height in water level and velocity points, if appropriate
    !
    if ( varsh ) then
       !
       if ( oned ) then
          !
          do m = mfu, ml
             !
             md = m - 1
             !
             nm  = kgrpnt(m ,1)
             nmd = kgrpnt(md,1)
             !
             if ( npors(nm) > 0.1 .and. npors(nm) < 0.99 ) hpors(nm) = 0.5 * ( hporf(nm) + hporf(nmd) )
             !
          enddo
          !
          do m = mf, ml
             !
             nm = kgrpnt(m,1)
             !
             if ( nporu(nm) > 0.1 .and. nporu(nm) < 0.99 ) hporu(nm) = hporf(nm)
             !
          enddo
          !
       else
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
                if ( npors(nm) > 0.1 .and. npors(nm) < 0.99 ) hpors(nm) = 0.25 * ( hporf(nm) + hporf(nmd) + hporf(ndm) + hporf(ndmd) )
                !
             enddo
          enddo
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
                if ( nporu(nm) > 0.1 .and. nporu(nm) < 0.99 ) hporu(nm) = 0.5 * ( hporf(nm) + hporf(ndm) )
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
                if ( nporv(nm) > 0.1 .and. nporv(nm) < 0.99 ) hporv(nm) = 0.5 * ( hporf(nm) + hporf(nmd) )
                !
             enddo
          enddo
          !
          ! set to default for permanently dry points
          !
          hpors(1) = ppor(2)
          hporu(1) = ppor(2)
          hporv(1) = ppor(2)
          !
       endif
       !
    endif
    !
    ! exchange structure heights with neighbouring subdomains
    !
    call SWEXCHG ( hpors, kgrpnt, 1, 1 )
    call SWEXCHG ( hporu, kgrpnt, 1, 1 )
    if ( .not.oned ) call SWEXCHG ( hporv, kgrpnt, 1, 1 )
    if (STPNOW()) return
    !
    ! synchronize structure heights at appropriate boundaries in case of repeating grid
    !
    call periodic ( hpors, kgrpnt, 1, 1 )
    call periodic ( hporu, kgrpnt, 1, 1 )
    if ( .not.oned ) call periodic ( hporv, kgrpnt, 1, 1 )
    !
    ! calculate the friction factors and added mass coefficient
    !
    if ( oned ) then
       !
       if ( vargs ) then
          !
          do m = mfu, ml
             !
             md = m - 1
             !
             nm  = kgrpnt(m ,1)
             nmd = kgrpnt(md,1)
             !
             dsiz = 0.5 * ( gsizf(nm) + gsizf(nmd) )
             !
             if ( dsiz /= 0. ) then
                !
                apors(nm) = alpha0 * kinvis * ( 1. - npors(nm) )**2. / ( npors(nm) * npors(nm) * dsiz * dsiz )
                bpors(nm) = beta0 * ( 1. - npors(nm) ) / ( npors(nm) * npors(nm) * dsiz )
                cpors(nm) = gamma * ( 1. - npors(nm) ) / npors(nm)
                !
             endif
             !
          enddo
          !
          do m = mf, ml
             !
             nm = kgrpnt(m,1)
             !
             dsiz = gsizf(nm)
             !
             if ( dsiz /= 0. ) then
                !
                aporu(nm) = alpha0 * kinvis * ( 1. - nporu(nm) )**2. / ( nporu(nm) * nporu(nm) * dsiz * dsiz )
                bporu(nm) = beta0 * ( 1. - nporu(nm) ) / ( nporu(nm) * nporu(nm) * dsiz )
                cporu(nm) = gamma * ( 1. - nporu(nm) ) / nporu(nm)
                !
             endif
             !
          enddo
          !
       else
          !
          do m = mfu, ml
             !
             nm = kgrpnt(m,1)
             !
             if ( dsiz /= 0. ) then
                !
                apors(nm) = alpha0 * kinvis * ( 1. - npors(nm) )**2. / ( npors(nm) * npors(nm) * dsiz * dsiz )
                bpors(nm) = beta0 * ( 1. - npors(nm) ) / ( npors(nm) * npors(nm) * dsiz )
                cpors(nm) = gamma * ( 1. - npors(nm) ) / npors(nm)
                !
             endif
             !
          enddo
          !
          do m = mf, ml
             !
             nm = kgrpnt(m,1)
             !
             if ( dsiz /= 0. ) then
                !
                aporu(nm) = alpha0 * kinvis * ( 1. - nporu(nm) )**2. / ( nporu(nm) * nporu(nm) * dsiz * dsiz )
                bporu(nm) = beta0 * ( 1. - nporu(nm) ) / ( nporu(nm) * nporu(nm) * dsiz )
                cporu(nm) = gamma * ( 1. - nporu(nm) ) / nporu(nm)
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
                dsiz = 0.25 * ( gsizf(nm) + gsizf(nmd) + gsizf(ndm) + gsizf(ndmd) )
                !
                if ( dsiz /= 0. ) then
                   !
                   apors(nm) = alpha0 * kinvis * ( 1. - npors(nm) )**2. / ( npors(nm) * npors(nm) * dsiz * dsiz )
                   bpors(nm) = beta0 * ( 1. - npors(nm) ) / ( npors(nm) * npors(nm) * dsiz )
                   cpors(nm) = gamma * ( 1. - npors(nm) ) / npors(nm)
                   !
                endif
                !
             enddo
          enddo
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
                dsiz = 0.5 * ( gsizf(nm) + gsizf(ndm) )
                !
                if ( dsiz /= 0. ) then
                   !
                   aporu(nm) = alpha0 * kinvis * ( 1. - nporu(nm) )**2. / ( nporu(nm) * nporu(nm) * dsiz * dsiz )
                   bporu(nm) = beta0 * ( 1. - nporu(nm) ) / ( nporu(nm) * nporu(nm) * dsiz )
                   cporu(nm) = gamma * ( 1. - nporu(nm) ) / nporu(nm)
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
                dsiz = 0.5 * ( gsizf(nm) + gsizf(nmd) )
                !
                if ( dsiz /= 0. ) then
                   !
                   aporv(nm) = alpha0 * kinvis * ( 1. - nporv(nm) )**2. / ( nporv(nm) * nporv(nm) * dsiz * dsiz )
                   bporv(nm) = beta0 * ( 1. - nporv(nm) ) / ( nporv(nm) * nporv(nm) * dsiz )
                   cporv(nm) = gamma * ( 1. - nporv(nm) ) / nporv(nm)
                   !
                endif
                !
             enddo
          enddo
          !
       else
          !
          ! loop over wl-points
          !
          do n = nfu, nl
             do m = mfu, ml
                !
                nm = kgrpnt(m,n)
                !
                if ( dsiz /= 0. ) then
                   !
                   apors(nm) = alpha0 * kinvis * ( 1. - npors(nm) )**2. / ( npors(nm) * npors(nm) * dsiz * dsiz )
                   bpors(nm) = beta0 * ( 1. - npors(nm) ) / ( npors(nm) * npors(nm) * dsiz )
                   cpors(nm) = gamma * ( 1. - npors(nm) ) / npors(nm)
                   !
                endif
                !
             enddo
          enddo
          !
          ! loop over u-points
          !
          do n = nfu, nl
             do m = mf, ml
                !
                nm = kgrpnt(m,n)
                !
                if ( dsiz /= 0. ) then
                   !
                   aporu(nm) = alpha0 * kinvis * ( 1. - nporu(nm) )**2. / ( nporu(nm) * nporu(nm) * dsiz * dsiz )
                   bporu(nm) = beta0 * ( 1. - nporu(nm) ) / ( nporu(nm) * nporu(nm) * dsiz )
                   cporu(nm) = gamma * ( 1. - nporu(nm) ) / nporu(nm)
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
                if ( dsiz /= 0. ) then
                   !
                   aporv(nm) = alpha0 * kinvis * ( 1. - nporv(nm) )**2. / ( nporv(nm) * nporv(nm) * dsiz * dsiz )
                   bporv(nm) = beta0 * ( 1. - nporv(nm) ) / ( nporv(nm) * nporv(nm) * dsiz )
                   cporv(nm) = gamma * ( 1. - nporv(nm) ) / nporv(nm)
                   !
                endif
                !
             enddo
          enddo
          !
       endif
       !
       ! synchronize friction factors and added mass coefficient at appropriate boundaries in case of repeating grid
       !
       call periodic ( aporu, kgrpnt, 1, 1 )
       call periodic ( bporu, kgrpnt, 1, 1 )
       call periodic ( cporu, kgrpnt, 1, 1 )
       call periodic ( aporv, kgrpnt, 1, 1 )
       call periodic ( bporv, kgrpnt, 1, 1 )
       call periodic ( cporv, kgrpnt, 1, 1 )
       call periodic ( apors, kgrpnt, 1, 1 )
       call periodic ( bpors, kgrpnt, 1, 1 )
       call periodic ( cpors, kgrpnt, 1, 1 )
       !
    endif
    !
end subroutine SwashPorousStruc
