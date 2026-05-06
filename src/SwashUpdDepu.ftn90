subroutine SwashUpdDepu ( u, v )
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
!    1.00, October 2021: New subroutine
!
!   Purpose
!
!   Update water depths in velocity points
!
!   Modules used
!
    use ocpcomm4
    use SwashCommdata3
    use m_genarr, only: kgrpnt
    use m_parall
    use SwashFlowdata
!
    implicit none
!
!   Argument variables
!
    real, dimension(mcgrd),           intent(in) :: u ! depth-averaged u-velocity at current time level
    real, dimension(mcgrd), optional, intent(in) :: v ! depth-averaged v-velocity at current time level
!
!   Local variables
!
    integer, save      :: ient = 0 ! number of entries in this subroutine
    integer            :: m        ! loop counter
    integer            :: md       ! index of point m-1
    integer            :: mend     ! end index of loop over u-points
    integer            :: msta     ! start index of loop over u-points
    integer            :: mu       ! index of point m+1
    integer            :: muu      ! index of point m+2
    integer            :: n        ! loop counter
    integer            :: nd       ! index of point n-1
    integer            :: ndm      ! pointer to m,n-1
    integer            :: nend     ! end index of loop over v-points
    integer            :: nm       ! pointer to m,n
    integer            :: nmd      ! pointer to m-1,n
    integer            :: nmu      ! pointer to m+1,n
    integer            :: nmuu     ! pointer to m+2,n
    integer            :: nsta     ! start index of loop over v-points
    integer            :: nu       ! index of point n+1
    integer            :: num      ! pointer to m,n+1
    integer            :: nuu      ! index of point n+2
    integer            :: nuum     ! pointer to m,n+2
    !
    real               :: depmin   ! local minimum of bottom depth
    real               :: fluxlim  ! flux limiter
    real               :: grad1    ! solution gradient
    real               :: grad2    ! another solution gradient
    real               :: htot     ! water depth in wl-point m,n
    real               :: htotu    ! water depth in wl-point m+1,n or m,n+1
    real               :: s1min    ! local minimum of water level
    !
    logical            :: STPNOW   ! indicates that program must stop
!
!   Structure
!
!   Description of the pseudo code
!
!   Source text
!
    if (ltrace) call strace (ient,'SwashUpdDepu')
    !
    if ( oned ) then
       !
       if ( .not.depcds ) then
          !
          ! compute the water depth in u-point based on upwinding
          !
          do m = mf, ml
             !
             mu = m + 1
             !
             nm  = kgrpnt(m ,1)
             nmu = kgrpnt(mu,1)
             !
             if ( u(nm) > 1.0e-5 ) then
                !
                hu(nm) = s1(nm) + dpu(nm)
                !
             else if ( u(nm) < -1.0e-5 ) then
                !
                hu(nm) = s1(nmu) + dpu(nm)
                !
             else
                !
                hu(nm) = max( s1(nm), s1(nmu) ) + dpu(nm)
                !
             endif
             !
          enddo
          !
       else
          !
          ! compute the water depth in u-point based on averaging
          !
          do m = mf, ml
             !
             mu = m + 1
             !
             nm  = kgrpnt(m ,1)
             nmu = kgrpnt(mu,1)
             !
             htot  = s1(nm ) + dps(nm )
             htotu = s1(nmu) + dps(nmu)
             !
             hu(nm) = 0.5 * ( htot + htotu )
             !
          enddo
          !
       endif
       !
       ! compute higher order correction to the water depth in internal u-point (if appropriate)
       !
       if ( corrdep ) then
          !
          propsc = nint(pnums(11))
          kappa  = pnums(12)
          mbound = pnums(13)
          phieby = pnums(14)
          !
          msta = mf + 1   ! first internal u-point
          mend = ml - 1   ! last  internal u-point
          !
          uloop: do m = msta, mend
             !
             md  = m  - 1
             mu  = m  + 1
             muu = mu + 1
             !
             nm   = kgrpnt(m  ,1)
             nmd  = kgrpnt(md ,1)
             nmu  = kgrpnt(mu ,1)
             nmuu = kgrpnt(muu,1)
             !
             if ( hu(nm) < ( dpu(nm) - flou(nm) ) ) then
                !
                if ( u(nm) > 1.0e-5 ) then
                   !
                   depmin = min( dps(nmd), dps(nm), dps(nmu) )
                   s1min  = min( s1 (nmd), s1 (nm), s1 (nmu) )
                   !
                   if ( s1min + depmin < 0. ) cycle uloop
                   !
                   grad1 = s1(nmu) - s1(nm )
                   grad2 = s1(nm ) - s1(nmd)
                   !
                   hu(nm) = hu(nm) + 0.5 * fluxlim(grad1,grad2)
                   !
                else if ( u(nm) < -1.0e-5 ) then
                   !
                   depmin = min( dps(nm), dps(nmu), dps(nmuu) )
                   s1min  = min( s1 (nm), s1 (nmu), s1 (nmuu) )
                   !
                   if ( s1min + depmin < 0. ) cycle uloop
                   !
                   grad1 = s1(nmu ) - s1(nm )
                   grad2 = s1(nmuu) - s1(nmu)
                   !
                   hu(nm) = hu(nm) - 0.5 * fluxlim(grad1,grad2)
                   !
                endif
                !
             endif
             !
          enddo uloop
          !
       endif
       !
       ! exchange water depth in u-point with neighbouring subdomains
       !
       call SWEXCHG ( hu, kgrpnt, 1, 1 )
       if (STPNOW()) return
       !
    else
       !
       if ( .not.depcds ) then
          !
          ! compute the water depth in u-point (including virtual ones) based on upwinding
          !
          do n = nf, nlu
             do m = mf, ml
                !
                mu = m + 1
                !
                nm  = kgrpnt(m ,n)
                nmu = kgrpnt(mu,n)
                !
                if ( nmu == 1 ) nmu = nm
                !
                if ( u(nm) > 1.0e-5 ) then
                   !
                   hu(nm) = s1(nm) + dpu(nm)
                   !
                else if ( u(nm) < -1.0e-5 ) then
                   !
                   hu(nm) = s1(nmu) + dpu(nm)
                   !
                else
                   !
                   hu(nm) = max( s1(nm), s1(nmu) ) + dpu(nm)
                   !
                endif
                !
             enddo
          enddo
          !
          ! compute the water depth in v-point (including virtual ones) based on upwinding
          !
          do m = mf, mlu
             do n = nf, nl
                !
                nu = n + 1
                !
                nm  = kgrpnt(m,n )
                num = kgrpnt(m,nu)
                !
                if ( num == 1 ) num = nm
                !
                if ( v(nm) > 1.0e-5 ) then
                   !
                   hv(nm) = s1(nm) + dpv(nm)
                   !
                else if ( v(nm) < -1.0e-5 ) then
                   !
                   hv(nm) = s1(num) + dpv(nm)
                   !
                else
                   !
                   hv(nm) = max( s1(nm), s1(num) ) + dpv(nm)
                   !
                endif
                !
             enddo
          enddo
          !
       else
          !
          ! compute the water depth in u-point (including virtual ones) based on averaging
          !
          do n = nf, nlu
             do m = mf, ml
                !
                mu = m + 1
                !
                nm  = kgrpnt(m ,n)
                nmu = kgrpnt(mu,n)
                !
                if ( nmu == 1 ) nmu = nm
                !
                htot  = s1(nm ) + dps(nm )
                htotu = s1(nmu) + dps(nmu)
                !
                hu(nm) = 0.5 * ( htot + htotu )
                !
             enddo
          enddo
          !
          ! compute the water depth in v-point (including virtual ones) based on averaging
          !
          do m = mf, mlu
             do n = nf, nl
                !
                nu = n + 1
                !
                nm  = kgrpnt(m,n )
                num = kgrpnt(m,nu)
                !
                if ( num == 1 ) num = nm
                !
                htot  = s1(nm ) + dps(nm )
                htotu = s1(num) + dps(num)
                !
                hv(nm) = 0.5 * ( htot + htotu )
                !
             enddo
          enddo
          !
       endif
       !
       ! compute higher order correction to the water depth in internal velocity points (if appropriate)
       !
       if ( corrdep ) then
          !
          propsc = nint(pnums(11))
          kappa  = pnums(12)
          mbound = pnums(13)
          phieby = pnums(14)
          !
          msta = mf + 1               ! first internal u-point
          if ( lreptx ) then
             mend = ml                ! last  internal u-point in case of repeating grid
          else
             mend = ml - 1            ! last  internal u-point
          endif
          !
          do n = nfu, nl
             !
             uuloop: do m = msta, mend
                !
                md  = m  - 1
                mu  = m  + 1
                muu = mu + 1
                !
                if ( lreptx .and. LMXL .and. muu > mlu ) muu = mlu
                !
                nm   = kgrpnt(m  ,n)
                nmd  = kgrpnt(md ,n)
                nmu  = kgrpnt(mu ,n)
                nmuu = kgrpnt(muu,n)
                !
                if ( nmd  == 1 ) nmd  = nm
                if ( nmu  == 1 ) nmu  = nm
                if ( nmuu == 1 ) nmuu = nmu
                !
                if ( hu(nm) < ( dpu(nm) - flou(nm) ) ) then
                   !
                   if ( u(nm) > 1.0e-5 ) then
                      !
                      depmin = min( dps(nmd), dps(nm), dps(nmu) )
                      s1min  = min( s1 (nmd), s1 (nm), s1 (nmu) )
                      !
                      if ( s1min + depmin < 0. ) cycle uuloop
                      !
                      grad1 = s1(nmu) - s1(nm )
                      grad2 = s1(nm ) - s1(nmd)
                      !
                      hu(nm) = hu(nm) + 0.5 * fluxlim(grad1,grad2)
                      !
                   else if ( u(nm) < -1.0e-5 ) then
                      !
                      depmin = min( dps(nm), dps(nmu), dps(nmuu) )
                      s1min  = min( s1 (nm), s1 (nmu), s1 (nmuu) )
                      !
                      if ( s1min + depmin < 0. ) cycle uuloop
                      !
                      grad1 = s1(nmu ) - s1(nm )
                      grad2 = s1(nmuu) - s1(nmu)
                      !
                      hu(nm) = hu(nm) - 0.5 * fluxlim(grad1,grad2)
                      !
                   endif
                   !
                endif
                !
             enddo uuloop
             !
          enddo
          !
          nsta = nf + 1               ! first internal v-point
          if ( lrepty ) then
             nend = nl                ! last  internal v-point in case of repeating grid
          else
             nend = nl - 1            ! last  internal v-point
          endif
          !
          do m = mfu, ml
             !
             vloop: do n = nsta, nend
                !
                nd  = n  - 1
                nu  = n  + 1
                nuu = nu + 1
                !
                if ( lrepty .and. LMYL .and. nuu > nlu ) nuu = nlu
                !
                nm   = kgrpnt(m,n  )
                ndm  = kgrpnt(m,nd )
                num  = kgrpnt(m,nu )
                nuum = kgrpnt(m,nuu)
                !
                if ( ndm  == 1 ) ndm  = nm
                if ( num  == 1 ) num  = nm
                if ( nuum == 1 ) nuum = num
                !
                if ( hv(nm) < ( dpv(nm) - flov(nm) ) ) then
                   !
                   if ( v(nm) > 1.0e-5 ) then
                      !
                      depmin = min( dps(ndm), dps(nm), dps(num) )
                      s1min  = min( s1 (ndm), s1 (nm), s1 (num) )
                      !
                      if ( s1min + depmin < 0. ) cycle vloop
                      !
                      grad1 = s1(num) - s1(nm )
                      grad2 = s1(nm ) - s1(ndm)
                      !
                      hv(nm) = hv(nm) + 0.5 * fluxlim(grad1,grad2)
                      !
                   else if ( v(nm) < -1.0e-5 ) then
                      !
                      depmin = min( dps(nm), dps(num), dps(nuum) )
                      s1min  = min( s1 (nm), s1 (num), s1 (nuum) )
                      !
                      if ( s1min + depmin < 0. ) cycle vloop
                      !
                      grad1 = s1(num ) - s1(nm )
                      grad2 = s1(nuum) - s1(num)
                      !
                      hv(nm) = hv(nm) - 0.5 * fluxlim(grad1,grad2)
                      !
                   endif
                   !
                endif
                !
             enddo vloop
             !
          enddo
          !
       endif
       !
       ! set to zero for permanently dry points
       !
       hu(1) = 0.
       hv(1) = 0.
       !
       ! exchange upwinded water depths in u- and v-points with neighbouring subdomains
       !
       call SWEXCHG ( hu, kgrpnt, 1, 1 )
       if (STPNOW()) return
       call SWEXCHG ( hv, kgrpnt, 1, 1 )
       if (STPNOW()) return
       !
       ! synchronize water depths in u- and v-points at appropriate boundaries in case of repeating grid
       !
       call periodic ( hu, kgrpnt, 1, 1 )
       call periodic ( hv, kgrpnt, 1, 1 )
       !
    endif
    !
end subroutine SwashUpdDepu
