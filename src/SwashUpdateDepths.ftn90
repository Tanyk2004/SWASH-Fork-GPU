subroutine SwashUpdateDepths ( u, v )
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
!    4.05: Dirk Rijnsdorp
!
!   Updates
!
!    1.00,   March 2010: New subroutine
!    4.05, January 2017: extension floating objects
!
!   Purpose
!
!   Initialize / update water depths in both water level and velocity points
!
!   Modules used
!
    use ocpcomm4
    use SwashCommdata3
    use m_genarr, only: kgrpnt
    use m_parall
    use SwashFlowdata
    use outp_data, only: hrunp
!
    implicit none
!
!   Argument variables
!
    real, dimension(mcgrd), intent(in) :: u ! depth-averaged u-velocity at current time level
    real, dimension(mcgrd), intent(in) :: v ! depth-averaged v-velocity at current time level
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
    real, dimension(2) :: rtmp     ! temporary array for communication
    real               :: s1min    ! local minimum of water level
    !
    logical            :: adapted  ! true if value of epsdry has been changed
    logical            :: STPNOW   ! indicates that program must stop
    !
    character(80)      :: msgstr   ! string to pass message
!
!   Structure
!
!   Description of the pseudo code
!
!   Source text
!
    if (ltrace) call strace (ient,'SwashUpdateDepths')
    !
    if ( oned ) then
       !
       ! compute and check water depth in wl-point (including virtual ones)
       !
       adapted = .false.
       !
       do m = mf, mlu
          !
          nm = kgrpnt(m,1)
          !
          hs(nm) = s1(nm) + dps(nm)
          !
          if ( hs(nm) < 0. .and. inewt == 0 ) then
             !
             if ( hs(nm) < -epsdry ) then
                !
                write (msgstr,'(a,i5,a,f8.2,a)') 'water depth is negative in m=',m+MXF-2,'; water depth = ',1000.*hs(nm),' mm'
                call msgerr (1, trim(msgstr) )
                !
                s1(nm)  = 0.99*epsdry - dps(nm)
                epsdry  = -1.01*hs(nm)
                adapted = .true.
                !
             else
                s1(nm) = 0.99*epsdry - dps(nm)
             endif
             !
             hs(nm) = s1(nm) + dps(nm)
             !
          endif
          !
          ! adapt water depth to include floating object
          !
          if ( ifloat /= 0 ) hs(nm) = min( dps(nm)-flos(nm), hs(nm) )
          !
          ! ... and in case of Newton iteration, water depth must be bounded below by bottom level
          !
          if ( inewt /= 0 ) hs(nm) = max( epsdry, hs(nm) )
          !
       enddo
       !
       ! exchange water depth in wl-point with neighbouring subdomains
       !
       call SWEXCHG ( hs, kgrpnt, 1, 1 )
       if (STPNOW()) return
       !
       ! check minimal depth in cells
       !
       if ( inewt == 0 ) then
          !
          rtmp(1) = epsdry
          if ( adapted ) then
             rtmp(2) = 1.
          else
             rtmp(2) = 0.
          endif
          !
          call SWREDUCE ( rtmp, 2, SWREAL, SWMAX )
          !
          epsdry = rtmp(1)
          if ( rtmp(2) /= 0. ) then
             adapted = .true.
          else
             adapted = .false.
          endif
          !
          if ( adapted ) then
             write (msgstr,'(a,e12.4)') 'new minimal depth for checking drying and flooding: DEPMIN = ',epsdry
             call msgerr (1, trim(msgstr) )
          endif
          !
          ! check minimal depth for drying and flooding (at most 1 cm)
          !
          if ( epsdry > 0.01 ) then
             !
             call msgerr ( 4, 'INSTABLE: water level is too far below the bottom level!' )
             call msgerr ( 0, '          Please reduce the time step!' )
             return
             !
          endif
          !
          ! also reset minimal depth at faces
          !
          epshu = epsdry
          !
       endif
       !
       ! compute the water depth in u-point based on averaging
       !
       humo = hum
       !
       do m = mf, ml
          !
          mu = m + 1
          !
          nm  = kgrpnt(m ,1)
          nmu = kgrpnt(mu,1)
          !
          hum(nm) = 0.5 * ( hs(nm) + hs(nmu) )
          !
       enddo
       !
       ! exchange averaged water depth in u-point with neighbouring subdomains
       !
       call SWEXCHG ( hum, kgrpnt, 1, 1 )
       if (STPNOW()) return
       !
       ! extrapolate water depth in u-point in time to improve accuracy of momentum-conservative time integration
       !
       humn = 1.5*hum - 0.5*humo
       !
       ! compute the water depth in u-point based on upwinding
       !
       if ( .not.depcds ) then
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
             if ( ifloat == 1 ) then
                !
                hu(nm) = min( dpu(nm)-flou(nm), hu(nm) )
                !
                if ( .not. hu(nm) < dpu(nm) - flou(nm) ) then
                   !
                   ! note: in free surface cell take averaging to prevent discontinuties
                   if ( hs(nm) < dps(nm) - flos(nm) ) then
                      !
                      hu(nm) = 0.5 * ( hs(nm) + hs(nmu) )
                      !
                   endif
                   !
                   if ( hs(nmu) < dps(nmu) - flos(nmu) ) then
                      !
                      hu(nm) = 0.5 * ( hs(nm) + hs(nmu) )
                      !
                   endif
                   !
                endif
                !
             else if ( ifloat == 2 ) then
                !
                hu(nm) = min( dpu(nm)-flou(nm), hu(nm) )
                !
             endif
             !
             ! ... and in case of Newton iteration, water depth must be bounded below by bottom level
             !
             if ( inewt /= 0 ) hu(nm) = max( epsdry, hu(nm) )
             !
          enddo
          !
       else
          !
          hu = hum
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
             if ( ifloat == 0 .or. hu(nm) < ( dpu(nm) - flou(nm) ) ) then
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
       ! exchange upwinded water depths in u-point with neighbouring subdomains
       !
       call SWEXCHG ( hu, kgrpnt, 1, 1 )
       if (STPNOW()) return
       !
    else
       !
       ! compute and check water depth in wl-point (including virtual ones)
       !
       adapted = .false.
       !
       do n = nf, nlu
          do m = mf, mlu
             !
             nm = kgrpnt(m,n)
             !
             hs(nm) = s1(nm) + dps(nm)
             !
             if ( hs(nm) < 0. .and. nm /= 1 .and. inewt == 0 ) then
                !
                if ( hs(nm) < -epsdry ) then
                   !
                   write (msgstr,'(a,i4,a,i4,a,f8.2,a)') 'water depth is negative in (m,n)=(',m+MXF-2,',',n+MYF-2,'); water depth = ',1000.*hs(nm),' mm'
                   call msgerr (1, trim(msgstr) )
                   !
                   s1(nm)  = 0.99*epsdry - dps(nm)
                   epsdry  = -1.01*hs(nm)
                   adapted = .true.
                   !
                else
                   s1(nm) = 0.99*epsdry - dps(nm)
                endif
                !
                hs(nm) = s1(nm) + dps(nm)
                !
             endif
             !
             ! adapt water depth to include floating object
             !
             if ( ifloat /= 0 ) hs(nm) = min( dps(nm)-flos(nm), hs(nm) )
             !
             ! ... and in case of Newton iteration, water depth must be bounded below by bottom level
             !
             if ( inewt /= 0 ) hs(nm) = max( epsdry, hs(nm) )
             !
          enddo
       enddo
       !
       ! set to zero for permanently dry points
       !
       hs(1) = 0.
       !
       ! exchange water depth in wl-point with neighbouring subdomains
       !
       call SWEXCHG ( hs, kgrpnt, 1, 1 )
       if (STPNOW()) return
       !
       ! check minimal depth in cells
       !
       if ( inewt == 0 ) then
          !
          rtmp(1) = epsdry
          if ( adapted ) then
             rtmp(2) = 1.
          else
             rtmp(2) = 0.
          endif
          !
          call SWREDUCE ( rtmp, 2, SWREAL, SWMAX )
          !
          epsdry = rtmp(1)
          if ( rtmp(2) /= 0. ) then
             adapted = .true.
          else
             adapted = .false.
          endif
          !
          if ( adapted ) then
             write (msgstr,'(a,e12.4)') 'new minimal depth for checking drying and flooding: DEPMIN = ',epsdry
             call msgerr (1, trim(msgstr) )
          endif
          !
          ! check minimal depth for drying and flooding (at most 1 cm)
          !
          if ( epsdry > 0.01 ) then
             !
             call msgerr ( 4, 'INSTABLE: water level is too far below the bottom level!' )
             call msgerr ( 0, '          Please reduce the time step!' )
             return
             !
          endif
          !
          ! also reset minimal depth at faces
          !
          epshu = epsdry
          !
       endif
       !
       ! compute the water depth in u-point (including virtual ones) based on averaging
       !
       humo = hum
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
             hum(nm) = 0.5 * ( hs(nm) + hs(nmu) )
             !
          enddo
       enddo
       !
       ! set to zero for permanently dry points
       !
       hum(1) = 0.
       !
       ! exchange averaged water depth in u-point with neighbouring subdomains
       !
       call SWEXCHG ( hum, kgrpnt, 1, 1 )
       if (STPNOW()) return
       !
       ! synchronize water depth in u-point at appropriate boundaries in case of repeating grid
       !
       call periodic ( hum, kgrpnt, 1, 1 )
       !
       ! extrapolate water depth in u-point in time to improve accuracy of momentum-conservative time integration
       !
       humn = 1.5*hum - 0.5*humo
       !
       ! compute the water depth in v-point (including virtual ones) based on averaging
       !
       hvmo = hvm
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
             hvm(nm) = 0.5 * ( hs(nm) + hs(num) )
             !
          enddo
       enddo
       !
       ! set to zero for permanently dry points
       !
       hvm(1) = 0.
       !
       ! exchange averaged water depth in v-point with neighbouring subdomains
       !
       call SWEXCHG ( hvm, kgrpnt, 1, 1 )
       if (STPNOW()) return
       !
       ! synchronize water depth in v-point at appropriate boundaries in case of repeating grid
       !
       call periodic ( hvm, kgrpnt, 1, 1 )
       !
       ! extrapolate water depth in v-point in time to improve accuracy of momentum-conservative time integration
       !
       hvmn = 1.5*hvm - 0.5*hvmo
       !
       ! compute the water depth in u-point (including virtual ones) based on upwinding
       !
       if ( .not.depcds ) then
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
                if ( ifloat == 1 ) then
                   !
                   hu(nm) = min( dpu(nm)-flou(nm), hu(nm) )
                   !
                   if ( .not. hu(nm) < dpu(nm) - flou(nm) ) then
                      !
                      ! note: in free surface cell take averaging to prevent discontinuties
                      if ( hs(nm) < dps(nm) - flos(nm) ) then
                         !
                         hu(nm) = 0.5 * ( hs(nm) + hs(nmu) )
                         !
                      endif
                      !
                      if ( hs(nmu) < dps(nmu) - flos(nmu) ) then
                         !
                         hu(nm) = 0.5 * ( hs(nm) + hs(nmu) )
                         !
                      endif
                      !
                   endif
                   !
                else if ( ifloat == 2 ) then
                   !
                   hu(nm) = min( dpu(nm)-flou(nm), hu(nm) )
                   !
                endif
                !
                ! ... and in case of Newton iteration, water depth must be bounded below by bottom level
                !
                if ( inewt /= 0 ) hu(nm) = max( epsdry, hu(nm) )
                !
             enddo
          enddo
          !
          ! synchronize water depth in u-point at appropriate boundaries in case of repeating grid
          !
          call periodic ( hu, kgrpnt, 1, 1 )
          !
       else
          !
          hu = hum
          !
       endif
       !
       ! compute the water depth in v-point (including virtual ones) based on upwinding
       !
       if ( .not.depcds ) then
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
                if ( ifloat == 1 ) then
                   !
                   hv(nm) = min( dpv(nm)-flov(nm), hv(nm) )
                   !
                   if ( .not. hv(nm) < dpv(nm) - flov(nm) ) then
                      !
                      ! note: in free surface cell take averaging to prevent discontinuties
                      if ( hs(nm) < dps(nm) - flos(nm) ) then
                         !
                         hv(nm) = 0.5 * ( hs(nm) + hs(num) )
                         !
                      endif
                      !
                      if ( hs(num) < dps(num) - flos(num) ) then
                         !
                         hv(nm) = 0.5 * ( hs(nm) + hs(num) )
                         !
                      endif
                      !
                   endif
                   !
                else if ( ifloat == 2 ) then
                   !
                   hv(nm) = min( dpv(nm)-flov(nm), hv(nm) )
                   !
                endif
                !
                ! ... and in case of Newton iteration, water depth must be bounded below by bottom level
                !
                if ( inewt /= 0 ) hv(nm) = max( epsdry, hv(nm) )
                !
             enddo
          enddo
          !
          ! synchronize water depth in v-point at appropriate boundaries in case of repeating grid
          !
          call periodic ( hv, kgrpnt, 1, 1 )
          !
       else
          !
          hv = hvm
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
                if ( ifloat == 0 .or. hu(nm) < ( dpu(nm) - flou(nm) ) ) then
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
          ! synchronize water depth in u-point at appropriate boundaries in case of repeating grid
          !
          call periodic ( hu, kgrpnt, 1, 1 )
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
                if ( ifloat == 0 .or. hv(nm) < ( dpv(nm) - flov(nm) ) ) then
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
          ! synchronize water depth in v-point at appropriate boundaries in case of repeating grid
          !
          call periodic ( hv, kgrpnt, 1, 1 )
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
    endif
    !
    ! compute inundation depth
    !
    do nm = 1, mcgrd
       if ( hindun(nm) /= 1. .and. hs(nm) > hrunp ) hindun(nm) = 1.
    enddo
    !
end subroutine SwashUpdateDepths
