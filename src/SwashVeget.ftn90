subroutine SwashVeget
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
!    1.00: Tomohiro Suzuki
!    5.01: Kenji Kumada, Tomohiro Suzuki
!
!   Updates
!
!    1.00, April 2012: New subroutine
!    5.01, March 2018: inclusion porosity effect
!
!   Purpose
!
!   Calculates vegetation coefficients
!
!   Method
!
!   The wave damping induced by aquatic vegetation is described by
!   the Morison type equation, modelling the plants as vertical,
!   noncompliant cylinders, neclecting swaying motions induced by
!   waves
!
!   The drag force consists of viscous effect and form drag around
!   the cylinders which is modelled as one process and is proportional
!   to the square of the characteristic velocity V, which is the
!   ambient flow. The inertia force due to local acceleration (the
!   time derivative of V) is optionally included
!
!   For the case of densely spaced cylinders (e.g. dense mangrove
!   fields, porous brushwood groins), the effect of (relatively
!   low) porosity should be included. In this case the characteristic
!   velocity is the pore velocity defined as V/n, with n the
!   porosity, i.e. the ratio between the fluid volume and the total
!   volume. The volume of the cylinders equals 1 - n, i.e. the
!   vegetation volume in the unit volume
!
!   Vegetation characteristics that are used as input are the
!   drag coefficient, vegetation height, plant density and stem
!   diameter. Based on these quantities the spatial occupation
!   of vegetation per unit volume is calculated in the velocity
!   points and temporarily stored in work arrays
!
!   The total force due to vegetation including inertia and porosity
!   effects consists of the drag force, the inertia force related to
!   the vegetation volume, and the inertia force related to the clear
!   fluid, and is given by Eq. (61) of Burcharth and Andersen (1995)
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
    integer       :: k        ! loop counter
    integer       :: l        ! loop counter
    integer       :: m        ! loop counter
    integer       :: md       ! index of point m-1
    integer       :: mu       ! index of point m+1
    integer       :: n        ! loop counter
    integer       :: nd       ! index of point n-1
    integer       :: ndm      ! pointer to m,n-1
    integer       :: nm       ! pointer to m,n / loop counter
    integer       :: nmd      ! pointer to m-1,n
    integer       :: nmu      ! pointer to m+1,n
    integer       :: nu       ! index of point n+1
    integer       :: num      ! pointer to m,n+1
    !
    real          :: fac      ! auxiliary factor
    real          :: fsum     ! auxiliary sum
    real          :: hvtot    ! total vegetation height
    real          :: np       ! porosity of vegetation field
    real          :: tdragf   ! total contribution to drag force over vegetation layers
    real          :: tinerf   ! total contribution to inertia force over vegetation layers
    real          :: zv       ! vegetation layer interface
    real          :: zvh      ! vegetation height level
    !
    logical       :: STPNOW   ! indicates that program must stop
!
!   Structure
!
!   Description of the pseudo code
!
!   Source text
!
    if (ltrace) call strace (ient,'SwashVeget')
    !
    ! initialize arrays to store spatial occupation of vegetation per unit volume in velocity points
    !
    work(:,1) = 0.     ! u-point
    work(:,2) = 0.     ! v-point
    !
    ! compute vegetation coefficient integrated over layers
    !
    hvtot = sum(hlayv)
    !
    if ( kmax == 1 ) then
       !
       if ( oned ) then
          !
          do m = mf, ml
             !
             nm = kgrpnt(m,1)
             !
             if ( wetu(nm) == 1 .and. hum(nm) > 0. ) then
                !
                tdragf = 0.
                tinerf = 0.
                !
                if ( hum(nm) > hvtot ) then
                   !
                   ! canopy is submerged
                   !
                   do l = 1, lmax
                      tdragf = tdragf + cdveg(l) * bveg(l) * nveg(l) * hlayv(l)
                      tinerf = tinerf +  nveg(l) * bveg(l) * bveg(l) * hlayv(l)
                   enddo
                   !
                else
                   !
                   ! canopy is emerged
                   !
                   zv = 0.
                   !
                   do l = 1, lmax
                      !
                      zvh = zv + hlayv(l)
                      !
                      if ( zv < hum(nm) .and. .not. zvh < hum(nm) ) then
                         !
                         ! upper part of layer or whole layer is vegetated
                         !
                         tdragf = tdragf + cdveg(l) * bveg(l) * nveg(l) * ( hum(nm) - zv )
                         tinerf = tinerf +  nveg(l) * bveg(l) * bveg(l) * ( hum(nm) - zv )
                         !
                      else if ( zvh < hum(nm) ) then
                         !
                         ! lower or middle part of layer is vegetated
                         !
                         tdragf = tdragf + cdveg(l) * bveg(l) * nveg(l) * hlayv(l)
                         tinerf = tinerf +  nveg(l) * bveg(l) * bveg(l) * hlayv(l)
                         !
                      endif
                      !
                      zv = zvh
                      !
                   enddo
                   !
                endif
                !
                work(nm,1) = 0.25 * pi * tinerf / hum(nm)
                !
                cvegu(nm,1,1) = 0.5 * tdragf / hum(nm)
                cvegu(nm,1,2) = cvm * work(nm,1)
                !
             endif
             !
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
                if ( wetu(nm) == 1 .and. hum(nm) > 0. ) then
                   !
                   tdragf = 0.
                   tinerf = 0.
                   !
                   if ( hum(nm) > hvtot ) then
                      !
                      ! canopy is submerged
                      !
                      do l = 1, lmax
                         tdragf = tdragf + cdveg(l) * bveg(l) * nveg(l) * hlayv(l)
                         tinerf = tinerf +  nveg(l) * bveg(l) * bveg(l) * hlayv(l)
                      enddo
                      !
                   else
                      !
                      ! canopy is emerged
                      !
                      zv = 0.
                      !
                      do l = 1, lmax
                         !
                         zvh = zv + hlayv(l)
                         !
                         if ( zv < hum(nm) .and. .not. zvh < hum(nm) ) then
                            !
                            ! upper part of layer or whole layer is vegetated
                            !
                            tdragf = tdragf + cdveg(l) * bveg(l) * nveg(l) * ( hum(nm) - zv )
                            tinerf = tinerf +  nveg(l) * bveg(l) * bveg(l) * ( hum(nm) - zv )
                            !
                         else if ( zvh < hum(nm) ) then
                            !
                            ! lower or middle part of layer is vegetated
                            !
                            tdragf = tdragf + cdveg(l) * bveg(l) * nveg(l) * hlayv(l)
                            tinerf = tinerf +  nveg(l) * bveg(l) * bveg(l) * hlayv(l)
                            !
                         endif
                         !
                         zv = zvh
                         !
                      enddo
                      !
                   endif
                   !
                   work(nm,1) = 0.25 * pi * tinerf / hum(nm)
                   !
                   cvegu(nm,1,1) = 0.5 * tdragf / hum(nm)
                   cvegu(nm,1,2) = cvm * work(nm,1)
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
                if ( wetv(nm) == 1 .and. hvm(nm) > 0. ) then
                   !
                   tdragf = 0.
                   tinerf = 0.
                   !
                   if ( hvm(nm) > hvtot ) then
                      !
                      ! canopy is submerged
                      !
                      do l = 1, lmax
                         tdragf = tdragf + cdveg(l) * bveg(l) * nveg(l) * hlayv(l)
                         tinerf = tinerf +  nveg(l) * bveg(l) * bveg(l) * hlayv(l)
                      enddo
                      !
                   else
                      !
                      ! canopy is emerged
                      !
                      zv = 0.
                      !
                      do l = 1, lmax
                         !
                         zvh = zv + hlayv(l)
                         !
                         if ( zv < hvm(nm) .and. .not. zvh < hvm(nm) ) then
                            !
                            ! upper part of layer or whole layer is vegetated
                            !
                            tdragf = tdragf + cdveg(l) * bveg(l) * nveg(l) * ( hvm(nm) - zv )
                            tinerf = tinerf +  nveg(l) * bveg(l) * bveg(l) * ( hvm(nm) - zv )
                            !
                         else if ( zvh < hvm(nm) ) then
                            !
                            ! lower or middle part of layer is vegetated
                            !
                            tdragf = tdragf + cdveg(l) * bveg(l) * nveg(l) * hlayv(l)
                            tinerf = tinerf +  nveg(l) * bveg(l) * bveg(l) * hlayv(l)
                            !
                         endif
                         !
                         zv = zvh
                         !
                      enddo
                      !
                   endif
                   !
                   work(nm,2) = 0.25 * pi * tinerf / hvm(nm)
                   !
                   cvegv(nm,1,1) = 0.5 * tdragf / hvm(nm)
                   cvegv(nm,1,2) = cvm * work(nm,2)
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
       if ( oned ) then
          !
          do m = mf, ml
             !
             nm = kgrpnt(m,1)
             !
             if ( wetu(nm) == 1 ) then
                !
                fsum = 0.
                !
                ! compute vegetation coefficient for each computational layer
                !
                do k = kmax, 1, -1
                   !
                   tdragf = 0.
                   tinerf = 0.
                   !
                   if ( zkum(nm,k) - zkum(nm,kmax) < hvtot ) then
                      !
                      zv = zkum(nm,kmax)
                      !
                      do l = 1, lmax
                         !
                         zvh = zv + hlayv(l)
                         !
                         if ( zkum(nm,k) < zv .and.  zv < zkum(nm,k-1) .and. .not. zvh < zkum(nm,k-1) ) then
                            !
                            ! upper part of layer is vegetated
                            !
                            tdragf = tdragf + cdveg(l) * bveg(l) * nveg(l) * ( zkum(nm,k-1) - zv )
                            tinerf = tinerf +  nveg(l) * bveg(l) * bveg(l) * ( zkum(nm,k-1) - zv )
                            !
                         else if ( .not. zkum(nm,k) < zv .and. zkum(nm,k) < zvh .and. zvh < zkum(nm,k-1) ) then
                            !
                            ! lower part of layer is vegetated
                            !
                            tdragf = tdragf + cdveg(l) * bveg(l) * nveg(l) * ( zvh - zkum(nm,k) )
                            tinerf = tinerf +  nveg(l) * bveg(l) * bveg(l) * ( zvh - zkum(nm,k) )
                            !
                         else if ( zkum(nm,k) < zv .and. zvh < zkum(nm,k-1) ) then
                            !
                            ! middle part of layer is vegetated
                            !
                            tdragf = tdragf + cdveg(l) * bveg(l) * nveg(l) * hlayv(l)
                            tinerf = tinerf +  nveg(l) * bveg(l) * bveg(l) * hlayv(l)
                            !
                         else if ( .not. zkum(nm,k) < zv .and. .not. zvh < zkum(nm,k-1) ) then
                            !
                            ! layer is completely vegetated
                            !
                            tdragf = tdragf + cdveg(l) * bveg(l) * nveg(l) * hkum(nm,k)
                            tinerf = tinerf +  nveg(l) * bveg(l) * bveg(l) * hkum(nm,k)
                            !
                         endif
                         !
                         zv = zvh
                         !
                      enddo
                      !
                   endif
                   !
                   cvegu(nm,k,1) = 0.5 * tdragf / hkum(nm,k)
                   cvegu(nm,k,2) = 0.25*pi * cvm * tinerf / hkum(nm,k)
                   !
                   fsum = fsum + tinerf
                   !
                enddo
                !
                work(nm,1) = 0.25 * pi * fsum / hum(nm)
                !
             endif
             !
          enddo
          !
       else
          !
          ! loop over u-points
          !
          do n = nfu, nl
             do m = mf, ml
                !
                mu = m + 1
                !
                nm  = kgrpnt(m ,n)
                nmu = kgrpnt(mu,n)
                !
                if ( wetu(nm) == 1 .and. nmu /= 1 ) then
                   !
                   fsum = 0.
                   !
                   ! compute vegetation coefficient for each computational layer
                   !
                   do k = kmax, 1, -1
                      !
                      tdragf = 0.
                      tinerf = 0.
                      !
                      if ( zkum(nm,k) - zkum(nm,kmax) < hvtot ) then
                         !
                         zv = zkum(nm,kmax)
                         !
                         do l = 1, lmax
                            !
                            zvh = zv + hlayv(l)
                            !
                            if ( zkum(nm,k) < zv .and.  zv < zkum(nm,k-1) .and. .not. zvh < zkum(nm,k-1) ) then
                               !
                               ! upper part of layer is vegetated
                               !
                               tdragf = tdragf + cdveg(l) * bveg(l) * nveg(l) * ( zkum(nm,k-1) - zv )
                               tinerf = tinerf +  nveg(l) * bveg(l) * bveg(l) * ( zkum(nm,k-1) - zv )
                               !
                            else if ( .not. zkum(nm,k) < zv .and. zkum(nm,k) < zvh .and. zvh < zkum(nm,k-1) ) then
                               !
                               ! lower part of layer is vegetated
                               !
                               tdragf = tdragf + cdveg(l) * bveg(l) * nveg(l) * ( zvh - zkum(nm,k) )
                               tinerf = tinerf +  nveg(l) * bveg(l) * bveg(l) * ( zvh - zkum(nm,k) )
                               !
                            else if ( zkum(nm,k) < zv .and. zvh < zkum(nm,k-1) ) then
                               !
                               ! middle part of layer is vegetated
                               !
                               tdragf = tdragf + cdveg(l) * bveg(l) * nveg(l) * hlayv(l)
                               tinerf = tinerf +  nveg(l) * bveg(l) * bveg(l) * hlayv(l)
                               !
                            else if ( .not. zkum(nm,k) < zv .and. .not. zvh < zkum(nm,k-1) ) then
                               !
                               ! layer is completely vegetated
                               !
                               tdragf = tdragf + cdveg(l) * bveg(l) * nveg(l) * hkum(nm,k)
                               tinerf = tinerf +  nveg(l) * bveg(l) * bveg(l) * hkum(nm,k)
                               !
                            endif
                            !
                            zv = zvh
                            !
                         enddo
                         !
                      endif
                      !
                      cvegu(nm,k,1) = 0.5 * tdragf / hkum(nm,k)
                      cvegu(nm,k,2) = 0.25*pi * cvm * tinerf / hkum(nm,k)
                      !
                      fsum = fsum + tinerf
                      !
                   enddo
                   !
                   work(nm,1) = 0.25 * pi * fsum / hum(nm)
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
                nu = n + 1
                !
                nm  = kgrpnt(m,n )
                num = kgrpnt(m,nu)
                !
                if ( wetv(nm) == 1 .and. num /= 1 ) then
                   !
                   fsum = 0.
                   !
                   ! compute vegetation coefficient for each computational layer
                   !
                   do k = kmax, 1, -1
                      !
                      tdragf = 0.
                      tinerf = 0.
                      !
                      if ( zkvm(nm,k) - zkvm(nm,kmax) < hvtot ) then
                         !
                         zv = zkvm(nm,kmax)
                         !
                         do l = 1, lmax
                            !
                            zvh = zv + hlayv(l)
                            !
                            if ( zkvm(nm,k) < zv .and.  zv < zkvm(nm,k-1) .and. .not. zvh < zkvm(nm,k-1) ) then
                               !
                               ! upper part of layer is vegetated
                               !
                               tdragf = tdragf + cdveg(l) * bveg(l) * nveg(l) * ( zkvm(nm,k-1) - zv )
                               tinerf = tinerf +  nveg(l) * bveg(l) * bveg(l) * ( zkvm(nm,k-1) - zv )
                               !
                            else if ( .not. zkvm(nm,k) < zv .and. zkvm(nm,k) < zvh .and. zvh < zkvm(nm,k-1) ) then
                               !
                               ! lower part of layer is vegetated
                               !
                               tdragf = tdragf + cdveg(l) * bveg(l) * nveg(l) * ( zvh - zkvm(nm,k) )
                               tinerf = tinerf +  nveg(l) * bveg(l) * bveg(l) * ( zvh - zkvm(nm,k) )
                               !
                            else if ( zkvm(nm,k) < zv .and. zvh < zkvm(nm,k-1) ) then
                               !
                               ! middle part of layer is vegetated
                               !
                               tdragf = tdragf + cdveg(l) * bveg(l) * nveg(l) * hlayv(l)
                               tinerf = tinerf +  nveg(l) * bveg(l) * bveg(l) * hlayv(l)
                               !
                            else if ( .not. zkvm(nm,k) < zv .and. .not. zvh < zkvm(nm,k-1) ) then
                               !
                               ! layer is completely vegetated
                               !
                               tdragf = tdragf + cdveg(l) * bveg(l) * nveg(l) * hkvm(nm,k)
                               tinerf = tinerf +  nveg(l) * bveg(l) * bveg(l) * hkvm(nm,k)
                               !
                            endif
                            !
                            zv = zvh
                            !
                         enddo
                         !
                      endif
                      !
                      cvegv(nm,k,1) = 0.5 * tdragf / hkvm(nm,k)
                      cvegv(nm,k,2) = 0.25*pi * cvm * tinerf / hkvm(nm,k)
                      !
                      fsum = fsum + tinerf
                      !
                   enddo
                   !
                   work(nm,2) = 0.25 * pi * fsum / hvm(nm)
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
    ! multiply vegetation coefficients with horizontally varying density, if appropriate
    !
    if ( varnpl ) then
       !
       if ( oned ) then
          !
          do m = mf, ml
             !
             nm = kgrpnt(m,1)
             !
             cvegu(nm,:,:) = cvegu(nm,:,:) * nplaf(nm)
             !
             work(nm,1) = work(nm,1) * nplaf(nm)
             !
          enddo
          !
       else
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
                fac = 0.5 * ( nplaf(nm) + nplaf(ndm) )
                !
                cvegu(nm,:,:) = cvegu(nm,:,:) * fac
                !
                work(nm,1) = work(nm,1) * fac
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
                fac = 0.5 * ( nplaf(nm) + nplaf(nmd) )
                !
                cvegv(nm,:,:) = cvegv(nm,:,:) * fac
                !
                work(nm,2) = work(nm,2) * fac
                !
             enddo
          enddo
          !
       endif
       !
    endif
    !
    ! adapt vegetation coefficients to include porosity effect, if appropriate
    ! (compare Eq. (61) with Eq. (54) of Burcharth and Andersen, 1995)
    !
    if ( iveg == 1 .and. iporos /= 0 ) then
       !
       if ( oned ) then
          !
          do m = mf, ml
             !
             nm = kgrpnt(m,1)
             !
             np = nporu(nm)
             !
             cvegu(nm,:,1) = cvegu(nm,:,1) / np / np / np
             cvegu(nm,:,2) = cvegu(nm,:,2) / np / np
             !
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
                np = nporu(nm)
                !
                cvegu(nm,:,1) = cvegu(nm,:,1) / np / np / np
                cvegu(nm,:,2) = cvegu(nm,:,2) / np / np
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
                np = nporv(nm)
                !
                cvegv(nm,:,1) = cvegv(nm,:,1) / np / np / np
                cvegv(nm,:,2) = cvegv(nm,:,2) / np / np
                !
             enddo
          enddo
          !
       endif
       !
    else if ( iveg == 2 ) then
       !
       if ( oned ) then
          !
          do m = mf, ml
             !
             nm = kgrpnt(m,1)
             !
             np = 1. - work(nm,1)
             !
             cvegu(nm,:,1) = cvegu(nm,:,1) / np / np / np
             cvegu(nm,:,2) = cvegu(nm,:,2) / np / np
             !
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
                np = 1. - work(nm,1)
                !
                cvegu(nm,:,1) = cvegu(nm,:,1) / np / np / np
                cvegu(nm,:,2) = cvegu(nm,:,2) / np / np
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
                np = 1. - work(nm,2)
                !
                cvegv(nm,:,1) = cvegv(nm,:,1) / np / np / np
                cvegv(nm,:,2) = cvegv(nm,:,2) / np / np
                !
             enddo
          enddo
          !
       endif
       !
    endif
    !
    ! exchange vegetation coefficients with neighbouring subdomains
    !
    call SWEXCHG ( cvegu(1,1,1), kgrpnt, 1, kmax )
    call SWEXCHG ( cvegu(1,1,2), kgrpnt, 1, kmax )
    if ( .not.oned ) then
       call SWEXCHG ( cvegv(1,1,1), kgrpnt, 1, kmax )
       call SWEXCHG ( cvegv(1,1,2), kgrpnt, 1, kmax )
    endif
    if (STPNOW()) return
    !
    ! synchronize vegetation coefficients at appropriate boundaries in case of repeating grid
    !
    call periodic ( cvegu(1,1,1), kgrpnt, 1, kmax )
    call periodic ( cvegu(1,1,2), kgrpnt, 1, kmax )
    if ( .not.oned ) then
       call periodic ( cvegv(1,1,1), kgrpnt, 1, kmax )
       call periodic ( cvegv(1,1,2), kgrpnt, 1, kmax )
    endif
    !
end subroutine SwashVeget
