subroutine SwashHorzVisc ( u, v )
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
!    1.00,    June 2010: New subroutine
!    1.05, January 2012: include water depth based mixing length model for wave breaking
!
!   Purpose
!
!   Calculates horizontal eddy viscosity coefficient
!
!   Method
!
!   The mixing length model is employed to account for wave breaking
!   The Smagorinsky model is utilized to account for subgrid turbulent mixing
!
!   Note: if eddy viscosity is larger than maximum, which is based on stability criterion, apply clipping
!
!   Modules used
!
    use ocpcomm4
    use SwashCommdata3
    use m_genarr, only: kgrpnt, guu, gvv
    use m_parall
    use SwashFlowdata, only: mf, mfu, ml, mlu, nf, nfu, nl, brks, hs, vnu2d
    use SwashTimeComm, only: dt
!
    implicit none
!
!   Argument variables
!
    real, dimension(mcgrd), intent(in) :: u ! depth-averaged u-velocity
    real, dimension(mcgrd), intent(in) :: v ! depth-averaged v-velocity
!
!   Local variables
!
    integer            :: icistb       ! counter for number of instable points
    integer, save      :: ient = 0     ! number of entries in this subroutine
    integer            :: m            ! loop counter
    integer            :: md           ! index of point m-1
    integer            :: mend         ! end index of loop over u-points
    integer            :: mu           ! index of point m+1
    integer            :: n            ! loop counter
    integer            :: nd           ! index of point n-1
    integer            :: ndm          ! pointer to m,n-1
    integer            :: ndmu         ! pointer to m+1,n-1
    integer            :: nend         ! end index of loop over v-points
    integer            :: nfm          ! pointer to m,nf
    integer            :: nfum         ! pointer to m,nfu
    integer            :: nldm         ! pointer to m,nl-1
    integer            :: nlm          ! pointer to m,nl
    integer            :: nm           ! pointer to m,n
    integer            :: nmd          ! pointer to m-1,n
    integer            :: nmf          ! pointer to mf,n
    integer            :: nmfu         ! pointer to mfu,n
    integer            :: nml          ! pointer to ml,n
    integer            :: nmld         ! pointer to ml-1,n
    integer            :: nmlu         ! pointer to mlu,n
    integer            :: nmu          ! pointer to m+1,n
    integer            :: nu           ! index of point n+1
    integer            :: num          ! pointer to m,n+1
    integer            :: numd         ! pointer to m-1,n+1
    integer            :: numu         ! pointer to m+1,n+1
    !
    real               :: dxl          ! local mesh size in x-direction
    real               :: dyl          ! local mesh size in y-direction
    real               :: rproc        ! auxiliary variable with percentage of instable points
    real               :: stabmx       ! auxiliary variable with maximum viscosity based stability criterion
    real               :: xshear       ! magnitude of shear squared in case of horizontal mixing
    !
    logical            :: STPNOW       ! indicates that program must stop
    !
    character(80)      :: msgstr       ! string to pass message
!
!   Structure
!
!   Description of the pseudo code
!
!   Source text
!
    if (ltrace) call strace (ient,'SwashHorzVisc')
    !
    icistb = 0
    !
    if ( oned .and. (ihvisc == 3 .or. ihvisc == 4) ) then
       !
       ! mixing length model for 1DH simulation
       !
       stabmx = 0.5 * dx * dx / dt
       !
       ! compute eddy viscosity coefficients at internal wl-points
       !
       do m = mfu, ml
          !
          md = m - 1
          !
          nm  = kgrpnt(m ,1)
          nmd = kgrpnt(md,1)
          !
          if ( ihvisc == 3 ) then
             !
             vnu2d(nm) = lmix * lmix * abs(u(nm) - u(nmd)) / dx
             !
          else if ( ihvisc == 4 ) then
             !
             if ( brks(nm) == 1 ) then
                !
                vnu2d(nm) = psurf(3) * hs(nm) * hs(nm) * abs(u(nm) - u(nmd)) / dx
                !
             else
                !
                vnu2d(nm) = 0.
                !
             endif
             !
          endif
          !
          if ( .not. vnu2d(nm) < stabmx ) then
             vnu2d(nm) = stabmx
             icistb    = icistb + 1
          endif
          !
       enddo
       !
       ! exchange eddy viscosity coefficients with neighbouring subdomains
       !
       call SWEXCHG ( vnu2d, kgrpnt, 1, 1 )
       if (STPNOW()) return
       !
       ! at real boundaries, copy them from internal points
       !
       nmf  = kgrpnt(mf ,1)
       nmfu = kgrpnt(mfu,1)
       nml  = kgrpnt(ml ,1)
       nmlu = kgrpnt(mlu,1)
       !
       if ( LMXF ) vnu2d(nmf ) = vnu2d(nmfu)
       if ( LMXL ) vnu2d(nmlu) = vnu2d(nml )
       !
    else if ( ihvisc > 1 ) then
       !
       ! Smagorinsky model or mixing length model
       !
       ! compute eddy viscosity coefficients at internal depth points only
       !
       if ( lreptx ) then
          mend = ml                ! last internal depth point in x-direction in case of repeating grid
       else
          mend = ml - 1            ! last internal depth point in x-direction
       endif
       !
       if ( lrepty ) then
          nend = nl                ! last internal depth point in y-direction in case of repeating grid
       else
          nend = nl - 1            ! last internal depth point in y-direction
       endif
       !
       do n = nf+1, nend
          !
          do m = mf+1, mend
             !
             md = m - 1
             mu = m + 1
             nd = n - 1
             nu = n + 1
             !
             nm   = kgrpnt(m ,n )
             nmd  = kgrpnt(md,n )
             nmu  = kgrpnt(mu,n )
             ndm  = kgrpnt(m ,nd)
             num  = kgrpnt(m ,nu)
             ndmu = kgrpnt(mu,nd)
             numd = kgrpnt(md,nu)
             numu = kgrpnt(mu,nu)
             !
             ! for permanently dry neighbours, corresponding values will be mirrored
             !
             if ( nmd  == 1 ) nmd  = nm
             if ( nmu  == 1 ) nmu  = nm
             if ( ndm  == 1 ) ndm  = nm
             if ( num  == 1 ) num  = nm
             if ( ndmu == 1 ) ndmu = ndm
             if ( numd == 1 ) numd = nmd
             if ( numu == 1 ) numu = nmu
             !
             dxl = 0.5 * ( gvv(nm) + gvv(nmu) )
             dyl = 0.5 * ( guu(nm) + guu(num) )
             !
             xshear = 2. * ( ( u(nmu) + u(numu) - u(nmd) - u(numd))/(4.*dxl) )**2 + &
                      2. * ( ( v(num) + v(numu) - v(ndm) - v(ndmu))/(4.*dyl) )**2 + &
                           ( ( v(nmu) - v(nm))/dxl + (u(num) - u(nm))/dyl )**2
             !
             if ( ihvisc == 2 ) then
                !
                vnu2d(nm) = csmag * csmag * dxl * dyl * sqrt(xshear)
                !
             else if ( ihvisc == 3 ) then
                !
                vnu2d(nm) = lmix * lmix * sqrt(xshear)
                !
             else if ( ihvisc == 4 ) then
                !
                if ( brks(nm) + brks(nmu) + brks(num) + brks(numu) /= 0 ) then
                   !
                   vnu2d(nm) = 0.0625 * psurf(3) * (hs(nm) + hs(nmu) + hs(num) + hs(numu))**2 * sqrt(xshear)
                   !
                else
                   !
                   vnu2d(nm) = 0.
                   !
                endif
                !
             endif
             !
             stabmx = 1./(2.*dt*(1./(dxl*dxl) + 1./(dyl*dyl)))
             !
             if ( .not. vnu2d(nm) < stabmx ) then
                vnu2d(nm) = stabmx
                icistb    = icistb + 1
             endif
             !
          enddo
          !
       enddo
       !
       ! set to zero for permanently dry points
       !
       vnu2d(1) = 0.
       !
       ! exchange eddy viscosity coefficients with neighbouring subdomains
       !
       call SWEXCHG ( vnu2d, kgrpnt, 1, 1 )
       if (STPNOW()) return
       !
       ! at real boundaries, copy them from internal points or apply periodicity
       !
       if ( .not.lreptx ) then
          !
          do n = nf, nl
             !
             nmf  = kgrpnt(mf  ,n)
             nmfu = kgrpnt(mfu ,n)
             nmld = kgrpnt(ml-1,n)
             nml  = kgrpnt(ml  ,n)
             !
             if ( LMXF ) then
                !
                vnu2d(nmf) = vnu2d(nmfu)
                !
             endif
             !
             if ( LMXL ) then
                !
                vnu2d(nml) = vnu2d(nmld)
                !
             endif
             !
          enddo
          !
       else
          !
          call periodic ( vnu2d, kgrpnt, 1, 1 )
          !
       endif
       !
       if ( .not.lrepty ) then
          !
          do m = mf, ml
             !
             nfm  = kgrpnt(m,nf  )
             nfum = kgrpnt(m,nfu )
             nldm = kgrpnt(m,nl-1)
             nlm  = kgrpnt(m,nl  )
             !
             if ( LMYF ) then
                !
                vnu2d(nfm) = vnu2d(nfum)
                !
             endif
             !
             if ( LMYL ) then
                !
                vnu2d(nlm) = vnu2d(nldm)
                !
             endif
             !
          enddo
          !
       else
          !
          call periodic ( vnu2d, kgrpnt, 1, 1 )
          !
       endif
       !
    endif
    !
    ! give warning for instable points
    !
    if ( icistb > 0 ) then
       !
       rproc = 100.*real(icistb)/real(mcgrd)
       !
       if ( .not. rproc < 1. ) then
          write (msgstr,'(a,f5.1)') 'percentage of instable points for computing horizontal eddy viscosity = ',rproc
          call msgerr (1, trim(msgstr) )
       endif
       !
    endif
    !
end subroutine SwashHorzVisc
