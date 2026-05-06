subroutine SwashElemOutp ( oqproc, mip, nvoqp, nvoqk, voqr, voq, voqk )
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
!    1.00,   March 2022: New subroutine
!
!   Purpose
!
!   Calculates requested output quantities that are stored in mesh elements
!
!   Modules used
!
    use ocpcomm4
    use SwashCommdata1
    use SwashCommdata3
    use SwashTimeComm
    use outp_data
    use m_genarr
    use SwashFlowdata
    use SwanGriddata
!
    implicit none
!
!   Argument variables
!
    integer, intent(in)                            :: mip      ! number of output points
    integer, intent(in)                            :: nvoqk    ! number of layer-dependent quantities per output point
    integer, intent(in)                            :: nvoqp    ! number of quantities per output point
    !
    integer, dimension(nmovar), intent(in)         :: voqr     ! place of each output quantity
    !
    real, dimension(mip,nvoqp), intent(out)        :: voq      ! interpolated output quantity at request
    real, dimension(mip,0:kmax,nvoqk), intent(out) :: voqk     ! interpolated layer-dependent output quantity at request
    !
    logical, dimension(nmovar), intent(in)         :: oqproc   ! indicates whether or not an output quantity must be processed
!
!   Local variables
!
    integer                               :: i        ! loop counter
    integer, save                         :: ient = 0 ! number of entries in this subroutine
    integer                               :: iface    ! face index
    integer                               :: jvx      ! index of u-velocity component
    integer                               :: jvy      ! index of v-velocity component
    integer                               :: k        ! loop counter over vertical layers
    !
    real                                  :: DEGCNV   ! indicates Cartesian or nautical degrees
    real                                  :: dirdeg   ! direction in degrees
    real                                  :: rval1    ! a value
    real                                  :: rval2    ! another value
!
!   Structure
!
!   Description of the pseudo code
!
!   Source text
!
    if (ltrace) call strace (ient,'SwashElemOutp')
    !
    if ( oqproc(95) .or. oqproc(96) ) call chkdiv
    !
    ! water depth
    !
    if ( oqproc(4) ) voq(:,voqr(4)) = hs(:)
    !
    ! bottom level
    !
    if ( oqproc(5) ) voq(:,voqr(5)) = dps(:)
    !
    ! water level
    !
    if ( oqproc(6) ) voq(:,voqr(6)) = s1(:)
    !
    ! non-hydrostatic pressure at bottom
    !
    if ( oqproc(12) ) then
       if ( ihydro /= 0 ) then
          voq(:,voqr(12)) = q(:,kmax)
       else
          voq(:,voqr(12)) = ovexcv(12)
       endif
    endif
    !
    ! pressure at bottom
    !
    if ( oqproc(13) ) then
       if ( ihydro == 0 ) then
          voq(:,voqr(13)) = 0.01 * rhow * grav * hs(:)
       else
          voq(:,voqr(13)) = 0.01 * rhow * ( grav * hs(:) + q(:,kmax) )
       endif
    endif
    !
    if ( oqproc(7) .or. oqproc(8) .or. oqproc(9) .or. oqproc(14) .or. oqproc(15) .or. oqproc(16) .or. oqproc(31) .or. oqproc(32) ) then
       !
       ! first, compute the depth-averaged flow velocity
       !
       if ( kmax == 1 ) then
          !
          uwrk(:) = u1(:,1)
          !
       else
          !
          uwrk = 0.
          !
          do iface = 1, nfaces
             !
             if ( hu(iface) > 0. ) then
                !
                do k = 1, kmax
                   !
                   uwrk(iface) = uwrk(iface) + hku(iface,k) * u1(iface,k)
                   !
                enddo
                !
                uwrk(iface) = uwrk(iface) / hu(iface)
                !
             endif
             !
          enddo
          !
       endif
       !
       ! next, compute cell-based velocity vector
       !
       call perot ( uwrk, 1, 1 )
       !
    endif
    !
    ! flow velocities in x- and y-directions
    !
    if ( oqproc(9) ) then
       jvx = voqr(9)
       jvy = jvx + 1
       voq(:,jvx) = uvc(:,1,1)
       voq(:,jvy) = uvc(:,1,2)
       do i = 1, mip
          rval1 = voq(i,jvx)
          rval2 = voq(i,jvy)
          voq(i,jvx) = coscq*rval1 - sincq*rval2
          voq(i,jvy) = sincq*rval1 + coscq*rval2
       enddo
    endif
    !
    ! velocity magnitude
    !
    if ( oqproc(7) ) then
       jvx = voqr(9)
       jvy = jvx + 1
       do i = 1, mip
          rval1 = voq(i,jvx)
          rval2 = voq(i,jvy)
          voq(i,voqr(7)) = sqrt(rval1*rval1 + rval2*rval2)
       enddo
    endif
    !
    ! velocity direction
    !
    if ( oqproc(8) ) then
       jvx = voqr(9)
       jvy = jvx + 1
       do i = 1, mip
          rval1 = voq(i,jvx)
          rval2 = voq(i,jvy)
          if ( rval1 /= 0. .or. rval2 /= 0. ) then
             if ( bnaut ) then
                dirdeg = atan2(rval2,rval1) * 180./pi
             else
                dirdeg = (alcq + atan2(rval2,rval1)) * 180./pi
             endif
             if ( dirdeg < 0. ) dirdeg = dirdeg + 360.
             !
             ! convert (if necessary) from nautical degrees to Cartesian degrees
             !
             voq(i,voqr(8)) = DEGCNV( dirdeg )
             !
          else
             voq(i,voqr(8)) = ovexcv(8)
          endif
       enddo
    endif
    !
    ! specific discharges in x- and y-directions
    !
    if ( oqproc(16) ) then
       jvx = voqr(16)
       jvy = jvx + 1
       voq(:,jvx) = uvc(:,1,1) * hs(:)
       voq(:,jvy) = uvc(:,1,2) * hs(:)
       do i = 1, mip
          rval1 = voq(i,jvx)
          rval2 = voq(i,jvy)
          voq(i,jvx) = coscq*rval1 - sincq*rval2
          voq(i,jvy) = sincq*rval1 + coscq*rval2
       enddo
    endif
    !
    ! magnitude of specific discharge
    !
    if ( oqproc(14) ) then
       jvx = voqr(16)
       jvy = jvx + 1
       do i = 1, mip
          rval1 = voq(i,jvx)
          rval2 = voq(i,jvy)
          voq(i,voqr(14)) = sqrt(rval1*rval1 + rval2*rval2)
       enddo
    endif
    !
    ! direction of specific discharge
    !
    if ( oqproc(15) ) then
       jvx = voqr(16)
       jvy = jvx + 1
       do i = 1, mip
          rval1 = voq(i,jvx)
          rval2 = voq(i,jvy)
          if ( rval1 /= 0. .or. rval2 /= 0. ) then
             if ( bnaut ) then
                dirdeg = atan2(rval2,rval1) * 180./pi
             else
                dirdeg = (alcq + atan2(rval2,rval1)) * 180./pi
             endif
             if ( dirdeg < 0. ) dirdeg = dirdeg + 360.
             !
             ! convert (if necessary) from nautical degrees to Cartesian degrees
             !
             voq(i,voqr(15)) = DEGCNV( dirdeg )
             !
          else
             voq(i,voqr(15)) = ovexcv(15)
          endif
       enddo
    endif
    !
    ! friction velocities in x- and y-directions
    !
    if ( oqproc(32) ) then
       !
       jvx = voqr(32)
       jvy = jvx + 1
       voq(:,jvx) = sqrt(cfricu(:)) * uvc(:,1,1)
       voq(:,jvy) = sqrt(cfricu(:)) * uvc(:,1,2)
       do i = 1, mip
          rval1 = voq(i,jvx)
          rval2 = voq(i,jvy)
          voq(i,jvx) = coscq*rval1 - sincq*rval2
          voq(i,jvy) = sincq*rval1 + coscq*rval2
       enddo
    endif
    !
    ! magnitude of friction velocity
    !
    if ( oqproc(31) ) then
       jvx = voqr(32)
       jvy = jvx + 1
       do i = 1, mip
          rval1 = voq(i,jvx)
          rval2 = voq(i,jvy)
          voq(i,voqr(31)) = sqrt(rval1*rval1 + rval2*rval2)
       enddo
    endif
    !
    ! salinity
    !
    if ( oqproc(19) ) voq(:,voqr(19)) = rp(:,1,lsal)
    !
    ! temperature
    !
    if ( oqproc(20) ) voq(:,voqr(20)) = rp(:,1,ltemp)
    !
    ! suspended sediment
    !
    if ( oqproc(39) ) voq(:,voqr(39)) = rp(:,1,lsed)
    !
    ! inundation depth / maximum horizontal runup
    !
    if ( oqproc(21) ) voq(:,voqr(21)) = hindun(:)
    !
    ! wave breaking points
    !
    if ( oqproc(38) ) voq(:,voqr(38)) = real(brks(:))
    !
    ! wave-induced setup
    !
    if ( oqproc(22) ) voq(:,voqr(22)) = setup(:)
    !
    ! significant wave height
    !
    if ( oqproc(23) ) voq(:,voqr(23)) = hsig(:)
    !
    ! RMS wave height
    !
    if ( oqproc(24) ) voq(:,voqr(24)) = 0.5 * sqrt(2.) * hsig(:)
    !
    if ( oqproc(33) .or. oqproc(34) .or. oqproc(35) ) then
       !
       ! first, compute depth-averaged mean velocity
       !
       if ( kmax == 1 ) then
          !
          uwrk(:) = mvelu(:,1)
          !
       else
          !
          uwrk = 0.
          !
          do iface = 1, nfaces
             !
             if ( hu(iface) > 0. ) then
                !
                do k = 1, kmax
                   !
                   uwrk(iface) = uwrk(iface) + hku(iface,k) * mvelu(iface,k)
                   !
                enddo
                !
                uwrk(iface) = uwrk(iface) / hu(iface)
                !
             endif
             !
          enddo
          !
       endif
       !
       ! next, compute cell-based velocity vector
       !
       call perot ( uwrk, 1, 1 )
       !
    endif
    !
    ! time-averaged or mean velocities in x- and y-directions
    !
    if ( oqproc(35) ) then
       jvx = voqr(35)
       jvy = jvx + 1
       voq(:,jvx) = uvc(:,1,1)
       voq(:,jvy) = uvc(:,1,2)
       do i = 1, mip
          rval1 = voq(i,jvx)
          rval2 = voq(i,jvy)
          voq(i,jvx) = coscq*rval1 - sincq*rval2
          voq(i,jvy) = sincq*rval1 + coscq*rval2
       enddo
    endif
    !
    ! time-averaged or mean velocity magnitude
    !
    if ( oqproc(33) ) then
       jvx = voqr(35)
       jvy = jvx + 1
       do i = 1, mip
          rval1 = voq(i,jvx)
          rval2 = voq(i,jvy)
          voq(i,voqr(33)) = sqrt(rval1*rval1 + rval2*rval2)
       enddo
    endif
    !
    ! time-averaged or mean velocity direction
    !
    if ( oqproc(34) ) then
       jvx = voqr(35)
       jvy = jvx + 1
       do i = 1, mip
          rval1 = voq(i,jvx)
          rval2 = voq(i,jvy)
          if ( rval1 /= 0. .or. rval2 /= 0. ) then
             if ( bnaut ) then
                dirdeg = atan2(rval2,rval1) * 180./pi
             else
                dirdeg = (alcq + atan2(rval2,rval1)) * 180./pi
             endif
             if ( dirdeg < 0. ) dirdeg = dirdeg + 360.
             !
             ! convert (if necessary) from nautical degrees to Cartesian degrees
             !
             voq(i,voqr(34)) = DEGCNV( dirdeg )
             !
          else
             voq(i,voqr(34)) = ovexcv(34)
          endif
       enddo
    endif
    !
    ! time-averaged or mean salinity
    !
    if ( oqproc(42) ) voq(:,voqr(42)) = mcons(:,1,lsal)
    !
    ! time-averaged or mean temperature
    !
    if ( oqproc(43) ) voq(:,voqr(43)) = mcons(:,1,ltemp)
    !
    ! time-averaged or mean suspended sediment
    !
    if ( oqproc(44) ) voq(:,voqr(44)) = mcons(:,1,lsed)
    !
    ! layer interfaces
    !
    if ( oqproc(51) ) then
       do k = 0, kmax
          voqk(:,k,voqr(51)) = zks(:,k)
       enddo
    endif
    !
    ! layer thicknesses
    !
    if ( oqproc(71) ) then
       do k = 1, kmax
          voqk(:,k,voqr(71)) = hks(:,k)
       enddo
    endif
    !
    ! physical vertical velocity
    !
    if ( oqproc(52) ) then
       if ( ihydro == 0 .or. ihydro == 3 ) then
          do k = 0, kmax
             voqk(:,k,voqr(52)) = wphy(:,k)
          enddo
       else if ( ihydro == 1 .or. ihydro == 2 ) then
          do k = 0, kmax
             voqk(:,k,voqr(52)) = w1(:,k)
          enddo
       endif
    endif
    !
    ! relative vertical velocity
    !
    if ( oqproc(53) ) then
       do k = 0, kmax
          voqk(:,k,voqr(53)) = wom(:,k)
       enddo
    endif
    !
    if ( oqproc(72) .or. oqproc(73) .or. oqproc(74) .or. oqproc(77) .or. oqproc(78) .or. oqproc(79) ) then
       !
       ! compute cell-based velocity vector
       !
       call perot ( u1, 1, kmax )
       !
    endif
    !
    ! layer-dependent flow velocities in x- and y-directions
    !
    if ( oqproc(74) ) then
       jvx = voqr(74)
       jvy = jvx + 1
       do k = 1, kmax
          voqk(:,k,jvx) = uvc(:,k,1)
          voqk(:,k,jvy) = uvc(:,k,2)
       enddo
       do i = 1, mip
          do k = 1, kmax
             rval1 = voqk(i,k,jvx)
             rval2 = voqk(i,k,jvy)
             voqk(i,k,jvx) = coscq*rval1 - sincq*rval2
             voqk(i,k,jvy) = sincq*rval1 + coscq*rval2
          enddo
       enddo
    endif
    !
    ! layer-dependent velocity magnitude
    !
    if ( oqproc(72) ) then
       jvx = voqr(74)
       jvy = jvx + 1
       do i = 1, mip
          do k = 1, kmax
             rval1 = voqk(i,k,jvx)
             rval2 = voqk(i,k,jvy)
             voqk(i,k,voqr(72)) = sqrt(rval1*rval1 + rval2*rval2)
          enddo
       enddo
    endif
    !
    ! layer-dependent velocity direction
    !
    if ( oqproc(73) ) then
       jvx = voqr(74)
       jvy = jvx + 1
       do i = 1, mip
          do k = 1, kmax
             rval1 = voqk(i,k,jvx)
             rval2 = voqk(i,k,jvy)
             if ( rval1 /= 0. .or. rval2 /= 0. ) then
                if ( bnaut ) then
                   dirdeg = atan2(rval2,rval1) * 180./pi
                else
                   dirdeg = (alcq + atan2(rval2,rval1)) * 180./pi
                endif
                if ( dirdeg < 0. ) dirdeg = dirdeg + 360.
                !
                ! convert (if necessary) from nautical degrees to Cartesian degrees
                !
                voqk(i,k,voqr(73)) = DEGCNV( dirdeg )
                !
             else
                voqk(i,k,voqr(73)) = ovexcv(73)
             endif
          enddo
       enddo
    endif
    !
    ! layer-dependent specific discharges in x- and y-directions
    !
    if ( oqproc(79) ) then
       jvx = voqr(79)
       jvy = jvx + 1
       do k = 1, kmax
          voqk(:,k,jvx) = uvc(:,k,1) * hks(:,k)
          voqk(:,k,jvy) = uvc(:,k,2) * hks(:,k)
       enddo
       do i = 1, mip
          do k = 1, kmax
             rval1 = voqk(i,k,jvx)
             rval2 = voqk(i,k,jvy)
             voqk(i,k,jvx) = coscq*rval1 - sincq*rval2
             voqk(i,k,jvy) = sincq*rval1 + coscq*rval2
          enddo
       enddo
    endif
    !
    ! layer-dependent magnitude of specific discharge
    !
    if ( oqproc(77) ) then
       jvx = voqr(79)
       jvy = jvx + 1
       do i = 1, mip
          do k = 1, kmax
             rval1 = voqk(i,k,jvx)
             rval2 = voqk(i,k,jvy)
             voqk(i,k,voqr(77)) = sqrt(rval1*rval1 + rval2*rval2)
          enddo
       enddo
    endif
    !
    ! layer-dependent direction of specific discharge
    !
    if ( oqproc(78) ) then
       jvx = voqr(79)
       jvy = jvx + 1
       do i = 1, mip
          do k = 1, kmax
             rval1 = voqk(i,k,jvx)
             rval2 = voqk(i,k,jvy)
             if ( rval1 /= 0. .or. rval2 /= 0. ) then
                if ( bnaut ) then
                   dirdeg = atan2(rval2,rval1) * 180./pi
                else
                   dirdeg = (alcq + atan2(rval2,rval1)) * 180./pi
                endif
                if ( dirdeg < 0. ) dirdeg = dirdeg + 360.
                !
                ! convert (if necessary) from nautical degrees to Cartesian degrees
                !
                voqk(i,k,voqr(78)) = DEGCNV( dirdeg )
                !
             else
                voqk(i,k,voqr(78)) = ovexcv(78)
             endif
          enddo
       enddo
    endif
    !
    ! layer-dependent non-hydrostatic pressure
    !
    if ( oqproc(82) ) then
       if ( ihydro == 1 .or. ihydro == 2 ) then
          do k = 1, kmax
             voqk(:,k,voqr(82)) = q(:,k)
          enddo
       else
          voqk(:,:,voqr(82)) = ovexcv(82)
       endif
    endif
    !
    ! layer-dependent pressure
    !
    if ( oqproc(83) ) then
       if ( ihydro == 0 ) then
          do k = 1, kmax
             voqk(:,k,voqr(83)) = 0.01 * rhow * grav * ( s1(:) - zks(:,k) )
          enddo
       else if ( ihydro /= 3 ) then
          do k = 1, kmax
             voqk(:,k,voqr(83)) = 0.01 * rhow * ( grav * ( s1(:) - zks(:,k) ) + q(:,k) )
          enddo
       endif
    endif
    !
    if ( oqproc(84) .or. oqproc(85) .or. oqproc(86) ) then
       !
       ! compute cell-based mean velocity vector
       !
       call perot ( mvelu, 1, kmax )
       !
    endif
    !
    ! time-averaged or mean layer-dependent flow velocities in x- and y-directions
    !
    if ( oqproc(86) ) then
       jvx = voqr(86)
       jvy = jvx + 1
       do k = 1, kmax
          voqk(:,k,jvx) = uvc(:,k,1)
          voqk(:,k,jvy) = uvc(:,k,2)
       enddo
       do i = 1, mip
          do k = 1, kmax
             rval1 = voqk(i,k,jvx)
             rval2 = voqk(i,k,jvy)
             voqk(i,k,jvx) = coscq*rval1 - sincq*rval2
             voqk(i,k,jvy) = sincq*rval1 + coscq*rval2
          enddo
       enddo
    endif
    !
    ! time-averaged or mean layer-dependent velocity magnitude
    !
    if ( oqproc(84) ) then
       jvx = voqr(86)
       jvy = jvx + 1
       do i = 1, mip
          do k = 1, kmax
             rval1 = voqk(i,k,jvx)
             rval2 = voqk(i,k,jvy)
             voqk(i,k,voqr(84)) = sqrt(rval1*rval1 + rval2*rval2)
          enddo
       enddo
    endif
    !
    ! time-averaged or mean layer-dependent velocity direction
    !
    if ( oqproc(85) ) then
       jvx = voqr(86)
       jvy = jvx + 1
       do i = 1, mip
          do k = 1, kmax
             rval1 = voqk(i,k,jvx)
             rval2 = voqk(i,k,jvy)
             if ( rval1 /= 0. .or. rval2 /= 0. ) then
                if ( bnaut ) then
                   dirdeg = atan2(rval2,rval1) * 180./pi
                else
                   dirdeg = (alcq + atan2(rval2,rval1)) * 180./pi
                endif
                if ( dirdeg < 0. ) dirdeg = dirdeg + 360.
                !
                ! convert (if necessary) from nautical degrees to Cartesian degrees
                !
                voqk(i,k,voqr(85)) = DEGCNV( dirdeg )
                !
             else
                voqk(i,k,voqr(85)) = ovexcv(85)
             endif
          enddo
       enddo
    endif
    !
    ! layer-dependent salinity
    !
    if ( oqproc(89) ) then
       do k = 1, kmax
          voqk(:,k,voqr(89)) = rp(:,k,lsal)
       enddo
    endif
    !
    ! layer-dependent temperature
    !
    if ( oqproc(90) ) then
       do k = 1, kmax
          voqk(:,k,voqr(90)) =  rp(:,k,ltemp)
       enddo
    endif
    !
    ! layer-dependent suspended sediment
    !
    if ( oqproc(91) ) then
       do k = 1, kmax
          voqk(:,k,voqr(91)) =  rp(:,k,lsed)
       enddo
    endif
    !
    ! time-averaged or mean layer-dependent salinity
    !
    if ( oqproc(92) ) then
       do k = 1, kmax
          voqk(:,k,voqr(92)) = mcons(:,k,lsal)
       enddo
    endif
    !
    ! time-averaged or mean layer-dependent temperature
    !
    if ( oqproc(93) ) then
       do k = 1, kmax
          voqk(:,k,voqr(93)) = mcons(:,k,ltemp)
       enddo
    endif
    !
    ! time-averaged or mean layer-dependent suspended sediment
    !
    if ( oqproc(94) ) then
       do k = 1, kmax
          voqk(:,k,voqr(94)) = mcons(:,k,lsed)
       enddo
    endif
    !
    ! horizontal divergence of flow velocity
    !
    if ( oqproc(95) ) then
       do k = 1, kmax
          voqk(:,k,voqr(95)) = divu(:,k)
       enddo
    endif
    !
    ! horizontal divergence of mass flux
    !
    if ( oqproc(96) ) then
       do k = 1, kmax
          voqk(:,k,voqr(96)) = divq(:,k)
       enddo
    endif
    !
end subroutine SwashElemOutp
