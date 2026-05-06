subroutine SwashUWindStress
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
!    1.00, March 2023: New subroutine
!
!   Purpose
!
!   Calculates wind stresses on triangular mesh
!
!   Method
!
!   The so-called Charnock drag coefficient formulation is proposed by Charnock (1955).
!   He assumed a logarithmic wind velocity profile in the turbulent layer above the free surface:
!
!     w10(z)     1        z
!     ------ = ----- ln (---)
!       u*     kappa      z0
!
!   where u* is the friction velocity, kappa is the Von Karman constant, z is the vertical height
!   above the free surface and z0 is the roughness height:
!
!     z0 = b * u*^2/g
!
!   with b the dimensionless Charnock coefficient and g the gravity acceleration.
!
!   Since, Cd = u*^2/w10^2 we have an implicit relation between drag coefficient Cd and wind speed w10.
!
!   Modules used
!
    use ocpcomm4
    use SwashCommdata2
    use SwashCommdata3
    use SwashTimecomm, only: dt
    use m_genarr
    use SwashFlowdata
    use SwanGriddata
    use SwanGridobjects
!
    implicit none
!
!   Parameter variables
!
    integer, parameter :: maxnit = 100  ! maximum number of iterations
    !
    real   , parameter :: eps  =  0.01  ! convergence criterion
    real   , parameter :: pp   =  0.55  ! first coefficient of 2nd order polynomial fit
    real   , parameter :: qq   =  2.97  ! second coefficient of 2nd order polynomial fit
    real   , parameter :: rr   = -1.49  ! third coefficient of 2nd order polynomial fit
    real   , parameter :: uref = 31.5   ! reference wind speed
!
!   Local variables
!
    integer            :: icelll        ! left cell of present face
    integer            :: icellr        ! right cell of present face
    integer, save      :: ient = 0      ! number of entries in this subroutine
    integer            :: iface         ! loop counter over faces
    integer            :: nit           ! number of iterations
    integer            :: vf1           ! first vertex of present face
    integer            :: vf2           ! second vertex of present face
    !
    real               :: alpha         ! correction factor for wind at 10 metres height to wind at surface
    real               :: cd            ! drag coefficient
    real               :: dsdt          ! vertical velocity of the free surface
    real               :: fac           ! factor
    real               :: finp          ! interpolation factor
    real               :: nx            ! x-component of normal to face
    real               :: ny            ! y-component of normal to face
    real               :: rdx           ! reciprocal of distance between circumcenters adjacent to face
    real               :: rwx           ! relative normal wind velocity
    real               :: rwy           ! relative tangential wind velocity
    real               :: s             ! magnitude u/ustar
    real               :: slope         ! slope of water level normal to face
    real               :: sold          ! ratio u/ustar at previous iteration
    real               :: utl           ! ratio u10/uref
    real               :: ux            ! x-component of wind velocity at midface
    real               :: uy            ! y-component of wind velocity at midface
    real               :: vt            ! tangential flow velocity component at midface
    real               :: w10           ! magnitude of local wind velocity
    real               :: wcel          ! local wave celerity
    real               :: wcrstp        ! percentage of height to apply wind stress on wave crest
    real               :: wfac          ! wind stress multiplication factor
    real               :: wlvl          ! local water level
    real               :: wspeed        ! local wind speed
    real               :: wxc           ! constant wind velocity component in x-direction
    real               :: wyc           ! constant wind velocity component in y-direction
    !
    type(facetype), dimension(:), pointer :: face ! datastructure for faces with their attributes
!
!   Structure
!
!   Description of the pseudo code
!
!   Source text
!
    if (ltrace) call strace (ient,'SwashUWindStress')
    !
    ! point to face object
    !
    face => gridobject%face_grid
    !
    alpha  = pwnd( 4)
    wcrstp = pwnd(10)
    !
    ! compute wind stress multiplication factor
    !
    wfac = rhoa / ( rhow * alpha * alpha )
    !
    ! compute constant wind velocity components
    !
    wxc = u10 * cos(wdic)
    wyc = u10 * sin(wdic)
    !
    ! determine wind velocity at faces
    !
    if ( varwi ) then
       !
       do iface = 1, nfaces
          !
          vf1 = face(iface)%atti(FACEV1)
          vf2 = face(iface)%atti(FACEV2)
          !
          nx = face(iface)%attr(FACENORMX)
          ny = face(iface)%attr(FACENORMY)
          !
          ux = 0.5 * ( wxf(vf1,2) + wxf(vf2,2) )
          uy = 0.5 * ( wyf(vf1,2) + wyf(vf2,2) )
          !
          windu(iface) = nx * ux + ny * uy
          !
          ! note: use array wndimp temporarily to store the tangential wind component
          !
          wndimp(iface) = nx * uy - ny * ux
          !
       enddo
       !
    else
       !
       do iface = 1, nfaces
          !
          nx = face(iface)%attr(FACENORMX)
          ny = face(iface)%attr(FACENORMY)
          !
          windu (iface) = nx * wxc + ny * wyc
          wndimp(iface) = nx * wyc - ny * wxc
          !
       enddo
       !
    endif
    !
    ! adjust normal wind velocity by wave celerity, if appropriate
    !
    if ( relwav ) then
       !
       do iface = 1, nfaces
          !
          if ( face(iface)%atti(FMARKER) == 0 .and. wetu(iface) == 1 ) then   ! internal wet face
             !
             finp = face(iface)%attr(FACELINPF)
             !
             rdx = face(iface)%attr(FACEDISTC)
             !
             ! consider left and right cells of current face
             !
             icelll = face(iface)%atti(FACECL)
             icellr = face(iface)%atti(FACECR)
             !
             slope = rdx * ( s0(icellr) - s0(icelll) )
             fac   = max( abs(slope), 1.e-6 )
             !
             dsdt = ( finp * s1(icelll) + (1.-finp) * s1(icellr) - finp * so(icelll) - (1.-finp) * so(icellr) ) / dt
             !
             wcel = min( abs(dsdt) / fac, sqrt( grav * hum(iface) ) )
             if ( slope < 0. ) wcel = -wcel
             !
             windu(iface) = windu(iface) - wcel
             !
          endif
          !
       enddo
       !
    endif
    !
    ! compute cell-based velocity vector
    !
    call perot ( u0, 1, 1 )
    !
    ! compute wind stress coefficient depending on drag coefficient and (relative) wind velocity
    !
    if ( iwind == 1 ) then
       !
       ! constant drag coefficient
       !
       cd = pwnd(1)
       cd = min( cdcap, cd )
       !
       do iface = 1, nfaces
          !
          if ( face(iface)%atti(FMARKER) == 0 .and. wetu(iface) == 1 ) then   ! internal wet face
             !
             if ( relwnd ) then
                !
                finp = face(iface)%attr(FACELINPF)
                !
                ! consider left and right cells of current face
                !
                icelll = face(iface)%atti(FACECL)
                icellr = face(iface)%atti(FACECR)
                !
                nx = face(iface)%attr(FACENORMX)
                ny = face(iface)%attr(FACENORMY)
                !
                vt = nx * ( finp * uvc(icelll,1,2) + (1.-finp) * uvc(icellr,1,2) ) - ny * ( finp * uvc(icelll,1,1) + (1.-finp) * uvc(icellr,1,1) )
                !
                rwx = alpha*windu (iface) - u0(iface,1)
                rwy = alpha*wndimp(iface) - vt
                !
             else
                !
                rwx = windu (iface)
                rwy = wndimp(iface)
                !
             endif
             !
             w10 = sqrt( rwx**2 + rwy**2 )
             !
             cwndu(iface) = wfac * cd * w10
             !
          endif
          !
       enddo
       !
    else if ( iwind == 2 ) then
       !
       ! drag coefficient based on Charnock drag formula
       !
       do iface = 1, nfaces
          !
          if ( face(iface)%atti(FMARKER) == 0 .and. wetu(iface) == 1 ) then   ! internal wet face
             !
             if ( relwnd ) then
                !
                finp = face(iface)%attr(FACELINPF)
                !
                ! consider left and right cells of current face
                !
                icelll = face(iface)%atti(FACECL)
                icellr = face(iface)%atti(FACECR)
                !
                nx = face(iface)%attr(FACENORMX)
                ny = face(iface)%attr(FACENORMY)
                !
                vt = nx * ( finp * uvc(icelll,1,2) + (1.-finp) * uvc(icellr,1,2) ) - ny * ( finp * uvc(icelll,1,1) + (1.-finp) * uvc(icellr,1,1) )
                !
                rwx = alpha*windu (iface) - u0(iface,1)
                rwy = alpha*wndimp(iface) - vt
                !
             else
                !
                rwx = windu (iface)
                rwy = wndimp(iface)
                !
             endif
             !
             w10 = sqrt( rwx**2 + rwy**2 )
             !
             nit = 0
             !
             ! initial value for s
             !
             s    = 22.4
             sold = 0.
             !
             ! Newton-Raphson iteration
             !
             do
                if ( abs(sold-s) < (eps*s) ) exit
                !
                nit  = nit + 1
                sold = s
                s    = sold*(log(pwnd(3)*grav*sold*sold / (max(0.001,pwnd(2)*w10*w10)))-2.) / (vonkar*sold-2.)
                !
                if ( .not. nit < maxnit ) then
                   call msgerr (1, 'no convergence in Charnock drag computation')
                   s = 22.4
                   exit
                endif
                !
             enddo
             !
             cd = 1./(s*s)
             cd = min( cdcap, cd )
             !
             cwndu(iface) = wfac * cd * w10
             !
          endif
          !
       enddo
       !
    else if ( iwind == 3 ) then
       !
       ! drag coefficient linear dependent on wind speed
       !
       do iface = 1, nfaces
          !
          if ( face(iface)%atti(FMARKER) == 0 .and. wetu(iface) == 1 ) then   ! internal wet face
             !
             rdx = face(iface)%attr(FACEDISTC)
             !
             ! consider left and right cells of current face
             !
             icelll = face(iface)%atti(FACECL)
             icellr = face(iface)%atti(FACECR)
             !
             if ( relwnd ) then
                !
                finp = face(iface)%attr(FACELINPF)
                !
                nx = face(iface)%attr(FACENORMX)
                ny = face(iface)%attr(FACENORMY)
                !
                vt = nx * ( finp * uvc(icelll,1,2) + (1.-finp) * uvc(icellr,1,2) ) - ny * ( finp * uvc(icelll,1,1) + (1.-finp) * uvc(icellr,1,1) )
                !
                rwx = alpha*windu (iface) - u0(iface,1)
                rwy = alpha*wndimp(iface) - vt
                !
             else
                !
                rwx = windu (iface)
                rwy = wndimp(iface)
                !
             endif
             !
             w10 = sqrt( rwx**2 + rwy**2 )
             !
             slope = rdx * ( s0(icellr) - s0(icelll) )
             wspeed = max(pwnd(8), min(w10,pwnd(9)))
             !
             cd = 0.001*( pwnd(5) + pwnd(6)*abs(slope) + pwnd(7)*wspeed )
             cd = min( cdcap, cd )
             !
             cwndu(iface) = wfac * cd * w10
             !
          endif
          !
       enddo
       !
    else if ( iwind == 4 ) then
       !
       ! drag coefficient based on 2nd order polynomial fit
       !
       do iface = 1, nfaces
          !
          if ( face(iface)%atti(FMARKER) == 0 .and. wetu(iface) == 1 ) then   ! internal wet face
             !
             if ( relwnd ) then
                !
                finp = face(iface)%attr(FACELINPF)
                !
                ! consider left and right cells of current face
                !
                icelll = face(iface)%atti(FACECL)
                icellr = face(iface)%atti(FACECR)
                !
                nx = face(iface)%attr(FACENORMX)
                ny = face(iface)%attr(FACENORMY)
                !
                vt = nx * ( finp * uvc(icelll,1,2) + (1.-finp) * uvc(icellr,1,2) ) - ny * ( finp * uvc(icelll,1,1) + (1.-finp) * uvc(icellr,1,1) )
                !
                rwx = alpha*windu (iface) - u0(iface,1)
                rwy = alpha*wndimp(iface) - vt
                !
             else
                !
                rwx = windu (iface)
                rwy = wndimp(iface)
                !
             endif
             !
             w10 = sqrt( rwx**2 + rwy**2 )
             !
             utl = w10/uref
             !
             cd = 0.001*( pp + qq*utl + rr*utl*utl )
             cd = min( cdcap, cd )
             !
             cwndu(iface) = wfac * cd * w10
             !
          endif
          !
       enddo
       !
    else
       !
       call msgerr ( 4, 'unknown drag formula for wind stress' )
       return
       !
    endif
    !
    ! compute wind stress term
    !
    do iface = 1, nfaces
       !
       windu(iface) = alpha * cwndu(iface) * windu(iface)
       !
    enddo
    !
    ! apply wind stress on wave crest only, if appropriate
    !
    if ( relwav ) then
       !
       do iface = 1, nfaces
          !
          if ( face(iface)%atti(FMARKER) == 0 .and. wetu(iface) == 1 ) then   ! internal wet face
             !
             finp = face(iface)%attr(FACELINPF)
             !
             ! consider left and right cells of current face
             !
             icelll = face(iface)%atti(FACECL)
             icellr = face(iface)%atti(FACECR)
             !
             wlvl = finp * s0(icelll) + (1.-finp) * s0(icellr)
             !
             if ( wlvl > smax(iface) ) smax(iface) = wlvl
             !
             if ( wlvl < (1.-wcrstp) * smax(iface) ) windu(iface) = 0.
             !
          endif
          !
       enddo
       !
    endif
    !
end subroutine SwashUWindStress
