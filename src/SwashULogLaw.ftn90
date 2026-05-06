subroutine SwashULogLaw
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
!    1.00,  February 2023: New subroutine
!
!   Purpose
!
!   Computes friction at bed based on the logarithmic or linear law of the wall in case of flexible mesh
!
!   Method
!
!   The following parameters will be used
!
!   z  = normal distance to the bottom
!   hr = bottom roughness height
!
!   z+ = cmu^(1/4) sqrt(tke) z / nu
!   roughness Reynolds number = cmu^(1/4) sqrt(tke) hr / nu
!
!   It is assumed that there is no transition layer between viscous and turbulent sublayers and so,
!   the edge of viscous sublayer is determined by the following implicit relation
!
!   ev = ln(E*ev) / kappa
!
!   from which it follows that ev = 11.6
!
!   when roughness Reynolds number < 11.6, bed is considered to be smooth, i.e. hr = 0
!
!   if bed is rough (hr /= 0)
!
!      if the location of first computational point is between the roughness elements then z = hr/30
!
!      tau_w = cmu^1/4 kappa sqrt(tke) u / ln(E_rough*z/hr)
!      u+    = ln( E_rough*z/hr )/kappa
!
!   if bed is smooth (hr = 0)
!
!      if the computational point lies in the viscous sublayer then
!
!         tau_w = nu * u / z
!         u+    = z+
!
!      else if the point lies in the logarithmic layer then
!
!         tau_w = cmu^1/4 kappa sqrt(tke) u / ln(E z+)
!         u+    = ln( E z+ )/kappa
!
!   B.E. Launder and D.B. Spalding, "The numerical computation of turbulent
!   flows", Comput. Meth. Appl. Mech. Engng., vol. 3, p. 269-289, 1974
!
!   In very shallow water (depth less than 5 cm) we assume a viscous sublayer
!
!   Various studies have shown that the sediment stratification effects in the bottom boundary layer (BBL)
!   can lead to a change in BBL dynamics. The velocity profile for a stable stratified logarithmic
!   boundary layer is given by Adams and Weatherly (1981). This actually amounts to the reduction of
!   the Von Karman constant, such that the velocity near the bed is enhanced.
!
!   C.E. Adams and G.L. Weatherly, "Some effects of suspended sediment
!   stratification on an oceanic bottom boundary layer", JGR, vol. 86,
!   no. C5, p. 4161-4172, 1981
!
!   Modules used
!
    use ocpcomm4
    use SwashCommdata2
    use SwashCommdata3
    use m_genarr, only: fricf
    use SwashFlowdata
    use SwanGriddata
    use SwanGridobjects
!
    implicit none
!
!   Parameter variables
!
    real,    parameter :: cmu    = 0.09  ! closure constant for standard k-eps model
    real,    parameter :: erough = 33.0  ! empirical constant for logarithmic log-law in case of rough beds
    real,    parameter :: esmoot = 9.0   ! empirical constant for logarithmic log-law in case of smooth beds
    real,    parameter :: ev     = 11.6  ! edge of viscous sublayer
    real,    parameter :: hmin   = 0.05  ! minimal water depth in meters
!
!   Local variables
!
    integer               :: icell    ! loop over cells
    integer, save         :: ient = 0 ! number of entries in this subroutine
    integer, dimension(3) :: v        ! vertices of present cell
    !
    real                  :: cf       ! local friction coefficient
    real                  :: drhodz   ! vertical gradient of density
    real                  :: dudz     ! vertical gradient of u-velocity
    real                  :: dvdz     ! vertical gradient of v-velocity
    real                  :: dz       ! local layer thickness
    real                  :: fac      ! stability factor
    real                  :: hr       ! local roughness height
    real                  :: r        ! roughness Reynolds number
    real                  :: richn    ! gradient Richardson number
    real                  :: shear    ! magnitude of shear squared
    real                  :: tke      ! turbulent kinetic energy at bed
    real                  :: uplus    ! dimensionless flow velocity
    real                  :: ustar    ! friction velocity
    real                  :: z0       ! roughness parameter
    real                  :: zplus    ! dimensionless normal distance to bed
    real                  :: zs       ! distance of half layer thickness to bed
    !
    type(celltype), dimension(:), pointer :: cell ! datastructure for cells with their attributes
!
!   Structure
!
!   Description of the pseudo code
!
!   Source text
!
    if (ltrace) call strace (ient,'SwashULogLaw')
    !
    ! point to cell object
    !
    cell => gridobject%cell_grid
    !
    if ( varfr ) then
       !
       do icell = 1, ncells
          !
          ! get vertices of present cell
          !
          v(1) = cell(icell)%atti(CELLV1)
          v(2) = cell(icell)%atti(CELLV2)
          v(3) = cell(icell)%atti(CELLV3)
          !
          cf = ( fricf(v(1),2) + fricf(v(2),2) + fricf(v(3),2) ) / 3.
          !
          zs = 0.5 * ( zks(icell,kmax-1) - zks(icell,kmax) )
          !
          if ( hs(icell) > hmin ) then
             !
             tke = rtur(icell,kmax,1)
             !
             ustar = ( cmu**.25 ) * sqrt(tke)
             !
             zplus = ustar * zs / kinvis
             r     = ustar * cf / kinvis
             !
             if ( r > ev ) then
                hr = cf
             else
                hr = 0.
             endif
             !
             z0 = hr / erough
             !
             uplus = zplus
             !
             if ( z0 /= 0. ) then
                !
                if ( .not. zs > z0 ) zs = zs + z0   ! the bottom is shifted to z = z0
                !
                uplus = log( zs/z0 ) / vonkar
                !
             else
                !
                if ( zplus > ev ) uplus = log( esmoot*zplus ) / vonkar
                !
             endif
             !
          else
             !
             ustar = 1.                   ! dummy value
             zplus = ustar * zs / kinvis
             uplus = zplus
             !
          endif
          !
          logfrc(icell,1) = ustar / uplus
          logfrc(icell,2) = uplus
          !
       enddo
       !
    else
       !
       do icell = 1, ncells
          !
          zs = 0.5 * ( zks(icell,kmax-1) - zks(icell,kmax) )
          !
          if ( hs(icell) > hmin ) then
             !
             tke = rtur(icell,kmax,1)
             !
             ustar = ( cmu**.25 ) * sqrt(tke)
             !
             zplus = ustar * zs      / kinvis
             r     = ustar * pbot(2) / kinvis
             !
             if ( r > ev ) then
                hr = pbot(2)
             else
                hr = 0.
             endif
             !
             z0 = hr / erough
             !
             uplus = zplus
             !
             if ( z0 /= 0. ) then
                !
                if ( .not. zs > z0 ) zs = zs + z0   ! the bottom is shifted to z = z0
                !
                uplus = log( zs/z0 ) / vonkar
                !
             else
                !
                if ( zplus > ev ) uplus = log( esmoot*zplus ) / vonkar
                !
             endif
             !
          else
             !
             ustar = 1.                   ! dummy value
             zplus = ustar * zs / kinvis
             uplus = zplus
             !
          endif
          !
          logfrc(icell,1) = ustar / uplus
          logfrc(icell,2) = uplus
          !
       enddo
       !
    endif
    !
    if ( lsed > 0 .and. psed(13) > 0. ) then
       !
       ! compute cell-based velocity vector
       !
       call perot ( u0, 1, kmax )
       !
       do icell = 1, ncells
          !
          if ( hs(icell) > hmin ) then
             !
             dz = 0.5 * ( hks(icell,kmax-1) + hks(icell,kmax) )
             !
             ! compute vertical shear squared near the bed
             !
             dudz = ( uvc(icell,kmax-1,1) - uvc(icell,kmax,1) ) / dz
             dudz = ( uvc(icell,kmax-1,2) - uvc(icell,kmax,2) ) / dz
             !
             shear = max ( dudz*dudz + dvdz*dvdz, 1.e-8 )
             !
             ! compute vertical gradient of density near the bed
             !
             drhodz = ( rho(icell,kmax-1) - rho(icell,kmax) ) / dz
             !
             ! compute gradient Richardson number near the bed
             !
             richn = -grav * drhodz / ( shear*rhow )
             !
             ! compute stability factor to reduce the bed shear stress due to the BBL stratification induced by resuspended load
             !
             fac = psed(13) * max( 0., min( 0.2, richn ) )
             !
             logfrc(icell,1) = logfrc(icell,1) / ( 1. + fac )
             logfrc(icell,2) = logfrc(icell,2) * ( 1. + fac )
             !
          endif
          !
       enddo
       !
    endif
    !
end subroutine SwashULogLaw
