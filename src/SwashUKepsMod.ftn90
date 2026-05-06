subroutine SwashUKepsMod
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
!   Builds and solves standard k-epsilon model for vertical mixing on triangular mesh
!
!   Method
!
!   The time integration is semi-implicit. First order timesplitting is used.
!
!   In the first stage the horizontal advection terms are integrated explictly.
!   In space, it is discretized with first order upwinding.
!
!   In the second stage the vertical transport terms are integrated implicitly in
!   conjunction with the source and sink terms. Sinks are treated implicitly, while
!   sources are taken explicitly.
!   The vertical advection term is approximated by means of first order upwinding.
!   In the vertical, the tri-diagonal systems of equations are solved by a double sweep.
!
!   Modules used
!
    use ocpcomm4
    use SwashCommdata3
    use SwashTimecomm
    use m_genarr, only: work
    use SwashFlowdata
    use SwanGriddata
    use SwanGridobjects
!
    implicit none
!
!   Parameter variables
!
    real,    parameter :: eps    = 1.e-7 ! small number
    real,    parameter :: hmin   = 1.e-3 ! minimal water depth in meters
    !
    real,    parameter :: ceps1  = 1.44  ! closure constant for standard k-eps model
    real,    parameter :: ceps2  = 1.92  ! closure constant for standard k-eps model
    real,    parameter :: cmu    = 0.09  ! closure constant for standard k-eps model
    real,    parameter :: sigmak = 1.0   ! closure constant for standard k-eps model
    real,    parameter :: sigmae = 1.3   ! closure constant for standard k-eps model
    real,    parameter :: sigrho = 0.5   ! closure constant for buoyancy term
!
!   Local variables
!
    integer               :: icell    ! loop over cells
    integer               :: icelll   ! left cell of present face
    integer               :: icellr   ! right cell of present face
    integer, save         :: ient = 0 ! number of entries in this subroutine
    integer               :: iface    ! face index
    integer               :: jf       ! loop counter over faces
    integer               :: k        ! loop counter over vertical layer interfaces
    integer               :: kend     ! end index of loop over vertical layer interfaces
    integer               :: ku       ! index of layer k+1
    integer               :: l        ! loop counter over turbulence quantities
    !
    real                  :: area     ! area of present cell
    real                  :: bi       ! inverse of main diagonal of the matrix
    real                  :: buotke   ! local buoyancy term representing exchange between turbulent energy and potential energy
    real                  :: contrib  ! total contribution of advective flux per cell
    real                  :: ctrkb    ! contribution of vertical advection term below considered point
    real                  :: ctrkt    ! contribution of vertical advection term above considered point
    real                  :: drhodz   ! vertical gradient of density
    real                  :: dudz     ! vertical gradient of u-velocity
    real                  :: dvdz     ! vertical gradient of v-velocity
    real                  :: dz       ! local layer thickness
    real                  :: fac      ! auxiliary factor
    real                  :: lf       ! length of present face
    real                  :: psm      ! Prandtl-Schmidt number
    real                  :: rsgn     ! sign for indicating face orientation
    real                  :: u        ! local flow velocity
    real                  :: ustar2   ! friction velocity squared
    real                  :: w        ! local relative vertical velocity
    real                  :: zs       ! distance of half layer thickness to bottom / free surface in waterlevel point
    !
    type(celltype), dimension(:), pointer :: cell ! datastructure for cells with their attributes
    type(facetype), dimension(:), pointer :: face ! datastructure for faces with their attributes
!
!   Structure
!
!   Description of the pseudo code
!
!   Source text
!
    if (ltrace) call strace (ient,'SwashUKepsMod')
    !
    ! point to cell and face objects
    !
    cell => gridobject%cell_grid
    face => gridobject%face_grid
    !
    ! first part of time splitting method (horizontal advection terms)
    !
    ! loop over turbulence quantities
    !
    do l = 1, ltur
       !
       if ( irough == 4 .and. l == 1 ) then
          !
          ! in case of logarithmic wall-law an equation for turbulent kinetic energy is build up at bottom
          !
          kend = kmax
       else
          kend = kmax-1
       endif
       !
       ! update turbulence quantity in cells due to advective transport
       !
       do k = 1, kend
          !
          ku = min(k+1,kmax)
          !
          do icell = 1, ncells
             !
             contrib = 0.
             !
             ! loop over faces of the cell
             !
             do jf = 1, cell(icell)%nof
                !
                ! face identifier
                !
                iface = cell(icell)%face(jf)%atti(FACEID)
                !
                ! get length of current face
                !
                lf = face(iface)%attr(FACELEN)
                !
                ! consider left and right cells of current face
                !
                icelll = face(iface)%atti(FACECL)
                icellr = face(iface)%atti(FACECR)
                !
                ! take into account orientation of the current face
                !
                if ( icell == icelll ) then
                   rsgn =  1.
                else if ( icell == icellr ) then
                   rsgn = -1.
                endif
                !
                u = ( hku(iface,ku)*u1(iface,k) + hku(iface,k)*u1(iface,ku) ) / ( hku(iface,k) + hku(iface,ku) )
                !
                ! compute the advective flux based on first order upwind and add to other faces of the cell
                !
                if ( u > 0. ) then
                   !
                   contrib = contrib + rsgn * lf * u * rtur(icelll,k,l)
                   !
                else
                   !
                   contrib = contrib + rsgn * lf * u * rtur(icellr,k,l)
                   !
                endif
                !
             enddo
             !
             ! update in wet cell
             !
             if ( hs(icell) > hmin ) then
                !
                area = cell(icell)%attr(CELLAREA)
                !
                rtur(icell,k,l) = rtur(icell,k,l) - dt * contrib / area
                !
             endif
             !
          enddo
          !
       enddo
       !
    enddo
    !
    ! second part of time splitting method (vertical advection and diffusion, production, buoyancy and dissipation rates)
    !
    ! initialize system of equations in dry cells
    !
    do icell = 1, ncells
       !
       if ( .not. hs(icell) > hmin ) then
          !
          amatt(icell,:,1) = 1.
          amatt(icell,:,2) = 0.
          amatt(icell,:,3) = 0.
          rhst (icell,:  ) = eps
          !
       endif
       !
    enddo
    !
    ! compute cell-based velocity vector
    !
    if ( irough /= 4 ) then
       call perot ( udep, 1, 1 )
       work(:,1) = uvc(:,1,1)
       work(:,2) = uvc(:,1,2)
    endif
    call perot ( u1, 1, kmax )
    !
    ! compute magnitude of shear squared in case of vertical mixing
    !
    do k = 1, kmax-1
       !
       do icell = 1, ncells
          !
          if ( hs(icell) > hmin ) then
             !
             dudz = 2. * ( uvc(icell,k,1) - uvc(icell,k+1,1) ) / ( hks(icell,k) + hks(icell,k+1) )
             !
             dvdz = 2. * ( uvc(icell,k,2) - uvc(icell,k+1,2) ) / ( hks(icell,k) + hks(icell,k+1) )
             !
             zshear(icell,k) = dudz*dudz + dvdz*dvdz
             !
          endif
          !
       enddo
       !
    enddo
    !
    ! loop over turbulence quantities
    !
    do l = 1, ltur
       !
       if ( l == 1 ) psm = sigmak
       if ( l == 2 ) psm = sigmae
       !
       if ( irough == 4 .and. l == 1 ) then
          kend = kmax
       else
          kend = kmax-1
       endif
       !
       ! build time-derivative
       !
       do k = 1, kend
          !
          do icell = 1, ncells
             !
             if ( hs(icell) > hmin ) then
                !
                ! time derivative
                !
                amatt(icell,k,1) = 1. / dt
                amatt(icell,k,2) = 0.
                amatt(icell,k,3) = 0.
                rhst (icell,k  ) = rtur(icell,k,l) / dt
                !
             else
                !
                rhst(icell,k) = rtur(icell,k,l)
                !
             endif
             !
          enddo
          !
       enddo
       !
       ! build vertical transport (implicit)
       !
       do k = 1, kmax-1
          !
          do icell = 1, ncells
             !
             if ( hs(icell) > hmin ) then
                !
                ! vertical advective term (first order upwinding)
                !
                w = wom(icell,k)
                !
                ctrkb = 0.5 * (w + abs(w)) / hks(icell,k+1)
                ctrkt = 0.5 * (w - abs(w)) / hks(icell,k  )
                !
                if ( k == kmax-1 ) ctrkb = 0.
                !
                amatt(icell,k,1) = amatt(icell,k,1) + ctrkb - ctrkt
                amatt(icell,k,2) =  ctrkt
                amatt(icell,k,3) = -ctrkb
                !
                ! vertical diffusive term
                !
                dz = 0.5 * ( hks(icell,k) + hks(icell,k+1) )
                !
                ctrkb = 0.5 * ( vnu3d(icell,k  ) + vnu3d(icell,k+1) ) / ( psm * dz * hks(icell,k+1) )
                ctrkt = 0.5 * ( vnu3d(icell,k-1) + vnu3d(icell,k  ) ) / ( psm * dz * hks(icell,k  ) )
                !
                amatt(icell,k,1) = amatt(icell,k,1) + ctrkb + ctrkt
                amatt(icell,k,2) = amatt(icell,k,2) - ctrkt
                amatt(icell,k,3) = amatt(icell,k,3) - ctrkb
                !
             endif
             !
          enddo
          !
       enddo
       !
       ! build source and sink terms
       !
       if ( l == 1 ) then
          !
          do k = 1, kmax-1
             !
             do icell = 1, ncells
                !
                if ( hs(icell) > hmin ) then
                   !
                   ! production rate term
                   !
                   rhst(icell,k) = rhst(icell,k) + vnu3d(icell,k) * zshear(icell,k)
                   !
                   ! dissipation rate term (Newton linearization)
                   !
                   amatt(icell,k,1) = amatt(icell,k,1) + 2. * rtur(icell,k,2) / ( rtur(icell,k,1)+eps )
                   rhst (icell,k  ) = rhst (icell,k  ) + rtur(icell,k,2)
                   !
                endif
                !
             enddo
             !
          enddo
          !
          if ( idens /= 0 ) then
             !
             ! buoyancy term
             !
             do k = 1, kmax-1
                !
                do icell = 1, ncells
                   !
                   if ( hs(icell) > hmin ) then
                      !
                      drhodz = 2. * ( rho(icell,k) - rho(icell,k+1) ) / ( hks(icell,k) + hks(icell,k+1) )
                      !
                      buotke = vnu3d(icell,k) * grav * drhodz / ( rhow*sigrho )
                      !
                      if ( buotke > 0. ) then
                         !
                         rhst(icell,k) = rhst(icell,k) + buotke
                         !
                      else
                         !
                         amatt(icell,k,1) = amatt(icell,k,1) - buotke / ( rtur(icell,k,1)+eps )
                         !
                      endif
                      !
                   endif
                   !
                enddo
                !
             enddo
             !
          endif
          !
       else if ( l == 2 ) then
          !
          do k = 1, kmax-1
             !
             do icell = 1, ncells
                !
                if ( hs(icell) > hmin ) then
                   !
                   ! production rate term
                   !
                   rhst(icell,k) = rhst(icell,k) + ceps1 * cmu * rtur(icell,k,1) * zshear(icell,k)
                   !
                   ! dissipation rate term (Newton linearization)
                   !
                   fac = ceps2 * rtur(icell,k,2) / ( rtur(icell,k,1)+eps )
                   !
                   amatt(icell,k,1) = amatt(icell,k,1) + 2. * fac
                   rhst (icell,k  ) = rhst (icell,k  ) + fac * rtur(icell,k,2)
                   !
                endif
                !
             enddo
             !
          enddo
          !
          if ( idens /= 0 ) then
             !
             ! buoyancy term
             !
             do k = 1, kmax-1
                !
                do icell = 1, ncells
                   !
                   if ( hs(icell) > hmin ) then
                      !
                      drhodz = 2. * ( rho(icell,k) - rho(icell,k+1) ) / ( hks(icell,k) + hks(icell,k+1) )
                      !
                      if ( drhodz > 0. ) then
                         !
                         rhst(icell,k) = rhst(icell,k) + ceps1 * cmu * rtur(icell,k,1) * grav * drhodz / ( rhow*sigrho )
                         !
                      endif
                      !
                   endif
                   !
                enddo
                !
             enddo
             !
          endif
          !
       endif
       !
       ! free surface
       !
       if ( l == 1 ) then
          !
          do icell = 1, ncells
             !
             if ( hs(icell) > hmin ) then
                !
                amatt(icell,0,1) =  1. + rtur(icell,2,1) / ( rtur(icell,0,1)+eps )
                amatt(icell,0,2) =  0.
                amatt(icell,0,3) = -2.
                rhst (icell,0  ) =  0.
                !
             else
                !
                rhst(icell,0) = rtur(icell,0,1)
                !
             endif
             !
          enddo
          !
       else if ( l == 2 ) then
          !
          do icell = 1, ncells
             !
             if ( hs(icell) > hmin ) then
                !
                zs = 0.5 * ( zks(icell,0) - zks(icell,1) )
                !
                amatt(icell,0,1) = 1.
                rhst (icell,0  ) = (cmu**.75) * (rtur(icell,0,1)**1.5) / ( vonkar*zs )
                !
             else
                !
                rhst(icell,0) = rtur(icell,0,2)
                !
             endif
             !
          enddo
          !
       endif
       !
       ! bottom
       !
       if ( l == 1 ) then
          !
          do icell = 1, ncells
             !
             if ( hs(icell) > hmin ) then
                !
                ! generation of turbulence due to bottom friction
                !
                if ( irough /= 4 ) then
                   !
                   ustar2 = cfricu(icell) * ( work(icell,1)*work(icell,1) + work(icell,2)*work(icell,2) )
                   !
                   amatt(icell,kmax,1) = 1.
                   rhst (icell,kmax  ) = ustar2 / sqrt(cmu)
                   !
                else
                   !
                   zs = 0.5 * ( zks(icell,kmax-1) - zks(icell,kmax) )
                   !
                   ! production rate term
                   !
                   rhst(icell,kmax) = rhst(icell,kmax) + logfrc(icell,1) * (uvc(icell,kmax,1)**2 + uvc(icell,kmax,2)**2) / zs
                   !
                   ! dissipation rate term
                   !
                   amatt(icell,kmax,1) = amatt(icell,kmax,1) + (cmu**.75) * sqrt(rtur(icell,kmax,1)) * logfrc(icell,2) / zs
                   !
                endif
                !
             else
                !
                rhst(icell,kmax) = rtur(icell,kmax,1)
                !
             endif
             !
          enddo
          !
       else if ( l == 2 ) then
          !
          do icell = 1, ncells
             !
             if ( hs(icell) > hmin ) then
                !
                zs = 0.5 * ( zks(icell,kmax-1) - zks(icell,kmax) )
                !
                amatt(icell,kmax,1) = 1.
                rhst (icell,kmax  ) = (cmu**.75) * (rtur(icell,kmax,1)**1.5) / ( vonkar*zs )
                !
             else
                !
                rhst(icell,kmax) = rtur(icell,kmax,2)
                !
             endif
             !
          enddo
          !
       endif
       !
       ! solve the equation
       !
       do icell = 1, ncells
          !
          bi = 1./amatt(icell,0,1)
          !
          amatt(icell,0,1) = bi
          amatt(icell,0,3) = amatt(icell,0,3)*bi
          rhst (icell,0  ) = rhst (icell,0  )*bi
          !
          do k = 1, kmax
             !
             bi = 1./(amatt(icell,k,1) - amatt(icell,k,2)*amatt(icell,k-1,3))
             amatt(icell,k,1) = bi
             amatt(icell,k,3) = amatt(icell,k,3)*bi
             rhst (icell,k  ) = (rhst(icell,k) - amatt(icell,k,2)*rhst(icell,k-1))*bi
             !
          enddo
          !
          rtur(icell,kmax,l) = rhst(icell,kmax)
          do k = kmax-1, 0, -1
             rtur(icell,k,l) = rhst(icell,k) - amatt(icell,k,3)*rtur(icell,k+1,l)
          enddo
          !
       enddo
       !
       ! if negative values occur, apply clipping or give warning
       !
       if ( ITEST < 30 ) then
          !
          do k = 0, kmax
             do icell = 1, ncells
                if ( rtur(icell,k,l) < 0. ) rtur(icell,k,l) = eps
             enddo
          enddo
          !
       else
          !
          do k = 0, kmax
             do icell = 1, ncells
                if ( rtur(icell,k,l) < 0. ) then
                   write (PRINTF,'(a,i1)') ' ++ keps: negative value found for l=',l
                   write (PRINTF,'(a,i7)') '          in cell   = ',icell
                   write (PRINTF,'(a,i3)') '          and layer = ',k
                   write (PRINTF,'(a,f14.8)') '          value     = ',rtur(icell,k,l)
                   rtur(icell,k,l) = eps
                endif
             enddo
          enddo
          !
       endif
       !
    enddo
    !
end subroutine SwashUKepsMod
