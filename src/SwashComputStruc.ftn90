subroutine SwashComputStruc
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
!    9.01: Dirk Rijnsdorp and Marcel Zijlema
!
!   Updates
!
!    1.00,   March 2010: New subroutine
!    9.01, October 2022: extension moving rigid bodies
!
!   Purpose
!
!   Performs one full simulation step with structured grid
!
!   Modules used
!
    use ocpcomm4
    use SwashCommdata3
    use m_parall
    use SwashFlowdata
    use SwashSolvedata, only: iamout
    use SwashRigBoddata
!
    implicit none
!
!   Local variables
!
    integer, save         :: ient = 0 ! number of entries in this subroutine
    integer               :: j        ! iteration counter
    integer               :: m        ! loop counter
    integer               :: maxit    ! maximum number of iterations
    integer, dimension(2) :: minf     ! m-th body with largest error
                                      ! minf(1): linear displacement
                                      ! minf(2): angular displacement
    integer               :: n        ! loop counter
    integer, dimension(2) :: ninf     ! n-th degree of freedom with largest error
                                      ! ninf(1): linear displacement
                                      ! ninf(2): angular displacement
    !
    real, dimension(2)    :: epslin   ! required accuracy. Meaning:
                                      ! epslin(1): linear displacement
                                      ! epslin(2): angular displacement
    real                  :: reps     ! accuracy of the final approximation
    real                  :: res      ! residual
    real, dimension(2)    :: resm     ! maximum error. Meaning:
                                      ! resm(1): linear displacement
                                      ! resm(2): angular displacement
    real                  :: xformx   ! maximum value of angular displacement of rigid bodies
    real                  :: xfotmx   ! maximum value of linear displacement of rigid bodies
    !
    logical               :: STPNOW   ! indicates that program must stop
!
!   Structure
!
!   Description of the pseudo code
!
!   Source text
!
    if (ltrace) call strace (ient,'SwashComputStruc')
    !
    ! determine maximum displacements and required accuracy
    !
    xfotmx    = maxval(abs(xfot0))
    xformx    = maxval(abs(xfor0))
    reps      = pship(8)
    epslin(1) = max(reps,reps*xfotmx)
    epslin(2) = max(reps,reps*xformx)
    !
    ! initialize
    !
    maxit   = nint(pship(9))
    j       = 0
    resm(1) = epslin(1) + 1.
    resm(2) = epslin(2) + 1.
    !
    ! start iteration process to complete the fluid-structure interaction (FSI)
    !
 10 if ( ( resm(1) > epslin(1) .or. resm(2) > epslin(2) ) .and. j < maxit ) then
       !
       j = j + 1
       !
       resm = 0.
       minf = 0
       ninf = 0
       !
       ! update kinematic boundary conditions at rigid bodies
       !
!TIMG       call SWTSTA(61)
       call SwashUpdKBCrigb
!TIMG       call SWTSTO(61)
       if (STPNOW()) return
       !
       ! compute flow
       !
!TIMG       call SWTSTA(51)
       if ( .not.momskip ) call SwashComputFlow
!TIMG       call SWTSTO(51)
       if (STPNOW()) return
       !
       ! compute forces and torques acting on rigid bodies
       !
!TIMG       call SWTSTA(62)
       call SwashForcesRigidBod
!TIMG       call SWTSTO(62)
       if (STPNOW()) return
       !
       ! set intermediate displacements for convergence check
       !
       xfoti = xfot1
       xfori = xfor1
       !
       ! compute motion of rigid bodies
       !
!TIMG       call SWTSTA(63)
       call SwashMotionRigidBod
!TIMG       call SWTSTO(63)
       if (STPNOW()) return
       !
       ! determine maximum error to check convergence of FSI coupling
       !
       do m = 1, mbod
          do n = 1, ndim
             !
             res = abs(xfot1(m,n) - xfoti(m,n))
             if ( res > resm(1) ) then
                resm(1) = res
                minf(1) = m
                ninf(1) = n
             endif
             !
             res = abs(xfor1(m,n) - xfori(m,n))
             if ( res > resm(2) ) then
                resm(2) = res
                minf(2) = m
                ninf(2) = n
             endif
             !
          enddo
       enddo
       !
       if ( iamout == 2 .and. INODE == MASTER .and. maxit > 1 ) then
          !
          if ( ninf(1) == 1 ) write (PRINTF,'(a,i4,a,e12.6,a,i2)') ' ++ FSI coupling: iter = ',j,' res = ',resm(1),' m of surged body no. ',minf(1)
          if ( ninf(1) == 2 ) write (PRINTF,'(a,i4,a,e12.6,a,i2)') ' ++ FSI coupling: iter = ',j,' res = ',resm(1),' m of swayed body no. ',minf(1)
          if ( ninf(1) == 3 ) write (PRINTF,'(a,i4,a,e12.6,a,i2)') ' ++ FSI coupling: iter = ',j,' res = ',resm(1),' m of heaved body no. ',minf(1)
          if ( ninf(2) == 1 ) write (PRINTF,'(a,i4,a,e12.6,a,i2)') ' ++ FSI coupling: iter = ',j,' res = ',resm(2)*180./pi,' deg of rolled body no. ' ,minf(2)
          if ( ninf(2) == 2 ) write (PRINTF,'(a,i4,a,e12.6,a,i2)') ' ++ FSI coupling: iter = ',j,' res = ',resm(2)*180./pi,' deg of pitched body no. ',minf(2)
          if ( ninf(2) == 3 ) write (PRINTF,'(a,i4,a,e12.6,a,i2)') ' ++ FSI coupling: iter = ',j,' res = ',resm(2)*180./pi,' deg of yawed body no. '  ,minf(2)
          !
       endif
       !
       goto 10
       !
    endif
    !
    ! check the reason for stopping, if appropriate
    !
    if ( ifloat == 2 .and. INODE == MASTER ) then
       !
       if ( resm(1) > epslin(1) .and. iamout > 0 ) then
          !
          write (PRINTF, '(a)') ' ++ FSI coupling: no convergence in linear displacement'
          write (PRINTF, '(a,i3)'   ) '                  total number of iterations = ',j
          write (PRINTF, '(a,e12.6)') '                  2-norm of the residual     = ',resm(1)
          write (PRINTF, '(a,e12.6)') '                  required accuracy          = ',epslin(1)
          !
       else if ( resm(2) > epslin(2) .and. iamout > 0 ) then
          !
          write (PRINTF, '(a)') ' ++ FSI coupling: no convergence in angular displacement'
          write (PRINTF, '(a,i3)'   ) '                  total number of iterations = ',j
          write (PRINTF, '(a,e12.6)') '                  2-norm of the residual     = ',resm(2)
          write (PRINTF, '(a,e12.6)') '                  required accuracy          = ',epslin(2)
          !
       else if ( iamout == 3 ) then
          !
          write (PRINTF, '(a,i3)'   ) ' ++ FSI coupling: total number of iterations            = ',j
          write (PRINTF, '(a,e12.6)') '                  maximum error in linear displacement  = ',resm(1)
          write (PRINTF, '(a,e12.6)') '                  maximum error in angular displacement = ',resm(2)
          !
       endif
       !
    endif
    !
    ! update water depths
    !
!TIMG    call SWTSTA(52)
    if ( kmax == 1 ) then
       call SwashUpdateDepths (   u1,   v1 )
    else
       call SwashUpdateDepths ( udep, vdep )
    endif
!TIMG    call SWTSTO(52)
    if (STPNOW()) return
    !
    ! update layer interfaces
    !
!TIMG    call SWTSTA(53)
    if ( kmax > 1 ) call SwashLayerIntfaces
!TIMG    call SWTSTO(53)
    !
    ! update mask arrays for wetting and drying
    !
!TIMG    call SWTSTA(54)
    call SwashDryWet
!TIMG    call SWTSTO(54)
    if (STPNOW()) return
    !
    ! update mask arrays for free surface / pressurized flow
    !
!TIMG    call SWTSTA(54)
    if ( ifloat /= 0 ) call SwashPresFlow
!TIMG    call SWTSTO(54)
    if (STPNOW()) return
    !
    ! update mask array for wave breaking
    !
!TIMG    call SWTSTA(54)
    if ( isurf /= 0 ) call SwashBreakPoint
!TIMG    call SWTSTO(54)
    if (STPNOW()) return
    !
    ! compute transport
    !
!TIMG    call SWTSTA(101)
    if ( itrans /= 0 ) call SwashComputTrans
!TIMG    call SWTSTO(101)
    if (STPNOW()) return
    !
    ! compute 3D turbulence
    !
!TIMG    call SWTSTA(102)
    if ( iturb /= 0 ) call SwashComputTurb
!TIMG    call SWTSTO(102)
    if (STPNOW()) return
    !
end subroutine SwashComputStruc
