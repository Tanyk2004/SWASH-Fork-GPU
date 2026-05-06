subroutine SwashDecOutQ ( outi, outr, ivtyp, rtype, psname, nvar, nvoqp, nvoqk, oqproc, voqr, oqt, logact )
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
!    1.00, March 2010: New subroutine
!
!   Purpose
!
!   Decodes output quantities
!
!   Modules used
!
    use ocpcomm4
    use SwashCommdata1
    use SwashCommdata3
    use SwashCommdata4
    use SwashTimeComm
    use m_parall
    use outp_data, only: etorq, ntvtk
    use SwashRigBoddata, only: mbod
!
    implicit none
!
!   Argument variables
!
    integer, intent(in)                     :: nvar     ! number of output quantities
    !
    integer, dimension(4), intent(in)       :: outi     ! integer coefficients
    integer, dimension(nvar), intent(in)    :: ivtyp    ! type of output variable
    integer, intent(inout)                  :: nvoqk    ! number of layer-dependent quantities per output point
    integer, intent(inout)                  :: nvoqp    ! number of quantities per output point
    integer, dimension(nmovar), intent(out) :: voqr     ! place of each output quantity
    !
    real*8, intent(out)                     :: oqt      ! indicates at what time requested output is processed
    real*8, dimension(2), intent(inout)     :: outr     ! code for one output request
    !
    logical, intent(out)                    :: logact   ! indicates whether output action should be taken or not
    logical, dimension(nmovar), intent(out) :: oqproc   ! indicates whether or not an output quantity must be processed
    !
    character(8), intent(in)                :: psname   ! name of output point set
    character(4), intent(inout)             :: rtype    ! type of output
!
!   Local variables
!
    integer                                 :: i        ! loop counter
    integer, save                           :: ient = 0 ! number of entries in this subroutine
    integer                                 :: ivtype   ! type of output variable (local variable)
    !
    real*8                                  :: dif      ! interval between end time and current time
    real*8                                  :: tnext    ! time of next requested output
!
!   Structure
!
!   --------------------------------------------------------------------------------------------
!   if dynamic mode
!   then determine tnext (time of next requested output)
!        determine Dif
!        if Dif is less than half time step and output interval is negative
!        then enable output (by making logact = true)
!        else if time of computation >= tnext and output interval is positive
!        then enable output
!        else disable output (by making logact = false)
!             return
!        ---------------------------------------------------------------------------------------
!   else enable output
!   --------------------------------------------------------------------------------------------
!   Set all oqproc = false (if oqproc is true corresponding quantity must be computed)
!   Make oqproc true for quantities Xp, Yp, Xc and Yc
!   Set all values of voqr = 0 (voqr indicates where value of a quantity is stored in array voq)
!   Make voqr nonzero for quantities Xp, Yp, Xc and Yc
!   Assign value to nvar depending on type of output request
!   --------------------------------------------------------------------------------------------
!
!   Remark: output interval negative means that output is made only at end of computation
!
!   Source text
!
    if (ltrace) call strace (ient,'SwashDecOutQ')
    !
    ! check time of output action
    !
    if ( nstatm == 1 ) then
       !
       dif = tfinc - timco
       if ( outr(1) < tinic ) then
          tnext = tinic
       else
          tnext = outr(1)
       endif
       if ( PARLL .and. oqt == -9999. ) oqt = outr(1)
       if ( PARLL ) then
          etorq(1) = timco
          if ( dt < etorq(2) ) etorq(2) = dt
       endif
       if ( ITEST >= 60 ) write (PRTEST, *) ' output times ', tnext, outr(2), dt, tfinc, timco
       if ( abs(dif) < 0.99*dt .and. outr(2) < 0. ) then
          outr(1) = timco
          logact  = .true.
       else if ( outr(2) > 0. .and. .not. real(timco) < real(tnext) ) then
          outr(1) = tnext + outr(2)
          logact  = .true.
       else
          logact = .false.
          return
       endif
       if ( logact ) ntvtk(outi(2)) = ntvtk(outi(2)) + 1
    else
       logact = .true.
    endif
    !
    ! action is taken, proceed with analysing output request
    !
    ! initialize voqr and oqproc
    !
    voqr   = 0
    oqproc = .false.
    !
    ! in case of multiple output quantities, reset ovsvty
    !
    do i = 1, nmovar
       if ( ovsvty(i) > 4 ) ovsvty(i) = 4
    enddo
    !
    ! ivtype 1 and 2 are Xp and Yp
    !
    oqproc(1) = .true.
    voqr(1)   = 1
    oqproc(2) = .true.
    voqr(2)   = 2
    !
    ! ivtype 25 and 26 are Xc and Yc
    !
    oqproc(25) = .true.
    voqr(25)   = 3
    oqproc(26) = .true.
    voqr(26)   = 4
    !
    nvoqp = 4
    !
    varloop: do i = 1, nvar
       !
       ivtype = ivtyp(i)
       !
       if ( ivtype < 1 .or. ivtype > nmovar ) then
          call msgerr (2, 'wrong value for ivtype')
          write (PRINTF, 201) rtype, psname, ivtype, nvar
          cycle varloop
       endif
       !
       if ( ovsvty(ivtype) <= 2 .and. .not.oqproc(ivtype) ) then
          ! output quantity is a scalar
          if ( ivtype < 51 .or. ivtype > 99 ) then
             nvoqp = nvoqp + 1
             voqr(ivtype) = nvoqp
          else
             nvoqk = nvoqk + 1
             voqr(ivtype) = nvoqk
          endif
       else if ( ovsvty(ivtype) == 3 .and. .not.oqproc(ivtype) ) then
          ! output quantity is a vector
          if ( ivtype < 51 .or. ivtype > 99 ) then
             nvoqp = nvoqp + 2
             voqr(ivtype) = nvoqp - 1
          else
             nvoqk = nvoqk + 2
             voqr(ivtype) = nvoqk - 1
          endif
       else if ( ovsvty(ivtype) == 4 .and. .not.oqproc(ivtype) ) then
          ! multiple output quantity
          if ( ivtype >= 101 .and. ivtype <= 113 ) then
             ! floating objects
             voqr(ivtype) = nvoqp + 1
             if ( mbod > 1 ) then
                nvoqp = nvoqp + mbod
                ovsvty(ivtype) = ovsvty(ivtype) + mbod
             else
                nvoqp = nvoqp + 1
                ovsvty(ivtype) = ovsvty(ivtype) + 1
             endif
          endif
       endif
       oqproc(ivtype) = .true.
       !
       if ( ITEST >= 80 .or. IOUTES >= 20 ) write (PRTEST, 202) i, ivtype, voqr(ivtype)
       !
       ! for velocity magnitude or direction add velocity components as output quantities
       !
       if ( ivtype == 7 .or. ivtype == 8 ) then
          if ( .not.oqproc(9) ) then
             nvoqp     = nvoqp + 2
             voqr(9)   = nvoqp - 1
             oqproc(9) = .true.
          endif
       endif
       !
       ! for magnitude or direction of specific discharge add components as output quantities
       !
       if ( ivtype == 14 .or. ivtype == 15 ) then
          if ( .not.oqproc(16) ) then
             nvoqp      = nvoqp + 2
             voqr(16)   = nvoqp - 1
             oqproc(16) = .true.
          endif
       endif
       !
       ! for wind velocity magnitude or direction add wind velocity components as output quantities
       !
       if ( ivtype == 27 .or. ivtype == 28 ) then
          if ( .not.oqproc(29) ) then
             nvoqp      = nvoqp + 2
             voqr(29)   = nvoqp - 1
             oqproc(29) = .true.
          endif
       endif
       !
       ! for time-averaged velocity magnitude or direction add velocity components as output quantities
       !
       if ( ivtype == 33 .or. ivtype == 34 ) then
          if ( .not.oqproc(35) ) then
             nvoqp      = nvoqp + 2
             voqr(35)   = nvoqp - 1
             oqproc(35) = .true.
          endif
       endif
       !
       ! for magnitude of friction velocity add components as output quantities
       !
       if ( ivtype == 31 ) then
          if ( .not.oqproc(32) ) then
             nvoqp      = nvoqp + 2
             voqr(32)   = nvoqp - 1
             oqproc(32) = .true.
          endif
       endif
       !
       ! for eddy viscosity add turbulent kinetic energy and dissipation rate as output quantities
       !
       if ( ivtype == 56 ) then
          if ( .not.oqproc(54) ) then
             nvoqk      = nvoqk + 1
             voqr(54)   = nvoqk
             oqproc(54) = .true.
          endif
          if ( .not.oqproc(55) ) then
             nvoqk      = nvoqk + 1
             voqr(55)   = nvoqk
             oqproc(55) = .true.
          endif
       endif
       !
       ! for mean eddy viscosity add mean turbulent kinetic energy and mean dissipation rate as output quantities
       !
       if ( ivtype == 59 ) then
          if ( .not.oqproc(57) ) then
             nvoqk      = nvoqk + 1
             voqr(57)   = nvoqk
             oqproc(57) = .true.
          endif
          if ( .not.oqproc(58) ) then
             nvoqk      = nvoqk + 1
             voqr(58)   = nvoqk
             oqproc(58) = .true.
          endif
       endif
       !
       ! for layer-dependent velocity magnitude or direction add velocity components as output quantities
       !
       if ( ivtype == 72 .or. ivtype == 73 ) then
          if ( .not.oqproc(74) ) then
             nvoqk      = nvoqk + 2
             voqr(74)   = nvoqk - 1
             oqproc(74) = .true.
          endif
       endif
       !
       ! for layer-dependent magnitude or direction of specific discharge add components as output quantities
       !
       if ( ivtype == 77 .or. ivtype == 78 ) then
          if ( .not.oqproc(79) ) then
             nvoqk      = nvoqk + 2
             voqr(79)   = nvoqk - 1
             oqproc(79) = .true.
          endif
       endif
       !
       ! for time-averaged, layer-dependent velocity magnitude or direction add velocity components as output quantities
       !
       if ( ivtype == 84 .or. ivtype == 85 ) then
          if ( .not.oqproc(86) ) then
             nvoqk      = nvoqk + 2
             voqr(86)   = nvoqk - 1
             oqproc(86) = .true.
          endif
       endif
       !
       ! change type of output in case of table format and layer-dependent quantities
       !
       if ( rtype(1:3) == 'TAB' .and. ivtype > 50 .and. ivtype < 100 ) then
          !
          if ( rtype == 'TABS' ) then
             rtype = 'TABU'
          else if ( rtype /= 'TABV' .and. rtype /= 'TABU' ) then
             if ( rtype /= 'TABP' ) call msgerr (1, 'type of output has been changed into STAB for layer-dependent quantities' )
             rtype = 'TABV'
          endif
          !
       endif
       !
    enddo varloop
    !
 201 format (' type, points, var: ', a4, 2x, a8, 2x, 2i8)
 202 format (' SwashDecOutQ: output quantity:', 3i6)
    !
end subroutine SwashDecOutQ
