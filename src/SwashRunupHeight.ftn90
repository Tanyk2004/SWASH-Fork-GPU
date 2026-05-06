subroutine SwashRunupHeight ( oqproc, nvoqp, voqr, voq )
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
!
!   Updates
!
!    1.00,   April 2018: New subroutine
!
!   Purpose
!
!   Calculates wave runup height for the purpose of output
!
!   Method
!
!   Based on the intersection between free surface and bottom level + threshold for runup
!
!   Note: we restrict ourselves to 1D only and wave propagation is pointing eastward!
!
!   Modules used
!
    use ocpcomm4
    use SwashCommdata1
    use SwashCommdata3, only: oned
    use SwashTimeComm, only: timco
    use outp_data, only: delrp
    use m_genarr, only: kgrpnt
    use m_parall
    use SwashFlowdata
!
    implicit none
!
!   Argument variables
!
    integer, intent(in)                     :: nvoqp  ! number of quantities
    integer, dimension(nmovar), intent(in)  :: voqr   ! place of each output quantity
    real   , dimension(nvoqp) , intent(out) :: voq    ! output quantity at request
    logical, dimension(nmovar), intent(in)  :: oqproc ! indicates whether or not an output quantity must be processed
!
!   Local variables
!
    integer, save :: ient = 0 ! number of entries in this subroutine
    integer, save :: jerr = 0 ! actual number of printing error statement
    integer       :: m        ! loop counter
    integer       :: mr       ! index of first wet point for runup calculation
    integer       :: nm       ! pointer to m
    integer       :: nmr      ! pointer to mr
    integer       :: nmru     ! pointer to mr+1
    !
    real          :: rh       ! runup height
!
!   Structure
!
!   Description of the pseudo code
!
!   Source text
!
    if (ltrace) call strace (ient,'SwashRunupHeight')
    !
    rh = ovexcv(115)
    !
    if ( oned ) then
       !
       ! check left boundary condition
       !
       if ( ibl(1) /= 2 .and. ibl(1) /= 3 .and. ibl(1) /= 7 .and. LMXF ) then
          !
          if ( jerr < 1 ) then
             call msgerr (1, 'wave runup height will not be computed as no waves are imposed on west side')
             jerr = jerr + 1
          endif
          !
       endif
       !
       mloop: do m = ml, mfu, -1
          !
          if ( jerr > 0 ) then
             mr = MXCGL+1
             exit mloop
          endif
          !
          nm = kgrpnt(m,1)
          !
          if ( hs(nm) > delrp ) then
             mr = m
             exit mloop
          endif
          !
       enddo mloop
       !
       ! take the largest point index in whole domain for runup calculation
       !
       m = mr + MXF   ! from local to global
       !
       call SWREDUCE ( m, 1, SWINT, SWMAX )
       !
       mr = m - MXF   ! from global to local
       !
       ! compute vertical runup level
       !
       if ( mr > mf .and. mr < mlu ) then
          !
          nmr  = kgrpnt(mr  ,1)
          nmru = kgrpnt(mr+1,1)
          !
          rh = s1(nmr) + ( delrp - hs(nmr) ) * ( s1(nmru) - s1(nmr) ) / ( hs(nmru) - hs(nmr) )
          !
       else
          !
          rh = ovexcv(115)
          !
       endif
       !
       call SWREDUCE ( rh, 1, SWREAL, SWMAX )
       !
    endif
    !
    ! store the requested output quantities
    !
    ! Tsec
    !
    if ( oqproc(41) ) voq(voqr(41)) = real(timco) - outpar(1)
    !
    ! runup height
    !
    voq(voqr(115)) = rh
    !
end subroutine SwashRunupHeight
