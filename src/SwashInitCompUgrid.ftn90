subroutine SwashInitCompUgrid ( logcom )
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering and Geosciences              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmer: Marcel Zijlema                                |
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
!    1.00, February 2010: New subroutine
!
!   Purpose
!
!   Initialises arrays for description of computational grid in case of unstructured grid
!
!   Modules used
!
    use ocpcomm4
    use SwashCommdata2
    use SwashCommdata3
    use m_genarr
    use m_parall
    use SwashFlowdata
    use SwanGriddata
!
    implicit none
!
!   Argument variables
!
    logical, dimension(6), intent(inout) :: logcom   ! indicates which commands have been given to know if all the
                                                     ! information for a certain command is available. Meaning:
                                                     ! (1) no meaning
                                                     ! (2) command CGRID has been carried out
                                                     ! (3) command READINP BOTTOM has been carried out
                                                     ! (4) command READ COOR has been carried out
                                                     ! (5) command READ UNSTRUC has been carried out
                                                     ! (6) arrays s1, u1 and v1 have been allocated
!
!   Local variables
!
    integer       :: i        ! loop counter
    integer, save :: ient = 0 ! number of entries in this subroutine
    integer       :: istat    ! indicate status of allocation
    !
    character(80) :: msgstr   ! string to pass message
!
!   Structure
!
!   Description of the pseudo code
!
!   Source text
!
    call strace (ient,'SwashInitCompUgrid')
    !
    ! compute coordinate offsets and reset grid coordinates
    !
    do i = 1, nverts
       if ( .not.lxoffs ) then
          xoffs = xcugrd(i)
          yoffs = ycugrd(i)
          lxoffs = .true.
          xcugrd(i) = 0.
          ycugrd(i) = 0.
       else
          xcugrd(i) = real(xcugrd(i) - dble(xoffs))
          ycugrd(i) = real(ycugrd(i) - dble(yoffs))
       endif
    enddo
    !
    ! check the grid
    !
    call SwashCheckGrid
    !
    ! compute xcgmin, xcgmax, ycgmin, ycgmax
    !
    xcgmin =  1.e9
    ycgmin =  1.e9
    xcgmax = -1.e9
    ycgmax = -1.e9
    do i = 1, nverts
       if ( xcugrd(i) < xcgmin ) xcgmin = xcugrd(i)
       if ( ycugrd(i) < ycgmin ) ycgmin = ycugrd(i)
       if ( xcugrd(i) > xcgmax ) xcgmax = xcugrd(i)
       if ( ycugrd(i) > ycgmax ) ycgmax = ycugrd(i)
    enddo
    !
    ! compute lengths of enclosure of computational domain
    !
    xclen = xcgmax - xcgmin
    yclen = ycgmax - ycgmin
    !
    istat = 0
    if ( .not.allocated(s1) ) allocate(s1(ncells), stat = istat)
    if ( istat == 0 ) then
       if ( .not.allocated(u1) ) allocate (u1(nfaces,kmax), stat = istat)
    endif
    if ( istat /= 0 ) then
       write (msgstr, '(a,i6)') 'allocation problem: array s1 or u1 and return code is ',istat
       call msgerr ( 4, trim(msgstr) )
       return
    endif
    s1 = 0.
    u1 = 0.
    logcom(6) = .true.
    !
    if ( kmax > 1 ) then
       !
       if ( .not.allocated(udep) ) allocate(udep(nfaces), stat = istat)
       if ( istat /= 0 ) then
          write (msgstr, '(a,i6)') 'allocation problem: array udep and return code is ',istat
          call msgerr ( 4, trim(msgstr) )
          return
       endif
       !
       udep = 0.
       !
    endif
    !
    ! set number of vertices and cells in global domain in case of serial run
    !
    if ( nvertsg == 0 ) then
       nvertsg = nverts
       ncellsg = ncells
    endif
    !
    ! the following arrays for regular grids (rectilinear and curvilinear) is allocated as empty one
    !
    if ( .not.allocated(KGRPGL) ) allocate(KGRPGL(0,0))
    if ( .not.allocated(kgrpnt) ) allocate(kgrpnt(0,0))
    if ( .not.allocated(XGRDGL) ) allocate(XGRDGL(0,0))
    if ( .not.allocated(YGRDGL) ) allocate(YGRDGL(0,0))
    !
    ! for sake of convenience, set mcgrd to nverts (for allocation of several arrays)
    !
    mcgrd   = nverts
    MCGRDGL = nvertsg
    !
end subroutine SwashInitCompUgrid
