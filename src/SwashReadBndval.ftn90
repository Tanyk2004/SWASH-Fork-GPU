subroutine SwashReadBndval ( bfiled, bctime, bcloc, bndval )
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
!    1.00, February 2010: New subroutine
!
!   Purpose
!
!   Reads boundary values from boundary file(s) and store them
!
!   Modules used
!
    use ocpcomm2
    use ocpcomm4
    use SwashCommdata3
    use SwashTimecomm
!
    implicit none
!
!   Argument variables
!
    integer, dimension(20), intent(inout)   :: bfiled   ! data concerning boundary condition files
    integer, dimension(nbval), intent(in)   :: bcloc    ! place in array bndval where to store interpolated boundary values
    !
    real*8, dimension(2), intent(inout)     :: bctime   ! data containing times of boundary values to be read
    real, dimension(nbval,2), intent(inout) :: bndval   ! array containing boundary values at two time levels
!
!   Local variables
!
    integer                                 :: binfo    ! information type for boundary values
                                                        ! =0; time series
                                                        ! =1; 1D spectrum
                                                        ! =2; 2D spectrum
    integer                                 :: i        ! loop counter
    integer                                 :: ibloc    ! counter for boundary point
    integer, save                           :: ient = 0 ! number of entries in this subroutine
    integer                                 :: iiopt    ! time coding option
    integer                                 :: j        ! loop counter
    integer                                 :: ndsd     ! unit reference number for data file
    integer                                 :: ndsl     ! unit reference number for namelist file
    !
    real                                    :: bval     ! actual boundary value
    real*8                                  :: timf1    ! time of reading previous boundary value
    real*8                                  :: timf2    ! time of reading current boundary value
    real*8                                  :: wf1      ! weighting factor for interpolation in time
    real*8                                  :: wf2      ! weighting factor for interpolation in time
    !
    character(18), dimension(2)             :: datitm   ! string containing data (e.g. time and boundary value)
    character(80)                           :: hedlin   ! heading line
    character(80)                           :: msgstr   ! string to pass message
!
!   Structure
!
!   Description of the pseudo code
!
!   Source text
!
    if (ltrace) call strace (ient,'SwashReadBndval')
    !
    ! if data file is exhausted, return
    !
    if ( bfiled(1) == -1 ) return
    !
    timf1 = bctime(1)
    timf2 = bctime(2)
    iiopt = bfiled(6)
    binfo = bfiled(7)
    !
    if ( ITEST >= 120 ) write (PRINTF, 201) bfiled(1), binfo, timf1, timf2, timco
    !
    bcloop: do
       !
       ! if current time > time of last reading value, read new boundary value
       !
       if ( timco > timf2 ) then
          !
          ! store current boundary values to previous time level
          !
          timf1 = timf2
          !
          do i = 1, bfiled(8)
             ibloc = bcloc(i)
             bndval(ibloc,1) = bndval(ibloc,2)
          enddo
          !
          ndsl = bfiled(4)
          ndsd = bfiled(5)
          !
 100      if ( .not. bfiled(1) > 0 ) timf2 = 0.
          !
          ! read heading lines per time step
          !
          do i = 1, bfiled(15)
             read (ndsd, '(a)', end=110, err=910) hedlin
             if ( ITEST >= 90 ) write (PRINTF, 202) hedlin
          enddo
          !
          if ( ITEST >= 60 ) write (PRINTF, 203) timf1, timf2, timco, bfiled(8), bfiled(15)
          !
          ! read data of all boundary locations
          !
          do i = 1, bfiled(8)
             !
             ibloc = bcloc(i)
             !
             ! read heading lines per boundary value
             !
             do j = 1, bfiled(16)
                read (ndsd, '(a)', end=110, err=910) hedlin
                if ( ITEST >= 100 ) write (PRINTF, 202) hedlin
             enddo
             !
             ! start reading boundary value and store it
             !
             if ( binfo == 0 ) then
                !
                if ( iiopt == 0 ) then
                   !
                   read (ndsd, *, end=110) timf2, bval
                   !
                else
                   !
                   read (ndsd, '(a)', end=110) hedlin
                   call LSPLIT (hedlin, datitm, 2)
                   call DTRETI (datitm(1), iiopt, timf2)
                   read (datitm(2), '(g20.0)') bval
                   !
                endif
                !
                if ( ITEST >= 60 ) write (PRTEST, *) ' boundary value ', timf2, bval
                bctime(2) = timf2
                bndval(ibloc,2) = bval
                !
             endif
             !
          enddo
          !
          if ( ITEST >= 60 ) write (PRINTF, 204) binfo, timf2
          !
          ! cycle back to next reading
          !
          cycle bcloop
          !
          ! if there are no more data on a boundary data file close this file, and see if there is a next one
          !
 110      close (ndsd)
          !
          ! read name of next boundary file and open it
          !
          if ( ndsl > 0 ) then
             !
             read (ndsl, '(a)', end=120, err=920) FILENM
             open (ndsd, file=FILENM, form='FORMATTED', status='OLD', err=920)
             !
             do i = 1, bfiled(14)
                ! read heading lines
                read (ndsd, '(a)', end=930, err=910) hedlin
                if ( ITEST >= 80 ) write (PRINTF, 202) hedlin
             enddo
             !
             ! go back to start again reading data from the next boundary file
             !
             goto 100
             !
          else
             !
             ! boundary values are read from a single file
             goto 130
             !
          endif
          !
          ! close file containing filenames
          !
 120      close (ndsl)
          !
          bfiled(5) = 0
          !
          ! write message and close file containing boundary values
          !
 130      call msgerr (1, 'data on boundary file exhausted')
          !
          bfiled(4) = 0
          bfiled(1) = -1
          timf2     = 999999999.
          cycle bcloop
          !
       else
          !
          ! all necessary data have been read from file(s), now interpolate in time
          !
          if ( timf1 /= timf2 ) then
             wf1 = ( timf2 - timco ) / ( timf2 - timf1 )
          else
             wf1 = 0.
          endif
          wf2 = 1. - wf1
          !
          do i = 1, bfiled(8)
             ibloc = bcloc(i)
             if ( i == 1 .and. ITEST >= 80 ) write (PRTEST, 205) timco, wf1, timf1, timf2, ibloc
             !
             ! note: interpolated value must be stored in bndval(:,1) !!!
             !
             bndval(ibloc,1) = wf1 * bndval(ibloc,1) + wf2 * bndval(ibloc,2)
             !
          enddo
          bctime(1) = timco
          bctime(2) = timf2
          !
          exit bcloop
          !
       endif
       !
    enddo bcloop
    !
    return
    !
 201 format (' boundary', i2, 2x, i1, ' times: ', 3f10.1)
 202 format (' heading line: ', a)
 203 format (' boundary times ', 3f12.0, 2x, 2i4)
 204 format (' boundary data type ', i2, ' processed, time: ', f10.1)
 205 format (' interpolation in time ', f14.1, f8.3, 2f14.1, i4)
    !
 910 inquire ( unit=ndsd, name=FILENM )
    call msgerr (4, 'error reading data from boundary file '//FILENM)
    return
 920 call msgerr (4, 'error opening boundary file '//FILENM)
    return
 930 inquire ( unit=ndsd, name=FILENM )
    call msgerr (4, 'unexpected end of file on boundary file '//FILENM)
    return
    !
end subroutine SwashReadBndval
