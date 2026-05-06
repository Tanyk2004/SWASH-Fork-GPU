subroutine SwashMain
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
!    1.00, December 2009: New subroutine
!
!   Purpose
!
!   Main subroutine, it carries out the following steps:
!
!   performs initialization of the simulation process
!   performs pre processing step
!   performs simulation run
!   performs post processing step
!
!   Modules used
!
    use ocpcomm2
    use ocpcomm4
    use SwashCommdata1
    use SwashCommdata3
    use SwashTimeComm
    use outp_data
    use m_bndspec
    use m_parall
!
    implicit none
!
!   Local variables
!
    integer                         :: i      ! loop counter
    integer, dimension(2)           :: if     ! first character in string
    integer, dimension(2)           :: il     ! last character in string
    integer                         :: inerr  ! number of the initialisation error
    integer                         :: iostat ! I/O status in call FOR
    integer                         :: iprog  ! indicates the progress of the simulation
    integer                         :: istep  ! counted number of integration steps since the last doubling/halving of the time step
    integer                         :: it     ! integration step counter
    integer                         :: it0    ! start counter for time integration
    integer                         :: nds    ! file reference number
    !
    real*8                          :: fac    ! auxiliary factor
    !
    real, dimension(:), allocatable :: blknd  ! array giving node number per subdomain
    real, dimension(:), allocatable :: blkndc ! auxiliary array for collecting the array blknd
    !
    logical                         :: lopen  ! indicates whether a file is open or not
    logical                         :: STPNOW ! indicates that program must stop
    !
    character(20), dimension(2)     :: chars  ! string containing an integer
    character(4)                    :: comput ! determines the sort of computation to be performed by Swash
                                              ! ='COMP'; computation requested
                                              ! ='NOCO'; no computation but output requested
                                              ! ='STOP'; make computation, output and stop
    character(18)                   :: DTTIWR ! to write a time string
    character(20)                   :: INTSTR ! converts integer to string
    character(80)                   :: msgstr ! string to pass message to call msgerr
!
!   Structure
!
!   ---------------------------------------------------------------------------
!   call SwashInit to initialize various common data
!   Repeat
!       call SwashReadInput to read and process user commands
!       if last command was STOP
!       then exit from repeat
!       -----------------------------------------------------------------------
!       call SwashCheckPrep to check input and prepare computation
!       if nonstationary computation is to be made
!       then start time step loop at it = 0 and
!       -----------------------------------------------------------------------
!       For requested number of time steps do
!           call SwashUpdateData to update boundary conditions and input fields
!           if it > 0
!           then call SwashCompute to solve the shallow water equations
!           -------------------------------------------------------------------
!           call SwashOutput to postprocess the results and create output
!           Update time
!   ---------------------------------------------------------------------------
!
!   Source text
    !
    ! initialize various data
!TIMG    !
!TIMG    dcumtm(:,1:2) = 0d0
!TIMG    ncumtm(:)     = 0
!TIMG    call SWTSTA(1)
    !
    LEVERR = 0
    MAXERR = 1
    ITRACE = 0
    inerr  = 0
    !
!TIMG    call SWTSTA(2)
    call SwashInit ( inerr )
!TIMG    call SWTSTO(2)
    if ( inerr > 0 ) return
    if (STPNOW()) return
    !
    comput = '    '
    !
    ! repeat
    !
    repeatloop: do
       !
       ! read and process user-defined commands
       !
!TIMG       call SWTSTA(3)
       call SwashReadInput ( comput )
!TIMG       call SWTSTO(3)
       if (STPNOW()) return
       !
       ! if last command was STOP then exit from repeat
       !
       if ( comput == 'STOP' ) then
          nds    = 0
          iostat = 0
          FILENM = 'norm_end'
          call FOR (nds, FILENM, 'UF', iostat)
          write (nds, *) ' Normal end of run ', PROJNR
          exit repeatloop
       endif
       !
       ! check consistencies and do some preparations before computation
       !
!TIMG       call SWTSTA(4)
       call SwashCheckPrep
!TIMG       call SWTSTO(4)
       if (STPNOW()) return
       !
       if ( LEVERR > MAXERR ) then
          !
          write (PRINTF, 201) LEVERR
          if ( LEVERR < 4 ) write (PRINTF, 202)
          !
       else
          !
          if ( ITEST >= 40 ) then
             if ( nstatc == 1 ) then
                write (PRINTF, '(" Type of computation: dynamic")')
             else
                if ( oned ) then
                   write (PRINTF, '(" Type of computation: static 1D")')
                else
                   write (PRINTF, '(" Type of computation: static 2D")')
                endif
             endif
          endif
          !
          ! show all the settings
          !
          if ( ITEST > 0 ) call SwashPrintSettings
          !
          ! synchronize nodes
          !
          call SWSYNC
          if (STPNOW()) return
          !
          ! loop over time steps
          !
          if ( nstatc == 1 ) then
             it0   = 0
          else
             it0   = 1
          endif
          !
          it    = it0
          istep = 0
          !
          timeloop: do
             !
             if ( real(timco) > real(tfinc) ) exit timeloop
             !
             if ( LEVERR > MAXERR ) then
                write (PRINTF, 203) LEVERR
                if ( LEVERR < 4 ) write (PRINTF, 202)
                exit timeloop
             endif
             !
             ! synchronize nodes
             !
             call SWSYNC
             if (STPNOW()) return
             !
             ! update boundary values and input fields
             !
!TIMG             call SWTSTA(7)
             if ( optg /= 5 ) then
                call SwashUpdateData ( it )
             else
                call SwashUpdateUData ( it )
             endif
!TIMG             call SWTSTO(7)
             if (STPNOW()) return
             !
             ! synchronize nodes
             !
             call SWSYNC
             if (STPNOW()) return
             !
             if ( comput /= 'NOCO' .and. it > 0 ) then
                !
                ! perform simulation at current time step
                !
!TIMG                call SWTSTA(8)
                if ( optg /= 5 ) then
                   !
                   ! structured grid
                   !
                   call SwashComputStruc
                   !
                else
                   !
                   ! unstructured grid
                   !
                   call SwashCompUnstruc
                   !
                endif
!TIMG                call SWTSTO(8)
                if (STPNOW()) return
                !
             endif
             !
             ! synchronize nodes
             !
             call SWSYNC
             if (STPNOW()) return
             !
             ! carry out the output requests
             !
!TIMG             call SWTSTA(9)
             call SwashOutput
!TIMG             call SWTSTO(9)
             if (STPNOW()) return
             !
             if ( errpts > 0 ) rewind ( errpts )
             !
             ! update time
             !
             if ( nstatc == 1 ) then
                !
                ! change time step in case of explicit time integration, if appropriate
                !
                if ( mtimei == 1 ) then
                   !
                   if ( cflmax > pnums(3) ) then
                      !
                      dt = 0.5*dt
                      write (msgstr,'(a,f12.5,a)') 'time step is halved! New time step: ',dt, ' sec'
                      call msgerr ( 1, trim(msgstr) )
                      !
                      istep = 0
                      !
                   else if ( cflmax < pnums(2) ) then
                      !
                      if ( .not. istep < 20 ) then
                         !
                         ! only increase actual time step if the time step has been constant for the last 20 integration steps
                         !
                         dt = 2.0*dt
                         write (msgstr,'(a,f12.5,a)') 'time step is doubled! New time step: ',dt,' sec'
                         call msgerr ( 1, trim(msgstr) )
                         !
                         istep = 0
                         !
                      endif
                      !
                   endif
                   !
                endif
                !
                timco = timco + dt
                chtime = DTTIWR(ITMOPT, timco)
                if ( .not. real(timco) > real(tfinc) ) then
                   if ( ITEST >= 1 ) write (PRINTF, 204) chtime, timco
                   if (SCREEN /= PRINTF .and. INODE == MASTER ) then
                      fac   = 100. * (timco-tinic) / (tfinc-tinic)
                      iprog = nint(fac)
                      chars(1) = INTSTR(iprog)
                      call TXPBLA(chars(1),if(1),il(1))
                      chars(2) = INTSTR(it+1)
                      call TXPBLA(chars(2),if(2),il(2))
                      write (SCREEN, 205) chars(1)(if(1):il(1)), chars(2)(if(2):il(2))
                   endif
                endif
                !
             else
                !
                timco = timco + dt
                !
             endif
             !
             it    = it    + 1
             istep = istep + 1
             !
          enddo timeloop
          !
          if ( LEVERR > MAXERR ) exit repeatloop
          !
       endif
       !
       ! delete linked list of Fourier series required for multiple COMPUTE
       !
       call delete ( fbfs )
       nullify( cubfs )
       lbfs = .false.
       !
    enddo repeatloop
    !
!TIMG    call SWTSTO(1)
    !
    ! finalize PVD collection files
    !
    do i = 1, max_outp_req
       nds = upvdf(i)
       if ( nds > 0 ) then
          write (nds,'(a)') trim(pvdlin3)
          write (nds,'(a)') trim(pvdlin4)
       endif
    enddo
    !
    ! close all open files except PRINT file
    !
    do i = 1, HIOPEN
       inquire ( unit=i, opened=lopen )
       if ( lopen .and. i /= PRINTF ) close (i)
    enddo
    !
    ! collect contents of individual process files for output requests in case of parallel computation
    !
    call SWSYNC
!TIMG    call SWTSTA(9)
    if ( PARLL .AND. IMRGE.EQ.1 ) then
       allocate (blknd(mxc*myc))
       blknd = real(INODE)
       if ( INODE == MASTER ) then
          allocate(blkndc(MXCGL*MYCGL))
          blkndc = 0.
       endif
       call SWCOLLECT ( blkndc, blknd )
       if (STPNOW()) return
       if ( INODE == MASTER ) then
          call SWCOLOUT ( blkndc )
          deallocate(blkndc)
       endif
    endif
!TIMG    call SWTSTO(9)
!
!TIMG    call SWPRTI
    !
    ! close PRINT file
    !
    inquire ( unit=PRINTF, opened=lopen )
    if ( lopen ) close (PRINTF)
    !
    ! deallocate all allocated arrays and linked lists
    !
    if ( allocated(blknd) ) deallocate(blknd)
    call SwashCleanMem
    !
 201 format (' ** no start of computation because of error level:', i3)
 202 format (' ** to ignore this error, change [maxerr] with the SET command')
 203 format (' ** no continuation of computation because of error level:',i3)
 204 format (' Time of simulation  ->  ',a18,' in sec:', f17.5)
 205 format ('[',a3,'%]   step ',a)
    !
end subroutine SwashMain
