!     SWAN - SERVICE ROUTINES
!
!  Contents of this file
!
!     READXY
!     REFIXY
!     CVCHEK                                                              30.60
!     CVMESH                                                              30.60
!     NEWTON                                                              30.60
!     EVALF                                                               30.60
!     SVALQI
!     SPRCON: Execution of some tests on the given model description
!     SINUPT
!     SINBTG
!     SINCMP
!     SIRAY : Searching the first point on a ray where the depth is DP
!     SWNMPS
!     SVARTP
!     BOUNPT
!     SWIPOL
!     CHGBAS                                                              40.00
!     GAMMAF
!     GAMMLN
!TIMG!     SWTSTA                                                              40.23
!TIMG!     SWTSTO                                                              40.23
!TIMG!     SWPRTI                                                              40.23
!     TXPBLA                                                              40.23
!     INTSTR                                                              40.23
!     NUMSTR                                                              40.23
!     SWCOPI                                                              40.23
!     SWCOPR                                                              40.23
!MatL4!     SWI2B                                                               40.30
!MatL4!     SWR2B                                                               40.30
!     MKPATH                                                              41.95
!
!  functions:
!  ----------
!  DEGCNV  (converts from cartesian convention to nautical and            32.01
!           vice versa)                                                   32.01
!  ANGRAD  (converts radians to degrees)                                  32.01
!  ANGDEG  (converts degrees to radians)                                  32.01
!
!***********************************************************************
!                                                                      *
      SUBROUTINE READXY (NAMX, NAMY, XX, YY, KONT, XSTA, YSTA)
!                                                                      *
!***********************************************************************
!
      USE OCPCOMM1                                                        40.41
      USE SwashCommdata3                                                  40.41
!
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
!  0. Authors
!     40.22: John Cazes and Tim Campbell
!     40.13: Nico Booij
!     40.51: Marcel Zijlema
!
!  1. UPDATE
!
!       Nov. 1996               offset values are added to standard values
!                               because they will be subtracted later
!     40.13, Nov. 01: a valid value for YY is required if a valid value
!                     for XX has been given; ocpcomm1.inc reactivated
!     40.51, Feb. 05: correction to location points equal to offset values
!
!  2. PURPOSE
!
!       Read x and y, initialize offset values XOFFS and YOFFS
!
!  3. METHOD
!
!       ---
!
!  4. PARAMETERLIST
!
!       NAMX, NAMY   inp char    names of the two coordinates as given in
!                                the user manual
!       XX, YY       out real    values of x and y taking into account offset
!       KONT         inp char    what to be done if values are missing
!                                see doc. of INDBLE (Ocean Pack doc.)
!       XSTA, YSTA   inp real    standard values of x and y
!
!  5. SUBROUTINES CALLING
!
!
!
!  6. SUBROUTINES USED
!
!       INDBLE (Ocean Pack)
      LOGICAL EQREAL
!
!  7. ERROR MESSAGES
!
!       ---
!
!  8. REMARKS
!
!       ---
!
!  9. STRUCTURE
!
!       ----------------------------------------------------------------
!       Read x and y in double prec.
!       If this is first couple of values
!       Then assign values to XOFFS and YOFFS
!            make LXOFFS True
!       ---------------------------------------------------------------
!       make XX and YY equal to x and y taking into account offset
!       ----------------------------------------------------------------
!
! 10. SOURCE TEXT
!
      DOUBLE PRECISION XTMP, YTMP
      CHARACTER  NAMX *(*), NAMY *(*), KONT *(*)
      SAVE  IENT
      DATA  IENT /0/
      CALL  STRACE (IENT,'READXY')
!
      CALL INDBLE (NAMX, XTMP, KONT, DBLE(XSTA)+DBLE(XOFFS))
      IF (CHGVAL) THEN                                                    40.13
!       a valid value was given for XX                                    40.13
        CALL INDBLE (NAMY, YTMP, 'REQ', DBLE(YSTA)+DBLE(YOFFS))           40.13
      ELSE                                                                40.13
        CALL INDBLE (NAMY, YTMP, KONT, DBLE(YSTA)+DBLE(YOFFS))
      ENDIF                                                               40.13
      IF (.NOT.LXOFFS) THEN
        XOFFS = REAL(XTMP)
        YOFFS = REAL(YTMP)
        LXOFFS = .TRUE.
      ENDIF
      IF (.NOT.EQREAL(XOFFS,REAL(XTMP))) THEN                             40.51
         XX = REAL(XTMP-DBLE(XOFFS))
      ELSE                                                                40.51
         XX = 0.                                                          40.51
      END IF                                                              40.51
      IF (.NOT.EQREAL(YOFFS,REAL(YTMP))) THEN                             40.51
         YY = REAL(YTMP-DBLE(YOFFS))
      ELSE                                                                40.51
         YY = 0.                                                          40.51
      END IF                                                              40.51
!
      RETURN
! * end of subroutine READXY  *
      END
!***********************************************************************
!                                                                      *
      SUBROUTINE REFIXY (NDS, XX, YY, IERR)
!                                                                      *
!***********************************************************************
!
      USE SwashCommdata3                                                  40.41
!
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
!  0. Authors
!     40.22: John Cazes and Tim Campbell
!     40.51: M. Zijlema
!
!  1. UPDATE
!
!       first version: 10.18 (Sept 1994)
!
!  2. PURPOSE
!
!       initialize offset values XOFFS and YOFFS, and shift XX and YY
!
!  3. METHOD
!
!       ---
!
!  4. PARAMETERLIST
!
!       NDS          in  int     file reference number
!       XX, YY       out real    values of x and y taking into account offset
!       IERR         out int     error indicator: IERR=0: no error, =-1: end-
!                                of-file, =-2: read error
!
!  5. SUBROUTINES CALLING
!
!
!
!  6. SUBROUTINES USED
!
      LOGICAL EQREAL
!
!  7. ERROR MESSAGES
!
!       ---
!
!  8. REMARKS
!
!       ---
!
!  9. STRUCTURE
!
!       ----------------------------------------------------------------
!       If this is first couple of values
!       Then assign values to XOFFS and YOFFS
!            make LXOFFS True
!       ---------------------------------------------------------------
!       make XX and YY equal to x and y taking into account offset
!       ----------------------------------------------------------------
!
! 10. SOURCE TEXT
!
      DOUBLE PRECISION XTMP, YTMP
      REAL             XX, YY
      SAVE  IENT
      DATA  IENT /0/
      CALL  STRACE (IENT,'REFIXY')
!
      READ (NDS, *, END=10, ERR=20) XTMP, YTMP
      IF (.NOT.LXOFFS) THEN
        XOFFS = REAL(XTMP)
        YOFFS = REAL(YTMP)
        LXOFFS = .TRUE.
      ENDIF
      IF (.NOT.EQREAL(XOFFS,REAL(XTMP))) THEN                             40.51
         XX = REAL(XTMP-DBLE(XOFFS))
      ELSE                                                                40.51
         XX = 0.                                                          40.51
      END IF                                                              40.51
      IF (.NOT.EQREAL(YOFFS,REAL(YTMP))) THEN                             40.51
         YY = REAL(YTMP-DBLE(YOFFS))
      ELSE                                                                40.51
         YY = 0.                                                          40.51
      END IF                                                              40.51
!
      IERR = 0
      RETURN
!     end of file
  10  IERR = -1
      RETURN
!     read error
  20  IERR = -2
      RETURN
! * end of subroutine REFIXY  *
      END
!****************************************************************
!
      SUBROUTINE CVCHEK (KGRPNT, XCGRID, YCGRID)                          30.72
!
!****************************************************************
!
      USE OCPCOMM4                                                        40.41
      USE SwashCommdata2                                                  40.41
      USE SwashCommdata3                                                  40.41
      USE SwashCommdata4                                                  40.41
!
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
!  0. Authors
!
!     30.72: IJsbrand Haagsma
!     40.13: Nico Booij
!     40.41: Marcel Zijlema
!
!  1. Updates
!
!            May  96: New subroutine
!     30.72, Feb. 98: Introduced generic names XCGRID, YCGRID and SPCSIG for SWAN
!     40.13, Mar. 01: messages corrected and extended
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!
!  2. Purpose
!
!     Checks whether the given curvilinear grid is correct
!     also set the value of CVLEFT.
!
!  3. Method
!
!     Going around a mesh in the same direction the interior
!     of the mesh must be always in the same side if the
!     coordinates are correct
!
!  4. Argument variables
!
!     KGRPNT: input  Array of indirect addressing
!
      INTEGER KGRPNT(MXC,MYC)                                             30.72
!
!     XCGRID: input  Coordinates of computational grid in x-direction     30.72
!     YCGRID: input  Coordinates of computational grid in y-direction     30.72
!
      REAL    XCGRID(MXC,MYC),    YCGRID(MXC,MYC)                         30.72
!
!
!     5. SUBROUTINES CALLING
!
!        SWRBC
!
!     6. SUBROUTINES USED
!
!        ---
!
!     7. ERROR MESSAGES
!
!        ---
!
!     8. REMARKS
!
!
!     9. STRUCTURE
!
!   ------------------------------------------------------------
!     FIRST = True
!     For ix=1 to MXC-1 do
!         For iy=1 to MYC-1 do
!             For iside=1 to 4 do
!                 Case iside=
!                 1: K1 = KGRPNT(ix,iy), K2 = KGRPNT(ix+1,iy),
!                    K3 = KGRPNT(ix+1,iy+1)
!                 2: K1 = KGRPNT(ix+1,iy), K2 = KGRPNT(ix+1,iy+1),
!                    K3 = KGRPNT(ix,iy+1)
!                 3: K1 = KGRPNT(ix+1,iy+1), K2 = KGRPNT(ix,iy+1),
!                    K3 = KGRPNT(ix,iy)
!                 4: K1 = KGRPNT(ix,iy+1), K2 = KGRPNT(ix,iy),
!                    K3 = KGRPNT(ix+1,iy)
!                 ---------------------------------------------------
!                 If K1>1 and K2>1 and K3>1
!                 Then Det = (xpg(K3)-xpg(K1))*(ypg(K2)-ypg(K1)) -
!                            (ypg(K3)-ypg(K1))*(xpg(K2)-xpg(K1))
!                      If FIRST
!                      Then Make FIRST = False
!                           If Det>0
!                           Then Make CVleft = False
!                           Else Make CVleft = True
!                      ----------------------------------------------
!                      If ((CVleft and Det<0) or (not CVleft and Det>0))
!                      Then Write error message with IX, IY, ISIDE
!   ------------------------------------------------------------
!
!     10. SOURCE
!
!****************************************************************
!
!
      LOGICAL  FIRST
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'CVCHEK')
!
!     test output
!
      IF (ITEST .GE. 150 .OR. INTES .GE. 30) THEN
        WRITE(PRINTF,186)
 186    FORMAT(/,' ... Subroutine CVCHEK...',
     &  /,2X,'POINT( IX, IY),  INDEX,      COORDX,       COORDY')
        ICON = 0
        DO 5 IIY = 1, MYC-1
          DO 6 IIX = 1, MXC-1
            ICON = ICON + 1
            WRITE(PRINTF,7)IIX-1,IIY-1,KGRPNT(IIX,IIY),
     &      XCGRID(IIX,IIY)+XOFFS, YCGRID(IIX,IIY)+YOFFS                  30.72 40.13
 6        CONTINUE
 5      CONTINUE
      ENDIF
 7    FORMAT(4X,I5,1X,I5,3X,I4,5X,F10.2,4X,F10.2)
!
      FIRST = .TRUE.
!
      DO 10 IX = 1,MXC-2
        DO 15 IY = 1,MYC-2
          DO 20 ISIDE = 1,4
            IF (ISIDE .EQ. 1) THEN
              IX1 = IX                                                    40.13
              IY1 = IY                                                    40.13
              IX2 = IX+1                                                  40.13
              IY2 = IY                                                    40.13
              IX3 = IX+1                                                  40.13
              IY3 = IY+1                                                  40.13
            ELSE IF (ISIDE .EQ. 2) THEN
              IX1 = IX+1                                                  40.13
              IY1 = IY                                                    40.13
              IX2 = IX+1                                                  40.13
              IY2 = IY+1                                                  40.13
              IX3 = IX                                                    40.13
              IY3 = IY+1                                                  40.13
            ELSE IF (ISIDE .EQ. 3) THEN
              IX1 = IX+1                                                  40.13
              IY1 = IY+1                                                  40.13
              IX2 = IX                                                    40.13
              IY2 = IY+1                                                  40.13
              IX3 = IX                                                    40.13
              IY3 = IY                                                    40.13
            ELSE IF (ISIDE .EQ. 4) THEN
              IX1 = IX                                                    40.13
              IY1 = IY+1                                                  40.13
              IX2 = IX                                                    40.13
              IY2 = IY                                                    40.13
              IX3 = IX+1                                                  40.13
              IY3 = IY                                                    40.13
            ENDIF
            K1  = KGRPNT(IX1,IY1)                                         40.13
            XC1 = XCGRID(IX1,IY1)                                         40.13 30.72
            YC1 = YCGRID(IX1,IY1)                                         40.13 30.72
            K2  = KGRPNT(IX2,IY2)                                         40.13
            XC2 = XCGRID(IX2,IY2)                                         40.13 30.72
            YC2 = YCGRID(IX2,IY2)                                         40.13 30.72
            K3  = KGRPNT(IX3,IY3)                                         40.13
            XC3 = XCGRID(IX3,IY3)                                         30.72
            YC3 = YCGRID(IX3,IY3)                                         30.72
            DET   = 0.
            IF (K1 .GE. 2 .AND. K2 .GE. 2 .AND. K3 .GE. 2) THEN
              DET = ((XC3 - XC1) * (YC2 - YC1)) -
     &              ((YC3 - YC1) * (XC2 - XC1))
              IF (DET .EQ. 0.) THEN
!               three grid points on one line                             40.13
                CALL MSGERR (2,'3 comp. grid points on one line')         40.13
                WRITE (PRINTF, 112)
     &               IX1-1, IY1-1, XC1+XOFFS, YC1+YOFFS,                  40.13
     &               IX2-1, IY2-1, XC2+XOFFS, YC2+YOFFS,                  40.13
     &               IX3-1, IY3-1, XC3+XOFFS, YC3+YOFFS                   40.13
 112            FORMAT (3(1X, 2I3, 2(1X, F14.4)))                         40.13
              ENDIF
!
              IF (FIRST) THEN
                FIRST = .FALSE.
                IF (DET .GT. 0.) THEN
                  CVLEFT = .FALSE.
                ELSE
                  CVLEFT = .TRUE.
                ENDIF
              ENDIF
              IF (     (      CVLEFT .AND. DET .GT. 0.)
     &            .OR. (.NOT. CVLEFT .AND. DET .LT. 0.)) THEN
!               crossing grid lines in a mesh                             40.13
                CALL MSGERR (2,'grid angle <0 or >180 degrees')           40.13
                WRITE (PRINTF, 112)
     &               IX1-1, IY1-1, XC1+XOFFS, YC1+YOFFS,                  40.13
     &               IX2-1, IY2-1, XC2+XOFFS, YC2+YOFFS,                  40.13
     &               IX3-1, IY3-1, XC3+XOFFS, YC3+YOFFS                   40.13
              ENDIF
            ENDIF
 20       CONTINUE
 15     CONTINUE
 10   CONTINUE
      RETURN
!     *** end of subroutine CVCHEK ***
      END
!***********************************************************************
!                                                                      *
      SUBROUTINE CVMESH (XP, YP, XC, YC, KGRPNT, XCGRID ,YCGRID)
!                                                                      *
!***********************************************************************
!
      USE OCPCOMM4                                                        40.41
      USE SwashCommdata2                                                  40.41
      USE SwashCommdata3                                                  40.41
      USE SwashCommdata4                                                  40.41
!
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
!  0. Authors
!
!     30.72: IJsbrand Haagsma
!     40.00, 40.13: Nico Booij
!     40.02: IJsbrand Haagsma
!     40.41: Marcel Zijlema
!
!  1. Updates
!
!     30.21, Jun. 96: New for curvilinear version
!     30.72, Feb. 98: Introduced generic names XCGRID, YCGRID and SPCSIG for SWAN
!     40.00, May  98: procedure for points outside grid accelerated
!     40.00, Feb  99: procedure extended for 1D case
!                     XOFFS and YOFFS added in write statements
!     40.02, Mar. 00: Fixed bug that placed dry testpoints outside computational grid
!     40.13, Mar. 01: message "CVMESH 2nd attempt .." suppressed
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!     40.41, Nov. 04: search for boundary points improved
!
!  2. Purpose
!
!     procedure to find location in curvilinear grid for a point
!     given in problem coordinates
!
!  3. Method
!
!     First attempt: use Newton-Raphson method to find XC and YC
!     (Note: in the program XC and YC indicate the mesh and position in
!     the mesh) in a few steps; this may be most efficient if a series of
!     points is processed, because the previous point provides a good
!     first estimate.
!     This procedure may fail if the number of iterations is larger than
!     a previously set limit (default=5).
!
!     If the first attempt fails then determine whether the points (XP,YP)
!     is inside the mesh.
!     A point is assumed to be inside a mesh if it is on the
!     interior side for all 4 sides of the mesh. Here we use common variable
!     CVLEFT (logical): if True interior of the mesh is always on the left
!     going along the mesh in the order: (ix,iy), (ix+1,iy), (ix+1,iy+1),
!     (ix,iy+1), (ix,iy). If it is False the interior is always on the right.
!     Whether a point  is on the left or on the right of a line can be
!     decided by looking at the sign of the determinant.
!     If the point is inside of any mesh of the computational grid then
!     the Newton-Raphson procedure is used
!     again with the pivoting point like first guess. Otherwise, scan the
!     boundaries whether the point is on the boundaries. If this fails, it
!     may be concluded that the point (XP,YP) is outside the grid.
!
!  4. Argument variables
!
!     XCGRID  input  Coordinates of computational grid in x-direction     30.72
!     YCGRID  input  Coordinates of computational grid in y-direction     30.72
!     XP, YP  input  a point given in problem coordinates
!     XC, YC  outp   same point in computational grid coordinates
!
      REAL     XCGRID(MXC,MYC),    YCGRID(MXC,MYC)                        30.72
      REAL     XP, YP, XC, YC
!
!     KGRPNT   input   array(MXC,MYC)  grid numbers
!                      if KGRPNT <= 1, point is not in comp. grid.
!
      INTEGER  KGRPNT(MXC,MYC)
!
!     Local variables
!
!     MXITNR   number of iterations in Newton-Raphson procedure
!     IX, IY   counter of computational grid point
!     K1       address of grid point
!     IXMIN    counter of grid point closest to (XP,YP)
!     IYMIN    counter of grid point closest to (XP,YP)
!     IBND     counter of boundary grid points
!
      INTEGER       :: IX, IY, K1, IXMIN, IYMIN, IBND
      INTEGER, SAVE :: MXITNR = 0
      INTEGER, SAVE :: IENT = 0
!
!     INMESH   if True, point (XP,YP) is inside the computational grid
!     FINDXY   if True, Newton-Raphson procedure succeeded
!     ONBND    if True, given point is on boundary                        40.41
!
      LOGICAL  INMESH ,FINDXY, ONBND
!
!     DISMIN   minimal distance found                                     40.41
!     XPC1     user coordinate of a computational grid point
!     YPC1     user coordinate of a computational grid point
!     XC0      grid coordinate of grid point closest to (XP,YP)
!     YC0      grid coordinate of grid point closest to (XP,YP)
!
      REAL       :: DISMIN                                                40.41
      REAL       :: XPC1, YPC1, XC0, YC0
      REAL       :: XC1, XC2, YC1, YC2, DET
!
!  5. SUBROUTINES CALLING
!
!     SINCMP
!
!  6. SUBROUTINES USED
!
!       NEWTON
!
!  7. ERROR MESSAGES
!
!       ---
!
!  8. REMARKS
!
!       ---
!
!  9. STRUCTURE
!
!     --------------------------------------------------------------
!     Determine XC and YC from XP and YP using Newton-Raphson iteration
!     process
!     If (XC and YC were found) then
!       Procedure is ready; Return values of XC and YC
!       return
!     else
!     ---------------------------------------------------------------------
!       Inmesh = True
!       For iside=1 to 4 do
!             Case iside=
!             1: K1 = KGRPNT(    2,    1), K2 = KGRPNT(MXC-1,    1)
!             2: K1 = KGRPNT(MXC-1,    2), K2 = KGRPNT(MXC-1,MYC-1)
!             3: K1 = KGRPNT(MXC-1,MYC-1), K2 = KGRPNT(    2,MYC-1)
!             4: K1 = KGRPNT(    1,MYC-1), K2 = KGRPNT(    1,    2)
!             ----------------------------------------------------------
!             Det = (xp-xpg(K1))*(ypg(K2)-ypg(K1)) -
!                   (yp-ypg(K1))*(xpg(K2)-xpg(K1))
!             If ((CVleft and Det>0) or (not CVleft and Det<0))
!             Then Make Inmesh = False
!             --------------------------------------------------------
!             If Inmesh
!             Then Determine XC and YC using Newton-Raphson iteration
!                  process
!                  Procedure is ready; Return values of XC and YC
!     ---------------------------------------------------------------------
!     No mesh is found: Make XC and YC = exception value
!     Return values of XC and YC
!     ---------------------------------------------------------------------
!
!****************************************************************
!
!
      IF (LTRACE) CALL STRACE (IENT,'CVMESH')
!
      IF (ONED) THEN
        CALL NEWT1D  (XP, YP, XCGRID, YCGRID, KGRPNT,                     40.00
     &                XC ,YC ,FINDXY)
        IF (.NOT.FINDXY) THEN
          XC = -99.
          YC = -99.
          IF (ITEST .GE. 150 .OR. INTES .GE. 20) THEN
            WRITE(PRINTF, 85) XP+XOFFS, YP+YOFFS                          40.00
          ENDIF
        ENDIF
        GOTO 99
      ELSE
!       two-dimensional computation
        XC = 1.
        YC = 1.
!       --- First attempt, to find XC,YC with Newton-Raphson method
        MXITNR = 5
        CALL NEWTON  (XP, YP, XCGRID, YCGRID,                             40.00
     &                MXITNR ,ITER, XC ,YC ,FINDXY)                       40.41 40.02
        IF ((ITEST .GE. 150 .OR. INTES .GE. 20) .AND. FINDXY) THEN        40.02
           WRITE(PRINTF,25) XP+XOFFS ,YP+YOFFS ,XC ,YC                    40.03
        ENDIF
 25     FORMAT (' CVMESH: (XP,YP)=','(',F12.4,',',F12.4,
     &          '), (XC,YC)=','(',F9.2,',',F9.2,')')
        IF (FINDXY) GOTO 80
!
        INMESH = .TRUE.
        DO ISIDE = 1,4
           IF (ISIDE .EQ. 1) THEN
              XC1 = XCGRID(    2,    1)
              YC1 = YCGRID(    2,    1)
              XC2 = XCGRID(MXC-1,    1)
              YC2 = YCGRID(MXC-1,    1)
           ELSE IF (ISIDE .EQ. 2) THEN
              XC1 = XCGRID(MXC-1,    2)
              YC1 = YCGRID(MXC-1,    2)
              XC2 = XCGRID(MXC-1,MYC-1)
              YC2 = YCGRID(MXC-1,MYC-1)
           ELSE IF (ISIDE .EQ. 3) THEN
              XC1 = XCGRID(MXC-1,MYC-1)
              YC1 = YCGRID(MXC-1,MYC-1)
              XC2 = XCGRID(    2,MYC-1)
              YC2 = YCGRID(    2,MYC-1)
           ELSE IF (ISIDE .EQ. 4) THEN
              XC1 = XCGRID(    1,MYC-1)
              YC1 = YCGRID(    1,MYC-1)
              XC2 = XCGRID(    1,    2)
              YC2 = YCGRID(    1,    2)
           ENDIF
!
           DET= (XP - XC1)*(YC2 - YC1) - (YP - YC1)*(XC2 - XC1)
           IF ((      CVLEFT .AND. DET .GE. 0) .OR.
     &         (.NOT. CVLEFT .AND. DET .LE. 0)) INMESH = .FALSE.

        ENDDO
!
 30     IF ( INMESH ) THEN
!         --- select grid point closest to (XP,YP)
          DISMIN = 1.E20
          DO 50 IX = 1,MXC-1
            DO 40 IY = 1,MYC-1
              K1  = KGRPNT(IX,IY)
              IF (K1.GT.1) THEN
                XPC1 = XCGRID(IX,IY)
                YPC1 = YCGRID(IX,IY)
                DISXY = SQRT ((XP-XPC1)**2 + (YP-YPC1)**2)
                IF (DISXY .LT. DISMIN) THEN
                  IXMIN  = IX
                  IYMIN  = IY
                  DISMIN = DISXY
                ENDIF
              ENDIF
  40        CONTINUE
  50      CONTINUE
!         second attempt using closest grid point as first guess
          MXITNR = 20
          XC0 = REAL(IXMIN)
          YC0 = REAL(IYMIN)
!         ITEST condition changed from 20 to 120                          40.13
          IF (ITEST.GE.120) WRITE (PRTEST, 55) XP+XOFFS ,YP+YOFFS ,       40.13
     &          XC0-1. ,YC0-1.
  55      FORMAT (' CVMESH 2nd attempt, (XP,YP)=','(',F12.4,',',F12.4,
     &          '), (XC,YC)=','(',F9.2,',',F9.2,')')
          DO KORNER = 1, 4
            IF (KORNER.EQ.1) THEN
              XC = XC0 + 0.2
              YC = YC0 + 0.2
            ELSE IF (KORNER.EQ.2) THEN
              XC = XC0 - 0.2
              YC = YC0 + 0.2
            ELSE IF (KORNER.EQ.3) THEN
              XC = XC0 - 0.2
              YC = YC0 - 0.2
            ELSE
              XC = XC0 + 0.2
              YC = YC0 - 0.2
            ENDIF
            CALL NEWTON  (XP, YP, XCGRID, YCGRID,                         40.00
     &                    MXITNR ,ITER, XC ,YC ,FINDXY)                   40.41 40.02
            IF (FINDXY) THEN
              IF (ITEST .GE. 150 .OR. INTES .GE. 20) THEN
                WRITE(PRINTF,25) XP+XOFFS ,YP+YOFFS ,XC ,YC               40.00
              ENDIF
              GOTO 80
            ENDIF
          ENDDO
          IF (ITER.GE.MXITNR) THEN                                        40.41
             WRITE (PRINTF, 75) XP+XOFFS, YP+YOFFS, MXITNR                40.00
  75         FORMAT (' search for point with location ', 2F12.4,          40.41
     &               ' fails in', I3, ' iterations')                      40.41
          END IF                                                          40.41
        ELSE
!         scan boundary to see whether the point is close to the boundary
          DISMIN=99999.                                                   40.41
          ONBND =.FALSE.                                                  40.41
          IX1 = 0                                                         40.41
          IY1 = 0                                                         40.51
          IX2 = 0
          DO ISIDE = 1, 4
            IF (ISIDE.EQ.1) THEN
               IXB = 1
               IYB = 1
               IXE = MXC-1
               IYE = 1
               MIP = MXC-1
            ELSE IF (ISIDE.EQ.2) THEN
               IXB = MXC-1
               IYB = 1
               IXE = MXC-1
               IYE = MYC-1
               MIP = MYC-1
            ELSE IF (ISIDE.EQ.3) THEN
               IXB = MXC-1
               IYB = MYC-1
               IXE = 1
               IYE = MYC-1
               MIP = MXC-1
            ELSE IF (ISIDE.EQ.4) THEN
               IXB = 1
               IYB = MYC-1
               IXE = 1
               IYE = 1
               MIP = MYC-1
            ENDIF
            DO IP = 1, MIP
              RR  = REAL(IP-1) / REAL(MIP-1)
              IF (IP.GT.1) THEN
                 IX1 = IX2
                 IY1 = IY2
                 XP1 = XP2
                 YP1 = YP2
              ENDIF
              IX2 = IXB + NINT(RR*REAL(IXE-IXB))
              IY2 = IYB + NINT(RR*REAL(IYE-IYB))
              IF (KGRPNT(IX2,IY2).GT.1) THEN
                XP2 = XCGRID(IX2,IY2)
                YP2 = YCGRID(IX2,IY2)
!               --- determine relative distance from boundary segment
!                   with respect to the length of that segment
                IF (IP.GT.1) THEN
                   SLEN2  = (XP2-XP1)**2 + (YP2-YP1)**2
                   RELDIS = ABS((XP-XP1)*(YP2-YP1)-(YP-YP1)*(XP2-XP1)) /
     &                      SLEN2
                ELSE
                   RELDIS = 1.
                ENDIF
                IF (RELDIS.LT.0.01) THEN                                  40.41
!                 --- determine location on the boundary section
                  IF (RELDIS-DISMIN.LE.0.01) THEN                         40.41
                     DISMIN = RELDIS                                      40.41
                     RELLOC = ((XP-XP1)*(XP2-XP1)+(YP-YP1)*(YP2-YP1)) /
     &                        SLEN2
                     IF (RELLOC.GE.-0.001 .AND. RELLOC.LE.1.001) THEN     40.41
                        RELLCM = RELLOC                                   40.41
                        IF (RELLCM.LT.0.01) RELLCM=0.                     40.41
                        IF (RELLCM.GT.0.99) RELLCM=1.                     40.41
                        IX1M  = IX1                                       40.41
                        IX2M  = IX2                                       40.41
                        IY1M  = IY1                                       40.41
                        IY2M  = IY2                                       40.41
                        ONBND = .TRUE.                                    40.41
                     ENDIF                                                40.41
                  ENDIF
                ENDIF
              ENDIF
            ENDDO
          ENDDO
          IF (ONBND) THEN                                                 40.41
             XC = FLOAT(IX1M) + RELLCM * FLOAT(IX2M-IX1M) - 1.            40.41
             YC = FLOAT(IY1M) + RELLCM * FLOAT(IY2M-IY1M) - 1.            40.41
             IF (ITEST .GE. 150 .OR. INTES .GE. 20) THEN
                WRITE(PRINTF, 65) XP+XOFFS, YP+YOFFS, XC, YC              40.00
  65            FORMAT (' CVMESH: (XP,YP)=','(',F12.4,',',F12.4,
     &                  ') is on the boundary, (XC,YC)=(',
     &                  F9.2,',',F9.2,')')
             ENDIF
             GOTO 80
          ENDIF
          XC = -99.
          YC = -99.
          IF (ITEST .GE. 150 .OR. INTES .GE. 20) THEN
            WRITE(PRINTF, 85) XP+XOFFS, YP+YOFFS                          40.00
  85        FORMAT (' CVMESH: (XP,YP)=','(',F12.4,',',F12.4,
     &              ') is outside grid')
          ENDIF
          GOTO 99
        ENDIF
      ENDIF                                                               40.00
  80  IF (KGRPNT(INT(XC+3.001)-2,INT(YC+3.001)-2).LE.1) THEN
         WRITE (PRINTF, 90) XP+XOFFS, YP+YOFFS                            40.41
  90     FORMAT (' point with location ',2F12.4,' is not active')         40.41
         XC = -99.
         YC = -99.
      ENDIF
  99  RETURN
      END
!***********************************************************************
!                                                                      *
      SUBROUTINE NEWTON (XP, YP, XCGRID, YCGRID,                          40.00
     &                   MXITNR, ITER, XC, YC, FIND)                      40.41 40.02
!                                                                      *
!***********************************************************************
!
      USE OCPCOMM4                                                        40.41
      USE SwashCommdata3                                                  40.41
      USE SwashCommdata4                                                  40.41
!
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
!  0. Authors
!
!     30.72: IJsbrand Haagsma
!     30.80: Nico Booij
!     30.82: IJsbrand Haagsma
!     40.41: Marcel Zijlema
!
!  1. Updates
!
!     30.21, Jun. 96: New for curvilinear version
!     30.72, Feb. 98: Introduced generic names XCGRID, YCGRID and SPCSIG for SWAN
!     30.82, Oct. 98: Updated description of several variables
!     30.80, Oct. 98: computation of update of XC,YC modified to avoid
!                     division by 0
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!
!  2. Purpose
!
!     Solve eqs. and find a point  (XC,YC) in a curvilinear grid (compt.
!     grid) for a given point (XP ,YP) in a cartesian grid (problem coord).
!
!  3. Method
!
!     In this subroutine the next equations are solved :
!
!                  @XP             @XP
!     XP(xc,yc) +  --- * @XC   +   --- * @YC  - XP(x,y) = 0
!                  @XC             @YC
!
!                  @YP             @YP
!     YP(xc,yc) +  --- * @XC   +    --- * @YC  - YP(x,y) = 0
!                  @XC             @YC
!
!     In the subroutine, next notation is used for the previous eqs.
!     XVC       + DXDXC * DXC   + DXDYC * DYC - XP  = 0.
!     YVC       + DYDXC * DXC   + DYDYC * DYC - YP  = 0.
!
!
!  4. Argument variables
!
! i   MXITNR: Maximum number of iterations                                30.82
!
      INTEGER MXITNR, ITER                                                40.41 40.02
!
!   o XC    : X-coordinate in computational coordinates                   30.82
! i   XCGRID: Coordinates of computational grid in x-direction            30.72
! i   XP    : X-coordinate in problem coordinates                         30.82
!   o YC    : Y-coordinate in computational coordinates                   30.82
! i   YCGRID: Coordinates of computational grid in y-direction            30.72
! i   YP    : Y-coordinate in problem coordinates                         30.82
!
      REAL    XC, XCGRID(MXC,MYC), XP                                     30.82
      REAL    YC, YCGRID(MXC,MYC), YP                                     30.82
!
!   o FIND  : Whether XC and YC are found                                 30.82
!
      LOGICAL FIND                                                        30.82
!
!  6. SUBROUTINES USED
!
!     STRACE
!
!
!  7. ERROR MESSAGES
!
!       ---
!
!  8. REMARKS
!
!       ---
!
!  9. STRUCTURE
!
!       -----------------------------------------------------------------
!       -----------------------------------------------------------------
!
! 13. Source text
!
      INTEGER, SAVE :: IENT = 0
      IF (LTRACE) CALL STRACE (IENT,'NEWTON')
!
      DXC    = 1000.
      DYC    = 1000.
      TOLDC  = 0.001
      FIND   = .FALSE.
!
      IF (ITEST .GE. 200) THEN
        WRITE(PRINTF,*) ' Coordinates in subroutine NEWTON '
        DO J = 1, MYC-1
          DO I = 1, MXC-1
            WRITE(PRINTF,30) I ,J ,XCGRID(I,J) ,YCGRID(I,J)               30.72
          ENDDO
        ENDDO
      ENDIF
 30   FORMAT(2(2X,I5),2(2X,E12.4))
!
      DO 14 K = 1 ,MXITNR
        ITER = K                                                          40.41
        I1   = INT(XC)                                                    40.00
        J1   = INT(YC)
        IF (I1 .EQ. MXC-1) I1 = I1 - 1
        IF (J1 .EQ. MYC-1) J1 = J1 - 1
        I2  = I1 + 1
        J2  = J1 + 1
        FJ1 = FLOAT(J1)
        FI1 = FLOAT(I1)
        FJ2 = FLOAT(J2)
        FI2 = FLOAT(I2)
!
        XVC   = (YC-FJ1)*((XC-FI1)*XCGRID(I2,J2)  +
     &                    (FI2-XC)*XCGRID(I1,J2)) +
     &          (FJ2-YC)*((XC-FI1)*XCGRID(I2,J1)  +
     &                    (FI2-XC)*XCGRID(I1,J1))
        YVC   = (YC-FJ1)*((XC-FI1)*YCGRID(I2,J2)  +
     &                    (FI2-XC)*YCGRID(I1,J2)) +
     &          (FJ2-YC)*((XC-FI1)*YCGRID(I2,J1)  +
     &                    (FI2-XC)*YCGRID(I1,J1))
        DXDXC = (YC -FJ1)*(XCGRID(I2,J2) - XCGRID(I1,J2)) +
     &          (FJ2-YC )*(XCGRID(I2,J1) - XCGRID(I1,J1))
        DXDYC = (XC -FI1)*(XCGRID(I2,J2) - XCGRID(I2,J1)) +
     &          (FI2-XC )*(XCGRID(I1,J2) - XCGRID(I1,J1))
        DYDXC = (YC -FJ1)*(YCGRID(I2,J2) - YCGRID(I1,J2)) +
     &          (FJ2-YC )*(YCGRID(I2,J1) - YCGRID(I1,J1))
        DYDYC = (XC -FI1)*(YCGRID(I2,J2) - YCGRID(I2,J1)) +
     &          (FI2-XC )*(YCGRID(I1,J2) - YCGRID(I1,J1))
!
        IF (ITEST .GE. 150)
     &    WRITE(PRINTF,35) K, XC-1., YC-1., XP, YP, XVC, YVC              40.00
 35     FORMAT(' NEWTON  iter=', I2, ' (XC,YC)=', 2(1X,F10.2),/,          40.00
     &         ' (XP,YP)=', 2(1X,F10.2),
     &         '  X,Y(XC,YC) = ', 2(1X,F10.2))
        IF (ITEST .GE. 180) WRITE(PRINTF,36)
     &     XCGRID(I1,J1), XCGRID(I1,J2), XCGRID(I2,J1), XCGRID(I2,J2),
     &     YCGRID(I1,J1), YCGRID(I1,J2), YCGRID(I2,J1), YCGRID(I2,J2),
     &                     DXDXC, DXDYC, DYDXC, DYDYC                     40.00
 36     FORMAT(' NEWTON grid coord:', 8(1x, F10.0), /
     &         '        deriv=', 4(1X,F10.2))                             40.00
!
!       *** the derivated terms of the eqs. are evaluated and  ***
!       *** the eqs. are solved                                ***
        DDEN = DXDXC*DYDYC - DYDXC*DXDYC                                  30.80
        DXP  = XP - XVC                                                   30.80
        DYP  = YP - YVC                                                   30.80
        IF ( DDEN.NE.0. ) THEN
           DXC = ( DYDYC*DXP - DXDYC*DYP) / DDEN                          30.80
           DYC = (-DYDXC*DXP + DXDXC*DYP) / DDEN                          30.80
        ENDIF
!
        XC = XC + DXC
        YC = YC + DYC
!
!       *** If the guess point (XC,YC) is outside of compt. ***
!       *** grid, put that point in the closest boundary    ***
        IF (XC .LT. 1.   ) XC = 1.
        IF (YC .LT. 1.   ) YC = 1.
        IF (XC .GT. MXC-1) XC = FLOAT(MXC-1)
        IF (YC .GT. MYC-1) YC = FLOAT(MYC-1)
!
        IF (ITEST .GE. 120 .OR. INTES .GE. 50 .OR. IOUTES .GE. 50)
     &    WRITE(PRINTF,42) DXC, DYC, XC-1., YC-1.                         40.00
 42     FORMAT(' (DXC,DYC)=', 2(1X,F10.2), ' (XC,YC)=', 2(1X,F10.2))      40.00
!
!       *** If the accuracy is reached stop the iteration,  ***
        IF (ABS(DXC) .LE. TOLDC .AND. ABS(DYC) .LE. TOLDC) THEN
!
          FIND = .TRUE.
          XC = XC -1.
          YC = YC -1.
          RETURN
        ENDIF
!
 14   CONTINUE
      RETURN
!     *** end of subroutine NEWTON ***
      END
!***********************************************************************
!                                                                      *
      SUBROUTINE NEWT1D (XP, YP, XCGRID, YCGRID, KGRPNT,                  40.00
     &                   XC, YC, FIND)
!                                                                      *
!***********************************************************************
!
      USE OCPCOMM4                                                        40.41
      USE SwashCommdata2                                                  40.41
      USE SwashCommdata3                                                  40.41
!
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
!  0. Authors
!
!     40.00, 40.13: Nico Booij
!     40.41: Marcel Zijlema
!
!  1. Updates
!
!     40.00, Feb. 99: New (adaptation from subr NEWTON for 1D case)
!     40.13, Feb. 01: DX and DY renamed to DELX and DELY (DX and DY are
!                     common var.); error in expression for RS corrected
!                     PRINTF replaced by PRTEST in test output
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!
!  2. Purpose
!
!     Finds broken coordinate XC for a given point XP in a rectilinear grid
!
!  3. Method
!
!     In this subroutine the step on the computational grid is selected
!     for which
!
!           (X-X1).(X2-X1)
!     0 <= --------------- <= 1
!          (X2-X1).(X2-X1)
!
!     where X, X1 and X2 are vectors; X corresponds to (Xp,Yp)
!     X1 and X2 are two neighbouring grid points
!
!  4. Argument variables
!
! i   KGRPNT: Grid adresses                                               40.00
!
      INTEGER KGRPNT(MXC,MYC)                                             40.00
!
!   o XC    : X-coordinate in computational coordinates                   30.82
! i   XCGRID: Coordinates of computational grid in x-direction            30.72
! i   XP    : X-coordinate in problem coordinates                         30.82
!   o YC    : Y-coordinate in computational coordinates                   30.82
! i   YCGRID: Coordinates of computational grid in y-direction            30.72
! i   YP    : Y-coordinate in problem coordinates                         30.82
!
      REAL    XC, XCGRID(MXC,MYC), XP                                     30.82
      REAL    YC, YCGRID(MXC,MYC), YP                                     30.82
!
!   o FIND  : Whether XC and YC are found                                 30.82
!
      LOGICAL FIND                                                        30.82
!
!     Local variables:

      REAL :: DELX, DELY   ! grid line                                    40.13
!
!  6. SUBROUTINES USED
!
!       ---
!
!  7. ERROR MESSAGES
!
!       ---
!
!  8. REMARKS
!
!       ---
!
!  9. STRUCTURE
!
!       -----------------------------------------------------------------
!       -----------------------------------------------------------------
!
! 13. Source text
!
      SAVE     IENT
      DATA     IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'NEWT1D')
!
      IF (ITEST .GE. 120) THEN                                            40.13
        WRITE(PRTEST,*) ' Coordinates in subroutine NEWT1D '              40.13
        DO I = 1, MXC
          WRITE(PRTEST,30) I, XCGRID(I,1)+XOFFS ,YCGRID(I,1)+YOFFS        40.13
        ENDDO
      ENDIF
 30   FORMAT(2X,I5,2(2X,E12.4))
!
      FIND = .FALSE.
      DO 40 IX = 2 ,MXC-1
        IF (KGRPNT(IX-1,1).GT.1) THEN
          X1 = XCGRID(IX-1,1)
          Y1 = YCGRID(IX-1,1)
        ELSE
          GOTO 40
        ENDIF
        IF (KGRPNT(IX,1).GT.1) THEN
          X2 = XCGRID(IX,1)
          Y2 = YCGRID(IX,1)
        ELSE
          GOTO 40
        ENDIF
!       both ends of the step are valid grid points
!       now verify whether projection of (Xp,Yp) is within the step
        DELX = X2 - X1                                                    40.13
        DELY = Y2 - Y1                                                    40.13
        RS = ((XP - X1) * DELX + (YP - Y1) * DELY) /                      40.13
     &              (DELX * DELX + DELY * DELY)                           40.13
        IF (RS.GE.0. .AND. RS.LE.1.) THEN
          FIND = .TRUE.
          XC = REAL(IX-2) + RS                                            40.00
          YC = 0.
          GOTO 50
        ENDIF
  40  CONTINUE
  50  RETURN
!     *** end of subroutine NEWT1D ***
      END
!***********************************************************************
!                                                                      *
      SUBROUTINE EVALF (XC ,YC ,XVC ,YVC ,XCGRID ,YCGRID)                 30.72
!                                                                      *
!***********************************************************************
!
      USE OCPCOMM4                                                        40.41
      USE SwashCommdata3                                                  40.41
!
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
!  0. Authors
!
!     30.72: IJsbrand Haagsma
!     40.41: Marcel Zijlema
!
!  1. Updates
!
!     30.21, Jun. 96: New for curvilinear version
!     30.72, Feb. 98: Introduced generic names XCGRID, YCGRID and SPCSIG for SWAN
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!
!  2. Purpose
!
!     Evaluate the coordinates (in problem coordinates) of point (XC,YC)
!     given in computational coordinates
!
!  3. Method
!
!     Bilinear interpolation
!
!  4. Argument variables
!
!     XCGRID: input  Coordinates of computational grid in x-direction     30.72
!     YCGRID: input  Coordinates of computational grid in y-direction     30.72
!
      REAL    XCGRID(MXC,MYC),    YCGRID(MXC,MYC)                         30.72
!
!       XC, YC      real, outp    point in computational grid coordinates
!       XVC, YCV    real, OUTP    same point  but in problem coordinates
!
!  6. SUBROUTINES USED
!
!       none
!
!  7. ERROR MESSAGES
!
!       ---
!
!  8. REMARKS
!
!
!  9. STRUCTURE
!
!       -----------------------------------------------------------------
!       -----------------------------------------------------------------
!
! 10. SOURCE TEXT
!
      SAVE     IENT
      DATA     IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'EVALF')
!
      I  = INT(XC)
      J  = INT(YC)
!
!     *** If the guess point (XC,YC) is in the boundary   ***
!     *** where I = MXC or/and J = MYC the interpolation  ***
!     *** is done in the mesh with pivoting point         ***
!     *** (MXC-1, J) or/and (I,MYC-1)                     ***
!
      IF (I .EQ. MXC-1) I = I - 1
      IF (J .EQ. MYC-1) J = J - 1
      T = XC - FLOAT(I)
      U = YC - FLOAT(J)
!     *** For x-coord. ***
      P1 = XCGRID(I,J)                                                    30.72
      P2 = XCGRID(I+1,J)                                                  30.72
      P3 = XCGRID(I+1,J+1)                                                30.72
      P4 = XCGRID(I,J+1)                                                  30.72
      XVC = (1.-T)*(1.-U)*P1+T*(1.-U)*P2+T*U*P3+(1.-T)*U*P4
!     *** For y-coord. ***
      P1 = YCGRID(I,J)                                                    30.72
      P2 = YCGRID(I+1,J)                                                  30.72
      P3 = YCGRID(I+1,J+1)                                                30.72
      P4 = YCGRID(I,J+1)                                                  30.72
      YVC = (1.-T)*(1.-U)*P1+T*(1.-U)*P2+T*U*P3+(1.-T)*U*P4
      RETURN
!     *** end of subroutine EVALF ***
      END
!***********************************************************************
!                                                                      *
      REAL FUNCTION SVALQI (XP, YP, IGRID, ARRINP, OUTVAL)                30.21
!                                                                      *
!***********************************************************************

      USE OCPCOMM4                                                        40.41
      USE SwashCommdata2                                                  40.41

      IMPLICIT NONE
!
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
!  0. Authors
!
!     30.60: Nico Booij
!     30.72: IJsbrand Haagsma
!     30.82: IJsbrand Haagsma
!     32.03: Nico Booij
!     40.30: Marcel Zijlema
!     40.41: Marcel Zijlema
!     40.80: Marcel Zijlema
!
!  1. Updates
!
!     30.72, Sept 97: INTEGER*4 replaced by INTEGER
!     30.60, Aug. 97: inequalities changed in view of bug reported by
!                     Ralf Kaiser (GT -> GE and LT -> LE)
!     32.03, Feb. 98: option for 1-D computation introduced
!                     real equality changed into inequality
!     30.82, Apr. 98: Replace statement with division through DYG to avoid division
!                     through zero in case of 1D.
!     30.82, Nov. 98: Now takes care of interpolation near points that
!                     contain exception values
!     40.30, Mar. 03: correcting indices IXC, IYC with offsets MXF, MYF
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!     40.80, Dec. 07: extension to unstructured grids
!
!  2. Purpose
!
!     Determining the value of a quantity from an input grid
!     such as depth and the current velocity components
!     for point given in problem coordinates
!
!  3. Method (updated...)
!
!     The required values are computed by bilinear interpolation. The
!     coordinates are given in the bottom grid as the number of meshes
!     in X- and Y-direction, IB and JB respectively (both real).
!
!           YB|
!             |
!             |--------------- *                *
!             | A
!             | |SYB1
!             | V
!         IYB-|-------------------- o
!             |
!             |                     |
!         JB1-|--------------- *    |           *
!             |                     |
!             |                |    |    SXB1   |
!             |                |    |<--------->|
!             +----------------------------------------------->
!                              |    |                       XB
!                             IB1  IXB
!
!                   *  bottom grid points
!                   o  point for interpolation
!
!  4. Argument variables
!
!     IGRID    Grid indicator
!
      INTEGER  IGRID
!
!     ARRINP   Array holding the values at the input grid locations
!     OUTVAL   If OUTVAL/=NEAREST, then value outside the input grid is
!              set to OUTVAL, otherwise the value is set to a value
!              nearest to boundary of the input grid
!     SVALQI   Value of quantity in (XP,YP)
!     XP       X-coordinate in computational gridpoint
!     YP       Y-coordinate in computational gridpoint
!
      REAL     ARRINP(*), OUTVAL, XP, YP
!
!  5. Parameter variables
!
!  6. Local variables
!
!     EQREAL   Boolean function which compares two REAL values
!     IB1      Grid counter in x-direction
!     IENT     Number of entries into this subroutine
!     II       Pointer number in ARRINP
!     INGRD    Boolean variable to determine whether point is in grid
!     IXB      Distance to origin in x-direction divided by meshsize
!              in x-direction
!     IYB      Distance to origin in y-direction divided by meshsize
!              in y-direction
!     JB1      Grid counter in y-direction
!     SUMWEXC  Sum of weight factors of points with exception value
!     SUMWREG  Sum of weight factors of points with regular value
!     SXB1     First weight factor for distance in x-direction
!     SXB2     Second weight factor for distance in x-direction
!     SYB1     First weight factor for distance in y-direction
!     SYB2     Second weight factor for distance in y-direction
!     WF1      Weight factor of point ARRINP(II)
!     WF2      Weight factor of point ARRINP(II+MXG(IGRID))
!     WF3      Weight factor of point ARRINP(II+1)
!     WF4      Weight factor of point ARRINP(II+1+MXG(IGRID))
!
      INTEGER  IB1, IENT, II, JB1
      REAL     IXB, IYB, SXB1, SXB2, SYB1, SYB2
      REAL     SUMWEXC, SUMWREG, WF1, WF2, WF3, WF4
      LOGICAL  INGRD
!
!  8. Subroutines used
!
!     LOGICAL FUNCTION EQREAL: Checks whether two reals are equal within certain margins
!     STRACE: Traces the entry into subroutines (test purposes)
!
      LOGICAL  EQREAL
!
!  9. Subroutines calling
!
!     SWDIM
!     INTEGER FUNCTION SIRAY
!     SWRBC
!     SNEXTI
!     FLFILE
!
! 10. Error messages
!
! 11. Remarks
!
! 12. Structure
!
!       ----------------------------------------------------------------
!       If the point is out of bottom grid in X-direction, then
!           Compute lines for interpolation and interpolation factors
!             such that the value at the side of the grid is taken
!       Else
!           Compute nearest line IX in the bottom grid and the interpo-
!             lation factor in X-direction
!       ----------------------------------------------------------------
!       If the point is out of bottom grid in Y-direction, then
!           Compute lines for interpolation and interpolation factors
!             such that the value at the side of the grid is taken
!       Else
!           Compute nearest line IY in the bottom grid and the interpo-
!             lation factor in Y-direction
!       ----------------------------------------------------------------
!       Compute pointer in arrays and interpolation factors in both
!       directions
!       Compute the depth to the reference level for the point
!       Add the water level to the depth
!       If depth > 0 and current is on, then
!           Interpolate X- and Y-component of current velocity
!       Else
!           Current components are zero
!       ----------------------------------------------------------------
!
! 13. Source text
!
      SAVE IENT                                                           30.72
      DATA IENT/0/
      CALL STRACE (IENT, 'SVALQI')
!
!     ***    Two different procedures in funcion of       ***
!     ***    grid type: rectilinear or curvilinear (staggered)*** ver 30.21
!
      IF (IGTYPE(IGRID) .EQ. 1) THEN                                      30.21
!
!     rectilinear input grid
!
        IXB = ( (XP-XPG(IGRID))*COSPG(IGRID) +
     &          (YP-YPG(IGRID))*SINPG(IGRID) ) / DXG(IGRID)
!
        INGRD = .TRUE.
        IF (IXB .LE. 0.) THEN                                             30.60
          IB1   = 1
          SXB2  = 0.
          IF (IXB.LT.-0.1) INGRD = .FALSE.                                20.3x
        ELSE IF (IXB .GE. FLOAT(MXG(IGRID)-1)) THEN                       30.60
          IB1   = MXG(IGRID)-1
          SXB2  = 1.
          IF (IXB.GT.FLOAT(MXG(IGRID))-0.9) INGRD = .FALSE.               20.3x
        ELSE
          IB1   = INT(IXB)
          SXB2  = IXB-REAL(IB1)
          IB1   = IB1+1
        ENDIF
        IF (MYG(IGRID).GT.1) THEN                                         32.03
          IYB = (-(XP-XPG(IGRID))*SINPG(IGRID) +                          30.82
     &            (YP-YPG(IGRID))*COSPG(IGRID) ) / DYG(IGRID)             30.82
          IF (IYB .LE. 0.) THEN                                           30.60
            JB1   = 1
            SYB2  = 0.
            IF (IYB.LT.-0.1) INGRD = .FALSE.
          ELSE IF (IYB .GE. FLOAT(MYG(IGRID)-1)) THEN                     30.60
            JB1   = MYG(IGRID)-1
            SYB2  = 1.
            IF (IYB.GT.FLOAT(MYG(IGRID))-0.9) INGRD = .FALSE.
          ELSE
            JB1   = INT(IYB)
            SYB2  = IYB-REAL(JB1)
            JB1   = JB1+1
          ENDIF
        ENDIF                                                             32.03
!
!       evaluate SVALQI (2D-mode):
!
        IF (.NOT.INGRD .AND. .NOT.EQREAL(OUTVAL,NEAREST)) THEN
          SVALQI = OUTVAL
        ELSE IF (MYG(IGRID).GT.1) THEN                                    32.03
          SXB1   = 1.- SXB2
          SYB1   = 1.- SYB2
          II     = IB1 + (JB1-1) * MXG(IGRID)
          WF1 = SXB1*SYB1
          WF2 = SXB1*SYB2
          WF3 = SXB2*SYB1
          WF4 = SXB2*SYB2
          SUMWEXC = 0.
          IF (EQREAL(ARRINP(II             ),EXCFLD(IGRID))) THEN
            SUMWEXC = SUMWEXC + WF1
            WF1 =0.
          ENDIF
          IF (EQREAL(ARRINP(II+  MXG(IGRID)),EXCFLD(IGRID))) THEN
            SUMWEXC = SUMWEXC + WF2
            WF2=0.
          ENDIF
          IF (EQREAL(ARRINP(II+1           ),EXCFLD(IGRID))) THEN
            SUMWEXC = SUMWEXC + WF3
            WF3=0.
          ENDIF
          IF (EQREAL(ARRINP(II+1+MXG(IGRID)),EXCFLD(IGRID))) THEN
            SUMWEXC = SUMWEXC + WF4
            WF4=0.
          ENDIF
          SUMWREG = 1. -SUMWEXC
!
          IF (SUMWEXC.GE.SUMWREG)   THEN
            SVALQI = EXCFLD(IGRID)
          ELSE
            SVALQI = ( WF1*ARRINP(II)   + WF2*ARRINP(II+MXG(IGRID))
     &               + WF3*ARRINP(II+1) + WF4*ARRINP(II+1+MXG(IGRID)) )
     &               / SUMWREG
          END IF
        ELSE
!
!       evaluate SVALQI (1D-mode):
!
          SXB1 = 1. - SXB2                                                32.03
          IF (EQREAL(ARRINP(IB1  ),EXCFLD(IGRID)).OR.                     30.82
     &        EQREAL(ARRINP(IB1+1),EXCFLD(IGRID))    ) THEN               30.82
!
!           One of the cornerpoints contains an exception value thus:     30.82
!
            SVALQI = EXCFLD(IGRID)                                        30.82
          ELSE                                                            30.82
            SVALQI = SXB1*ARRINP(IB1)                                     32.03
     &             + SXB2*ARRINP(IB1+1)                                   32.03
          ENDIF                                                           30.82
        ENDIF
      ELSEIF ( IGTYPE(IGRID).EQ.3 ) THEN                                  40.80
!
!     unstructured input grid
!
        CALL SwanInterpolatePoint(SVALQI, XP, YP, ARRINP, EXCFLD(IGRID))
!
      ELSE
!
!     curvilinear input grid                                              32.03
!
        WRITE (PRINTF,*) 'No interpolation in curvilinear input grid!'
!
      ENDIF
!
!     ***** test *****
      IF (ITEST .GE. 280)
!     &   WRITE(PRINTF, 6010) SVALQI,IGRID,XP,YP,IXB,IYB,II,ARRINP(II)
! 6010 FORMAT(' SVALQI  IGRID       XP      YP        IXB',
!     &       '       IYB  II  ARRINP(II)', /
!     &      ,E10.3,I3,1X,4E10.3,I4,E10.3)
     &   WRITE(PRINTF, 6010) XP, YP, IXB, IYB, SVALQI
 6010 FORMAT(' Test SVALQI:',5F10.3)
!
      RETURN
!     end of subroutine SVALQI
      END
!***********************************************************************
!                                                                      *
      SUBROUTINE SPRCON (XCGRID, YCGRID, KGRPNT)                          40.31 40.00
!                                                                      *
!***********************************************************************
!
      USE OCPCOMM4                                                        40.41
      USE SwashCommdata1                                                  40.41
      USE SwashCommdata2                                                  40.41
      USE SwashCommdata3                                                  40.41
      USE outp_data                                                       40.31
      USE SwanGriddata                                                    40.80
!
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
!  0. AUTHORS
!
!     30.72: IJsbrand Haagsma
!     32.02: Roeland Ris & Cor van der Schelde (1D-version)
!     32.03: Roeland Ris & Cor van der Schelde (1D-version)
!     40.02: IJsbrand Haagsma
!     40.31: Marcel Zijlema
!     40.41: Marcel Zijlema
!     40.80: Marcel Zijlema
!
!  1. Updates
!
!     30.72, Sept 97: INTEGER*4 replaced by INTEGER
!     30.72, Sept 97: Replaced DO-block with one CONTINUE to DO-block with
!                     two CONTINUE's
!     32.02, Jan. 98: Introduced 1D-version
!     32.03  Feb. 98: corrections processed
!     30.72, Feb. 98: Introduced generic names XCGRID, YCGRID and SPCSIG for SWAN
!     40.02, Feb. 00: Removed obsolescent DO-construct
!     40.02, Oct. 00: Avoided scalar/array conflict
!     40.31, Dec. 03: removing POOL construction
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!     40.80, Sep. 07: extension to unstructured grids
!
!  2. Purpose
!
!     Execution of some tests on the given model description
!
!  3. Method
!
!     This subroutine carries out the following tests:
!     - check if bottom and computational grid are defined (if MODIF=0)
!     - check on the location of the corner points of the computational
!       grid
!     - check on the location of output point sets
!
!  4. Argument variables
!
!     XCGRID: input  Coordinates of computational grid in x-direction     30.72
!     YCGRID: input  Coordinates of computational grid in y-direction     30.72
!
      REAL    :: XCGRID(MXC,MYC), YCGRID(MXC,MYC)                         30.72
!
!  6. Local variables
!
!     I, J    counters
!
      INTEGER I, J
!
!     XR      x (comp. grid coord.)
!     YR      y (comp. grid coord.)
!
      REAL    XR, YR
!
!  8. Subroutines used
!
!     COPYCH
!     MSGERR (both Ocean Pack)
!     SINBTG
!     SINUPT (both SWAN/SER)
!
!  9. Subroutines calling
!
!     SWREAD
!
! 10. Error messages
!
!     ---
!
! 11. Remarks
!
!     ---
!
! 12. Structure
!
!     ----------------------------------------------------------------
!     Check whether bottom grid and computational grid are defined
!     and if the bottom and current are read
!     For every corner point of the computational grid do
!         Compute problem grid coordinates
!         If corner point is outside bottom grid (SINBTG = false), then
!             Call MSGERR to generate a warning
!     If computational grid is rotating (ICOMP = 0), then
!         If the rotation point is out of the bottom grid, then
!             Call MSGERR to generate an error message
!     Else
!         Compute coordinates of the center of the computational grid
!         If the center is outside the bottom grid, then
!             Call MSGERR to generate an error message
!     Read number of output pointsets in array IOUTD
!     If number of pointsets is not 0, then
!         For every pointset do
!             Call COPYCH to read the name from IOUTD
!             If the pointset is not the bottom grid, comp. grid or set
!               of lines or places and recordlength > 0, then
!                 Read type of pointset
!                 If pointset is of type F (frame), then
!                     Check location of the cornerpoints in bottom grid
!                       and computational grid
!                 If pointset is of type C (curve), then
!                     Check locations of all end points of the curves
!                 If pointset is of type P (points), then
!                     Check location of all the points
!                 If pointset is of typr R (rays), then
!                     Check end points of first and last ray
!                 If pointset is of type G (grid), then
!                     Check location of the cornerpoints in bottom grid
!                       and computational grid
!     ----------------------------------------------------------------
!
! 13. Source text
!
      INTEGER  KGRPNT(MXC,MYC)                                            40.02
      LOGICAL  SINBTG
      CHARACTER STYPE *1
      TYPE(OPSDAT), POINTER :: CUOPS                                      40.31
      SAVE IENT
      DATA IENT /0/
      CALL STRACE (IENT,'SPRCON')
!
!     ***** test location of computational grid *****
!
      IF (OPTG .EQ. 1) THEN                                               30.21
!
!       rectilinear grid
!
        IF (ONED) THEN                                                    32.02
!
!         For 1D-version:
!         *** Check angles of bottom grid and computational grid ***      32.02
!
          IF ( ABS((ALPG(1) - ALPC)) .GT. 0.0017 ) THEN                   32.02
            CALL MSGERR( 2, ' Difference between angle of bottom grid'//  32.02
     &                      ' (alpinp) and computational grid (alpc)')    32.02
            CALL MSGERR( 2, ' greater than 0.1 degrees.')                 32.02
          ENDIF                                                           32.02
!
!         *** Check location of computational grid ***                    32.02
!
          DO  I=0,1                                                       30.72
            XR = I*XCLEN
            XP = XPC + XR*COSPC
            YP = YPC + XR*SINPC
            IF (.NOT.SINBTG(XP,YP) ) THEN
              CALL MSGERR(1,'Corner of comp grid outside bottom grid')
              WRITE (PRINTF, 6010) XP+XOFFS, YP+YOFFS
            ENDIF                                                         30.72
          ENDDO
        ELSE                                                              32.02
!
!         two-dimensional case
!
          DO 11 I=0,1                                                     30.72
            DO 10 J=0,1
              XR = I*XCLEN
              YR = J*YCLEN
              XP = XPC + XR*COSPC - YR*SINPC
              YP = YPC + XR*SINPC + YR*COSPC
              IF (.NOT.SINBTG(XP,YP) ) THEN
                CALL MSGERR(1,'Corner of comp grid outside bottom grid')
                WRITE (PRINTF, 6010) XP+XOFFS, YP+YOFFS
 6010           FORMAT (' Coordinates :',2F10.2)
              ENDIF                                                       30.72
   10       CONTINUE                                                      30.72
   11     CONTINUE
        ENDIF
!
        XR = 0.5*XCLEN
        YR = 0.5*YCLEN
        XP = XPC + XR*COSPC - YR*SINPC
        YP = YPC + XR*SINPC + YR*COSPC
        IF (.NOT. SINBTG(XP,YP) ) THEN
          CALL MSGERR (2,' Centre of comp. grid outside bottom grid')
        ENDIF
      ELSEIF (OPTG.EQ.5) THEN                                             40.80
!
!       --- check location of computational grid                          40.80
!
        DO I = 1, nverts                                                  40.80
           IF ( vmark(I) /= 0 ) THEN                                      40.80
              XP = xcugrd(I)                                              40.80
              YP = ycugrd(I)                                              40.80
              IF (.NOT.SINBTG(XP,YP) ) THEN                               40.80
                CALL MSGERR(1,'Corner of comp grid outside bottom grid')  40.80
                WRITE (PRINTF, 6010) XP+XOFFS, YP+YOFFS                   40.80
              ENDIF                                                       40.80
           ENDIF                                                          40.80
        ENDDO                                                             40.80
      ENDIF                                                               30.21
!
!     ***** test location of output pointsets *****
      IF (LOPS) THEN                                                      40.31
        CUOPS => FOPS                                                     40.31
        DO                                                                40.31
!         --- get name of point set                                       40.31
          SNAME = CUOPS%PSNAME                                            40.31
          IF (ITEST.GE.80) WRITE (PRTEST, 12) SNAME                       40.31
  12      FORMAT (' test SPRCON ', A8)                                    40.31
!
!         bottom grid is excluded from the test
          IF (SNAME.EQ.'BOTTGRID') GOTO 100
!
!         computational grid is excluded from the test
          IF (SNAME.EQ.'COMPGRID') GOTO 100
!
!         wind grid is excluded from the test
          IF (SNAME.EQ.'WXGRID' .OR. SNAME.EQ.'WYGRID') GOTO 100
!
!         velocity grid is excluded from the test
          IF (SNAME.EQ.'VXGRID' .OR. SNAME.EQ.'VYGRID') GOTO 100
!
!         waterlevel grid is excluded from the test
          IF (SNAME.EQ.'WLEVGRID') GOTO 100
!
!         friction grid is excluded from the test
          IF (SNAME.EQ.'FRICGRID') GOTO 100
!
!         ship grid is excluded from the test
          IF (SNAME.EQ.'SHIPGRID') GOTO 100
!
!         no grid is excluded from the test
          IF (SNAME.EQ.'NOGRID') GOTO 100
!
!         element mesh is excluded from the test
          IF (SNAME.EQ.'ELEMESH') GOTO 100
!
          STYPE = CUOPS%PSTYPE                                            40.31
!                                                                         32.02
!         *** Check other output locations ***                            32.02
!                                                                         32.02
          IF (STYPE.EQ.'F' .AND. OPTG .EQ. 1) THEN
!
!           check the four corners of the frame
!
            XQLEN = CUOPS%OPR(3)                                          40.31
            YQLEN = CUOPS%OPR(4)                                          40.31
            XPQ   = CUOPS%OPR(1)                                          40.31
            YPQ   = CUOPS%OPR(2)                                          40.31
            ALPQ  = CUOPS%OPR(5)                                          40.31
            COSPQ = COS(ALPQ)
            SINPQ = SIN(ALPQ)
            IF (ONED) THEN                                                32.02
              DO  I=0,1                                                   30.72
                XQ = I*XQLEN
                XP = XPQ + XQ*COSPQ
                YP = YPQ + XQ*SINPQ
                CALL SINUPT (SNAME, XP, YP, XCGRID, YCGRID, KGRPNT)       30.72
              ENDDO                                                       30.72
            ELSE                                                          32.02
              DO 21 I=0,1                                                 30.72
                DO 20 J=0,1
                  XQ = I*XQLEN
                  YQ = J*YQLEN
                  XP = XPQ + XQ*COSPQ - YQ*SINPQ
                  YP = YPQ + XQ*SINPQ + YQ*COSPQ
                  CALL SINUPT (SNAME, XP, YP, XCGRID, YCGRID, KGRPNT)     30.72
   20           CONTINUE                                                  30.72
   21         CONTINUE                                                    30.72
            ENDIF                                                         32.02
          ENDIF
!         --------------------------------------------------------------
          IF (STYPE .EQ. 'C') THEN
!
!           check first and last point of a curve
!
            MXK = CUOPS%MIP                                               40.31
            IF ( MXK/=0 ) THEN
              XP = CUOPS%XP(1)                                            40.31
              YP = CUOPS%YP(1)                                            40.31
              CALL SINUPT (SNAME, XP, YP, XCGRID, YCGRID, KGRPNT)         30.72
              XP = CUOPS%XP(MXK)
              YP = CUOPS%YP(MXK)
              CALL SINUPT (SNAME, XP, YP, XCGRID, YCGRID, KGRPNT)
            ENDIF
          ENDIF
!         --------------------------------------------------------------
          IF (STYPE .EQ. 'P' .OR. STYPE .EQ. 'U') THEN                    40.80
!
!           check all individual output points
!
            MIP = CUOPS%MIP                                               40.31
            DO 40 IP=1,MIP
              XP = CUOPS%XP(IP)                                           40.31
              YP = CUOPS%YP(IP)                                           40.31
              CALL SINUPT (SNAME, XP, YP, XCGRID, YCGRID, KGRPNT)         30.72
   40       CONTINUE
          ENDIF
!         --------------------------------------------------------------
          IF (STYPE .EQ. 'R') THEN
            MIP = CUOPS%MIP                                               40.31
            DO 50 IP=1,MIP, MIP
              XP = CUOPS%XP(IP)                                           40.31
              YP = CUOPS%YP(IP)                                           40.31
              XQ = CUOPS%XQ(IP)                                           40.31
              YQ = CUOPS%YQ(IP)                                           40.31
              CALL SINUPT (SNAME, XP, YP, XCGRID, YCGRID, KGRPNT)         30.72
              CALL SINUPT (SNAME, XQ, YQ, XCGRID, YCGRID, KGRPNT)         30.72
   50       CONTINUE
          ENDIF
!         --------------------------------------------------------------
!         stype = 'N'     VERSION 20.63
!
          IF (STYPE.EQ.'N') THEN
            MIP   = CUOPS%MIP                                             40.31
            XQLEN = CUOPS%OPR(1)                                          40.31
            IF (XQLEN.NE.-999.) THEN                                      40.80
!              nested grid is regular, check corners
               YQLEN = CUOPS%OPR(2)                                       40.31
               XPQ   = CUOPS%OPR(3)                                       40.31
               YPQ   = CUOPS%OPR(4)                                       40.31
               ALPQ  = CUOPS%OPR(5)                                       40.31
               COSPQ = COS(ALPQ)
               SINPQ = SIN(ALPQ)
               DO I=0,1                                                   40.02
                 DO J=0,1                                                 40.02
                   XQ = I*XQLEN
                   YQ = J*YQLEN
                   XP = XPQ + XQ*COSPQ - YQ*SINPQ
                   YP = YPQ + XQ*SINPQ + YQ*COSPQ
                   CALL SINUPT (SNAME, XP, YP, XCGRID, YCGRID, KGRPNT)    30.72
                 ENDDO                                                    40.02
               ENDDO                                                      40.02
            ELSE                                                          40.80
!              nested grid is unstructured, check whole outline
               DO IP=1,MIP                                                40.80
                  XP = CUOPS%XP(IP)                                       40.80
                  YP = CUOPS%YP(IP)                                       40.80
                  CALL SINUPT (SNAME, XP, YP, XCGRID, YCGRID, KGRPNT)     40.80
               ENDDO                                                      40.80
            ENDIF                                                         40.80
          ENDIF
!
  100     IF (.NOT.ASSOCIATED(CUOPS%NEXTOPS)) EXIT                        40.31
          CUOPS => CUOPS%NEXTOPS                                          40.31
        END DO                                                            40.31
      ENDIF
!
      RETURN
!   * end of subroutine SPRCON *
      END
!***********************************************************************
!                                                                      *
      SUBROUTINE SINUPT (PSNAME, XP, YP, XCGRID, YCGRID, KGRPNT)
!                                                                      *
!***********************************************************************
!
      USE OCPCOMM4                                                        40.41
      USE SwashCommdata2                                                  40.41
      USE SwashCommdata3                                                  40.41
!
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
!  0. Authors
!
!     30.72: IJsbrand Haagsma
!     40.41: Marcel Zijlema
!
!  1. Updates
!
!      0.0 , Mar. 87: Heading added, IF..GOTO.. changed into IF..THEN..
!     30.72, Feb. 98: Introduced generic names XCGRID, YCGRID and SPCSIG for SWAN
!     40.00, Feb. 99: test skipped for irregular bottom grid
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!
!  2. Purpose
!
!     Checking whether the point XP, YP (given in problem coordinates)
!     of the output pointset SNAME is located in the computational grid
!     and bottom grid or not. If not, a warning is generated.
!
!  3. Method
!
!     ---
!
!  4. Argument variables
!
!     KGRPNT: input  Adresses of the computational grid points
!
      INTEGER KGRPNT(MXC,MYC)
!
!     XCGRID: input  Coordinates of computational grid in x-direction     30.72
!     XP    : input  X-coordinate of the point (problem coordinates)
!     YCGRID: input  Coordinates of computational grid in y-direction     30.72
!     YP    : input  Y-coordinate of the point (problem coordinates)
!
      REAL    XP, YP
      REAL    XCGRID(MXC,MYC),    YCGRID(MXC,MYC)                         30.72
!
!     PSNAME: input  Name of the output pointset (any type)
!
      CHARACTER PSNAME *(*)
!
!  5. SUBROUTINES CALLING
!
!     SPRCON (SWAN/SWREAD)
!
!  6. SUBROUTINES USED
!
!     SINBTG, SINCMP (both SWAN/SER) and MSGERR (Ocean Pack)
!
      LOGICAL SINBTG, SINCMP
!
!  7. ERROR MESSAGES
!
!     ---
!
!  8. REMARKS
!
!     ---
!
!  9. STRUCTURE
!
!     ----------------------------------------------------------------
!     If point (XP,YP) is not in the bottom grid (SINBTG = FALSE), then
!         Call MSGERR to generate a warning
!     If point (XP,YP) is not in the comp. grid (SINCMP = FALSE), then
!         Call MSGERR to generate a warning
!     ----------------------------------------------------------------
!
! 10. SOURCE TEXT
!
!
      SAVE IENT
      DATA IENT/0/
      CALL STRACE(IENT,'SINUPT')
!
      IF (.NOT. SINBTG (XP,YP) ) THEN
         CALL MSGERR(1,'(corner)point outside bottom grid')
         WRITE (PRINTF, 6010) PSNAME, XP+XOFFS, YP+YOFFS
      ENDIF
      IF (.NOT.SINCMP (XP, YP, XCGRID, YCGRID, KGRPNT)) THEN
        CALL MSGERR(1,'(corner)point outside comp. grid')
        WRITE (PRINTF, 6010) PSNAME, XP+XOFFS, YP+YOFFS
      ENDIF
 6010 FORMAT('       Set of output locations: ',A8,
     &       '  coordinates:', 2F12.2)
!
      RETURN
!     end of subroutine SINUPT *
      END
!
!***********************************************************************
!                                                                      *
      LOGICAL FUNCTION SINBTG (XP, YP)
!                                                                      *
!***********************************************************************
!
      USE SwashCommdata2                                                  40.41
!
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
!  0. Authors
!
!     32.02: Roeland Ris & Cor van der Schelde (1D-version)
!     40.41: Marcel Zijlema
!
!  1. Updates
!
!      0.0 , Mar. 87: name of function changed from INBODP into SINBTG
!     32.02, Jan. 98: Introduced 1D-version
!     40.00, Feb. 99: 1D procedure simplified, tolerance introduced
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!
!  2. Purpose
!
!     Checking whether a point given in problem coordinates is in the
!     bottom grid (SINBTG = true) or not (SINBTG = false).
!
!  3. Method
!
!     ---
!
!  4. Argument variables
!
!     XP      REAL   input    X-coordinate (problem grid) of the point
!     YP      REAL   input    Y-coordinate (problem grid) of the point
!
!  6. Local variables
!
!     XB      x-coordinate (of bottom grid)
!     YB      y-coordinate (of bottom grid)
!     XLENB   length of bottom grid in x-direction (of bottom grid)
!     YLENB   length of bottom grid in y-direction (of bottom grid)
!     BTOL    tolerance length
!
      REAL    XB, YB, XLENB, YLENB, BTOL
!
!  8. Subroutines used
!
!     ---
!
!  9. Subroutines calling
!
!     SPRCON (SWAN/MAIN)
!     SINUPT (SWAN/SER)
!
! 10. Error messages
!
!     ---
!
! 11. Remarks
!
!     ---
!
! 12. Structure
!
!     ----------------------------------------------------------------
!     If the bottom grid is defined (DXB>0 and DYB>0), then
!         Compute coordinates XB,YB in the bottom grid
!         Give SINBTG initial value TRUE
!         If XB < 0, XB > X-length of grid, YB < 0 or YB > .. then
!            SINBTG is FALSE
!     ----------------------------------------------------------------
!
! 13. Source text
!
      DATA  IENT /0/
      CALL  STRACE (IENT,'SINBTG')
!
      SINBTG = .TRUE.
      IF ( IGTYPE(1).NE.1 ) RETURN
!
      XLENB = (MXG(1)-1)*DXG(1)
      YLENB = (MYG(1)-1)*DYG(1)
      BTOL  = 0.01 * (XLENB+YLENB)
!
!     ***** compute bottom grid coordinates from problem coordinates ****
!
      XB =  (XP-XPG(1))*COSPG(1) + (YP-YPG(1))*SINPG(1)
      YB = -(XP-XPG(1))*SINPG(1) + (YP-YPG(1))*COSPG(1)
!
!     ***** check location of point *****
      IF (XB .LT. -BTOL) SINBTG = .FALSE.
      IF (XB .GT. XLENB+BTOL) SINBTG = .FALSE.
      IF (YB .LT. -BTOL) SINBTG = .FALSE.
      IF (YB .GT. YLENB+BTOL) SINBTG = .FALSE.
!
      RETURN
!   * end of subroutine SINBTG *
      END
!***********************************************************************
!                                                                      *
      LOGICAL FUNCTION SINCMP (XP, YP ,XCGRID ,YCGRID ,KGRPNT)
!                                                                      *
!***********************************************************************
!
      USE SwashCommdata2                                                  40.41
      USE SwashCommdata3                                                  40.41
      USE M_PARALL                                                        40.31
!
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
!  0. Authors
!
!     30.72: IJsbrand Haagsma
!     32.02: Roeland Ris & Cor van der Schelde
!     30.60, 40.00: Nico Booij
!     40.41: Marcel Zijlema
!     40.80: Marcel Zijlema
!
!  1. Updates
!
!     00.00, Mar. 87: name changed from INREKP into SINCMP, heading added
!     30.60, Aug. 97: assignment of SINCMP moved
!     30.72, Sept 97: INTEGER*4 replaced by INTEGER
!     32.02, Jan. 98: Introduced 1D-version
!     30.72, Feb. 98: Introduced generic names XCGRID, YCGRID and SPCSIG for SWAN
!     40.00, June 98: argument KGRBND added, call CVMESH modified
!            Febr 99: separate 1D code removed, margin introduced
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!     40.80, Sep. 07: extension to unstructured grids
!
!  2. Purpose
!
!     Checking whether a point given in problem coordinates is in the
!     computational grid (SINCMP = true) or not (SINCMP = false).
!
!  3. Method
!
!     ---
!
!  4. Argument variables
!
!     KGRPNT  input  grid point addresses
!
      INTEGER KGRPNT(MXC,MYC)
!
!     XCGRID: input  Coordinates of computational grid in x-direction     30.72
!     XP      REAL   input    X-coordinate (problem grid) of the point
!     YCGRID: input  Coordinates of computational grid in y-direction     30.72
!     YP      REAL   input    Y-coordinate (problem grid) of the point
!
      REAL    XCGRID(MXC,MYC),    YCGRID(MXC,MYC)                         30.72
      REAL    XP,     YP
!
!  6. Local variables
!
!     CTOL    tolerance value (margin around comput. grid)                40.00
!
      REAL    CTOL
!
!  8. Subroutines used
!
!     ---
!
!  9. Subroutines calling
!
!     SINUPT
!
! 10. Error messages
!
!     ---
!
! 11. Remarks
!
!     ---
!
! 12. Structure
!
!     ----------------------------------------------------------------
!     Compute coordinates XC,YC in the computational grid
!     Give SINCMP initial value TRUE
!     If XC < 0, XC > XCLEN, YC < 0 or YC > YCLEN, then
!         SINCMP = FALSE
!     ----------------------------------------------------------------
!
! 13. Source text
!
      SAVE IENT
      DATA IENT/0/
      CALL STRACE(IENT,'SINCMP')
!
!     *** Different procedure depending on grid type **
      IF (OPTG .EQ. 1) THEN                                               30.21
!
!       rectilinear grid: compute comp. coordinates from problem coordinates
!
        XC   =  (XP-XPC)*COSPC+(YP-YPC)*SINPC
        YC   = -(XP-XPC)*SINPC+(YP-YPC)*COSPC
!       XC and YC are in m
!
!       ***** check for location *****
        SINCMP = .TRUE.
        CTOL   = 0.01 * (XCLEN+YCLEN)                                     40.00
        IF (XC .LT. -CTOL) SINCMP = .FALSE.                               40.00
        IF (XC .GT. XCLEN+CTOL) SINCMP = .FALSE.                          40.00
        IF (YC .LT. -CTOL) SINCMP = .FALSE.                               40.00
        IF (YC .GT. YCLEN+CTOL) SINCMP = .FALSE.                          40.00
      ELSE IF (OPTG .EQ. 3) THEN                                          40.80
!
!       curvilinear grid
!
        CALL CVMESH (XP, YP, XC, YC, KGRPNT, XCGRID ,YCGRID)
!       XC and YC are nondimensional; equivalent to grid index
!
!       ***** check for location *****
        SINCMP = .TRUE.
        IF (XC .LT. -0.01) SINCMP = .FALSE.                               40.00
        IF (XC .GT. REAL(MXC-2)+0.01) SINCMP = .FALSE.                    40.00
        IF (YC .LT. -0.01) SINCMP = .FALSE.                               40.00
        IF (YC .GT. REAL(MYC-2)+0.01) SINCMP = .FALSE.                    40.00
      ELSE IF (OPTG.EQ.5) THEN                                            40.80
!
!       unstructured grid
!
        SINCMP = .TRUE.                                                   40.80
        CALL SwanFindPoint ( XP, YP, K )                                  40.80
        IF ( K.LT.0 ) SINCMP = .FALSE.                                    40.80
      ENDIF
!
!     --- check if output location is in global subdomain                 40.31
!
      IF ( PARLL .AND. .NOT.SINCMP ) THEN                                 40.31
         IF ( XP.GE.XCGMIN .AND. XP.LE.XCGMAX .AND.                       40.31
     &        YP.GE.YCGMIN .AND. YP.LE.YCGMAX ) SINCMP = .TRUE.           40.31
      END IF                                                              40.31
!
      RETURN
!   * end of subroutine SINCMP *
      END
!***********************************************************************
!                                                                      *
      INTEGER FUNCTION SIRAY (DP, XP1, YP1, XP2, YP2, XX, YY, BOTDEP,     30.70
     &                        BOTLEV, WATLEV)                             30.70
!                                                                      *
!***********************************************************************
!
      USE OCPCOMM4                                                        40.41
      USE SwashCommdata2                                                  40.41
      USE SwashCommdata3                                                  40.41
!
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
!  0. AUTHORS
!
!     30.72: IJsbrand Haagsma
!     40.41: Marcel Zijlema
!
!  1. UPDATE
!
!     00.00, Mar. 87: heading added, name of routine changed from
!                     IRAAI in SIRAY
!     30.72, Oct. 97: logical function EQREAL introduced for floating point
!                     comparisons
!     30.70, Nov. 97: changed into INTEGER function
!                     test output added
!                     arguments BOTDEP, BOTLEV, WATLEV added
!     40.03, Nov. 99: X2= etc. moved out of IF-ENDIF group
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!
!  2. PURPOSE
!
!     Searching the first point on a ray where the depth is DP
!
!  3. METHOD
!
!     ---
!
!  4. PARAMETERLIST
!
!     DP      REAL   input    depth
!     XP1     REAL   input    X-coordinate start point of ray
!     YP1     REAL   input    Y-coordinate start point of ray
!     XP2     REAL   input    X-coordinate end point of ray
!     YP2     REAL   input    Y-coordinate end point of ray
!     XX      REAL   input    X-coordinate point with depth DP
!     YY      REAL   input    Y-coordinate point with depth DP
!
!  5. SUBROUTINES CALLING
!
!     SWREPS (SWAN/READ)
!
!  6. SUBROUTINES USED
!
!     SVALQI (SWAN/READ)
!
!  7. ERROR MESSAGES
!
!     ---
!
!  8. REMARKS
!
!     ---
!
!  9. STRUCTURE
!
!     ----------------------------------------------------------------
!     Give SIRAY initial value 0
!     Compute stepsize, raylength and number of steps along the ray
!     Compute bottom coordinates  of startpoint as number of meshes
!     Call SVALQI to interpolate depth in startpoint of ray
!     For every step along the ray do
!         Compute coordinates of the intermediate point in problem
!           grid and bottom grid
!         Call SVALQI to interpolate the depth for this point
!         If the required depth is in the interval, then
!             Compute coordinates of the point with depth DP
!             Set SIRAY 1
!         Else
!             Coordinates and depth at start of new interval are values
!               at end of old interval
!     ----------------------------------------------------------------
!  10. SOURCE TEXT
!
      LOGICAL   EQREAL, BOTDEP                                            30.72
      REAL      BOTLEV(*), WATLEV(*)                                      30.70
      SAVE IENT
      DATA IENT/0/
      CALL STRACE (IENT,'SIRAY')
!
      SIRAY   = 0
      DIFDEP = 1e+10
      DSTEP  = MIN(DXG(1), DYG(1))
      RAYLEN = SQRT ((XP2-XP1)*(XP2-XP1) + (YP2-YP1)*(YP2-YP1))
      NSTEP  = 1 + INT(1.5*RAYLEN/DSTEP + 0.5)                            30.70
!
      DO 10 JJ = 0, NSTEP                                                 30.70
        X3  = XP1 + REAL(JJ)*(XP2-XP1)/REAL(NSTEP)
        Y3  = YP1 + REAL(JJ)*(YP2-YP1)/REAL(NSTEP)
        IF (BOTDEP) THEN
          D3   = SVALQI (X3, Y3, 1, BOTLEV, NEAREST)                      30.70
        ELSE
          D3   = SVALQI (X3, Y3, 1, BOTLEV, NEAREST)                      30.70
          IF (LEDS(7).GE.2) THEN                                          30.70
             D3 = D3 + SVALQI (X3, Y3, 7, WATLEV, NEAREST)                30.70
          ELSE
             D3 = D3 + SWL
          ENDIF
        ENDIF
        IF (ITEST.GE.160) WRITE (PRTEST, 14) X3+XOFFS, Y3+YOFFS, D3       30.70
  14    FORMAT (' SIRAY, scan point', 2(1X,F8.0), 1X, F8.2)               30.70
        IF (ABS(D3-DP).LT.DIFDEP) THEN                                    10.20
          DIFDEP=ABS(D3-DP)                                               10.20
          JDMINMAX=JJ                                                     10.20
        ENDIF
        IF (JJ.GT.0) THEN
          IF ((DP-D2)*(DP-D3).LE.0) THEN                                  40.03
            IF (EQREAL(D2,D3)) THEN                                       30.72
              XX = X2
              YY = Y2
            ELSE
              XX = X2+(X3-X2)*(D2-DP)/(D2-D3)
              YY = Y2+(Y3-Y2)*(D2-DP)/(D2-D3)
            ENDIF
            SIRAY = 1
            GOTO 20
          ENDIF                                                           40.03
        ENDIF                                                             40.03
        X2 = X3
        Y2 = Y3
        D2 = D3
   10 CONTINUE
!
!     exact depth not found, take closest value:
!
      X3 = XP1 + REAL(JDMINMAX)*(XP2-XP1)/REAL(NSTEP)                     10.20
      Y3 = YP1 + REAL(JDMINMAX)*(YP2-YP1)/REAL(NSTEP)                     10.20
      XX = X3                                                             10.20
      YY = Y3                                                             10.20
!
  20  IF (ITEST.GE.140) WRITE (PRTEST, 24) XX+XOFFS, YY+YOFFS             30.70
  24  FORMAT (' SIRAY, result ', 2(1X,F8.0))                              30.70
      RETURN
! * end of function SIRAY *
      END
!************************************************************************
!                                                                       *
      SUBROUTINE SWNMPS (PSNAME, PSTYPE, MIP, IERR)                       40.31
!                                                                       *
!************************************************************************
!
      USE OCPCOMM1                                                        40.41
      USE OCPCOMM4                                                        40.41
      USE SwashCommdata1                                                  40.41
      USE outp_data                                                       40.31
!
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
!  0. Authors
!     40.41: Marcel Zijlema
!
!  1. UPDATE
!
!       Oct. 1996, ver. 30.50: new subr.
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!
!  2. PURPOSE
!
!       Read name of set of output points; get type and number of
!       points in the set
!
!  3. METHOD
!
!
!  4. PARAMETERLIST
!
!       PSNAME   char   output   name
!       PSTYPE   char   output   type
!       MIP      int    output   number of points
!
!  5. SUBROUTINES CALLING
!
!       SPREOQ
!
!  6. SUBROUTINES USED
!
!       INCSTR (Ocean Pack)
!
!  7. ERROR MESSAGES
!
!       ---
!
!  8. REMARKS
!
!       ---
!
!  9. STRUCTURE
!
!
! 10. SOURCE TEXT
!
      INTEGER   MIP                                                       40.31
      CHARACTER PSNAME *(*), PSTYPE *1                                    40.31
      TYPE(OPSDAT), POINTER :: CUOPS                                       40.31
      SAVE IENT
      DATA IENT/0/
      CALL STRACE (IENT,'SWNMPS')
!
      IERR = 0
      CALL INCSTR ('SNAME', PSNAME, 'STA', 'BOTTGRID')
      IF (LENCST.GT.8) CALL MSGERR (2, 'SNAME is too long')
      CUOPS => FOPS                                                       40.31
      DO                                                                  40.31
        IF (CUOPS%PSNAME.EQ.PSNAME) EXIT                                  40.31
        IF (.NOT.ASSOCIATED(CUOPS%NEXTOPS)) THEN                          40.31
           CALL MSGERR(2, 'Set of output locations is not known')         40.31
           GOTO 900                                                       40.31
        END IF                                                            40.31
        CUOPS => CUOPS%NEXTOPS                                            40.31
      END DO                                                              40.31
      PSTYPE = CUOPS%PSTYPE                                               40.31
      IF (PSTYPE.EQ.'F' .OR. PSTYPE.EQ.'H') THEN                          40.31
         MIP = CUOPS%OPI(1) * CUOPS%OPI(2)                                40.31
!        get direction of frame in case of coordinates plotting
         ALPQ = CUOPS%OPR(5)                                              40.31
      ELSE
         MIP  = CUOPS%MIP                                                 40.31
         ALPQ = 0.                                                        40.31
      ENDIF
 800  IF (ITEST.GE.100) WRITE (PRTEST, 802) PSNAME, PSTYPE, MIP
 802  FORMAT (' exit SWNMPS, name:', A8, ' type:', A1,
     &        '  num of p:', I5)
      RETURN
 900  PSTYPE = ' '
      MIP    = 0
      IERR   = 1
      RETURN
      END
!************************************************************************
!                                                                       *
      SUBROUTINE SVARTP (IVTYPE)
!                                                                       *
!************************************************************************
!
      USE SwashCommdata1                                                  40.41
!
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
!  0. Authors
!
!     32.02: Roeland Ris & Cor van der Schelde
!     40.03: Nico Booij
!
!  1. Updates
!
!     10.09, Aug. 94: output quantity RPER added
!     20.61, Sep. 95: quantities TM02 and FWID added
!     20.67, Dec. 95: FWID renamed FSPR (freq. spread)
!     32.02, Feb. 98: 1D-version introduced
!     40.00, Apr. 98: subr simplified using new array OVKEYW
!     40.03, Sep. 00: inconsistency with manual corrected
!
!  2. Purpose
!
!     Converting keyword into integer
!
!  3. Method
!
!     This subroutine determines an integer value indicating the
!     required output variable from the keyword denoting the same
!     for storage in array with output requests.
!
!  4. PARAMETERLIST
!
!     IVTYPE  INT    output   type number output variable
!
!  5. Subroutines calling
!
!     SPROUT
!
!  6. Subroutines used
!
!     KEYWIS (Ocean Pack)
!
! 10. Error messages
!
!     ---
!
! 11. Remarks
!
!     ---
!
! 12. Structure
!
!     -----------------------------------------------------------------
!     If the keyword is equal to given string, then
!         IVTYPE is given integer value
!     -----------------------------------------------------------------
!
! 13. Source text
!
      LOGICAL KEYWIS
      SAVE IENT
      DATA IENT/0/
      CALL STRACE (IENT,'SVARTP')
!
      IVTYPE  =  0
!
      CALL INKEYW ('STA', 'ZZZZ')
!     check if given keyword corresponds to output quantity
      DO IVT = NMOVAR, 1, -1
!       loop in reverse order to check more specific names first
!       e.g. VELK before VEL
        IF (KEYWIS (OVKEYW(IVT))) THEN                                    40.00
          IVTYPE = IVT                                                    40.00
          GOTO 40
        ENDIF
      ENDDO
!     keyword OUTPUT means that output times will be entered              40.00
      IF (KEYWIS ('OUT')) IVTYPE = 98                                     40.03
!     keyword ZZZZ means end of list of output quantities                 40.00
      IF (KEYWIS ('ZZZZ')) IVTYPE = 999
!
      IF (IVTYPE .EQ. 0) CALL WRNKEY
!
  40  RETURN
!     end of subroutine SVARTP *
      END
!*********************************************************************
!                                                                    *
      LOGICAL FUNCTION BOUNPT (IX,IY,KGRPNT)
!                                                                    *
!*********************************************************************
!
      USE SwashCommdata3                                                  40.41
!
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
!  0. Authors
!
!     40.00, 40.03: Nico Booij
!
!  1. UPDATE
!
!       Feb. 1998, ver. 40.00: new subroutine
!       40.03, Sep. 00: inconsistency with manual corrected
!
!  2. PURPOSE
!
!       determine whether a grid point is a point where a boundary condition
!       can be applied
!
!  3. METHOD
!
!
!  4. PARAMETERLIST
!
!       IX, IY  int   inp    grid point indices
!       KGRPNT  int   inp    indirect addresses of grid points
!
!  5. SUBROUTINES CALLING
!
!       BCFILE
!
!  6. SUBROUTINES USED
!
!       ---
!
!  7. ERROR MESSAGES
!
!       ---
!
!  8. REMARKS
!
!     KGRPNT(IX,IY)=1 means that (IX,IY) is not an active grid point
!
!  9. STRUCTURE
!
!     -----------------------------------------------------------------
!     Make BOUNPT = False
!     If the grid point is not active
!     Then return
!     -----------------------------------------------------------------
!     If grid point is on the outer boundary
!     Then make BOUNPT = True
!          return
!     -----------------------------------------------------------------
!     If a neighbouring grid point is inactive
!     Then make BOUNPT = True
!          return
!     -----------------------------------------------------------------
!
! 10. SOURCE TEXT
!
      INTEGER IX,IY,KGRPNT(MXC,MYC)
      SAVE IENT
      DATA IENT/0/
      CALL STRACE (IENT, 'BOUNPT')
!
      BOUNPT = .FALSE.
!
      IF (IX.LE.0)   RETURN
      IF (IY.LE.0)   RETURN
      IF (IX.GE.MXC) RETURN
      IF (IY.GE.MYC) RETURN
!
!     If the grid point is not active
!     Then return
!
      IF (KGRPNT(IX,IY).LE.1) RETURN
!
!     If grid point is on the outer boundary
!     Then make BOUNPT = True
!          return
!
      IF (IX.EQ.1) THEN
        BOUNPT = .TRUE.
        RETURN
      ENDIF
      IF (IX.EQ.MXC-1) THEN
        BOUNPT = .TRUE.
        RETURN
      ENDIF
      IF (IY.EQ.1) THEN
        BOUNPT = .TRUE.
        RETURN
      ENDIF
      IF (IY.EQ.MYC-1) THEN
        BOUNPT = .TRUE.
        RETURN
      ENDIF
!
!     If a neighbouring grid point is inactive
!     Then make BOUNPT = True
!          return
!
      IF (KGRPNT(IX-1,IY).LE.1) THEN
        BOUNPT = .TRUE.
        RETURN
      ENDIF
      IF (KGRPNT(IX+1,IY).LE.1) THEN
        BOUNPT = .TRUE.
        RETURN
      ENDIF
      IF (KGRPNT(IX,IY-1).LE.1) THEN
        BOUNPT = .TRUE.
        RETURN
      ENDIF
      IF (KGRPNT(IX,IY+1).LE.1) THEN
        BOUNPT = .TRUE.
        RETURN
      ENDIF
      RETURN
      END
!***********************************************************************
!                                                                      *
      SUBROUTINE SWIPOL (FINP, EXCVAL, XC, YC, MIP, FOUTP,                40.86
     &                   KGRPNT, DEP2)                                    40.86 40.00
!                                                                      *
!***********************************************************************
!
      USE OCPCOMM4                                                        40.41
      USE SwashCommdata3                                                  40.41
      USE SwashCommdata4                                                  40.41
!
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
!  0. Authors
!
!     40.13: Nico Booij
!     40.41: Marcel Zijlema
!     40.86: Nico Booij
!
!  1. UPDATE
!
!     40.00, July 98: no interpolation if one or more corners are dry
!                     argument DEP2 added
!                     margin around comp. grid introduced
!     40.13, Aug. 01: provision for repeating grid
!                     swcomm4.inc reactivated
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!     40.86, Feb. 08: interpolation over an obstacle prevented
!      8.01, Mar. 22: interpolation is allowed in case of temporarily dry points
!                     argument DEP2 is now superfluous
!
!  2. PURPOSE
!
!       Interpolate the function FINP to the point given by computational
!       grid coordinates XC and YC; result appears in array FOUTP
!
!  3. METHOD
!
!       This subroutine computes the contributions from surrounding
!       points to the function value in an output point. The points used
!       are indicated in the sketch below.
!
!                                                            Y
!           +-------------------------------------------------->
!           |
!           |       .       .       .       .       .       .
!           |
!           |
!           |
!           |       .       .       *       *       .       .
!           |
!           |
!           |                         o
!           |       .       .       *       *       .       .
!           |
!           |
!           |
!           |       .       .       .       .       .       .
!           |
!           |
!           |
!         X |       .       .       .       .       .       .
!           |
!           V
!
!                 *   point of the computational grid contributing
!                       to output point (o)
!                 .   other grid points
!
!  4. PARAMETERLIST
!
!       FINP    real a input    array of function values defined on the
!                               computational grid
!       EXCVAL  real   input    exception value (assigned if point is outside
!                               computational grid)
!       XC, YC  real a input    array containing computational grid coordinates
!                               of output points
!       MIP     INT    input    number of output points
!       FOUTP   real a output   array of interpolated values for the output
!                               points
!
!  5. SUBROUTINES CALLING
!
!       SwashQuanOutp
!
!  6. SUBROUTINES USED
!
!       ---
!
!  7. ERROR MESSAGES
!
!       ---
!
!  8. REMARKS
!
!       IINTPC=1: bilinear interpolation
!       IINTPC=2: higher order interpolation using functions G1 and G2
!
!  9. STRUCTURE
!
!       ----------------------------------------------------------------
!       For every output point do
!           If the output point is near line XCL then
!               Determine points contributing to the output point
!               Compute contribution to the projection of the output
!                 point on line XCL
!               Compute multiplication factor for interpolation in X-
!                 direction
!               Compute contribution for the output point
!               Add result to value of variable for the output point in
!                 array IFOP
!       ----------------------------------------------------------------
!
! 10. SOURCE TEXT
!
      REAL FINP(MCGRD), FOUTP(MIP), XC(MIP), YC(MIP), DEP2(MCGRD)         40.00
      LOGICAL OUTSID
      INTEGER  KGRPNT(MXC,MYC)                                            30.21
!
      REAL*8 :: WW(1:4)  ! Interpolation weights for the 4 corners        40.86
      REAL*8 :: SUMWW    ! sum of the weights                             40.86
      INTEGER :: JX(1:4), JY(1:4) ! grid counters for the 4 corners       40.86
      INTEGER :: INDX(1:4)     ! grid counters for the 4 corners          40.86
      INTEGER :: JC            ! corner counter                           40.86
!
      SAVE IENT
      DATA IENT /0/
      IF (LTRACE) CALL  STRACE (IENT, 'SWIPOL')
!
        IF (ITEST.GE.150) WRITE (PRTEST, 61)                             060997
  61    FORMAT ('   XC    , YC  ,',
     &  '   JX1, JY1, JX2,  JY2  SX1,  SY1, FOUTP(IP),',
     &  '    INDX1  INDX2  INDX3  INDX4')
!
      DO 100 IP=1,MIP
        IF (XC(IP) .LE. -0.5 .OR. YC(IP) .LE. -0.5) THEN                  40.00
          FOUTP(IP) = EXCVAL
          JX1   = 0
          JY1   = 0
          JX2   = 0
          JY2   = 0
          SX1   = 0.
          SX2   = 0.
          INDX(1:4) = 1                                                   40.86
          WW  (1:4) = 0.
          GOTO 80
        ENDIF                                                             30.21
        OUTSID = .FALSE.
        FOUTP(IP) = 0.
        JX1 = INT(XC(IP)+3.001) - 2
        JX2 = JX1 + 1
        SX2 = XC(IP) + 1. - FLOAT(JX1)
        SX1 = 1. - SX2
        IF (JX1.LT.0) OUTSID = .TRUE.
        IF (JX1.EQ.0) JX1 = 1
        IF (.NOT.LREPTX) THEN                                             40.13
          IF (JX1.GT.MXC-1) OUTSID = .TRUE.
          IF (JX1.EQ.MXC-1) JX2 = MXC-1
        ELSE
          IF (JX1.GT.MXC) OUTSID = .TRUE.
          IF (JX1.EQ.MXC) JX2 = MXC
        ENDIF
        IF (ONED) THEN
          JY1 = 1                                                         40.86
          JY2 = 1                                                         40.86
          SY1 = 0.5                                                       40.86
          SY2 = 0.5                                                       40.86
        ELSE
          JY1 = INT(YC(IP)+3.001) - 2
          JY2 = JY1 + 1
          SY2 = YC(IP) + 1. - FLOAT(JY1)
          SY1 = 1. - SY2
          IF (JY1.LT.0) OUTSID = .TRUE.
          IF (JY1.EQ.0) JY1 = 1
          IF (.NOT.LREPTY) THEN
             IF (JY1.GT.MYC-1) OUTSID = .TRUE.
             IF (JY1.EQ.MYC-1) JY2 = MYC-1
          ELSE
             IF (JY1.GT.MYC) OUTSID = .TRUE.
             IF (JY1.EQ.MYC) JY2 = MYC
          ENDIF
        ENDIF
        IF (OUTSID) THEN
          FOUTP(IP) = EXCVAL
        ELSE
          JX(1) = JX1                                                     40.86 30.21
          JY(1) = JY1                                                     40.86 30.21
          WW(1) = SX1*SY1                                                 40.86
          JX(2) = JX2                                                     40.86 30.21
          JY(2) = JY1                                                     40.86 30.21
          WW(2) = SX2*SY1                                                 40.86
          JX(3) = JX1                                                     40.86 30.21
          JY(3) = JY2                                                     40.86 30.21
          WW(3) = SX1*SY2                                                 40.86
          JX(4) = JX2                                                     40.86 30.21
          JY(4) = JY2                                                     40.86 30.21
          WW(4) = SX2*SY2                                                 40.86
          DO JC = 1, 4                                                    40.86
            INDX(JC) = KGRPNT(JX(JC),JY(JC))                              40.86 30.21
            IF (WW(JC).LT.0.01) THEN
              WW(JC) = 0.                                                 40.86
            ELSE
              IF (INDX(JC).LE.1) THEN                                     40.86
                WW(JC) = 0.                                               40.86
              ENDIF
            ENDIF
          ENDDO
          SUMWW = SUM(WW(1:4))                                            40.86
          IF (OUTSID) THEN
            FOUTP(IP) = EXCVAL
          ELSE
            IF (SUMWW.GT.0.1) THEN                                        40.86
              FOUTP(IP) = SUM(WW*FINP(INDX)) / SUMWW                      40.86
            ELSE
              FOUTP(IP) = EXCVAL                                          40.86
            ENDIF
          ENDIF
        ENDIF
  80    IF (ITEST.GE.150) WRITE (PRTEST, 82)
     &  XC(IP) , YC(IP) ,JX1, JY1, JX2,JY2, (WW(JC),JC=1,4),              40.86
     &  (INDX(JC), JC=1,4), (FINP(INDX(JC)), JC=1,4)                      40.86
  82    FORMAT (2(F7.1,1X),4I5, 4(1X,F5.2), 3X,4(2X,I5), 3X,              40.86
     &  4(1X,E9.3))                                                       40.86
 100  CONTINUE
!
      RETURN
! * end of subroutine SWIPOL *
      END
!*****************************************************************
!                                                                *
      SUBROUTINE CHGBAS (X1, X2, PERIOD, Y1, Y2, N1, N2,
     &                   ITEST, PRTEST)
!                                                                *
!*****************************************************************
!
!   --|-----------------------------------------------------------|--
!     |            Delft University of Technology                 |
!     | Faculty of Civil Engineering, Fluid Mechanics Group       |
!     | P.O. Box 5048,  2600 GA  Delft, the Netherlands           |
!     |                                                           |
!     | Authors :  G. van Vledder, N. Booij                       |
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
!  0. Update history
!
!       ver 20.48: also accomodates periodic variables such as directions
!
!  1. Purpose
!
!       change x-basis of a discretized y-function
!
!  2. Method
!
!     A piecewise constant representation of the functions is assumed
!
!     first boundaries of a cell in X1 are determined
!     then it is determined whether there are overlaps with cells
!     in X2. if so Y1*common length is added to Y2
!     Finally Y2 values are divided by cell lengths
!
!  3. Parameter list
!
!     Name    I/O  Type  Description
!
!     X1       i    ra   x-coordinates of input grid
!     X2       i    ra   x-coordinates of output grid
!     PERIOD   i    r    period, i.e. x-axis is periodic if period>0
!                        e.g. spectral directions
!     Y1       i    ra   function values of input grid
!     Y2       o    ra   function values of output grid
!     N1       i    i    number of x-values of input grid
!     N2       i    i    number of x-values of output grid
!
!  4. Subroutines used
!
!     ---
!
!  5. Error messages
!
!  6. Remarks
!
!       Cell boundaries in X1 are: X1A and X1B
!       X2 is assumed to be monotonically increasing; this is checked
!       X1 is assumed to be monotonous but not necessarily increasing
!
!  7. Structure
!
!       ------------------------------------------------------------------
!       Make all values of Y2 = 0
!       For each cell in X1 do
!           determine boundaries of cell in X1
!           --------------------------------------------------------------
!           For each cell in X2 do
!               determine overlap with cell in X1; limits: RLOW and RUPP
!               add to Y2: Y1 * length of overlapping interval
!       ------------------------------------------------------------------
!       For each cell in X2 do
!           divide Y2 value by cell length
!       ------------------------------------------------------------------
!
!  8. Source text
!
      INTEGER  I1, I2, N1, N2, ITEST, PRTEST
      REAL     X1(N1), Y1(N1), X2(N2), Y2(N2), PERIOD
      REAL     X1A, X1B, X2A, X2B, RLOW, RUPP
      LOGICAL  TWICE
      SAVE IENT
      DATA IENT /0/
      CALL STRACE (IENT, 'CHGBAS')
!
!     initialize output data
!
      DO I2 = 1, N2
        Y2(I2) = 0.
      ENDDO
      DO I2 = 2, N2
        IF (X2(I2).LE.X2(I2-1))
     &    CALL MSGERR (2, 'subr. CHGBAS: values of X2 not increasing')
      ENDDO
!     boundaries of the range in X2
      X2LO  = 1.5 * X2(1)  - 0.5 * X2(2)
      X2HI  = 1.5 * X2(N2) - 0.5 * X2(N2-1)
      TWICE = .FALSE.
!
!     loop over cells in X1
!
      DO 300 I1 = 1, N1
        IF (ABS(Y1(I1)) .LT. 1.E-20) GOTO 300
!
!       determine cell boundaries in X1
!
        IF (I1.EQ.1) THEN
          X1A = 1.5 * X1(1) - 0.5 * X1(2)
        ELSE
          X1A = 0.5 * (X1(I1) + X1(I1-1))
        ENDIF

        IF (I1.EQ.N1) THEN
          X1B = 1.5 * X1(N1) - 0.5 * X1(N1-1)
        ELSE
          X1B = 0.5 * (X1(I1) + X1(I1+1))
        ENDIF
!
!       swap X1A and X1B if X1A > X1B
!
        IF (X1A.GT.X1B) THEN
          RR  = X1A
          X1A = X1B
          X1B = RR
        ENDIF

        IF (PERIOD.LE.0.) THEN
          IF (X1A.GT.X2HI) GOTO 300
          IF (X1B.LT.X2LO) GOTO 300
        ELSE
!         X is periodic; move interval in X1 if necessary
          TWICE = .FALSE.
          IADD = 0
  60      IF (X1B.GT.X2HI) THEN
            X1A = X1A - PERIOD
            X1B = X1B - PERIOD
            IADD = IADD + 1
            IF (IADD.GT.99)
     &         CALL MSGERR (2, 'endless loop in CHGBAS')
            GOTO 60
          ENDIF
  70      IF (X1A.LT.X2LO) THEN
            X1A = X1A + PERIOD
            X1B = X1B + PERIOD
            IADD = IADD + 1
            IF (IADD.GT.99)
     &           CALL MSGERR (2, 'endless loop in CHGBAS')
            GOTO 70
          ENDIF
          IF (X1A.GT.X2HI) GOTO 300
          IF (X1B.LT.X2LO) GOTO 300
          IF (X1A.LT.X2LO .AND. X1A+PERIOD.LT.X2HI) TWICE = .TRUE.
          IF (X1B.GT.X2HI .AND. X1B-PERIOD.GT.X2LO) TWICE = .TRUE.
        ENDIF
!
!       loop over cells in X2
!
 100    DO 200 I2 = 1, N2

          IF (I2.EQ.1) THEN
            X2A = X2LO
          ELSE
            X2A = 0.5 * (X2(I2) + X2(I2-1))
          ENDIF

          IF (I2.EQ.N2) THEN
            X2B = X2HI
          ELSE
            X2B = 0.5 * (X2(I2) + X2(I2+1))
          ENDIF
!
!         (RLOW,RUPP) is overlapping interval of (X1A,X1B) and (X2A,X2B)
!
          IF (X1A.LT.X2B) THEN
            RLOW = MAX (X1A, X2A)
          ELSE
            GOTO 200
          ENDIF

          IF (X1B.GT.X2A) THEN
            RUPP = MIN (X1B, X2B)
          ELSE
            GOTO 200
          ENDIF

          IF (RUPP.LT.RLOW) THEN
            CALL MSGERR (3, 'interpolation error')
            WRITE (PRTEST, 140) I1, X1A, X1B, I2, X2A, X2B
 140        FORMAT (' I, XA, XB ', 2(I3, 2(1X,E12.4)))
          ELSE
            Y2(I2) = Y2(I2) + Y1(I1) * (RUPP-RLOW)
          ENDIF
 200    CONTINUE
!
!       Cell in X1 covers both ends of sector boundary
        IF (TWICE) THEN
          IF (X1A.LT.X2LO) THEN
             X1A = X1A + PERIOD
             X1B = X1B + PERIOD
          ENDIF
          IF (X1B.GT.X2HI) THEN
             X1A = X1A - PERIOD
             X1B = X1B - PERIOD
          ENDIF
          TWICE = .FALSE.
          GOTO 100
        ENDIF
 300  CONTINUE
!
      DO I2 = 1, N2
        IF (I2.EQ.1) THEN
          CELLEN = X2(2) - X2(1)
        ELSE IF (I2.EQ.N2) THEN
          CELLEN = X2(N2) - X2(N2-1)
        ELSE
          CELLEN = 0.5 * (X2(I2+1) - X2(I2-1))
        ENDIF
!       divide Y2 by cell length
        Y2(I2) = Y2(I2) / CELLEN
      ENDDO
      IF (ITEST.GE.160) THEN
        WRITE (PRTEST, 84) N1, N2
  84    FORMAT (' test CHGBAS ', 2I5)
        WRITE (PRTEST, 85) (X1(II), II = 1, N1)
        WRITE (PRTEST, 85) (Y1(II), II = 1, N1)
        WRITE (PRTEST, 85) (X2(II), II = 1, N2)
        WRITE (PRTEST, 85) (Y2(II), II = 1, N2)
  85    FORMAT (10 (1X,E10.3))
      ENDIF
!
      RETURN
      END
!***********************************************************************
!                                                                      *
      REAL FUNCTION DEGCNV (DEGREE)
!                                                                      *
!***********************************************************************
!
      USE SwashCommdata3                                                  40.41
!
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
!  0. Authors
!
!  1. UPDATE
!
!       SEP 1997: New for SWAN 32.01
!                 Cor van der Schelde - Delft Hydraulics
!       30.70, Feb. 98: test output suppressed (causes problem if subr
!                       is used during output
!
!  2. PURPOSE
!
!       Transform degrees from nautical to cartesian or vice versa.
!
!  3. METHOD
!
!       DEGCNV = 180 + dnorth - degree
!
!  4. PARAMETERLIST
!
!       DEGCNV      direction in cartesian or nautical degrees.
!       DEGREE      direction in nautical or cartesian degrees.
!
!  5. SUBROUTINES CALLING
!
!       ---
!
!  6. SUBROUTINES USED
!
!       NONE
!
!  7. ERROR MESSAGES
!
!       NONE
!
!  8. REMARKS
!
!           Nautical convention           Cartesian convention
!
!                    0                             90
!                    |                              |
!                    |                              |
!                    |                              |
!                    |                              |
!        270 --------+-------- 90       180 --------+-------- 0
!                    |                              |
!                    |                              |
!                    |                              |
!                    |                              |
!                   180                            270
!
!  9. STRUCTURE
!
!     ---------------------------------
!     IF (NAUTICAL DEGREES) THEN
!       CONVERT DEGREES
!     IF (DEGREES > 360 OR < 0) THEN
!       CORRECT DEGREES WITHIN 0 - 360
!     ---------------------------------
!
! 10. SOURCE TEXT
!
!***********************************************************************
!
      SAVE IENT
      DATA IENT /0/
      CALL STRACE(IENT,'DEGCNV')
!
      IF ( BNAUT ) THEN
        DEGCNV = 180. + DNORTH - DEGREE
      ELSE
        DEGCNV = DEGREE
      ENDIF
!
      IF (DEGCNV .GE. 360.) THEN
        DEGCNV = MOD (DEGCNV, 360.)
      ELSE IF (DEGCNV .LT. 0.) THEN
        DEGCNV = MOD (DEGCNV, 360.) + 360.
      ELSE
!       DEGCNV between 0 and 360; do nothing
      ENDIF
!
!
!     *** end of subroutine DEGCNV ***
!
      RETURN
      END
!
!***********************************************************************
!                                                                      *
      REAL FUNCTION ANGRAD (DEGREE)
!                                                                      *
!***********************************************************************
!
      USE SwashCommdata3                                                  40.41
!
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
!  0. Authors
!
!  1. UPDATE
!
!       SEP 1997: New for SWAN 32.01
!                 Cor van der Schelde - Delft Hydraulics
!       30.70, Feb. 98: test output suppressed (causes problem if subr
!                       is used during output
!
!  2. PURPOSE
!
!       Transform degrees to radians
!
!  3. METHOD
!
!       ANGRAD = DEGREE * PI / 180
!
!  4. PARAMETERLIST
!
!       ANGRAD      radians
!       DEGREE      degrees
!
!  5. SUBROUTINES CALLING
!
!       ---
!
!  6. SUBROUTINES USED
!
!       NONE
!
!  7. ERROR MESSAGES
!
!       NONE
!
!  8. REMARKS
!
!       NONE
!
!  9. STRUCTURE
!
!     ---------------------------------
!     ANGLE[radian] = ANGLE[degrees} * PI / 180
!     ---------------------------------
!
! 10. SOURCE TEXT
!
!***********************************************************************
!
      SAVE IENT
      DATA IENT /0/
      CALL STRACE(IENT,'ANGRAD')
!
      ANGRAD = DEGREE * PI / 180.
!
!
!     *** end of subroutine ANGRAD ***
!
      RETURN
      END
!
!***********************************************************************
!                                                                      *
      REAL FUNCTION ANGDEG (RADIAN)
!                                                                      *
!***********************************************************************
!
      USE SwashCommdata3                                                  40.41
!
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
!  0. Authors
!
!  1. UPDATE
!
!       SEP 1997: New for SWAN 32.01
!                 Cor van der Schelde - Delft Hydraulics
!       30.70, Feb. 98: test output suppressed (causes problem if subr
!                       is used during output
!
!  2. PURPOSE
!
!       Transform radians to degrees
!
!  3. METHOD
!
!       ANGDEG = RADIAN * 180 / PI
!
!  4. PARAMETERLIST
!
!       RADIAN      radians
!       ANGDEG      degrees
!
!  5. SUBROUTINES CALLING
!
!       ---
!
!  6. SUBROUTINES USED
!
!       NONE
!
!  7. ERROR MESSAGES
!
!       NONE
!
!  8. REMARKS
!
!       NONE
!
!  9. STRUCTURE
!
!     ---------------------------------
!     ANGLE[degrees] = ANGLE[radians} * 180 / PI
!     ---------------------------------
!
! 10. SOURCE TEXT
!
!***********************************************************************
!
      SAVE IENT
      DATA IENT /0/
      CALL STRACE(IENT,'ANGDEG')
!
      ANGDEG = RADIAN * 180. / PI
!
!
!     *** end of subroutine ANGDEG ***
!
      RETURN
      END
!********************************************************************
!                                                                   *
      REAL FUNCTION GAMMAF(XX)
!                                                                   *
!********************************************************************
!
!   Updates
!     ver 30.70, Oct 1997 by N.Booij: new subroutine
!
!   Purpose
!     Compute the transcendental function Gamma
!
!   Subroutines used
!     GAMMLN  (Numerical Recipes)
!
      REAL XX, YY, ABIG                                                   40.00
      SAVE IENT, ABIG
      DATA IENT /0/, ABIG /30./
      CALL STRACE (IENT, 'GAMMAF')
      YY = GAMMLN(XX)
      IF (YY.GT.ABIG) YY = ABIG
      IF (YY.LT.-ABIG) YY = -ABIG
      GAMMAF = EXP(YY)
      RETURN
      END
!********************************************************************
!                                                                   *
      FUNCTION GAMMLN(XX)
!                                                                   *
!********************************************************************
!
!   Method:
!     function is copied from: Press et al., "Numerical Recipes"
!
      DOUBLE PRECISION  COF(6),STP,HALF,ONE,FPF,X,TMP,SER
      DATA COF,STP/76.18009173D0,-86.50532033D0,24.01409822D0,
     &    -1.231739516D0,.120858003D-2,-.536382D-5,2.50662827465D0/
      DATA HALF,ONE,FPF/0.5D0,1.0D0,5.5D0/
      X=XX-ONE
      TMP=X+FPF
      TMP=(X+HALF)*LOG(TMP)-TMP
      SER=ONE
      DO 11 J=1,6
        X=X+ONE
        SER=SER+COF(J)/X
11    CONTINUE
      GAMMLN=TMP+LOG(STP*SER)
      RETURN
      END
!TIMG!****************************************************************
!TIMG!
!TIMG      SUBROUTINE SWTSTA (ITIMER)
!TIMG!
!TIMG!****************************************************************
!TIMG!
!TIMG      USE SwashTimecomm                                                   40.41
!TIMG      USE OCPCOMM4                                                        40.41
!TIMG!
!TIMG      IMPLICIT NONE
!TIMG!
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
!TIMG!
!TIMG!  0. Authors
!TIMG!
!TIMG!     40.23: Marcel Zijlema
!TIMG!     40.41: Marcel Zijlema
!TIMG!
!TIMG!  1. Updates
!TIMG!
!TIMG!     40.23, Aug. 02: New subroutine
!TIMG!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!TIMG!
!TIMG!  2. Purpose
!TIMG!
!TIMG!     Start timing
!TIMG!
!TIMG!  3. Method
!TIMG!
!TIMG!     Get cpu and wall-clock times and store
!TIMG!
!TIMG!  4. Argument variables
!TIMG!
!TIMG!     ITIMER      number of timer to be used
!TIMG!
!TIMG      INTEGER :: ITIMER
!TIMG!
!TIMG!  6. Local variables
!TIMG!
!TIMG!     C     :     clock count of processor
!TIMG!     I     :     index in LISTTM, loop variable
!TIMG!     IFOUND:     index in LISTTM, location of ITIMER
!TIMG!     IFREE :     index in LISTTM, first free position
!TIMG!     M     :     maximum clock count
!TIMG!     R     :     number of clock counts per second
!TIMG!F95!     TIMER :     current real cpu-time
!TIMG!     TIMER1:     current cpu-time used
!TIMG!     TIMER2:     current wall-clock time used
!TIMG!
!TIMG      INTEGER          :: I, IFOUND, IFREE
!TIMG      INTEGER          :: C, R, M
!TIMG!F95      REAL             :: TIMER
!TIMG      DOUBLE PRECISION :: TIMER1, TIMER2
!TIMG!
!TIMG!  7. Common blocks used
!TIMG!
!TIMG!
!TIMG!  8. Subroutines used
!TIMG!
!TIMG!F95!     CPU_TIME         Returns real value from cpu-time clock
!TIMG!     SYSTEM_CLOCK     Returns integer values from a real-time clock
!TIMG!
!TIMG!  9. Subroutines calling
!TIMG!
!TIMG!     SWMAIN, SWCOMP, SWOMPU
!TIMG!
!TIMG! 12. Structure
!TIMG!
!TIMG!     Get and store the cpu and wall-clock times
!TIMG!
!TIMG! 13. Source text
!TIMG!
!TIMG
!TIMG!
!TIMG!     --- check whether a valid timer number is given
!TIMG!
!TIMG      IF (ITIMER.LE.0 .OR. ITIMER.GT.NSECTM) THEN
!TIMG         WRITE(PRINTF,*) 'SWTSTA: ITIMER out of range: ',
!TIMG     &                   ITIMER, 1, NSECTM
!TIMG         STOP
!TIMG      END IF
!TIMG!
!TIMG!     --- check whether timing for ITIMER was started already,
!TIMG!         also determine first free location in LISTTM
!TIMG!
!TIMG      IFOUND=0
!TIMG      IFREE =0
!TIMG      I     =0
!TIMG 100  IF (I.LT.LASTTM .AND. (IFOUND.EQ.0 .OR. IFREE.EQ.0)) THEN
!TIMG         I=I+1
!TIMG         IF (LISTTM(I).EQ.ITIMER) THEN
!TIMG            IFOUND=I
!TIMG         END IF
!TIMG         IF (IFREE.EQ.0 .AND. LISTTM(I).EQ.-1) THEN
!TIMG            IFREE =I
!TIMG         END IF
!TIMG         GOTO 100
!TIMG      END IF
!TIMG
!TIMG      IF (IFOUND.EQ.0 .AND. IFREE.EQ.0 .AND. LASTTM.LT.MXTIMR) THEN
!TIMG         LASTTM=LASTTM+1
!TIMG         IFREE =LASTTM
!TIMG      END IF
!TIMG!
!TIMG!     --- produce warning if found in the list
!TIMG!
!TIMG      IF (IFOUND.GT.0) THEN
!TIMG         WRITE(PRINTF,*)
!TIMG     &      'SWTSTA: warning: previous timing for section ',
!TIMG     &      ITIMER,' not closed properly/will be ignored.'
!TIMG      END IF
!TIMG!
!TIMG!     --- produce error if not found and no free position available
!TIMG!
!TIMG      IF (IFOUND.EQ.0 .AND. IFREE.EQ.0) THEN
!TIMG         WRITE(PRINTF,*)
!TIMG     &      'SWTSTA: maximum number of simultaneous timers',
!TIMG     &      ' exceeded:',MXTIMR
!TIMG         STOP
!TIMG      END IF
!TIMG!
!TIMG!     --- register ITIMER in appropriate location of LISTTM
!TIMG!
!TIMG      IF (IFOUND.EQ.0) THEN
!TIMG         IFOUND=IFREE
!TIMG      END IF
!TIMG      LISTTM(IFOUND)=ITIMER
!TIMG!
!TIMG!     --- get current cpu/wall-clock time and store in TIMERS
!TIMG!
!TIMG      TIMER1=0D0
!TIMG!F95      CALL CPU_TIME (TIMER)
!TIMG!F95      TIMER1=DBLE(TIMER)
!TIMG      CALL SYSTEM_CLOCK (C,R,M)
!TIMG      TIMER2=DBLE(C)/DBLE(R)
!TIMG
!TIMG      TIMERS(IFOUND,1)=TIMER1
!TIMG      TIMERS(IFOUND,2)=TIMER2
!TIMG
!TIMG      RETURN
!TIMG      END
!TIMG!****************************************************************
!TIMG!
!TIMG      SUBROUTINE SWTSTO (ITIMER)
!TIMG!
!TIMG!****************************************************************
!TIMG!
!TIMG      USE SwashTimecomm                                                   40.41
!TIMG      USE OCPCOMM4                                                        40.41
!TIMG!
!TIMG      IMPLICIT NONE
!TIMG!
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
!TIMG!
!TIMG!  0. Authors
!TIMG!
!TIMG!     40.23: Marcel Zijlema
!TIMG!     40.41: Marcel Zijlema
!TIMG!
!TIMG!  1. Updates
!TIMG!
!TIMG!     40.23, Aug. 02: New subroutine
!TIMG!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!TIMG!
!TIMG!  2. Purpose
!TIMG!
!TIMG!     Stop timing
!TIMG!
!TIMG!  3. Method
!TIMG!
!TIMG!     Get cpu and wall-clock times and store
!TIMG!
!TIMG!  4. Argument variables
!TIMG!
!TIMG!     ITIMER      number of timer to be used
!TIMG!
!TIMG      INTEGER :: ITIMER
!TIMG!
!TIMG!  6. Local variables
!TIMG!
!TIMG!     C     :     clock count of processor
!TIMG!     I     :     index in LISTTM, loop variable
!TIMG!     IFOUND:     index in LISTTM, location of ITIMER
!TIMG!     M     :     maximum clock count
!TIMG!     R     :     number of clock counts per second
!TIMG!F95!     TIMER :     current real cpu-time
!TIMG!     TIMER1:     current cpu-time used
!TIMG!     TIMER2:     current wall-clock time used
!TIMG!
!TIMG      INTEGER          :: I, IFOUND
!TIMG      INTEGER          :: C, R, M
!TIMG!F95      REAL             :: TIMER
!TIMG      DOUBLE PRECISION :: TIMER1, TIMER2
!TIMG!
!TIMG!  7. Common blocks used
!TIMG!
!TIMG!
!TIMG!  8. Subroutines used
!TIMG!
!TIMG!F95!     CPU_TIME         Returns real value from cpu-time clock
!TIMG!     SYSTEM_CLOCK     Returns integer values from a real-time clock
!TIMG!
!TIMG!  9. Subroutines calling
!TIMG!
!TIMG!     SWMAIN, SWCOMP, SWOMPU
!TIMG!
!TIMG! 12. Structure
!TIMG!
!TIMG!     Get and store the cpu and wall-clock times
!TIMG!
!TIMG! 13. Source text
!TIMG!
!TIMG
!TIMG!
!TIMG!     --- check whether a valid timer number is given
!TIMG!
!TIMG      IF (ITIMER.LE.0 .OR. ITIMER.GT.NSECTM) THEN
!TIMG         WRITE(PRINTF,*) 'SWTSTO: ITIMER out of range: ',
!TIMG     &                   ITIMER, 1, NSECTM
!TIMG         STOP
!TIMG      END IF
!TIMG!
!TIMG!     --- check whether timing for ITIMER was started already,
!TIMG!         also determine first free location in LISTTM
!TIMG!
!TIMG      IFOUND=0
!TIMG      I     =0
!TIMG 100  IF (I.LT.LASTTM .AND. IFOUND.EQ.0) THEN
!TIMG         I=I+1
!TIMG         IF (LISTTM(I).EQ.ITIMER) THEN
!TIMG            IFOUND=I
!TIMG         END IF
!TIMG         GOTO 100
!TIMG      END IF
!TIMG!
!TIMG!     --- produce error if not found
!TIMG!
!TIMG      IF (IFOUND.EQ.0) THEN
!TIMG         WRITE(PRINTF,*)
!TIMG     &      'SWTSTO: section ',ITIMER,' not found',
!TIMG     &      ' in list of active timings'
!TIMG         STOP
!TIMG      END IF
!TIMG!
!TIMG!     --- get current cpu/wall-clock time
!TIMG!
!TIMG      TIMER1=0D0
!TIMG!F95      CALL CPU_TIME (TIMER)
!TIMG!F95      TIMER1=DBLE(TIMER)
!TIMG      CALL SYSTEM_CLOCK (C,R,M)
!TIMG      TIMER2=DBLE(C)/DBLE(R)
!TIMG!
!TIMG!     --- calculate elapsed time since start of timing,
!TIMG!         store in appropriate location in DCUMTM,
!TIMG!         increment number of timings for current section
!TIMG!
!TIMG      DCUMTM(ITIMER,1)=DCUMTM(ITIMER,1)+(TIMER1-TIMERS(IFOUND,1))
!TIMG      DCUMTM(ITIMER,2)=DCUMTM(ITIMER,2)+(TIMER2-TIMERS(IFOUND,2))
!TIMG      NCUMTM(ITIMER)  =NCUMTM(ITIMER)+1
!TIMG!
!TIMG!     --- free appropriate location of LISTTM,
!TIMG!         adjust last occupied position of LISTTM
!TIMG!
!TIMG      IF (IFOUND.GT.0) THEN
!TIMG         LISTTM(IFOUND)=-1
!TIMG      END IF
!TIMG 200  IF (LASTTM.GT.1 .AND. LISTTM(LASTTM).EQ.-1) THEN
!TIMG         LASTTM=LASTTM-1
!TIMG         GOTO 200
!TIMG      END IF
!TIMG      IF (LISTTM(LASTTM).EQ.-1) LASTTM=0
!TIMG
!TIMG      RETURN
!TIMG      END
!TIMG!****************************************************************
!TIMG!
!TIMG      SUBROUTINE SWPRTI
!TIMG!
!TIMG!****************************************************************
!TIMG!
!TIMG      USE OCPCOMM4                                                        40.41
!TIMG      USE SwashTimecomm                                                   40.41
!TIMG      USE M_PARALL                                                        40.31
!TIMG
!TIMG      IMPLICIT NONE
!TIMG!
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
!TIMG!
!TIMG!  0. Authors
!TIMG!
!TIMG!     40.23: Marcel Zijlema
!TIMG!     40.30: Marcel Zijlema
!TIMG!     40.41: Marcel Zijlema
!TIMG!
!TIMG!  1. Updates
!TIMG!
!TIMG!     40.23, Aug. 02: New subroutine
!TIMG!     40.30, Jan. 03: introduction distributed-memory approach using MPI
!TIMG!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!TIMG!
!TIMG!  2. Purpose
!TIMG!
!TIMG!     Print timings info
!TIMG!
!TIMG!  6. Local variables
!TIMG!
!TIMG!     IDEBUG:     level of timing output requested:
!TIMG!                 0 - no output for detailed timings
!TIMG!                 1 - aggregate output for detailed timings
!TIMG!                 2 - complete output for all detailed timings
!TIMG!     IENT  :     number of entries
!TIMG!     J     :     loop counter
!TIMG!     K     :     loop counter
!TIMG!     TABLE :     array for computing aggregate cpu- and wallclock-times
!TIMG!
!TIMG      INTEGER          :: IENT, J, K, IDEBUG
!TIMG      DOUBLE PRECISION :: TABLE(30,2)
!TIMG      PARAMETER (IDEBUG=0)
!TIMG!
!TIMG!  7. Common blocks used
!TIMG!
!TIMG!
!TIMG!  8. Subroutines used
!TIMG!
!TIMG!     STRACE           Tracing routine for debugging
!TIMG!
!TIMG!  9. Subroutines calling
!TIMG!
!TIMG!     SWMAIN (in SWANMAIN)
!TIMG!
!TIMG! 12. Structure
!TIMG!
!TIMG!     Compile table with overview of cpu/wall clock time used in
!TIMG!     important parts of SWAN and write to PRINT file
!TIMG!
!TIMG! 13. Source text
!TIMG!
!TIMG      SAVE IENT
!TIMG      DATA IENT/0/
!TIMG      IF (LTRACE) CALL STRACE (IENT,'SWPRTI')
!TIMG!
!TIMG!     --- compile table with overview of cpu/wall clock time used in
!TIMG!         important parts of SWAN and write to PRINT file
!TIMG!
!TIMG      IF ( ITEST.GE.1 .OR. IDEBUG.GE.1 ) THEN
!TIMG!
!TIMG!        --- initialise table to zero
!TIMG!
!TIMG         DO K = 1, 30
!TIMG            DO J = 1, 2
!TIMG               TABLE(K,J) = 0D0
!TIMG            END DO
!TIMG         END DO
!TIMG!
!TIMG!        --- compute times for basic blocks
!TIMG!
!TIMG         DO J = 1, 2
!TIMG!
!TIMG!           --- total run-time
!TIMG!
!TIMG            TABLE(1,J) = DCUMTM(1,J)
!TIMG!
!TIMG!           --- initialisation, reading, preparation:
!TIMG!
!TIMG            DO K = 2, 7
!TIMG               TABLE(2,J) = TABLE(2,J) + DCUMTM(K,J)
!TIMG            END DO
!TIMG!
!TIMG!           --- domain decomposition:
!TIMG!
!TIMG            TABLE(2,J) = TABLE(2,J) + DCUMTM(211,J)
!TIMG            TABLE(2,J) = TABLE(2,J) + DCUMTM(212,J)
!TIMG            TABLE(2,J) = TABLE(2,J) + DCUMTM(201,J)
!TIMG!
!TIMG!           --- total calculation including communication:
!TIMG!
!TIMG            TABLE(3,J) = TABLE(3,J) + DCUMTM(8,J)
!TIMG!
!TIMG!           --- output:
!TIMG!
!TIMG            TABLE(5,J) = TABLE(5,J) + DCUMTM(9,J)
!TIMG!
!TIMG!           --- exchanging data:
!TIMG!
!TIMG            TABLE(7,J) = TABLE(7,J) + DCUMTM(213,J)
!TIMG!
!TIMG!           --- solving system:
!TIMG!
!TIMG            TABLE(9,J) = TABLE(9,J) + DCUMTM(119,J)
!TIMG            TABLE(9,J) = TABLE(9,J) + DCUMTM(120,J)
!TIMG!
!TIMG!           --- global reductions:
!TIMG!
!TIMG            TABLE(10,J) = TABLE(10,J) + DCUMTM(202,J)
!TIMG!
!TIMG!           --- collecting data:
!TIMG!
!TIMG            TABLE(11,J) = TABLE(11,J) + DCUMTM(214,J)
!TIMG!
!TIMG!           --- setup:
!TIMG!
!TIMG            TABLE(12,J) = TABLE(12,J) + DCUMTM(106,J)
!TIMG!
!TIMG!           --- propagation velocities:
!TIMG!
!TIMG            TABLE(14,J) = TABLE(14,J) + DCUMTM(111,J)
!TIMG            TABLE(14,J) = TABLE(14,J) + DCUMTM(113,J)
!TIMG            TABLE(14,J) = TABLE(14,J) + DCUMTM(114,J)
!TIMG!
!TIMG!           --- x-y advection:
!TIMG!
!TIMG            TABLE(15,J) = TABLE(15,J) + DCUMTM(140,J)
!TIMG!
!TIMG!           --- sigma advection:
!TIMG!
!TIMG            TABLE(16,J) = TABLE(16,J) + DCUMTM(141,J)
!TIMG!
!TIMG!           --- theta advection:
!TIMG!
!TIMG            TABLE(17,J) = TABLE(17,J) + DCUMTM(142,J)
!TIMG!
!TIMG!           --- wind:
!TIMG!
!TIMG            TABLE(18,J) = TABLE(18,J) + DCUMTM(132,J)
!TIMG!
!TIMG!           --- whitecapping:
!TIMG!
!TIMG            TABLE(19,J) = TABLE(19,J) + DCUMTM(133,J)
!TIMG!
!TIMG!           --- bottom friction:
!TIMG!
!TIMG            TABLE(20,J) = TABLE(20,J) + DCUMTM(130,J)
!TIMG!
!TIMG!           --- wave breaking:
!TIMG!
!TIMG            TABLE(21,J) = TABLE(21,J) + DCUMTM(131,J)
!TIMG!
!TIMG!           --- quadruplets:
!TIMG!
!TIMG            TABLE(22,J) = TABLE(22,J) + DCUMTM(135,J)
!TIMG!
!TIMG!           --- triads:
!TIMG!
!TIMG            TABLE(23,J) = TABLE(23,J) + DCUMTM(134,J)
!TIMG!
!TIMG!           --- limiter:
!TIMG!
!TIMG            TABLE(24,J) = TABLE(24,J) + DCUMTM(122,J)
!TIMG!
!TIMG!           --- rescaling:
!TIMG!
!TIMG            TABLE(25,J) = TABLE(25,J) + DCUMTM(121,J)
!TIMG!
!TIMG!           --- reflections:
!TIMG!
!TIMG            TABLE(26,J) = TABLE(26,J) + DCUMTM(136,J)
!DIFFR!TIMG!
!DIFFR!TIMG!           --- diffraction:
!DIFFR!TIMG!
!DIFFR!TIMG            TABLE(27,J) = TABLE(27,J) + DCUMTM(137,J)
!MUD!TIMG!
!MUD!TIMG!           --- fluid mud:
!MUD!TIMG!
!MUD!TIMG            TABLE(28,J) = TABLE(28,J) + DCUMTM(138,J)
!VEGET!TIMG!
!VEGET!TIMG!           --- vegetation:
!VEGET!TIMG!
!VEGET!TIMG            TABLE(29,J) = TABLE(29,J) + DCUMTM(139,J)
!TIMG
!TIMG         END DO
!TIMG!
!TIMG!        --- add up times for some basic blocks
!TIMG!
!TIMG         DO J = 1, 2
!TIMG!
!TIMG!           --- total calculation:
!TIMG!
!TIMG            TABLE(3,J) = TABLE(3,J) - TABLE( 7,J)
!TIMG            TABLE(3,J) = TABLE(3,J) - TABLE(10,J)
!TIMG            IF ( TABLE(3,J).LT.0D0 ) TABLE(3,J) = 0D0
!TIMG!
!TIMG!           --- total communication:
!TIMG!                * exchanging data
!TIMG!                * global reductions
!TIMG!                * collecting data
!TIMG!
!TIMG            TABLE(4,J) = TABLE(4,J) + TABLE( 7,J)
!TIMG            TABLE(4,J) = TABLE(4,J) + TABLE(10,J)
!TIMG            TABLE(4,J) = TABLE(4,J) + TABLE(11,J)
!TIMG!
!TIMG!           --- total propagation:
!TIMG!                * velocities and derivatives
!TIMG!
!TIMG            TABLE(6,J) = TABLE(6,J) + TABLE(14,J)
!TIMG            TABLE(6,J) = TABLE(6,J) + TABLE(15,J)
!TIMG            TABLE(6,J) = TABLE(6,J) + TABLE(16,J)
!TIMG            TABLE(6,J) = TABLE(6,J) + TABLE(17,J)
!TIMG!
!TIMG!           --- sources:
!TIMG!                * wind, whitecapping, friction, breaking,
!TIMG!                * quadruplets, triads, limiter, rescaling,
!TIMG!                * reflections
!TIMG!
!TIMG            DO K = 18, 26
!TIMG               TABLE(8,J) = TABLE(8,J) + TABLE(K,J)
!TIMG            END DO
!DIFFR!TIMG!
!DIFFR!TIMG!                * diffraction
!DIFFR!TIMG!
!DIFFR!TIMG            TABLE(8,J) = TABLE(8,J) + TABLE(27,J)
!MUD!TIMG!
!MUD!TIMG!                * fluid mud
!MUD!TIMG!
!MUD!TIMG            TABLE(8,J) = TABLE(8,J) + TABLE(28,J)
!VEGET!TIMG!
!VEGET!TIMG!                * vegetation
!VEGET!TIMG!
!VEGET!TIMG            TABLE(8,J) = TABLE(8,J) + TABLE(29,J)
!TIMG!
!TIMG!           --- other computing:
!TIMG!
!TIMG            TABLE(13,J) = TABLE(13,J) + TABLE( 3,J)
!TIMG            TABLE(13,J) = TABLE(13,J) - TABLE( 6,J)
!TIMG            TABLE(13,J) = TABLE(13,J) - TABLE( 8,J)
!TIMG            TABLE(13,J) = TABLE(13,J) - TABLE( 9,J)
!TIMG            TABLE(13,J) = TABLE(13,J) - TABLE(12,J)
!TIMG            IF ( TABLE(13,J).LT.0D0 ) TABLE(13,J) = 0D0
!TIMG
!TIMG         END DO
!TIMG!
!TIMG!        --- print CPU-times used in important parts of SWAN
!TIMG!
!TIMG         WRITE(PRINTF,'(/)')
!TIMG         WRITE(PRINTF,110) INODE
!TIMG         WRITE(PRINTF,111) INODE
!TIMG         WRITE(PRINTF,110) INODE
!TIMG         WRITE(PRINTF,112) INODE
!TIMG         WRITE(PRINTF,110) INODE
!TIMG         WRITE(PRINTF,115) INODE,'total time:'       ,(TABLE(1,J),J=1,2)
!TIMG         WRITE(PRINTF,110) INODE
!TIMG         WRITE(PRINTF,115) INODE,'total pre-processing:',
!TIMG     &                                                (TABLE(2,j),j=1,2)
!TIMG         WRITE(PRINTF,115) INODE,'total calculation:',(TABLE(3,j),j=1,2)
!TIMG         WRITE(PRINTF,115) INODE,'total communication:',
!TIMG     &                                                (TABLE(4,j),j=1,2)
!TIMG         WRITE(PRINTF,115) INODE,'total post-processing:',
!TIMG     &                                                (TABLE(5,j),j=1,2)
!TIMG         WRITE(PRINTF,110) INODE
!TIMG         WRITE(PRINTF,113) INODE
!TIMG         WRITE(PRINTF,110) INODE
!TIMG         WRITE(PRINTF,115) INODE,'calc. propagation:',(TABLE(6,j),j=1,2)
!TIMG         WRITE(PRINTF,115) INODE,'exchanging data:'  ,(TABLE(7,j),j=1,2)
!TIMG         WRITE(PRINTF,115) INODE,'calc. sources:'    ,(TABLE(8,j),j=1,2)
!TIMG         WRITE(PRINTF,115) INODE,'solving system:'   ,(TABLE(9,j),j=1,2)
!TIMG         WRITE(PRINTF,115) INODE,'reductions:'      ,(TABLE(10,j),j=1,2)
!TIMG         WRITE(PRINTF,115) INODE,'collecting data:' ,(TABLE(11,j),j=1,2)
!TIMG         WRITE(PRINTF,115) INODE,'calc. setup:'     ,(TABLE(12,j),j=1,2)
!TIMG         WRITE(PRINTF,115) INODE,'other computing:' ,(TABLE(13,j),j=1,2)
!TIMG         WRITE(PRINTF,110) INODE
!TIMG         WRITE(PRINTF,114) INODE
!TIMG         WRITE(PRINTF,110) INODE
!TIMG         WRITE(PRINTF,115) INODE,'prop. velocities:',(TABLE(14,j),j=1,2)
!TIMG         WRITE(PRINTF,115) INODE,'x-y advection:'   ,(TABLE(15,j),j=1,2)
!TIMG         WRITE(PRINTF,115) INODE,'sigma advection:' ,(TABLE(16,j),j=1,2)
!TIMG         WRITE(PRINTF,115) INODE,'theta advection:' ,(TABLE(17,j),j=1,2)
!TIMG         WRITE(PRINTF,115) INODE,'wind:'            ,(TABLE(18,j),j=1,2)
!TIMG         WRITE(PRINTF,115) INODE,'whitecapping:'    ,(TABLE(19,j),j=1,2)
!TIMG         WRITE(PRINTF,115) INODE,'bottom friction:' ,(TABLE(20,j),j=1,2)
!MUD!TIMG         WRITE(PRINTF,115) INODE,'fluid mud:'       ,(TABLE(28,j),j=1,2)
!VEGET!TIMG         WRITE(PRINTF,115) INODE,'vegetation:'      ,(TABLE(29,j),j=1,2)
!TIMG         WRITE(PRINTF,115) INODE,'wave breaking:'   ,(TABLE(21,j),j=1,2)
!TIMG         WRITE(PRINTF,115) INODE,'quadruplets:'     ,(TABLE(22,j),j=1,2)
!TIMG         WRITE(PRINTF,115) INODE,'triads:'          ,(TABLE(23,j),j=1,2)
!TIMG         WRITE(PRINTF,115) INODE,'limiter:'         ,(TABLE(24,j),j=1,2)
!TIMG         WRITE(PRINTF,115) INODE,'rescaling:'       ,(TABLE(25,j),j=1,2)
!TIMG         WRITE(PRINTF,115) INODE,'reflections:'     ,(TABLE(26,j),j=1,2)
!DIFFR!TIMG         WRITE(PRINTF,115) INODE,'diffraction:'     ,(TABLE(27,j),j=1,2)
!TIMG
!TIMG      END IF
!TIMG
!TIMG      IF ( IDEBUG.GE.2 ) THEN
!TIMG         WRITE(PRINTF,120) INODE
!TIMG         DO J = 1, NSECTM
!TIMG            IF (NCUMTM(J).GT.0)
!TIMG     &         WRITE(PRINTF,121) INODE,J,DCUMTM(J,1),DCUMTM(J,2),
!TIMG     &                           NCUMTM(J)
!TIMG         END DO
!TIMG      END IF
!TIMG
!TIMG  110 FORMAT(i3,1x,'#')
!TIMG  111 FORMAT(i3,' # Details on timings of the simulation:')
!TIMG  112 FORMAT(i3,1x,'#',26x,'cpu-time',1x,'wall-clock')
!TIMG  113 FORMAT(i3,' # Splitting up calc. + comm. times:')
!TIMG  114 FORMAT(i3,' # Overview source contributions:')
!TIMG  115 FORMAT(i3,1x,'#',1x,a22,2f11.2)
!TIMG
!TIMG  120 FORMAT(/,i3,' #    item     cpu-time    real time     count')
!TIMG  121 FORMAT(i3,1x,'#',4x,i4,2f13.4,i10)
!TIMG
!TIMG      RETURN
!TIMG      END
!****************************************************************
!
      SUBROUTINE TXPBLA(TEXT,IF,IL)
!
!****************************************************************
!
      IMPLICIT NONE
!
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
!  0. Authors
!
!     40.23: Marcel Zijlema
!
!  1. Updates
!
!     40.23, Feb. 03: New subroutine
!
!  2. Purpose
!
!     determines the position of the first and the last non-blank
!     (or non-tabulator) character in the text-string
!
!  4. Argument variables
!
!     IF          position of the first non-blank character in TEXT
!     IL          position of the last non-blank character in TEXT
!     TEXT        text string
!
      INTEGER IF, IL
      CHARACTER*(*) TEXT
!
!  6. Local variables
!
!     FOUND :     TEXT is found or not
!     ITABVL:     integer value of tabulator character
!     LENTXT:     length of TEXT
!
      INTEGER LENTXT, ITABVL
      LOGICAL FOUND
!
! 12. Structure
!
!     Trivial.
!
! 13. Source text
!
!DOS      ITABVL = ICHAR('	')
!UNIX      ITABVL = ICHAR('	')
      LENTXT = LEN (TEXT)
      IF = 1
      FOUND = .FALSE.
  100 IF (IF .LE. LENTXT .AND. .NOT. FOUND) THEN
         IF (.NOT. (TEXT(IF:IF) .EQ. ' ' .OR.
     &              ICHAR(TEXT(IF:IF)) .EQ. ITABVL)) THEN
            FOUND = .TRUE.
         ELSE
            IF = IF + 1
         ENDIF
         GOTO 100
      ENDIF
      IL = LENTXT + 1
      FOUND = .FALSE.
  200 IF (IL .GT. 1 .AND. .NOT. FOUND) THEN
         IL = IL - 1
         IF (.NOT. (TEXT(IL:IL) .EQ. ' ' .OR.
     &              ICHAR(TEXT(IL:IL)) .EQ. ITABVL)) THEN
            FOUND = .TRUE.
         ENDIF
         GOTO 200
      ENDIF

      RETURN
      END
!****************************************************************
!
      CHARACTER*20 FUNCTION INTSTR ( IVAL )
!
!****************************************************************
!
      IMPLICIT NONE
!
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
!  0. Authors
!
!     40.23: Marcel Zijlema
!
!  1. Updates
!
!     40.23, Feb. 03: New subroutine
!
!  2. Purpose
!
!     Convert integer to string
!
!  4. Argument variables
!
!     IVAL        integer to be converted
!
      INTEGER IVAL
!
!  6. Local variables
!
!     CVAL  :     character represented an integer of mantisse
!     I     :     counter
!     IPOS  :     position in mantisse
!     IQUO  :     whole quotient
!
      INTEGER I, IPOS, IQUO
      CHARACTER*1, ALLOCATABLE :: CVAL(:)
!
! 12. Structure
!
!     Trivial.
!
! 13. Source text
!
      IPOS = 1
 100  CONTINUE
      IF (IVAL/10**IPOS.GE.1.) THEN
         IPOS = IPOS + 1
         GO TO 100
      END IF
      ALLOCATE(CVAL(IPOS))

      DO I=IPOS,1,-1
         IQUO=IVAL/10**(I-1)
         CVAL(IPOS-I+1)=CHAR(INT(IQUO)+48)
         IVAL=IVAL-IQUO*10**(I-1)
      END DO

      WRITE (INTSTR,*) (CVAL(I), I=1,IPOS)

      RETURN
      END
!****************************************************************
!
      CHARACTER*20 FUNCTION NUMSTR ( IVAL, RVAL, FORM )
!
!****************************************************************
!
      USE OCPCOMM4                                                        40.41
!
      IMPLICIT NONE
!
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
!  0. Authors
!
!     40.23: Marcel Zijlema
!     40.41: Marcel Zijlema
!
!  1. Updates
!
!     40.23, Feb. 03: New subroutine
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!
!  2. Purpose
!
!     Convert integer or real to string with given format
!
!  4. Argument variables
!
!     IVAL        integer to be converted
!     FORM        given format
!     RVAL        real to be converted
!
      INTEGER   IVAL
      REAL      RVAL
      CHARACTER FORM*20
!
!  6. Local variables
!
! 12. Structure
!
!     Trivial.
!
! 13. Source text
!
      IF ( IVAL.NE.INAN ) THEN
         WRITE (NUMSTR,FORM) IVAL
      ELSE IF ( RVAL.NE.RNAN ) THEN
         WRITE (NUMSTR,FORM) RVAL
      ELSE
         NUMSTR = ''
      END IF

      RETURN
      END
!****************************************************************
!
      SUBROUTINE SWCOPI ( IARR1, IARR2, LENGTH )
!
!****************************************************************
!
      USE OCPCOMM4                                                        40.41
!
      IMPLICIT NONE
!
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
!  0. Authors
!
!     40.23: Marcel Zijlema
!     40.41: Marcel Zijlema
!
!  1. Updates
!
!     40.23, Feb. 03: New subroutine
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!
!  2. Purpose
!
!     Copies integer array IARR1 to IARR2
!
!  3. Method
!
!     ---
!
!  4. Argument variables
!
!     IARR1       source array
!     IARR2       target array
!     LENGTH      array length
!
      INTEGER LENGTH
      INTEGER IARR1(LENGTH), IARR2(LENGTH)
!
!  6. Local variables
!
!     I     :     loop counter
!     IENT  :     number of entries
!
      INTEGER I, IENT
!
!  8. Subroutines used
!
!     MSGERR           Writes error message
!     STRACE           Tracing routine for debugging
!
!  9. Subroutines calling
!
!     ---
!
! 10. Error messages
!
!     ---
!
! 11. Remarks
!
!     ---
!
! 12. Structure
!
!     Trivial.
!
! 13. Source text
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'SWCOPI')

!     --- check array length

      IF ( LENGTH.LE.0 ) THEN
         CALL MSGERR( 3, 'Array length should be positive' )
      END IF

!     --- copy elements of array IARR1 to IARR2

      DO 10 I = 1, LENGTH
         IARR2(I) = IARR1(I)
  10  CONTINUE

      RETURN
      END
!****************************************************************
!
      SUBROUTINE SWCOPR ( ARR1, ARR2, LENGTH )
!
!****************************************************************
!
      USE OCPCOMM4                                                        40.41
!
      IMPLICIT NONE
!
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
!  0. Authors
!
!     40.23: Marcel Zijlema
!     40.41: Marcel Zijlema
!
!  1. Updates
!
!     40.23, Feb. 03: New subroutine
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!
!  2. Purpose
!
!     Copies real array ARR1 to ARR2
!
!  3. Method
!
!     ---
!
!  4. Argument variables
!
!     ARR1        source array
!     ARR2        target array
!     LENGTH      array length
!
      INTEGER LENGTH
      REAL    ARR1(LENGTH), ARR2(LENGTH)
!
!  6. Local variables
!
!     I     :     loop counter
!     IENT  :     number of entries
!
      INTEGER I, IENT
!
!  8. Subroutines used
!
!     MSGERR           Writes error message
!     STRACE           Tracing routine for debugging
!
!  9. Subroutines calling
!
!     ---
!
! 10. Error messages
!
!     ---
!
! 11. Remarks
!
!     ---
!
! 12. Structure
!
!     Trivial.
!
! 13. Source text
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'SWCOPR')

!     --- check array length

      IF ( LENGTH.LE.0 ) THEN
         CALL MSGERR( 3, 'Array length should be positive' )
      END IF

!     --- copy elements of array ARR1 to ARR2

      DO 10 I = 1, LENGTH
         ARR2(I) = ARR1(I)
  10  CONTINUE

      RETURN
      END
!MatL4!****************************************************************
!MatL4!
!MatL4      SUBROUTINE SWI2B ( IVAL, BVAL )
!MatL4!
!MatL4!****************************************************************
!MatL4!
!MatL4      USE OCPCOMM4                                                        40.41
!MatL4!
!MatL4      IMPLICIT NONE
!MatL4!
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
!MatL4!
!MatL4!  0. Authors
!MatL4!
!MatL4!     40.30: Marcel Zijlema
!MatL4!     40.41: Marcel Zijlema
!MatL4!
!MatL4!  1. Updates
!MatL4!
!MatL4!     40.30, May 03: New subroutine
!MatL4!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!MatL4!
!MatL4!  2. Purpose
!MatL4!
!MatL4!     Calculates 32-bit representation of an integer number
!MatL4!
!MatL4!  3. Method
!MatL4!
!MatL4!     The representation of an integer number is divided into 4 parts
!MatL4!     of 8 bits each, resulting in 32-bit word in memory. Generally,
!MatL4!     storage words are represented with bits counted from the right,
!MatL4!     making bit 0 the lower-order bit and bit 31 the high-order bit,
!MatL4!     which is also the sign bit.
!MatL4!
!MatL4!     The integer number is always an exact representation of an
!MatL4!     integer of value positive, negative, or zero. Each bit, except
!MatL4!     the leftmost bit, corresponds to the actual exponent as power
!MatL4!     of two.
!MatL4!
!MatL4!     For representing negative numbers, the method called
!MatL4!     "excess 2**(m - 1)" is used, which represents an m-bit number by
!MatL4!     storing it as the sum of itself and 2**(m - 1). For a 32-bit
!MatL4!     machine, m = 32. This results in a positive number, so the
!MatL4!     leftmost bit need to be reversed. This method is identical to the
!MatL4!     two's complement method.
!MatL4!
!MatL4!     An example:
!MatL4!
!MatL4!        the 32-bit representation of 5693 is
!MatL4!
!MatL4!        decimal    :     0        0       22       61
!MatL4!        hexidecimal:     0        0       16       3D
!MatL4!        binary     : 00000000 00000000 00010110 00111101
!MatL4!
!MatL4!        since,
!MatL4!
!MatL4!        5693 = 2^12 + 2^10 + 2^9 + 2^5 + 2^4 + 2^3 + 2^2 + 2^0
!MatL4!
!MatL4!  4. Argument variables
!MatL4!
!MatL4!     BVAL        a byte value as a part of the representation of
!MatL4!                 integer number
!MatL4!     IVAL        integer number
!MatL4!
!MatL4      INTEGER BVAL(4), IVAL
!MatL4!
!MatL4!  6. Local variables
!MatL4!
!MatL4!     I     :     loop counter
!MatL4!     IENT  :     number of entries
!MatL4!     IQUOT :     auxiliary integer with quotient
!MatL4!     M     :     maximal exponent number possible (for 32-bit machine, m=32)
!MatL4!
!MatL4      INTEGER I, IENT, IQUOT, M
!MatL4!
!MatL4! 12. Structure
!MatL4!
!MatL4!     initialise 4 parts of the representation
!MatL4!     if integer < 0, increased it by 2**(m-1)
!MatL4!     compute the actual part of the representation
!MatL4!     determine the sign bit
!MatL4!
!MatL4! 13. Source text
!MatL4!
!MatL4      SAVE IENT
!MatL4      DATA IENT/0/, M/32/
!MatL4      IF (LTRACE) CALL STRACE (IENT,'SWI2B')
!MatL4
!MatL4!     --- initialise 4 parts of the representation
!MatL4
!MatL4      DO 10 I = 1, 4
!MatL4         BVAL(I) = 0
!MatL4 10   CONTINUE
!MatL4
!MatL4      IQUOT = IVAL
!MatL4
!MatL4!     --- if integer < 0, increased it by 2*(m-1)
!MatL4
!MatL4      IF ( IVAL.LT.0 ) IQUOT = IQUOT + 2**(M-1)
!MatL4
!MatL4!     --- compute the actual part of the representation
!MatL4
!MatL4      DO 20 I = 4, 1, -1
!MatL4         BVAL(I) = MOD(IQUOT,256)
!MatL4         IQUOT = INT(IQUOT/256)
!MatL4 20   CONTINUE
!MatL4
!MatL4!     --- determine the sign bit
!MatL4
!MatL4      IF ( IVAL.LT.0 ) BVAL(1) = BVAL(1) + 128
!MatL4
!MatL4      RETURN
!MatL4      END
!MatL4!****************************************************************
!MatL4!
!MatL4      SUBROUTINE SWR2B ( RVAL, BVAL )
!MatL4!
!MatL4!****************************************************************
!MatL4!
!MatL4      USE OCPCOMM4                                                        40.41
!MatL4!
!MatL4      IMPLICIT NONE
!MatL4!
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
!MatL4!
!MatL4!  0. Authors
!MatL4!
!MatL4!     40.30: Marcel Zijlema
!MatL4!     40.41: Marcel Zijlema
!MatL4!
!MatL4!  1. Updates
!MatL4!
!MatL4!     40.30, May 03: New subroutine
!MatL4!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!MatL4!
!MatL4!  2. Purpose
!MatL4!
!MatL4!     Calculates 32-bit representation of a floating-point number
!MatL4!
!MatL4!  3. Method
!MatL4!
!MatL4!     The representation of a floating-point number is divided into 4
!MatL4!     parts of 8 bits each, resulting in 32-bit word in memory.
!MatL4!     Generally, storage words are represented with bits counted from
!MatL4!     the right, making bit 0 the lower-order bit and bit 31 the
!MatL4!     high-order bit, which is also the sign bit.
!MatL4!
!MatL4!     The floating-point number is a processor approximation. Its format
!MatL4!     has an 8-bit biased exponent and a 23-bit fraction or mantissa. The
!MatL4!     leftmost bit is the sign bit which is zero for plus and 1 for
!MatL4!     minus. The biased exponent equals the bias and the actual exponent
!MatL4!     (power of two) of the number. For a 32-bit machine, bias=127.
!MatL4!
!MatL4!     Furthermore, the floating-point number is usually stored in the
!MatL4!     normalized form, i.e. it has a binary point to the left of the
!MatL4!     mantissa and an implied leading 1 to the left of the binary point.
!MatL4!     Thus, if X is a floating-point number, then it is calculated as
!MatL4!     follows:
!MatL4!
!MatL4!         X = (-1)**sign bit 1.fraction * 2**(biased exponent-bias)
!MatL4!
!MatL4!     There are several exceptions. Let a fraction, biased exponent
!MatL4!     and sign bit be denoted as F, E and S, respectively. The following
!MatL4!     formats adhere to IEEE standard:
!MatL4!
!MatL4!     S = 0, E = 00000000 and F  = 00 ... 0 : X = 0
!MatL4!     S = 0, E = 00000000 and F <> 00 ... 0 : X = +0.fraction * 2**(1-bias)
!MatL4!     S = 1, E = 00000000 and F <> 00 ... 0 : X = -0.fraction * 2**(1-bias)
!MatL4!     S = 0, E = 11111111 and F  = 00 ... 0 : X = +Inf
!MatL4!     S = 1, E = 11111111 and F  = 00 ... 0 : X = -Inf
!MatL4!     S = 0, E = 11111111 and F <> 00 ... 0 : X = NaN
!MatL4!
!MatL4!     A NaN (Not a Number) is a value reserved for signalling an
!MatL4!     attempted invalid operation, like 0/0. Its representation
!MatL4!     equals the representation of +Inf plus 1, i.e. 2**31 - 2**23 + 1
!MatL4!
!MatL4!     An example:
!MatL4!
!MatL4!        the 32-bit representation of 23.1 is
!MatL4!
!MatL4!        decimal    :    65      184      204      205
!MatL4!        hexidecimal:    41       B8       CC       CD
!MatL4!        binary     : 01000001 10111000 11001100 11001101
!MatL4!
!MatL4!        since,
!MatL4!
!MatL4!        23.1 = 2^4 + 2^2 + 2^1 + 2^0 + 2^-4 + 2^-5 + 2^-8 + 2^-9 +
!MatL4!               2^-12 + 2^-13 + 2^-16 + 2^-17 + 2^-19
!MatL4!
!MatL4!        so that the biased exponent = 4 + 127 = 131 = 10000011 = E
!MatL4!        and the sign bit = 0 = S. The remaining of the 32-bit word is
!MatL4!        the fraction, which is
!MatL4!
!MatL4!    3 2 1 0 -1 -2 -3 -4 -5 -6 -7 -8 -9 -10 -11 -12 -13 -14 -15 -16 -17 -18 -19
!MatL4!
!MatL4! F= 0 1 1 1  0  0  0  1  1  0  0  1  1   0   0   1   1   0   0   1   1   0   1
!MatL4!
!MatL4!  4. Argument variables
!MatL4!
!MatL4!     BVAL        a byte value as a part of the representation of
!MatL4!                 floating-point number
!MatL4!     RVAL        floating-point number
!MatL4!
!MatL4      INTEGER BVAL(4)
!MatL4      REAL    RVAL
!MatL4!
!MatL4!  6. Local variables
!MatL4!
!MatL4!     ACTEXP:     actual exponent in the representation
!MatL4!     BEXPO :     biased exponent
!MatL4!     BIAS  :     bias (for 32-bit machine, bias=127)
!MatL4!     EXPO  :     calculated exponent of floating-point number
!MatL4!     FRAC  :     fraction of floating-point number
!MatL4!     I     :     loop counter
!MatL4!     IENT  :     number of entries
!MatL4!     IPART :     i-the part of the representation
!MatL4!     IQUOT :     auxiliary integer with quotient
!MatL4!     LEADNR:     leading number of floating-point number
!MatL4!     LFRAC :     length of fraction in representation
!MatL4!           :     (for 32-bit machine, lfrac=23)
!MatL4!     RFRAC :     auxiliary real with fraction
!MatL4!
!MatL4      INTEGER ACTEXP, EXPO, I, IENT, IPART, LFRAC, BIAS, BEXPO
!MatL4      INTEGER*8 LEADNR, IQUOT
!MatL4      REAL FRAC, RFRAC, MAXREAL
!MatL4!
!MatL4! 12. Structure
!MatL4!
!MatL4!     initialise 4 parts of the representation and biased exponent
!MatL4!     determine leading number and fraction
!MatL4!     do while leading number >= 1
!MatL4!        calculate positive exponent as power of two
!MatL4!     or while fraction > 0
!MatL4!        calculate negative exponent as power of two
!MatL4!     end do
!MatL4!     compute the actual part of the representation
!MatL4!
!MatL4! 13. Source text
!MatL4!
!MatL4      SAVE IENT
!MatL4      DATA IENT/0/, LFRAC/23/, BIAS/127/
!MatL4      IF (LTRACE) CALL STRACE (IENT,'SWR2B')
!MatL4
!MatL4!     --- initialise 4 parts of the representation and biased exponent
!MatL4
!MatL4      DO 10 I = 1, 4
!MatL4         BVAL(I) = 0
!MatL4 10   CONTINUE
!MatL4      BEXPO  = -1
!MatL4
!MatL4!     --- clamp reals to the max integer*8, to prevent overflow
!MatL4
!MatL4      MAXREAL=9.0E+18
!MatL4      IF ( RVAL .GT. MAXREAL ) THEN
!MatL4         RVAL=MAXREAL
!MatL4      ELSEIF ( RVAL .LT. -MAXREAL ) THEN
!MatL4         RVAL=-MAXREAL
!MatL4      END IF
!MatL4
!MatL4!     --- determine leading number and fraction
!MatL4
!MatL4      IF ( ABS(RVAL).LT.1.E-7 ) THEN
!MatL4         LEADNR = 0
!MatL4         FRAC   = 0.
!MatL4      ELSE
!MatL4         LEADNR = INT(ABS(RVAL),KIND=8)
!MatL4         FRAC   = ABS(RVAL) - REAL(LEADNR)
!MatL4      END IF
!MatL4
!MatL4 20   IF ( LEADNR.GE.1 ) THEN
!MatL4
!MatL4!        --- calculate positive exponent as power of two
!MatL4
!MatL4         EXPO  = 0
!MatL4         IQUOT = LEADNR
!MatL4 30      IF ( IQUOT.GE.2 ) THEN
!MatL4
!MatL4              IQUOT = INT(IQUOT/2,KIND=8)
!MatL4              EXPO  = EXPO + 1
!MatL4
!MatL4         GOTO 30
!MatL4         END IF
!MatL4
!MatL4      ELSE IF ( FRAC.GT.0. ) THEN
!MatL4
!MatL4!        --- calculate negative exponent as power of two
!MatL4
!MatL4         EXPO = 0
!MatL4         RFRAC = FRAC
!MatL4 40      IF ( RFRAC.LT.1. ) THEN
!MatL4
!MatL4              RFRAC = RFRAC * 2.
!MatL4              EXPO  = EXPO - 1
!MatL4
!MatL4         GOTO 40
!MatL4         END IF
!MatL4
!MatL4      ELSE
!MatL4
!MatL4         GOTO 50
!MatL4
!MatL4      END IF
!MatL4
!MatL4!     --- compute the actual part of the representation
!MatL4
!MatL4      IF ( BEXPO.EQ.-1 ) THEN
!MatL4
!MatL4!        --- determine biased exponent
!MatL4
!MatL4         BEXPO = EXPO + BIAS
!MatL4
!MatL4!        --- the first seven bits of biased exponent belong
!MatL4!            to first part of the representation
!MatL4
!MatL4         BVAL(1) = INT(BEXPO/2)
!MatL4
!MatL4!        --- determine the sign bit
!MatL4
!MatL4         IF ( RVAL.LT.0. ) BVAL(1) = BVAL(1) + 128
!MatL4
!MatL4!        --- the eighth bit of biased component is the leftmost
!MatL4!            bit of second part of the representation
!MatL4
!MatL4         BVAL(2) = MOD(BEXPO,2)*2**7
!MatL4         IPART = 2
!MatL4
!MatL4      ELSE
!MatL4
!MatL4!        --- compute the actual exponent of bit 1 in i-th part of
!MatL4!            the representation
!MatL4
!MatL4         ACTEXP = (IPART-2)*8 + 7 - BEXPO + BIAS + EXPO
!MatL4         IF ( ACTEXP.LT.0 ) THEN
!MatL4            ACTEXP = ACTEXP + 8
!MatL4            IPART = IPART + 1
!MatL4            IF ( IPART.GT.4 ) GOTO 50
!MatL4         END IF
!MatL4         BVAL(IPART) = BVAL(IPART) + 2**ACTEXP
!MatL4
!MatL4      END IF
!MatL4
!MatL4      IF ( EXPO.LT.(BEXPO-BIAS-LFRAC) ) GOTO 50
!MatL4      LEADNR = LEADNR - 2.**EXPO
!MatL4      IF ( EXPO.LT.0 ) FRAC = FRAC - 2.**EXPO
!MatL4
!MatL4      GOTO 20
!MatL4 50   CONTINUE
!MatL4
!MatL4      RETURN
!MatL4      END
!****************************************************************
!
      SUBROUTINE MKPATH ( PATH, IERR )
!
!****************************************************************
!
      USE OCPCOMM4
!
      IMPLICIT NONE
!
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
!  0. Authors
!
!     41.95: Marcel Zijlema
!
!  1. Updates
!
!     41.95, Jul. 22: New subroutine
!
!  2. Purpose
!
!     Creates a directory on OS (e.g. Windows, Linux and macOS)
!
!  3. Method
!
!     Use of a Fortran 2008 standard EXECUTE_COMMAND_LINE
!
!  4. Argument variables
!
!     IERR  :     status error
!                 =0 : creating path successful
!                 /=0: creating path failed
!     PATH  :     string to pass path
!
      INTEGER          :: IERR
      CHARACTER(LEN=*) :: PATH
!
!  6. Local variables
!
!     CSTAT :     command status
!     CMSG  :     command error message
!     ESTAT :     exit status
!     IENT  :     number of entries
!     MSGSTR:     string to pass message
!
      INTEGER            :: CSTAT, ESTAT, IENT
!
      CHARACTER(LEN=100) :: CMSG
      CHARACTER(LEN=140) :: MSGSTR
!
! 13. Source text
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'MKPATH')
!
      IERR = 0
!
      CALL EXECUTE_COMMAND_LINE('mkdir '//TRIM(PATH), EXITSTAT=ESTAT,
     &                          CMDSTAT=CSTAT, CMDMSG=CMSG)
      IF (CSTAT.GT.0) THEN
         WRITE (MSGSTR,'(A)') 'Command execution failed with error '//
     &                        TRIM(CMSG)
         CALL MSGERR( 1, TRIM(MSGSTR) )
         IERR = 1
      ELSE IF (CSTAT.LT.0) THEN
         CALL MSGERR( 2, ' Command execution not supported' )
         IERR = 2
      ELSE IF (ESTAT.NE.0) THEN
         WRITE (MSGSTR, '(A,I5)')
     &                        'Error while creating path '//TRIM(PATH)//
     &                        ' - exit status number is ', ESTAT
         CALL MSGERR( 2, TRIM(MSGSTR) )
         IERR = 3
      END IF
!
      RETURN
      END
