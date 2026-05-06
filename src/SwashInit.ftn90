subroutine SwashInit ( inerr )
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
!    1.00, January 2010: New subroutine
!
!   Purpose
!
!   Initializes several variables and arrays
!
!   Modules used
!
    use ocpcomm1
    use ocpcomm2
    use ocpcomm3
    use ocpcomm4
    use SwashCommdata1
    use SwashCommdata2
    use SwashCommdata3
    use SwashCommdata4
    use SwashTimecomm
    use SwashSolvedata
    use SwashRigBoddata
    use m_bndspec
    use outp_data
    use SwanGriddata
    use m_parall
!
    implicit none
!
!   Argument variables
!
    integer, intent(in) :: inerr   ! number of the initialisation error
!
!   Local variables
!
    integer             :: i       ! loop counter
    integer             :: ivtype  ! type of output quantity
    logical             :: STPNOW  ! indicates that program must stop
!
!   Structure
!
!   Description of the pseudo code
!
!   Source text
!
    VERTXT = BLANK
    VERNUM = 12.01
    write (VERTXT, '(f5.2)') VERNUM
    !call BUGFIX ('A')
    !call BUGFIX ('B')
    !call BUGFIX ('C')
    !
    call OCPINI ( 'swashinit', .true., inerr )
    if ( inerr > 0 ) return
    if (STPNOW()) return
    !
    write (PRINTF, 101) VERTXT
    if (SCREEN /= PRINTF .and. INODE == MASTER ) write (SCREEN,102)
    !
    ! initial values for common variables
    !
    PROJID = 'SWASH'
    PROJNR = BLANK
    PROJT1 = BLANK
    PROJT2 = BLANK
    PROJT3 = BLANK
    !
    ! units
    !
    uh     = 'm'
    uv     = 'm/s'
    ut     = 'sec'
    ul     = 'm'
    uq     = 'm2/s'
    udi    = 'degr'
    uc     = 'ppt'
    ud     = 'kg/m3'
    uk     = 'o C'
    uf     = 'N/m2'
    up     = 'hPa'
    ue     = 'm2/s'
    !
    ! some constants
    !
    pi     = 4.*atan(1.)
    pi2    = 2.*pi
    degrad = pi/180.
    !
    epsdry = 5.e-5
    epshu  = 1.e-3
    epsuf  = 1.e-8
    !
    swl    = 0.
    !
    dnorth = 90.
    ylatc  = 52.
    bnaut  = .false.
    !
    myseed = 12345678
    !
    ! physical parameters
    !
    grav   = 9.813
    rhow   = -999.        ! reference density (subject to change)
    rhoa   = 1.205
    dynvis = 1.e-3
    vonkar = 0.4
    !
    ! density parameters
    !
    idens  = 0
    tempw  = 14.
    salw   = 31.
    rhos   = 2650.
    lmixt  = .true.
    !
    ! viscosity and diffusivity parameters
    !
    ihvisc = 0
    bvisc  = -999.
    hvisc  = 0.
    hdiff  = 0.
    !
    ! spherical coordinates
    !
    mproj  = 0
    kspher = 0
    rearth = 2.e7/pi
    lendeg = 2.e7/180.
    !
    ! number of computational grid points
    !
    mcgrd   = 1
    mcgrdgl = 1
    nverts  = 0
    nvertsg = 0
    !
    ! computational grid
    !
    optg   = 1
    mxc    = 0
    myc    = 0
    mxcgl  = 0
    mycgl  = 0
    mxf    = 1
    mxl    = 0
    myf    = 1
    myl    = 0
    alpc   = 0.
    !
    oned   = .false.
    momskip= .false.
    linswe = .false.
    !
    kmax   = 1
    kpmax  = 1
    lsubg  = .false.
    !
    ! grid offsets
    !
    xoffs  = 0.
    yoffs  = 0.
    lxoffs = .false.
    !
    ! repeating grid / periodicity
    !
    lreptx = .false.
    lrepty = .false.
    !
    ! time of computation
    !
    timco  = -1.e10
    tinic  = 0.
    dt     = 0.
    chtime = '    '
    mtc    = 1
    nstatm = 1
    nstatc = -1
    ncompt = 0
    !
    ! initial conditions
    !
    tkeini  = 1.e-7
    epsini  = 1.e-7
    !
    restrt  = .false.
    instead = .false.
    !
    ! boundary conditions
    !
    nbfils    = 0
    nbval     = 0
    nbgrpt    = 0
    nbggl     = 0
    fbfs%nbfs = -999
    lrampf    = .false.
    alobnd    = .false.
    !
    spwidb = 0.
    spwidl = 0.
    spwidr = 0.
    spwidt = 0.
    !
    nspl   = 0
    !
    ! internal wave generation
    !
    iwvgen  = 0
    piwg    = 0.
    piwg(4) = 0.5
    !
    lsrcfx = .false.
    lsrcfy = .false.
    !
    ! input grids
    !
    do i = 1, numgrd
       !
       igtype(i) = 0
       xpg(i)    = 0.
       ypg(i)    = 0.
       alpg(i)   = 0.
       cospg(i)  = 1.
       sinpg(i)  = 0.
       dxg(i)    = 0.
       dyg(i)    = 0.
       mxg(i)    = 0
       myg(i)    = 0
       leds(i)   = 0
       lflgrd(i) = .true.
       lstag(i)  = .false.
       excfld(i) = -1.e20
       ifldyn(i) = 0
       iflndf(i) = 0
       ifltim(i) = -1.e20
       ifllay(i) = 1
       !
    enddo
    !
    initsf = .false.
    inituf = .false.
    initvf = .false.
    !
    varwi  = .false.
    varfr  = .false.
    vargs  = .false.
    varsh  = .false.
    varnpl = .false.
    varva  = .false.
    varea  = .false.
    !
    ! wind speed and direction
    !
    u10  = 0.
    wdip = 0.
    !
    ! wind stress coefficients
    !
    iwind   = 0
    pwnd    = 0.
    pwnd(1) = 0.002
    pwnd(4) = 1.
    pwnd(10)= 0.4
    relwnd  = .false.
    relwav  = .false.
    cdcap   = 99999.
    !
    ! space varying wind and pressure
    !
    svwp   = .false.
    prmean = -1.
    !
    ! parameters for Coriolis force
    !
    coriolis = .false.
    pcor(1) = -999.
    pcor(2) = 0.1
    !
    ! bottom friction coefficients
    !
    irough = 0
    pbot   = 0.
    !
    ! vegetation
    !
    iveg  = 0
    ivegw = 0
    lmax  = 0
    cvm   = 0.
    !
    ! control wave breaking
    !
    isurf    = 0
    psurf    = 0.
    psurf(1) = 0.6
    psurf(2) = -1.
    psurf(3) = 1.0
    !
    ! porosity
    !
    iporos  = 0
    ppor    = 0.
    ppor(2) = 99999.
    ppor(3) = 200.
    ppor(4) = 1.1
    !
    ! floating rigid bodies
    !
    ifloat    = 0
    mbod      = 0
    mfoti     = 2
    mxmli     = 0
    mxfen     = 0
    pship     = 0.
    pship( 1) = 0.
    pship( 2) = 1.
    pship( 3) = 0.25
    pship( 4) = 0.5
    pship( 5) = 1.
    pship( 8) = 0.0001
    pship( 9) = 50.
    pship(10) = 0.3
    pship(11) = 0.5
    !
    ! ambient current
    !
    iamb = 0
    uamb = 0.
    vamb = 0.
    eamb = 0.
    !
    ! transport of constituent
    !
    itrans = 0
    lsal   = 0
    ltemp  = 0
    lsed   = 0
    ltrans = 0
    icreep = 0
    tcret  = 0.
    !
    ! suspended sediment transport
    !
    psed    = 0.
    psed(3) = 0.7
    psed(4) = 5.5
    psed(5) = 0.05
    psed(6) = 0.00033
    !
    ! turbulence model
    !
    iturb    = 0
    ltur     = 0
    pturb    = 0.
    pturb(2) = 0.07
    pturb(3) = 0.16
    !
    ! discretization parameters
    !
    mtimei     = -1
    jhk        =  1     ! RK3
    dpsopt     =  2
    corrdep    = .true.
    depcds     = .false.
    stricthead = .false.
    strictmom  = .false.
    horwinc    = .false.
    verwinc    = .false.
    mimetic    = .false.
    ccn        = .true.
    !
    pnums     = 0.
    pnums( 1) = 0.5
    pnums( 2) = 0.4
    pnums( 3) = 0.8
    pnums( 4) = 0.5
    pnums( 6) = -1.
    pnums( 7) = 0.
    pnums(11) = 6.
    pnums(12) = 0.
    pnums(13) = 2.
    pnums(16) = 3.
    pnums(17) = -1.
    pnums(31) = 0.5
    pnums(32) = 0.5
    pnums(33) = 0.5
    pnums(36) = 6.
    pnums(37) = 0.
    pnums(38) = 2.
    pnums(41) = 6.
    pnums(42) = 0.
    pnums(43) = 2.
    pnums(46) = 5.
    pnums(47) = 0.
    pnums(51) = 5.
    pnums(52) = 0.
    pnums(61) = 4.
    pnums(64) = 1.
    pnums(66) = 4.
    pnums(69) = 1.
    pnums(71) = 1.
    !
    ! parameters for non-hydrostatic computation
    !
    ihydro   = 0
    iproj    = -999
    pnums(5) = -1.
    !
    ! parameters for pressure projection method
    !
    lpproj = .false.
    pnums(58) = 0.0001
    pnums(59) = 50.
    !
    qlay = 0
    qmax = 1
    !
    ! parameters for linear solvers
    !
    iamout    = 0
    icond     = -999
    pnums(21) = -1.
    pnums(22) = 0.01
    pnums(23) = 0.
    pnums(24) = 100.
    pnums(25) = 500.
    pnums(26) = -1.
    pnums(27) = -1.
    pnums(28) = 1.
    !
    inewt     = -1
    !
    ! test output control
    !
    ITEST  = 1
    intes  = 0
    ioutes = 0
    ltrace = .false.
    testfl = .false.
    nptst  = 0
    nptsta = 1
    maxmes = 200
    ifpar  = 0
    !
    upvdf  = 0
    !
    ! parameters for wave, mean current, mean constituent and mean turbulence output
    !
    lcuroutp = .false.
    ltraoutp = .false.
    lturoutp = .false.
    lwavoutp = .false.
    ncuroutp =  0.
    ntraoutp =  0.
    nturoutp =  0.
    nwavoutp =  0.
    tcuroutp = -1.
    ttraoutp = -1.
    tturoutp = -1.
    twavoutp = -1.
    !
    ! inundation depth for horizontal runup
    !
    hrunp = -999.
    !
    ! threshold depth for runup level
    !
    delrp = -999.
    !
    ! center of gravity of floating object for computing moments
    !
    cogx = 0.
    cogy = 0.
    cogz = 0.
    !
    ! rotation angle of floating object for computing forces and moments
    !
    alpobj = 0.
    !
    ! plot output
    !
    do i = 1, nmovar
      ovkeyw(i) = 'XXXX'
    enddo
    !
    ! properties of output variables
    !
    ivtype = 1
    !
    ! keyword used in command file
    !
    ovkeyw(ivtype) = 'XP'
    !
    ! short name
    !
    ovsnam(ivtype) = 'Xp'
    !
    ! long name
    !
    ovlnam(ivtype) = 'X user coordinate'
    !
    ! unit name
    !
    ovunit(ivtype) = ul
    !
    ! type (scalar/vector etc.)
    !
    ovsvty(ivtype) = 1
    !
    ! lower and upper limit
    !
    ovllim(ivtype) = -1.e10
    ovulim(ivtype) = 1.e10
    !
    ! lowest and highest expected value
    !
    ovlexp(ivtype) = -1.e10
    ovhexp(ivtype) = 1.e10
    !
    ! exception value
    !
    ovexcv(ivtype) = -1.e10
    !
    ivtype = 2
    ovkeyw(ivtype) = 'YP'
    ovsnam(ivtype) = 'Yp'
    ovlnam(ivtype) = 'Y user coordinate'
    ovunit(ivtype) = ul
    ovsvty(ivtype) = 1
    ovllim(ivtype) = -1.e10
    ovulim(ivtype) = 1.e10
    ovlexp(ivtype) = -1.e10
    ovhexp(ivtype) = 1.e10
    ovexcv(ivtype) = -1.e10
    !
    ivtype = 3
    ovkeyw(ivtype) = 'DIST'
    ovsnam(ivtype) = 'Dist'
    ovlnam(ivtype) = 'distance along output curve'
    ovunit(ivtype) = ul
    ovsvty(ivtype) = 1
    ovllim(ivtype) = 0.
    ovulim(ivtype) = 1.e10
    ovlexp(ivtype) = 0.
    ovhexp(ivtype) = 1.e10
    ovexcv(ivtype) = -99.
    !
    ivtype = 4
    ovkeyw(ivtype) = 'DEP'
    ovsnam(ivtype) = 'Depth'
    ovlnam(ivtype) = 'Depth'
    ovunit(ivtype) = uh
    ovsvty(ivtype) = 1
    ovllim(ivtype) = -1.e4
    ovulim(ivtype) = 1.e4
    ovlexp(ivtype) = -100.
    ovhexp(ivtype) = 100.
    ovexcv(ivtype) = -99.
    !
    ivtype = 5
    ovkeyw(ivtype) = 'BOTL'
    ovsnam(ivtype) = 'Botlev'
    ovlnam(ivtype) = 'Bottom level'
    ovunit(ivtype) = uh
    ovsvty(ivtype) = 1
    ovllim(ivtype) = -1.e4
    ovulim(ivtype) = 1.e4
    ovlexp(ivtype) = -100.
    ovhexp(ivtype) = 100.
    ovexcv(ivtype) = -99.
    !
    ivtype = 6
    ovkeyw(ivtype) = 'WATL'
    ovsnam(ivtype) = 'Watlev'
    ovlnam(ivtype) = 'Water level'
    ovunit(ivtype) = uh
    ovsvty(ivtype) = 1
    ovllim(ivtype) = -1.e4
    ovulim(ivtype) = 1.e4
    ovlexp(ivtype) = -100.
    ovhexp(ivtype) = 100.
    ovexcv(ivtype) = -99.
    !
    ivtype = 7
    ovkeyw(ivtype) = 'VMAG'
    ovsnam(ivtype) = 'Vmag'
    ovlnam(ivtype) = 'Velocity magnitude'
    ovunit(ivtype) = uv
    ovsvty(ivtype) = 1
    ovllim(ivtype) = -100.
    ovulim(ivtype) = 100.
    ovlexp(ivtype) = -2.
    ovhexp(ivtype) = 2.
    ovexcv(ivtype) = 0.
    !
    ivtype = 8
    ovkeyw(ivtype) = 'VDIR'
    ovsnam(ivtype) = 'Vdir'
    ovlnam(ivtype) = 'Velocity direction'
    ovunit(ivtype) = udi
    ovsvty(ivtype) = 2
    ovllim(ivtype) = 0.
    ovulim(ivtype) = 360.
    ovlexp(ivtype) = 0.
    ovhexp(ivtype) = 360.
    ovexcv(ivtype) = -999.
    !
    ivtype = 9
    ovkeyw(ivtype) = 'VEL'
    ovsnam(ivtype) = 'vel'
    ovlnam(ivtype) = 'Flow velocity'
    ovunit(ivtype) = uv
    ovsvty(ivtype) = 3
    ovllim(ivtype) = -100.
    ovulim(ivtype) = 100.
    ovlexp(ivtype) = -2.
    ovhexp(ivtype) = 2.
    ovexcv(ivtype) = 0.
    !
    ivtype = 10
    ovkeyw(ivtype) = 'VKSI'
    ovsnam(ivtype) = 'Vksi'
    ovlnam(ivtype) = 'Grid-oriented U-velocity'
    ovunit(ivtype) = uv
    ovsvty(ivtype) = 1
    ovllim(ivtype) = -100.
    ovulim(ivtype) = 100.
    ovlexp(ivtype) = -2.
    ovhexp(ivtype) = 2.
    ovexcv(ivtype) = 0.
    !
    ivtype = 11
    ovkeyw(ivtype) = 'VETA'
    ovsnam(ivtype) = 'Veta'
    ovlnam(ivtype) = 'Grid-oriented V-velocity'
    ovunit(ivtype) = uv
    ovsvty(ivtype) = 1
    ovllim(ivtype) = -100.
    ovulim(ivtype) = 100.
    ovlexp(ivtype) = -2.
    ovhexp(ivtype) = 2.
    ovexcv(ivtype) = 0.
    !
    ivtype = 12
    ovkeyw(ivtype) = 'NHPRES'
    ovsnam(ivtype) = 'Nhpres'
    ovlnam(ivtype) = 'Non-hydrostatic pressure at bottom'
    ovunit(ivtype) = uf
    ovsvty(ivtype) = 1
    ovllim(ivtype) = -1.e-4
    ovulim(ivtype) = 1.e4
    ovlexp(ivtype) = -100.
    ovhexp(ivtype) = 100.
    ovexcv(ivtype) = -99.
    !
    ivtype = 13
    ovkeyw(ivtype) = 'PRESS'
    ovsnam(ivtype) = 'Press'
    ovlnam(ivtype) = 'Pressure at bottom'
    ovunit(ivtype) = up
    ovsvty(ivtype) = 1
    ovllim(ivtype) = -1.e-4
    ovulim(ivtype) = 1.e4
    ovlexp(ivtype) = -2.
    ovhexp(ivtype) = 2.
    ovexcv(ivtype) = -99.
    !
    ivtype = 14
    ovkeyw(ivtype) = 'QMAG'
    ovsnam(ivtype) = 'Qmag'
    ovlnam(ivtype) = 'Magnitude of discharge per unit width'
    ovunit(ivtype) = uq
    ovsvty(ivtype) = 1
    ovllim(ivtype) = -100.
    ovulim(ivtype) = 100.
    ovlexp(ivtype) = -2.
    ovhexp(ivtype) = 2.
    ovexcv(ivtype) = 0.
    !
    ivtype = 15
    ovkeyw(ivtype) = 'QDIR'
    ovsnam(ivtype) = 'Qdir'
    ovlnam(ivtype) = 'Direction of discharge per unit width'
    ovunit(ivtype) = udi
    ovsvty(ivtype) = 2
    ovllim(ivtype) = 0.
    ovulim(ivtype) = 360.
    ovlexp(ivtype) = 0.
    ovhexp(ivtype) = 360.
    ovexcv(ivtype) = -999.
    !
    ivtype = 16
    ovkeyw(ivtype) = 'DISCH'
    ovsnam(ivtype) = 'disch'
    ovlnam(ivtype) = 'Discharge per unit width'
    ovunit(ivtype) = uq
    ovsvty(ivtype) = 3
    ovllim(ivtype) = -100.
    ovulim(ivtype) = 100.
    ovlexp(ivtype) = -2.
    ovhexp(ivtype) = 2.
    ovexcv(ivtype) = 0.
    !
    ivtype = 17
    ovkeyw(ivtype) = 'QKSI'
    ovsnam(ivtype) = 'Qksi'
    ovlnam(ivtype) = 'Grid-oriented U-discharge per unit width'
    ovunit(ivtype) = uq
    ovsvty(ivtype) = 1
    ovllim(ivtype) = -100.
    ovulim(ivtype) = 100.
    ovlexp(ivtype) = -2.
    ovhexp(ivtype) = 2.
    ovexcv(ivtype) = 0.
    !
    ivtype = 18
    ovkeyw(ivtype) = 'QETA'
    ovsnam(ivtype) = 'Qeta'
    ovlnam(ivtype) = 'Grid-oriented V-discharge per unit width'
    ovunit(ivtype) = uq
    ovsvty(ivtype) = 1
    ovllim(ivtype) = -100.
    ovulim(ivtype) = 100.
    ovlexp(ivtype) = -2.
    ovhexp(ivtype) = 2.
    ovexcv(ivtype) = 0.
    !
    ivtype = 19
    ovkeyw(ivtype) = 'SAL'
    ovsnam(ivtype) = 'Sal'
    ovlnam(ivtype) = 'Salinity'
    ovunit(ivtype) = uc
    ovsvty(ivtype) = 1
    ovllim(ivtype) = 0.
    ovulim(ivtype) = 50.
    ovlexp(ivtype) = 0.
    ovhexp(ivtype) = 50.
    ovexcv(ivtype) = -99.
    !
    ivtype = 20
    ovkeyw(ivtype) = 'TEMP'
    ovsnam(ivtype) = 'Temp'
    ovlnam(ivtype) = 'Temperature'
    ovunit(ivtype) = uk
    ovsvty(ivtype) = 1
    ovllim(ivtype) = -100.
    ovulim(ivtype) = 100.
    ovlexp(ivtype) = -100.
    ovhexp(ivtype) = 100.
    ovexcv(ivtype) = -99.
    !
    ivtype = 21
    ovkeyw(ivtype) = 'HRUN'
    ovsnam(ivtype) = 'Hrunup'
    ovlnam(ivtype) = 'maximum horizontal runup'
    ovunit(ivtype) = ' '
    ovsvty(ivtype) = 1
    ovllim(ivtype) = 0.
    ovulim(ivtype) = 1.
    ovlexp(ivtype) = 0.
    ovhexp(ivtype) = 1.e10
    ovexcv(ivtype) = -99.
    !
    ivtype = 22
    ovkeyw(ivtype) = 'SETUP'
    ovsnam(ivtype) = 'Setup'
    ovlnam(ivtype) = 'Setup due to waves'
    ovunit(ivtype) = uh
    ovsvty(ivtype) = 1
    ovllim(ivtype) = -1.
    ovulim(ivtype) = 1.
    ovlexp(ivtype) = -1.
    ovhexp(ivtype) = 1.
    ovexcv(ivtype) = -9.
    !
    ivtype = 23
    ovkeyw(ivtype) = 'HS'
    ovsnam(ivtype) = 'Hsig'
    ovlnam(ivtype) = 'Significant wave height'
    ovunit(ivtype) = uh
    ovsvty(ivtype) = 1
    ovllim(ivtype) = 0.
    ovulim(ivtype) = 100.
    ovlexp(ivtype) = 0.
    ovhexp(ivtype) = 10.
    ovexcv(ivtype) = -9.
    !
    ivtype = 24
    ovkeyw(ivtype) = 'HRMS'
    ovsnam(ivtype) = 'Hrms'
    ovlnam(ivtype) = 'RMS wave height'
    ovunit(ivtype) = uh
    ovsvty(ivtype) = 1
    ovllim(ivtype) = 0.
    ovulim(ivtype) = 100.
    ovlexp(ivtype) = 0.
    ovhexp(ivtype) = 10.
    ovexcv(ivtype) = -9.
    !
    ivtype = 25
    ovkeyw(ivtype) = 'XC'
    ovsnam(ivtype) = 'Xc'
    ovlnam(ivtype) = 'X computational grid coordinate'
    ovunit(ivtype) = ' '
    ovsvty(ivtype) = 1
    ovllim(ivtype) = 0.
    ovulim(ivtype) = 1000.
    ovlexp(ivtype) = 0.
    ovhexp(ivtype) = 100.
    ovexcv(ivtype) = -9.
    !
    ivtype = 26
    ovkeyw(ivtype) = 'YC'
    ovsnam(ivtype) = 'Yc'
    ovlnam(ivtype) = 'Y computational grid coordinate'
    ovunit(ivtype) = ' '
    ovsvty(ivtype) = 1
    ovllim(ivtype) = 0.
    ovulim(ivtype) = 1000.
    ovlexp(ivtype) = 0.
    ovhexp(ivtype) = 100.
    ovexcv(ivtype) = -9.
    !
    ivtype = 27
    ovkeyw(ivtype) = 'WMAG'
    ovsnam(ivtype) = 'Wndvel'
    ovlnam(ivtype) = 'Wind velocity at 10 m above sea level'
    ovunit(ivtype) = uv
    ovsvty(ivtype) = 1
    ovllim(ivtype) = -100.
    ovulim(ivtype) = 100.
    ovlexp(ivtype) = -50.
    ovhexp(ivtype) = 50.
    ovexcv(ivtype) = 0.
    !
    ivtype = 28
    ovkeyw(ivtype) = 'WDIR'
    ovsnam(ivtype) = 'Windir'
    ovlnam(ivtype) = 'Wind direction at 10 m above sea level'
    ovunit(ivtype) = udi
    ovsvty(ivtype) = 2
    ovllim(ivtype) = 0.
    ovulim(ivtype) = 360.
    ovlexp(ivtype) = 0.
    ovhexp(ivtype) = 360.
    ovexcv(ivtype) = -999.
    !
    ivtype = 29
    ovkeyw(ivtype) = 'WIND'
    ovsnam(ivtype) = 'wind'
    ovlnam(ivtype) = 'Wind velocity at 10 m above sea level'
    ovunit(ivtype) = uv
    ovsvty(ivtype) = 3
    ovllim(ivtype) = -100.
    ovulim(ivtype) = 100.
    ovlexp(ivtype) = -50.
    ovhexp(ivtype) = 50.
    ovexcv(ivtype) = 0.
    !
    ivtype = 30
    ovkeyw(ivtype) = 'FRC'
    ovsnam(ivtype) = 'FrCoef'
    ovlnam(ivtype) = 'Bottom friction coefficient'
    ovunit(ivtype) = ' '
    ovsvty(ivtype) = 1
    ovllim(ivtype) = 0.
    ovulim(ivtype) = 1.
    ovlexp(ivtype) = 0.
    ovhexp(ivtype) = 1.
    ovexcv(ivtype) = -9.
    !
    ivtype = 31
    ovkeyw(ivtype) = 'USTAR'
    ovsnam(ivtype) = 'Ustar'
    ovlnam(ivtype) = 'Magnitude of friction velocity'
    ovunit(ivtype) = uv
    ovsvty(ivtype) = 1
    ovllim(ivtype) = -100.
    ovulim(ivtype) = 100.
    ovlexp(ivtype) = -2.
    ovhexp(ivtype) = 2.
    ovexcv(ivtype) = 0.
    !
    ivtype = 32
    ovkeyw(ivtype) = 'UFRIC'
    ovsnam(ivtype) = 'Ufric'
    ovlnam(ivtype) = 'Friction velocity'
    ovunit(ivtype) = uv
    ovsvty(ivtype) = 3
    ovllim(ivtype) = -100.
    ovulim(ivtype) = 100.
    ovlexp(ivtype) = -2.
    ovhexp(ivtype) = 2.
    ovexcv(ivtype) = 0.
    !
    ivtype = 33
    ovkeyw(ivtype) = 'MVMAG'
    ovsnam(ivtype) = 'MVmag'
    ovlnam(ivtype) = 'Mean velocity magnitude'
    ovunit(ivtype) = uv
    ovsvty(ivtype) = 1
    ovllim(ivtype) = -100.
    ovulim(ivtype) = 100.
    ovlexp(ivtype) = -2.
    ovhexp(ivtype) = 2.
    ovexcv(ivtype) = 0.
    !
    ivtype = 34
    ovkeyw(ivtype) = 'MVDIR'
    ovsnam(ivtype) = 'MVdir'
    ovlnam(ivtype) = 'Mean velocity direction'
    ovunit(ivtype) = udi
    ovsvty(ivtype) = 2
    ovllim(ivtype) = 0.
    ovulim(ivtype) = 360.
    ovlexp(ivtype) = 0.
    ovhexp(ivtype) = 360.
    ovexcv(ivtype) = -999.
    !
    ivtype = 35
    ovkeyw(ivtype) = 'MVEL'
    ovsnam(ivtype) = 'Mvel'
    ovlnam(ivtype) = 'Mean velocity'
    ovunit(ivtype) = uv
    ovsvty(ivtype) = 3
    ovllim(ivtype) = -100.
    ovulim(ivtype) = 100.
    ovlexp(ivtype) = -2.
    ovhexp(ivtype) = 2.
    ovexcv(ivtype) = 0.
    !
    ivtype = 36
    ovkeyw(ivtype) = 'MVKSI'
    ovsnam(ivtype) = 'MVksi'
    ovlnam(ivtype) = 'Mean grid-oriented U-velocity'
    ovunit(ivtype) = uv
    ovsvty(ivtype) = 1
    ovllim(ivtype) = -100.
    ovulim(ivtype) = 100.
    ovlexp(ivtype) = -2.
    ovhexp(ivtype) = 2.
    ovexcv(ivtype) = 0.
    !
    ivtype = 37
    ovkeyw(ivtype) = 'MVETA'
    ovsnam(ivtype) = 'MVeta'
    ovlnam(ivtype) = 'Mean grid-oriented V-velocity'
    ovunit(ivtype) = uv
    ovsvty(ivtype) = 1
    ovllim(ivtype) = -100.
    ovulim(ivtype) = 100.
    ovlexp(ivtype) = -2.
    ovhexp(ivtype) = 2.
    ovexcv(ivtype) = 0.
    !
    ivtype = 38
    ovkeyw(ivtype) = 'BRKP'
    ovsnam(ivtype) = 'Brkpnt'
    ovlnam(ivtype) = 'wave breaking point'
    ovunit(ivtype) = ' '
    ovsvty(ivtype) = 1
    ovllim(ivtype) = 0.
    ovulim(ivtype) = 1.
    ovlexp(ivtype) = 0.
    ovhexp(ivtype) = 1.e10
    ovexcv(ivtype) = -99.
    !
    ivtype = 39
    ovkeyw(ivtype) = 'SED'
    ovsnam(ivtype) = 'Sed'
    ovlnam(ivtype) = 'Suspended sediment'
    ovunit(ivtype) = ud
    ovsvty(ivtype) = 1
    ovllim(ivtype) = 0.
    ovulim(ivtype) = 1.e4
    ovlexp(ivtype) = 0.
    ovhexp(ivtype) = 1.e4
    ovexcv(ivtype) = -99.
    !
    ivtype = 40
    ovkeyw(ivtype) = 'TIME'
    ovsnam(ivtype) = 'Time'
    ovlnam(ivtype) = 'Date-time'
    ovunit(ivtype) = ' '
    ovsvty(ivtype) = 1
    ovllim(ivtype) = 0.
    ovulim(ivtype) = 1.
    ovlexp(ivtype) = 0.
    ovhexp(ivtype) = 1.
    ovexcv(ivtype) = -99999.
    !
    ivtype = 41
    ovkeyw(ivtype) = 'TSEC'
    ovsnam(ivtype) = 'Tsec'
    ovlnam(ivtype) = 'Time in seconds from reference time'
    ovunit(ivtype) = 's'
    ovsvty(ivtype) = 1
    ovllim(ivtype) = 0.
    ovulim(ivtype) = 100000.
    ovlexp(ivtype) = -100000.
    ovhexp(ivtype) = 10.
    ovexcv(ivtype) = -99999.
    !
    ivtype = 42
    ovkeyw(ivtype) = 'MSAL'
    ovsnam(ivtype) = 'MSal'
    ovlnam(ivtype) = 'Mean salinity'
    ovunit(ivtype) = uc
    ovsvty(ivtype) = 1
    ovllim(ivtype) = 0.
    ovulim(ivtype) = 50.
    ovlexp(ivtype) = 0.
    ovhexp(ivtype) = 50.
    ovexcv(ivtype) = -99.
    !
    ivtype = 43
    ovkeyw(ivtype) = 'MTEMP'
    ovsnam(ivtype) = 'MTemp'
    ovlnam(ivtype) = 'Mean temperature'
    ovunit(ivtype) = uk
    ovsvty(ivtype) = 1
    ovllim(ivtype) = -100.
    ovulim(ivtype) = 100.
    ovlexp(ivtype) = -100.
    ovhexp(ivtype) = 100.
    ovexcv(ivtype) = -99.
    !
    ivtype = 44
    ovkeyw(ivtype) = 'MSED'
    ovsnam(ivtype) = 'MSed'
    ovlnam(ivtype) = 'Mean suspended sediment'
    ovunit(ivtype) = ud
    ovsvty(ivtype) = 1
    ovllim(ivtype) = 0.
    ovulim(ivtype) = 1.e4
    ovlexp(ivtype) = 0.
    ovhexp(ivtype) = 1.e4
    ovexcv(ivtype) = -99.
    !
    ivtype = 45
    ovkeyw(ivtype) = 'VORT'
    ovsnam(ivtype) = 'vort'
    ovlnam(ivtype) = 'Vorticity'
    ovunit(ivtype) = '1/s'
    ovsvty(ivtype) = 1
    ovllim(ivtype) = -100.
    ovulim(ivtype) = 100.
    ovlexp(ivtype) = -2.
    ovhexp(ivtype) = 2.
    ovexcv(ivtype) = -99.
    !
    ivtype = 46
    ovkeyw(ivtype) = 'DRAF'
    ovsnam(ivtype) = 'Draft'
    ovlnam(ivtype) = 'Draft of floating object'
    ovunit(ivtype) = uh
    ovsvty(ivtype) = 1
    ovllim(ivtype) = -1.e4
    ovulim(ivtype) = 1.e4
    ovlexp(ivtype) = -100.
    ovhexp(ivtype) = 100.
    ovexcv(ivtype) = -99.
    !
    ivtype = 50
    ovkeyw(ivtype) = 'AAAA'
    ovsnam(ivtype) = 'Aux'
    ovlnam(ivtype) = 'auxiliary variable'
    ovunit(ivtype) = ' '
    ovsvty(ivtype) = 1
    ovllim(ivtype) = -1.e10
    ovulim(ivtype) = 1.e10
    ovlexp(ivtype) = -1.e10
    ovhexp(ivtype) = 1.e10
    ovexcv(ivtype) = -1.e10
    !
    ! layer-dependent quantities (k = 0:kmax)
    !
    ivtype = 51
    ovkeyw(ivtype) = 'ZK'
    ovsnam(ivtype) = 'zk'
    ovlnam(ivtype) = 'Layer interface'
    ovunit(ivtype) = uh
    ovsvty(ivtype) = 1
    ovllim(ivtype) = -1.e4
    ovulim(ivtype) = 1.e4
    ovlexp(ivtype) = -100.
    ovhexp(ivtype) = 100.
    ovexcv(ivtype) = -99.
    !
    ivtype = 52
    ovkeyw(ivtype) = 'VZ'
    ovsnam(ivtype) = 'w'
    ovlnam(ivtype) = 'Velocity in z-direction'
    ovunit(ivtype) = uv
    ovsvty(ivtype) = 1
    ovllim(ivtype) = -100.
    ovulim(ivtype) = 100.
    ovlexp(ivtype) = -2.
    ovhexp(ivtype) = 2.
    ovexcv(ivtype) = -99.
    !
    ivtype = 53
    ovkeyw(ivtype) = 'VOMEGA'
    ovsnam(ivtype) = 'omega'
    ovlnam(ivtype) = 'Vertical velocity relative to sigma plane'
    ovunit(ivtype) = uv
    ovsvty(ivtype) = 1
    ovllim(ivtype) = -100.
    ovulim(ivtype) = 100.
    ovlexp(ivtype) = -2.
    ovhexp(ivtype) = 2.
    ovexcv(ivtype) = -99.
    !
    ivtype = 54
    ovkeyw(ivtype) = 'TKE'
    ovsnam(ivtype) = 'Tke'
    ovlnam(ivtype) = 'Turbulent kinetic energy'
    ovunit(ivtype) = 'm2/s2'
    ovsvty(ivtype) = 1
    ovllim(ivtype) = 0.
    ovulim(ivtype) = 1000.
    ovlexp(ivtype) = 0.
    ovhexp(ivtype) = 1000.
    ovexcv(ivtype) = -99.
    !
    ivtype = 55
    ovkeyw(ivtype) = 'EPS'
    ovsnam(ivtype) = 'Eps'
    ovlnam(ivtype) = 'Dissipation rate of turbulent kinetic energy'
    ovunit(ivtype) = 'm2/s3'
    ovsvty(ivtype) = 1
    ovllim(ivtype) = 0.
    ovulim(ivtype) = 1000.
    ovlexp(ivtype) = 0.
    ovhexp(ivtype) = 1000.
    ovexcv(ivtype) = -99.
    !
    ivtype = 56
    ovkeyw(ivtype) = 'VISC'
    ovsnam(ivtype) = 'Visc'
    ovlnam(ivtype) = 'Vertical eddy viscosity'
    ovunit(ivtype) = ue
    ovsvty(ivtype) = 1
    ovllim(ivtype) = 0.
    ovulim(ivtype) = 1000.
    ovlexp(ivtype) = 0.
    ovhexp(ivtype) = 1000.
    ovexcv(ivtype) = -99.
    !
    ivtype = 57
    ovkeyw(ivtype) = 'MTKE'
    ovsnam(ivtype) = 'MTke'
    ovlnam(ivtype) = 'Mean turbulent kinetic energy'
    ovunit(ivtype) = 'm2/s2'
    ovsvty(ivtype) = 1
    ovllim(ivtype) = 0.
    ovulim(ivtype) = 1000.
    ovlexp(ivtype) = 0.
    ovhexp(ivtype) = 1000.
    ovexcv(ivtype) = -99.
    !
    ivtype = 58
    ovkeyw(ivtype) = 'MEPS'
    ovsnam(ivtype) = 'MEps'
    ovlnam(ivtype) = 'Mean dissipation rate of turbulent kinetic energy'
    ovunit(ivtype) = 'm2/s3'
    ovsvty(ivtype) = 1
    ovllim(ivtype) = 0.
    ovulim(ivtype) = 1000.
    ovlexp(ivtype) = 0.
    ovhexp(ivtype) = 1000.
    ovexcv(ivtype) = -99.
    !
    ivtype = 59
    ovkeyw(ivtype) = 'MVISC'
    ovsnam(ivtype) = 'MVisc'
    ovlnam(ivtype) = 'Mean vertical eddy viscosity'
    ovunit(ivtype) = ue
    ovsvty(ivtype) = 1
    ovllim(ivtype) = 0.
    ovulim(ivtype) = 1000.
    ovlexp(ivtype) = 0.
    ovhexp(ivtype) = 1000.
    ovexcv(ivtype) = -99.
    !
    ! layer-dependent quantities (k = 1:kmax)
    !
    ivtype = 71
    ovkeyw(ivtype) = 'HK'
    ovsnam(ivtype) = 'hk'
    ovlnam(ivtype) = 'Layer thickness'
    ovunit(ivtype) = uh
    ovsvty(ivtype) = 1
    ovllim(ivtype) = -1.e4
    ovulim(ivtype) = 1.e4
    ovlexp(ivtype) = -100.
    ovhexp(ivtype) = 100.
    ovexcv(ivtype) = -99.
    !
    ivtype = 72
    ovkeyw(ivtype) = 'VMAGK'
    ovsnam(ivtype) = 'Vmag_k'
    ovlnam(ivtype) = 'Velocity magnitude per layer'
    ovunit(ivtype) = uv
    ovsvty(ivtype) = 1
    ovllim(ivtype) = -100.
    ovulim(ivtype) = 100.
    ovlexp(ivtype) = -2.
    ovhexp(ivtype) = 2.
    ovexcv(ivtype) = 0.
    !
    ivtype = 73
    ovkeyw(ivtype) = 'VDIRK'
    ovsnam(ivtype) = 'Vdir_k'
    ovlnam(ivtype) = 'Velocity direction per layer'
    ovunit(ivtype) = udi
    ovsvty(ivtype) = 2
    ovllim(ivtype) = 0.
    ovulim(ivtype) = 360.
    ovlexp(ivtype) = 0.
    ovhexp(ivtype) = 360.
    ovexcv(ivtype) = -999.
    !
    ivtype = 74
    ovkeyw(ivtype) = 'VELK'
    ovsnam(ivtype) = 'vel_k'
    ovlnam(ivtype) = 'Flow velocity per layer'
    ovunit(ivtype) = uv
    ovsvty(ivtype) = 3
    ovllim(ivtype) = -100.
    ovulim(ivtype) = 100.
    ovlexp(ivtype) = -2.
    ovhexp(ivtype) = 2.
    ovexcv(ivtype) = 0.
    !
    ivtype = 75
    ovkeyw(ivtype) = 'VKSIK'
    ovsnam(ivtype) = 'Vksi_k'
    ovlnam(ivtype) = 'Grid-oriented U-velocity per layer'
    ovunit(ivtype) = uv
    ovsvty(ivtype) = 1
    ovllim(ivtype) = -100.
    ovulim(ivtype) = 100.
    ovlexp(ivtype) = -2.
    ovhexp(ivtype) = 2.
    ovexcv(ivtype) = 0.
    !
    ivtype = 76
    ovkeyw(ivtype) = 'VETAK'
    ovsnam(ivtype) = 'Veta_k'
    ovlnam(ivtype) = 'Grid-oriented V-velocity per layer'
    ovunit(ivtype) = uv
    ovsvty(ivtype) = 1
    ovllim(ivtype) = -100.
    ovulim(ivtype) = 100.
    ovlexp(ivtype) = -2.
    ovhexp(ivtype) = 2.
    ovexcv(ivtype) = 0.
    !
    ivtype = 77
    ovkeyw(ivtype) = 'QMAGK'
    ovsnam(ivtype) = 'Qmag_k'
    ovlnam(ivtype) = 'Magnitude of discharge per unit width per layer'
    ovunit(ivtype) = uq
    ovsvty(ivtype) = 1
    ovllim(ivtype) = -100.
    ovulim(ivtype) = 100.
    ovlexp(ivtype) = -2.
    ovhexp(ivtype) = 2.
    ovexcv(ivtype) = 0.
    !
    ivtype = 78
    ovkeyw(ivtype) = 'QDIRK'
    ovsnam(ivtype) = 'Qdir_k'
    ovlnam(ivtype) = 'Direction of discharge per unit width per layer'
    ovunit(ivtype) = udi
    ovsvty(ivtype) = 2
    ovllim(ivtype) = 0.
    ovulim(ivtype) = 360.
    ovlexp(ivtype) = 0.
    ovhexp(ivtype) = 360.
    ovexcv(ivtype) = -999.
    !
    ivtype = 79
    ovkeyw(ivtype) = 'DISCHK'
    ovsnam(ivtype) = 'disc_k'
    ovlnam(ivtype) = 'Discharge per unit width per layer'
    ovunit(ivtype) = uq
    ovsvty(ivtype) = 3
    ovllim(ivtype) = -100.
    ovulim(ivtype) = 100.
    ovlexp(ivtype) = -2.
    ovhexp(ivtype) = 2.
    ovexcv(ivtype) = 0.
    !
    ivtype = 80
    ovkeyw(ivtype) = 'QKSIK'
    ovsnam(ivtype) = 'Qksi_k'
    ovlnam(ivtype) = 'Grid-oriented U-discharge per unit width per layer'
    ovunit(ivtype) = uq
    ovsvty(ivtype) = 1
    ovllim(ivtype) = -100.
    ovulim(ivtype) = 100.
    ovlexp(ivtype) = -2.
    ovhexp(ivtype) = 2.
    ovexcv(ivtype) = 0.
    !
    ivtype = 81
    ovkeyw(ivtype) = 'QETAK'
    ovsnam(ivtype) = 'Qeta_k'
    ovlnam(ivtype) = 'Grid-oriented V-discharge per unit width per layer'
    ovunit(ivtype) = uq
    ovsvty(ivtype) = 1
    ovllim(ivtype) = -100.
    ovulim(ivtype) = 100.
    ovlexp(ivtype) = -2.
    ovhexp(ivtype) = 2.
    ovexcv(ivtype) = 0.
    !
    ivtype = 82
    ovkeyw(ivtype) = 'NHPRSK'
    ovsnam(ivtype) = 'Nprs_k'
    ovlnam(ivtype) = 'Non-hydrostatic pressure per layer'
    ovunit(ivtype) = uf
    ovsvty(ivtype) = 1
    ovllim(ivtype) = -1.e-4
    ovulim(ivtype) = 1.e4
    ovlexp(ivtype) = -100.
    ovhexp(ivtype) = 100.
    ovexcv(ivtype) = -99.
    !
    ivtype = 83
    ovkeyw(ivtype) = 'PRESSK'
    ovsnam(ivtype) = 'Pres_k'
    ovlnam(ivtype) = 'Pressure per layer'
    ovunit(ivtype) = up
    ovsvty(ivtype) = 1
    ovllim(ivtype) = -1.e-4
    ovulim(ivtype) = 1.e4
    ovlexp(ivtype) = -2.
    ovhexp(ivtype) = 2.
    ovexcv(ivtype) = -99.
    !
    ivtype = 84
    ovkeyw(ivtype) = 'MVMAGK'
    ovsnam(ivtype) = 'MVmagk'
    ovlnam(ivtype) = 'Mean velocity magnitude per layer'
    ovunit(ivtype) = uv
    ovsvty(ivtype) = 1
    ovllim(ivtype) = -100.
    ovulim(ivtype) = 100.
    ovlexp(ivtype) = -2.
    ovhexp(ivtype) = 2.
    ovexcv(ivtype) = 0.
    !
    ivtype = 85
    ovkeyw(ivtype) = 'MVDIRK'
    ovsnam(ivtype) = 'MVdirk'
    ovlnam(ivtype) = 'Mean velocity direction per layer'
    ovunit(ivtype) = udi
    ovsvty(ivtype) = 2
    ovllim(ivtype) = 0.
    ovulim(ivtype) = 360.
    ovlexp(ivtype) = 0.
    ovhexp(ivtype) = 360.
    ovexcv(ivtype) = -999.
    !
    ivtype = 86
    ovkeyw(ivtype) = 'MVELK'
    ovsnam(ivtype) = 'Mvel_k'
    ovlnam(ivtype) = 'Mean velocity per layer'
    ovunit(ivtype) = uv
    ovsvty(ivtype) = 3
    ovllim(ivtype) = -100.
    ovulim(ivtype) = 100.
    ovlexp(ivtype) = -2.
    ovhexp(ivtype) = 2.
    ovexcv(ivtype) = 0.
    !
    ivtype = 87
    ovkeyw(ivtype) = 'MVKSIK'
    ovsnam(ivtype) = 'MVksik'
    ovlnam(ivtype) = 'Mean grid-oriented U-velocity per layer'
    ovunit(ivtype) = uv
    ovsvty(ivtype) = 1
    ovllim(ivtype) = -100.
    ovulim(ivtype) = 100.
    ovlexp(ivtype) = -2.
    ovhexp(ivtype) = 2.
    ovexcv(ivtype) = 0.
    !
    ivtype = 88
    ovkeyw(ivtype) = 'MVETAK'
    ovsnam(ivtype) = 'MVetak'
    ovlnam(ivtype) = 'Mean grid-oriented V-velocity per layer'
    ovunit(ivtype) = uv
    ovsvty(ivtype) = 1
    ovllim(ivtype) = -100.
    ovulim(ivtype) = 100.
    ovlexp(ivtype) = -2.
    ovhexp(ivtype) = 2.
    ovexcv(ivtype) = 0.
    !
    ivtype = 89
    ovkeyw(ivtype) = 'SALK'
    ovsnam(ivtype) = 'Salk'
    ovlnam(ivtype) = 'Salinity per layer'
    ovunit(ivtype) = uc
    ovsvty(ivtype) = 1
    ovllim(ivtype) = 0.
    ovulim(ivtype) = 50.
    ovlexp(ivtype) = 0.
    ovhexp(ivtype) = 50.
    ovexcv(ivtype) = -99.
    !
    ivtype = 90
    ovkeyw(ivtype) = 'TEMPK'
    ovsnam(ivtype) = 'Tempk'
    ovlnam(ivtype) = 'Temperature per layer'
    ovunit(ivtype) = uk
    ovsvty(ivtype) = 1
    ovllim(ivtype) = -100.
    ovulim(ivtype) = 100.
    ovlexp(ivtype) = -100.
    ovhexp(ivtype) = 100.
    ovexcv(ivtype) = -99.
    !
    ivtype = 91
    ovkeyw(ivtype) = 'SEDK'
    ovsnam(ivtype) = 'Sedk'
    ovlnam(ivtype) = 'Suspended sediment per layer'
    ovunit(ivtype) = ud
    ovsvty(ivtype) = 1
    ovllim(ivtype) = 0.
    ovulim(ivtype) = 1.e4
    ovlexp(ivtype) = 0.
    ovhexp(ivtype) = 1.e4
    ovexcv(ivtype) = -99.
    !
    ivtype = 92
    ovkeyw(ivtype) = 'MSALK'
    ovsnam(ivtype) = 'MSalk'
    ovlnam(ivtype) = 'Mean salinity per layer'
    ovunit(ivtype) = uc
    ovsvty(ivtype) = 1
    ovllim(ivtype) = 0.
    ovulim(ivtype) = 50.
    ovlexp(ivtype) = 0.
    ovhexp(ivtype) = 50.
    ovexcv(ivtype) = -99.
    !
    ivtype = 93
    ovkeyw(ivtype) = 'MTEMPK'
    ovsnam(ivtype) = 'MTempk'
    ovlnam(ivtype) = 'Mean temperature per layer'
    ovunit(ivtype) = uk
    ovsvty(ivtype) = 1
    ovllim(ivtype) = -100.
    ovulim(ivtype) = 100.
    ovlexp(ivtype) = -100.
    ovhexp(ivtype) = 100.
    ovexcv(ivtype) = -99.
    !
    ivtype = 94
    ovkeyw(ivtype) = 'MSEDK'
    ovsnam(ivtype) = 'MSedk'
    ovlnam(ivtype) = 'Mean suspended sediment per layer'
    ovunit(ivtype) = ud
    ovsvty(ivtype) = 1
    ovllim(ivtype) = 0.
    ovulim(ivtype) = 1.e4
    ovlexp(ivtype) = 0.
    ovhexp(ivtype) = 1.e4
    ovexcv(ivtype) = -99.
    !
    ivtype = 95
    ovkeyw(ivtype) = 'DIVU'
    ovsnam(ivtype) = 'divU'
    ovlnam(ivtype) = 'Horizontal divergence of flow velocity'
    ovunit(ivtype) = '1/s'
    ovsvty(ivtype) = 1
    ovllim(ivtype) = -100.
    ovulim(ivtype) = 100.
    ovlexp(ivtype) = -2.
    ovhexp(ivtype) = 2.
    ovexcv(ivtype) = 0.
    !
    ivtype = 96
    ovkeyw(ivtype) = 'DIVQ'
    ovsnam(ivtype) = 'divQ'
    ovlnam(ivtype) = 'Horizontal divergence of mass flux'
    ovunit(ivtype) = '1/s'
    ovsvty(ivtype) = 1
    ovllim(ivtype) = -100.
    ovulim(ivtype) = 100.
    ovlexp(ivtype) = -2.
    ovhexp(ivtype) = 2.
    ovexcv(ivtype) = 0.
    !
    ! hydrodynamic forces and moments
    !
    ivtype = 101
    ovkeyw(ivtype) = 'FORCEX'
    ovsnam(ivtype) = 'Fx'
    ovlnam(ivtype) = 'Force in x-direction'
    ovunit(ivtype) = 'kN'
    ovsvty(ivtype) = 4
    ovllim(ivtype) = -1.e-4
    ovulim(ivtype) = 1.e4
    ovlexp(ivtype) = -2.
    ovhexp(ivtype) = 2.
    ovexcv(ivtype) = -99.
    !
    ivtype = 102
    ovkeyw(ivtype) = 'FORCEY'
    ovsnam(ivtype) = 'Fy'
    ovlnam(ivtype) = 'Force in y-direction'
    ovunit(ivtype) = 'kN'
    ovsvty(ivtype) = 4
    ovllim(ivtype) = -1.e-4
    ovulim(ivtype) = 1.e4
    ovlexp(ivtype) = -2.
    ovhexp(ivtype) = 2.
    ovexcv(ivtype) = -99.
    !
    ivtype = 103
    ovkeyw(ivtype) = 'FORCEZ'
    ovsnam(ivtype) = 'Fz'
    ovlnam(ivtype) = 'Force in z-direction'
    ovunit(ivtype) = 'kN'
    ovsvty(ivtype) = 4
    ovllim(ivtype) = -1.e-4
    ovulim(ivtype) = 1.e4
    ovlexp(ivtype) = -1000.
    ovhexp(ivtype) = 1000.
    ovexcv(ivtype) = -99.
    !
    ivtype = 104
    ovkeyw(ivtype) = 'MOMX'
    ovsnam(ivtype) = 'Mx'
    ovlnam(ivtype) = 'Moment in x-direction'
    ovunit(ivtype) = 'kNm'
    ovsvty(ivtype) = 4
    ovllim(ivtype) = -1.e-4
    ovulim(ivtype) = 1.e4
    ovlexp(ivtype) = -2.
    ovhexp(ivtype) = 2.
    ovexcv(ivtype) = -99.
    !
    ivtype = 105
    ovkeyw(ivtype) = 'MOMY'
    ovsnam(ivtype) = 'My'
    ovlnam(ivtype) = 'Moment in y-direction'
    ovunit(ivtype) = 'kNm'
    ovsvty(ivtype) = 4
    ovllim(ivtype) = -1.e-4
    ovulim(ivtype) = 1.e4
    ovlexp(ivtype) = -2.
    ovhexp(ivtype) = 2.
    ovexcv(ivtype) = -99.
    !
    ivtype = 106
    ovkeyw(ivtype) = 'MOMZ'
    ovsnam(ivtype) = 'Mz'
    ovlnam(ivtype) = 'Moment in z-direction'
    ovunit(ivtype) = 'kNm'
    ovsvty(ivtype) = 4
    ovllim(ivtype) = -1.e-4
    ovulim(ivtype) = 1.e4
    ovlexp(ivtype) = -2.
    ovhexp(ivtype) = 2.
    ovexcv(ivtype) = -99.
    !
    ! rigid body motions
    !
    ivtype = 107
    ovkeyw(ivtype) = 'TRAX'
    ovsnam(ivtype) = 'Tx'
    ovlnam(ivtype) = 'Translation in x-direction'
    ovunit(ivtype) = 'm'
    ovsvty(ivtype) = 4
    ovllim(ivtype) = -1.e4
    ovulim(ivtype) = 1.e4
    ovlexp(ivtype) = -100.
    ovhexp(ivtype) = 100.
    ovexcv(ivtype) = -9999.
    !
    ivtype = 108
    ovkeyw(ivtype) = 'TRAY'
    ovsnam(ivtype) = 'Ty'
    ovlnam(ivtype) = 'Translation in y-direction'
    ovunit(ivtype) = 'm'
    ovsvty(ivtype) = 4
    ovllim(ivtype) = -1.e4
    ovulim(ivtype) = 1.e4
    ovlexp(ivtype) = -100.
    ovhexp(ivtype) = 100.
    ovexcv(ivtype) = -9999.
    !
    ivtype = 109
    ovkeyw(ivtype) = 'TRAZ'
    ovsnam(ivtype) = 'Tz'
    ovlnam(ivtype) = 'Translation in z-direction'
    ovunit(ivtype) = 'm'
    ovsvty(ivtype) = 4
    ovllim(ivtype) = -1.e4
    ovulim(ivtype) = 1.e4
    ovlexp(ivtype) = -100.
    ovhexp(ivtype) = 100.
    ovexcv(ivtype) = -9999.
    !
    ivtype = 110
    ovkeyw(ivtype) = 'ROTX'
    ovsnam(ivtype) = 'Rx'
    ovlnam(ivtype) = 'Rotation around x-axis'
    ovunit(ivtype) = 'degr'
    ovsvty(ivtype) = 4
    ovllim(ivtype) = -1.e4
    ovulim(ivtype) = 1.e4
    ovlexp(ivtype) = -100.
    ovhexp(ivtype) = 100.
    ovexcv(ivtype) = -9999.
    !
    ivtype = 111
    ovkeyw(ivtype) = 'ROTY'
    ovsnam(ivtype) = 'Ry'
    ovlnam(ivtype) = 'Rotation around y-axis'
    ovunit(ivtype) = 'degr'
    ovsvty(ivtype) = 4
    ovllim(ivtype) = -1.e4
    ovulim(ivtype) = 1.e4
    ovlexp(ivtype) = -100.
    ovhexp(ivtype) = 100.
    ovexcv(ivtype) = -9999.
    !
    ivtype = 112
    ovkeyw(ivtype) = 'ROTZ'
    ovsnam(ivtype) = 'Rz'
    ovlnam(ivtype) = 'Rotation around z-axis'
    ovunit(ivtype) = 'degr'
    ovsvty(ivtype) = 4
    ovllim(ivtype) = -1.e4
    ovulim(ivtype) = 1.e4
    ovlexp(ivtype) = -100.
    ovhexp(ivtype) = 100.
    ovexcv(ivtype) = -9999.
    !
    ivtype = 113
    ovkeyw(ivtype) = 'PTO'
    ovsnam(ivtype) = 'PTOpow'
    ovlnam(ivtype) = 'PTO power'
    ovunit(ivtype) = 'kW'
    ovsvty(ivtype) = 4
    ovllim(ivtype) = 0.
    ovulim(ivtype) = 1000.
    ovlexp(ivtype) = 0.
    ovhexp(ivtype) = 1000.
    ovexcv(ivtype) = -99.
    !
    ivtype = 115
    ovkeyw(ivtype) = 'RUNUP'
    ovsnam(ivtype) = 'Runlev'
    ovlnam(ivtype) = 'Runup height'
    ovunit(ivtype) = uh
    ovsvty(ivtype) = 1
    ovllim(ivtype) = -1.e4
    ovulim(ivtype) = 1.e4
    ovlexp(ivtype) = -100.
    ovhexp(ivtype) = 100.
    ovexcv(ivtype) = -9999.
    !
    ! various parameters for computation of output quantities
    !
    ! reference time for TSEC
    !
    outpar(1) = 0.
    !
    ! indicator for direction
    ! =0: direction always w.r.t. user coordinates; =1: direction w.r.t. frame
    !
    outpar(2) = 0.
    !
    ! indicator for interpolating output quantities in output frame
    ! =0: interpolation; =1: no interpolation
    !
    outpar(3) = 0.
    !
 101 format (/,20X,'---------------------------------------',   &
             /,20X,'                 SWASH',                    &
             /,20X,'     WAVE MOTION IN COASTAL WATERS     ',   &
             /,20X,'         VERSION NUMBER ', A,               &
             /,20X,'---------------------------------------',//)
    !
 102 format (/, ' SWASH is preparing computation',/)
    !
end subroutine SwashInit
