!
!  Routines for the solution of system of equations
!
!  Contents of this file
!
!     avm
!     ave
!     avmp
!     ivl
!     ivu
!     pcg  (single precision)
!     pcg2 (double precision)
!     sip
!     ilu
!     iludr
!     iluds
!     bicgstab
!     tridiag  (single precision)
!     tridiag2 (double precision)
!     dac  (single precision)
!     dac2 (double precision)
!     newton1D
!     newton2D
!
subroutine avm ( amat, vi, vo )
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
!    1.00, July 2010: New subroutine
!
!   Purpose
!
!   Performs the matrix-vector multiplication for the unpreconditioned system of equations
!
!   Modules used
!
    use ocpcomm4
    use SwashCommdata3, only: mcgrd, qmax, oned
    use SwashFlowdata, only: mfu, ml, nfu, nl, nconct, ishif
    use m_genarr, only: kgrpnt
!
    implicit none
!
!   Argument variables
!
    real, dimension(mcgrd,qmax,nconct), intent(in)  :: amat ! coefficient matrix of the system of equations
    real, dimension(mcgrd,qmax       ), intent(in)  :: vi   ! input vector to be multiplied
    real, dimension(mcgrd,qmax       ), intent(out) :: vo   ! output vector as result of multiplication
!
!   Local variables
!
    integer, save :: ient = 0 ! number of entries in this subroutine
    integer       :: j        ! loop counter
    integer       :: k        ! loop counter over vertical layers
    integer       :: kd       ! index of layer k-1
    integer       :: kdd      ! index of layer k-2
    integer       :: ku       ! index of layer k+1
    integer       :: kuu      ! index of layer k+2
    integer       :: m        ! loop counter
    integer       :: md       ! index of point m-1
    integer       :: mu       ! index of point m+1
    integer       :: n        ! loop counter
    integer       :: nd       ! index of point n-1
    integer       :: ndm      ! pointer to m,n-1
    integer       :: nm       ! pointer to m,n
    integer       :: nmd      ! pointer to m-1,n
    integer       :: nmu      ! pointer to m+1,n
    integer       :: nu       ! index of point n+1
    integer       :: num      ! pointer to m,n+1
!
!   Structure
!
!   perform matrix-vector multiplication vo = A*vi
!
!   Source text
!
    if (ltrace) call strace (ient,'avm')
    !
    if ( oned ) then
       !
       do m = mfu, ml
          !
          md = m - 1
          mu = m + 1
          !
          nm  = kgrpnt(m ,1)
          nmd = kgrpnt(md,1)
          nmu = kgrpnt(mu,1)
          !
          do k = 1, qmax
             !
             kd  = max(k-1,1   )
             kdd = max(k-2,1   )
             ku  = min(k+1,qmax)
             kuu = min(k+2,qmax)
             !
             vo(nm,k) = amat(nm,k, 2)*vi(nmd,k  ) + amat(nm,k, 3)*vi(nmu,k  ) + amat(nm,k, 4)*vi(nm ,kd ) +  &
                        amat(nm,k, 8)*vi(nmu,kd ) + amat(nm,k,10)*vi(nmd,kd ) + amat(nm,k,12)*vi(nmu,ku ) +  &
                        amat(nm,k,14)*vi(nmd,ku ) + amat(nm,k,16)*vi(nm ,kdd) + amat(nm,k,18)*vi(nmu,kdd) +  &
                        amat(nm,k,20)*vi(nmd,kdd) + amat(nm,k,22)*vi(nmu,kuu) + amat(nm,k,24)*vi(nmd,kuu)
             !
             ! note: nconct-23 is equivalent with kmax-1 and 2 for Keller-box and standard pressure gradient approximations, respectively
             !
             do j = 0, nconct-23
                !
                vo(nm,k) = vo(nm,k) + amat(nm,k,ishif(j))*vi(nm,min(k+j,qmax))
                !
             enddo
             !
          enddo
          !
       enddo
       !
    else
       !
       do n = nfu, nl
          do m = mfu, ml
             !
             md = m - 1
             mu = m + 1
             nd = n - 1
             nu = n + 1
             !
             nm  = kgrpnt(m ,n )
             nmd = kgrpnt(md,n )
             nmu = kgrpnt(mu,n )
             ndm = kgrpnt(m ,nd)
             num = kgrpnt(m ,nu)
             !
             do k = 1, qmax
                !
                kd  = max(k-1,1   )
                kdd = max(k-2,1   )
                ku  = min(k+1,qmax)
                kuu = min(k+2,qmax)
                !
                vo(nm,k) = amat(nm,k, 2)*vi(nmd,k  ) + amat(nm,k, 3)*vi(nmu,k  ) + amat(nm,k, 4)*vi(nm ,kd ) + amat(nm,k, 6)*vi(ndm,k  ) +  &
                           amat(nm,k, 7)*vi(num,k  ) + amat(nm,k, 8)*vi(nmu,kd ) + amat(nm,k, 9)*vi(num,kd ) + amat(nm,k,10)*vi(nmd,kd ) +  &
                           amat(nm,k,11)*vi(ndm,kd ) + amat(nm,k,12)*vi(nmu,ku ) + amat(nm,k,13)*vi(num,ku ) + amat(nm,k,14)*vi(nmd,ku ) +  &
                           amat(nm,k,15)*vi(ndm,ku ) + amat(nm,k,16)*vi(nm ,kdd) + amat(nm,k,18)*vi(nmu,kdd) + amat(nm,k,19)*vi(num,kdd) +  &
                           amat(nm,k,20)*vi(nmd,kdd) + amat(nm,k,21)*vi(ndm,kdd) + amat(nm,k,22)*vi(nmu,kuu) + amat(nm,k,23)*vi(num,kuu) +  &
                           amat(nm,k,24)*vi(nmd,kuu) + amat(nm,k,25)*vi(ndm,kuu)
                !
                do j = 0, nconct-23
                   !
                   vo(nm,k) = vo(nm,k) + amat(nm,k,ishif(j))*vi(nm,min(k+j,qmax))
                   !
                enddo
                !
             enddo
             !
          enddo
       enddo
       !
    endif
    !
end subroutine avm
!
subroutine ave ( amat, vi, vo, vt )
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
!    1.00, February 2012: New subroutine
!
!   Purpose
!
!   Performs the matrix-vector multiplication for the preconditioned system of equations by means of the Eisenstat trick
!
!   Method
!
!   The amount of work of matrix-vector multiplication for the preconditioned system is approximately two times
!   as much as for the unpreconditioned system. In the paper of Eisenstat it is shown that much of the extra work
!   can be avoided which is implemented in this subroutine.
!
!   S.C. Eisenstat
!   Efficient implementation of a class of preconditioned conjugate gradient methods
!   SIAM J. Sci. Stat. Comput., vol. 2, 1-4, 1981
!
!   Modules used
!
    use ocpcomm4
    use SwashCommdata3, only: mcgrd, qmax, oned
    use SwashFlowdata, only: mfu, ml, nfu, nl, nconct, ishif
    use m_genarr, only: kgrpnt
    use m_parall
!
    implicit none
!
!   Argument variables
!
    real, dimension(mcgrd,qmax,nconct), intent(in)  :: amat ! coefficient matrix of the system of equations
    real, dimension(mcgrd,qmax       ), intent(in)  :: vi   ! input vector to be multiplied
    real, dimension(mcgrd,qmax       ), intent(out) :: vo   ! output vector as result of multiplication
    real, dimension(mcgrd,qmax       ), intent(out) :: vt   ! auxiliary variable to store a vector temporarily
!
!   Local variables
!
    integer, save :: ient = 0 ! number of entries in this subroutine
    integer       :: j        ! loop counter
    integer       :: k        ! loop counter over vertical layers
    integer       :: kd       ! index of layer k-1
    integer       :: kdd      ! index of layer k-2
    integer       :: ku       ! index of layer k+1
    integer       :: kuu      ! index of layer k+2
    integer       :: m        ! loop counter
    integer       :: md       ! index of point m-1
    integer       :: mend     ! end index of loop over wl-points in x-direction
    integer       :: mu       ! index of point m+1
    integer       :: n        ! loop counter
    integer       :: nd       ! index of point n-1
    integer       :: ndm      ! pointer to m,n-1
    integer       :: nend     ! end index of loop over wl-points in y-direction
    integer       :: nm       ! pointer to m,n
    integer       :: nmd      ! pointer to m-1,n
    integer       :: nmu      ! pointer to m+1,n
    integer       :: nu       ! index of point n+1
    integer       :: num      ! pointer to m,n+1
!
!   Structure
!                                              -1  -1
!   perform matrix-vector multiplication vo = L  AU  vi  by means of the Eisenstat implementation, as follows:
!
!         -1      -1                     -1
!   vo = U  vi + L  ( vi + (diag(A)-2I) U  vi )
!
!   Hence, the following steps will be carried out:
!
!            -1
!   1) vt = U  vi
!
!            -1
!   2) vo = L  ( vi + (amat(1)-2)vt )
!
!   3) vo = vo + vt
!
!   Source text
!
    if (ltrace) call strace (ient,'ave')
    !
    ! not computed for points ml and nl at subdomain interfaces since, these end points are owned by the neighbouring subdomains
    !
    mend = ml - 1
    if ( LMXL ) mend = ml
    !
    nend = nl - 1
    if ( LMYL ) nend = nl
    !
    if ( oned ) then
       !
       do m = mend, mfu, -1
          !
          mu = m + 1
          !
          nm  = kgrpnt(m ,1)
          nmu = kgrpnt(mu,1)
          !
          do k = qmax, 1, -1
             !
             kd  = max(k-1,1   )
             kdd = max(k-2,1   )
             ku  = min(k+1,qmax)
             kuu = min(k+2,qmax)
             !
             vt(nm,k) = vi(nm,k) - amat(nm,k,3)*vt(nmu,k) - amat(nm,k,8)*vt(nmu,kd) - amat(nm,k,12)*vt(nmu,ku) - amat(nm,k,18)*vt(nmu,kdd) - amat(nm,k,22)*vt(nmu,kuu)
             !
             ! note: nconct-23 is equivalent with kmax-1 and 2 for Keller-box and standard pressure gradient approximations, respectively
             !
             do j = 1, nconct-23
                !
                vt(nm,k) = vt(nm,k) - amat(nm,k,ishif(j))*vt(nm,min(k+j,qmax))
                !
             enddo
             !
          enddo
          !
       enddo
       !
       do m = mfu, mend
          !
          md = m - 1
          !
          nm  = kgrpnt(m ,1)
          nmd = kgrpnt(md,1)
          !
          do k = 1, qmax
             !
             kd  = max(k-1,1   )
             kdd = max(k-2,1   )
             ku  = min(k+1,qmax)
             kuu = min(k+2,qmax)
             !
             vo(nm,k) = vi(nm,k) + (amat(nm,k,1) - 2.)*vt(nm,k)                                                                         &
                        - amat(nm,k, 2)*vo(nmd,k  ) - amat(nm,k, 4)*vo(nm ,kd ) - amat(nm,k,10)*vo(nmd,kd ) - amat(nm,k,14)*vo(nmd,ku)  &
                        - amat(nm,k,16)*vo(nm ,kdd) - amat(nm,k,20)*vo(nmd,kdd) - amat(nm,k,24)*vo(nmd,kuu)
             !
          enddo
          !
       enddo
       !
       do m = mfu, mend
          !
          nm = kgrpnt(m,1)
          !
          do k = 1, qmax
             !
             vo(nm,k) = vo(nm,k) + vt(nm,k)
             !
          enddo
          !
       enddo
       !
    else
       !
       do n = nend, nfu, -1
          do m = mend, mfu, -1
             !
             mu = m + 1
             nu = n + 1
             !
             nm  = kgrpnt(m ,n )
             nmu = kgrpnt(mu,n )
             num = kgrpnt(m ,nu)
             !
             do k = qmax, 1, -1
                !
                kd  = max(k-1,1   )
                kdd = max(k-2,1   )
                ku  = min(k+1,qmax)
                kuu = min(k+2,qmax)
                !
                vt(nm,k) = vi(nm,k) - amat(nm,k, 3)*vt(nmu,k  ) - amat(nm,k, 7)*vt(num,k  ) - amat(nm,k, 8)*vt(nmu,kd )  &
                                    - amat(nm,k, 9)*vt(num,kd ) - amat(nm,k,12)*vt(nmu,ku ) - amat(nm,k,13)*vt(num,ku )  &
                                    - amat(nm,k,18)*vt(nmu,kdd) - amat(nm,k,19)*vt(num,kdd) - amat(nm,k,22)*vt(nmu,kuu)  &
                                    - amat(nm,k,23)*vt(num,kuu)
                !
                do j = 1, nconct-23
                   !
                   vt(nm,k) = vt(nm,k) - amat(nm,k,ishif(j))*vt(nm,min(k+j,qmax))
                   !
                enddo
                !
             enddo
             !
          enddo
       enddo
       !
       do n = nfu, nend
          do m = mfu, mend
             !
             md = m - 1
             nd = n - 1
             !
             nm  = kgrpnt(m ,n )
             nmd = kgrpnt(md,n )
             ndm = kgrpnt(m ,nd)
             !
             do k = 1, qmax
                !
                kd  = max(k-1,1   )
                kdd = max(k-2,1   )
                ku  = min(k+1,qmax)
                kuu = min(k+2,qmax)
                !
                vo(nm,k) = vi(nm,k) + (amat(nm,k,1) - 2.)*vt(nm,k)                                                                          &
                           - amat(nm,k, 2)*vo(nmd,k  ) - amat(nm,k, 4)*vo(nm ,kd ) - amat(nm,k, 6)*vo(ndm,k  ) - amat(nm,k,10)*vo(nmd,kd )  &
                           - amat(nm,k,11)*vo(ndm,kd ) - amat(nm,k,14)*vo(nmd,ku ) - amat(nm,k,15)*vo(ndm,ku ) - amat(nm,k,16)*vo(nm ,kdd)  &
                           - amat(nm,k,20)*vo(nmd,kdd) - amat(nm,k,21)*vo(ndm,kdd) - amat(nm,k,24)*vo(nmd,kuu) - amat(nm,k,25)*vo(ndm,kuu)
                !
             enddo
             !
          enddo
       enddo
       !
       do n = nfu, nend
          do m = mfu, mend
             !
             nm = kgrpnt(m,n)
             !
             do k = 1, qmax
                !
                vo(nm,k) = vo(nm,k) + vt(nm,k)
                !
             enddo
             !
          enddo
       enddo
       !
    endif
    !
end subroutine ave
!
subroutine avmp ( amat, vi, vo, vt )
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
!    1.00, April 2010: New subroutine
!
!   Purpose
!
!   Performs the matrix-vector multiplication for the preconditioned system of equations by means of the Eisenstat trick
!
!   Method
!
!   The amount of work of matrix-vector multiplication for the preconditioned system is approximately two times
!   as much as for the unpreconditioned system. In the paper of Eisenstat it is shown that much of the extra work
!   can be avoided which is implemented in this subroutine.
!
!   S.C. Eisenstat
!   Efficient implementation of a class of preconditioned conjugate gradient methods
!   SIAM J. Sci. Stat. Comput., vol. 2, 1-4, 1981
!
!   Modules used
!
    use ocpcomm4
    use SwashCommdata3, only: mcgrd
    use SwashFlowdata, only: mfu, ml, nfu, nl
    use m_genarr, only: kgrpnt
    use m_parall
!
    implicit none
!
!   Argument variables
!
    real, dimension(mcgrd,5), intent(in)  :: amat ! coefficient matrix of the system of equations
    real, dimension(mcgrd  ), intent(in)  :: vi   ! input vector to be multiplied
    real, dimension(mcgrd  ), intent(out) :: vo   ! output vector as result of multiplication
    real, dimension(mcgrd  ), intent(out) :: vt   ! auxiliary variable to store a vector temporarily
!
!   Local variables
!
    integer, save :: ient = 0 ! number of entries in this subroutine
    integer       :: m        ! loop counter
    integer       :: md       ! index of point m-1
    integer       :: mend     ! end index of loop over wl-points in x-direction
    integer       :: mu       ! index of point m+1
    integer       :: n        ! loop counter
    integer       :: nd       ! index of point n-1
    integer       :: ndm      ! pointer to m,n-1
    integer       :: nend     ! end index of loop over wl-points in y-direction
    integer       :: nm       ! pointer to m,n
    integer       :: nmd      ! pointer to m-1,n
    integer       :: nmu      ! pointer to m+1,n
    integer       :: nu       ! index of point n+1
    integer       :: num      ! pointer to m,n+1
!
!   Structure
!                                              -1  -1
!   perform matrix-vector multiplication vo = L  AU  vi  by means of the Eisenstat implementation, as follows:
!
!         -1      -1                     -1
!   vo = U  vi + L  ( vi + (diag(A)-2I) U  vi )
!
!   Hence, the following steps will be carried out:
!
!            -1
!   1) vt = U  vi
!
!            -1
!   2) vo = L  ( vi + (amat(1)-2)vt )
!
!   3) vo = vo + vt
!
!   Source text
!
    if (ltrace) call strace (ient,'avmp')
    !
    ! not computed for points ml and nl at subdomain interfaces since, these end points are owned by the neighbouring subdomains
    !
    mend = ml - 1
    if ( LMXL ) mend = ml
    !
    nend = nl - 1
    if ( LMYL ) nend = nl
    !
    do n = nend, nfu, -1
       do m = mend, mfu, -1
          !
          mu = m + 1
          nu = n + 1
          !
          nm  = kgrpnt(m ,n )
          nmu = kgrpnt(mu,n )
          num = kgrpnt(m ,nu)
          !
          vt(nm) = vi(nm) - amat(nm,3)*vt(nmu) - amat(nm,5)*vt(num)
          !
       enddo
    enddo
    !
    do n = nfu, nend
       do m = mfu, mend
          !
          md = m - 1
          nd = n - 1
          !
          nm  = kgrpnt(m ,n )
          nmd = kgrpnt(md,n )
          ndm = kgrpnt(m ,nd)
          !
          vo(nm) = vi(nm) + (amat(nm,1) - 2.)*vt(nm) - amat(nm,2)*vo(nmd) - amat(nm,4)*vo(ndm)
          !
       enddo
    enddo
    !
    do n = nfu, nend
       do m = mfu, mend
          !
          nm  = kgrpnt(m,n)
          !
          vo(nm) = vo(nm) + vt(nm)
          !
       enddo
    enddo
    !
end subroutine avmp
!
subroutine ivl ( prec, v )
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
!    1.00, July 2010: New subroutine
!
!   Purpose
!
!   Performs the multiplication of a vector by the inverse of a lower triangular matrix
!
!   Modules used
!
    use ocpcomm4
    use SwashCommdata3, only: mcgrd, qmax, oned
    use SwashFlowdata, only: mfu, ml, nfu, nl, nconct
    use m_genarr, only: kgrpnt
!
    implicit none
!
!   Argument variables
!
    real, dimension(mcgrd,-1:qmax+2,nconct ), intent(in)    :: prec ! preconditioning matrix
    real, dimension(mcgrd,-1:qmax+nconct-23), intent(inout) :: v    ! vector to be multiplied
!
!   Local variables
!
    integer, save :: ient = 0 ! number of entries in this subroutine
    integer       :: k        ! loop counter over vertical layers
    integer       :: m        ! loop counter
    integer       :: md       ! index of point m-1
    integer       :: n        ! loop counter
    integer       :: nd       ! index of point n-1
    integer       :: ndm      ! pointer to m,n-1
    integer       :: nm       ! pointer to m,n
    integer       :: nmd      ! pointer to m-1,n
!
!   Structure
!                                         -1
!   perform matrix-vector multiplication L  * v by means of forward substitution
!
!   Note: v will be overwritten!
!
!   Source text
!
    if (ltrace) call strace (ient,'ivl')
    !
    if ( oned ) then
       !
       do m = mfu, ml
          !
          md = m - 1
          !
          nm  = kgrpnt(m ,1)
          nmd = kgrpnt(md,1)
          !
          do k = 1, qmax
             !
             v(nm,k) = v(nm,k) - prec(nm,k, 2)*v(nmd,k  ) - prec(nm,k, 4)*v(nm ,k-1) - prec(nm,k,10)*v(nmd,k-1) - prec(nm,k,14)*v(nmd,k+1)  &
                               - prec(nm,k,16)*v(nm ,k-2) - prec(nm,k,20)*v(nmd,k-2) - prec(nm,k,24)*v(nmd,k+2)
             !
          enddo
          !
       enddo
       !
    else
       !
       do n = nfu, nl
          do m = mfu, ml
             !
             md = m - 1
             nd = n - 1
             !
             nm  = kgrpnt(m ,n )
             nmd = kgrpnt(md,n )
             ndm = kgrpnt(m ,nd)
             !
             do k = 1, qmax
                !
                v(nm,k) = v(nm,k) - prec(nm,k, 2)*v(nmd,k  ) - prec(nm,k, 4)*v(nm ,k-1) - prec(nm,k, 6)*v(ndm,k  ) - prec(nm,k,10)*v(nmd,k-1)  &
                                  - prec(nm,k,11)*v(ndm,k-1) - prec(nm,k,14)*v(nmd,k+1) - prec(nm,k,15)*v(ndm,k+1) - prec(nm,k,16)*v(nm ,k-2)  &
                                  - prec(nm,k,20)*v(nmd,k-2) - prec(nm,k,21)*v(ndm,k-2) - prec(nm,k,24)*v(nmd,k+2) - prec(nm,k,25)*v(ndm,k+2)
                !
             enddo
             !
          enddo
       enddo
       !
    endif
    !
end subroutine ivl
!
subroutine ivu ( prec, v )
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
!    1.00, July 2010: New subroutine
!
!   Purpose
!
!   Performs the multiplication of a vector by the inverse of an upper triangular matrix
!
!   Modules used
!
    use ocpcomm4
    use SwashCommdata3, only: mcgrd, qmax, oned
    use SwashFlowdata, only: mfu, ml, nfu, nl, nconct, ishif
    use m_genarr, only: kgrpnt
    use SwashSolvedata, only: icond
!
    implicit none
!
!   Argument variables
!
    real, dimension(mcgrd,-1:qmax+2,nconct ), intent(in)    :: prec ! preconditioning matrix
    real, dimension(mcgrd,-1:qmax+nconct-23), intent(inout) :: v    ! vector to be multiplied
!
!   Local variables
!
    integer, save :: ient = 0 ! number of entries in this subroutine
    integer       :: j        ! loop counter
    integer       :: k        ! loop counter over vertical layers
    integer       :: m        ! loop counter
    integer       :: mu       ! index of point m+1
    integer       :: n        ! loop counter
    integer       :: nm       ! pointer to m,n
    integer       :: nmu      ! pointer to m+1,n
    integer       :: nu       ! index of point n+1
    integer       :: num      ! pointer to m,n+1
!
!   Structure
!                                         -1
!   perform matrix-vector multiplication U  * v by means of backward substitution
!
!   Note: v will be overwritten!
!
!   Source text
!
    if (ltrace) call strace (ient,'ivu')
    !
    if ( oned ) then
       !
       if ( icond == 2 ) then
          !
          do m = ml, mfu, -1
             !
             mu = m + 1
             !
             nm  = kgrpnt(m ,1)
             nmu = kgrpnt(mu,1)
             !
             do k = qmax, 1, -1
                !
                v(nm,k) = v(nm,k) - prec(nm,k,3)*v(nmu,k) - prec(nm,k,8)*v(nmu,k-1) - prec(nm,k,12)*v(nmu,k+1) - prec(nm,k,18)*v(nmu,k-2) - prec(nm,k,22)*v(nmu,k+2)
                !
                ! note: nconct-23 is equivalent with kmax-1 and 2 for Keller-box and standard pressure gradient approximations, respectively
                !
                do j = 1, nconct-23
                   !
                   v(nm,k) = v(nm,k) - prec(nm,k,ishif(j))*v(nm,k+j)
                   !
                enddo
                !
             enddo
             !
          enddo
          !
       else
          !
          do m = ml, mfu, -1
             !
             mu = m + 1
             !
             nm  = kgrpnt(m ,1)
             nmu = kgrpnt(mu,1)
             !
             do k = qmax, 1, -1
                !
                v(nm,k) = v(nm,k) - prec(nm,k,3)*v(nmu,k) - prec(nm,k,8)*v(nmu,k-1) - prec(nm,k,12)*v(nmu,k+1) - prec(nm,k,18)*v(nmu,k-2) - prec(nm,k,22)*v(nmu,k+2)
                !
                ! note: nconct-23 is equivalent with kmax-1 and 2 for Keller-box and standard pressure gradient approximations, respectively
                !
                do j = 1, nconct-23
                   !
                   v(nm,k) = v(nm,k) - prec(nm,k,ishif(j))*v(nm,k+j)
                   !
                enddo
                !
                v(nm,k) = v(nm,k)*prec(nm,k,1)
                !
             enddo
             !
          enddo
          !
       endif
       !
    else
       !
       if ( icond == 2 ) then
          !
          do n = nl, nfu, -1
             do m = ml, mfu, -1
                !
                mu = m + 1
                nu = n + 1
                !
                nm  = kgrpnt(m ,n )
                nmu = kgrpnt(mu,n )
                num = kgrpnt(m ,nu)
                !
                do k = qmax, 1, -1
                   !
                   v(nm,k) = v(nm,k) - prec(nm,k, 3)*v(nmu,k  ) - prec(nm,k, 7)*v(num,k  ) - prec(nm,k, 8)*v(nmu,k-1)  &
                                     - prec(nm,k, 9)*v(num,k-1) - prec(nm,k,12)*v(nmu,k+1) - prec(nm,k,13)*v(num,k+1)  &
                                     - prec(nm,k,18)*v(nmu,k-2) - prec(nm,k,19)*v(num,k-2) - prec(nm,k,22)*v(nmu,k+2)  &
                                     - prec(nm,k,23)*v(num,k+2)
                   !
                   do j = 1, nconct-23
                      !
                      v(nm,k) = v(nm,k) - prec(nm,k,ishif(j))*v(nm,k+j)
                      !
                   enddo
                   !
                enddo
                !
             enddo
          enddo
          !
       else
          !
          do n = nl, nfu, -1
             do m = ml, mfu, -1
                !
                mu = m + 1
                nu = n + 1
                !
                nm  = kgrpnt(m ,n )
                nmu = kgrpnt(mu,n )
                num = kgrpnt(m ,nu)
                !
                do k = qmax, 1, -1
                   !
                   v(nm,k) = v(nm,k) - prec(nm,k, 3)*v(nmu,k  ) - prec(nm,k, 7)*v(num,k  ) - prec(nm,k, 8)*v(nmu,k-1)  &
                                     - prec(nm,k, 9)*v(num,k-1) - prec(nm,k,12)*v(nmu,k+1) - prec(nm,k,13)*v(num,k+1)  &
                                     - prec(nm,k,18)*v(nmu,k-2) - prec(nm,k,19)*v(num,k-2) - prec(nm,k,22)*v(nmu,k+2)  &
                                     - prec(nm,k,23)*v(num,k+2)
                   !
                   do j = 1, nconct-23
                      !
                      v(nm,k) = v(nm,k) - prec(nm,k,ishif(j))*v(nm,k+j)
                      !
                   enddo
                   !
                   v(nm,k) = v(nm,k)*prec(nm,k,1)
                   !
                enddo
                !
             enddo
          enddo
          !
       endif
       !
    endif
    !
end subroutine ivu
!
subroutine pcg ( amat, rhs, x )
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
!    1.00, April 2010: New subroutine
!
!   Purpose
!
!   Solves system of equations by means of preconditioned conjugate gradient algorithm
!
!   Method
!
!   System of equations in single precision
!
!   The system is supposed to be symmetric and positive definite
!
!   A diagonal matrix D is computed in order to scale the system with the inverse of D
!
!   Based on the scaled system an incomplete decomposition is obtained, using the following rules:
!
!   A = LU - R
!
!   (a)  diag(L) = diag(U) = I,
!   (b)  the off-diagonal parts of L and U are equal to the corresponding parts of A
!   (c)  diag(LU) = diag(A)
!
!   This incomplete decomposition is used as the ILUD split preconditioner. If the last rule (c) is
!   replaced by
!
!   rowsum(LU) = rowsum(A)
!
!   the MILUD preconditioner is obtained. We also used a RILUD preconditioner, which is the
!   average of ILUD and MILUD. The contributions on the diagonal are weighted by the scalar amod,
!   amod = 0 corresponds with ILUD, whereas amod = 1 corresponds with MILUD.
!
!   Using this incomplete LU decomposition, the system Ax = b is split preconditioned as follows:
!
!    -1  -1     -1
!   L  AU  y = L  b    with  y = Ux
!
!   The well-known CG method is applied for the solution of preconditioned Ax = b
!
!   Modules used
!
    use ocpcomm4
    use SwashCommdata3, only: mcgrd, pnums
    use SwashFlowdata, only: mfu, ml, nfu, nl
    use m_genarr, only: kgrpnt
    use m_parall
    use SwashSolvedata
!
    implicit none
!
!   Argument variables
!
    real, dimension(mcgrd,5), intent(inout) :: amat ! coefficient matrix of the system of equations
    real, dimension(mcgrd  ), intent(inout) :: rhs  ! right-hand side of the system of equations
    real, dimension(mcgrd  ), intent(inout) :: x    ! solution of the system of equations
!
!   Parameter variables
!
    double precision, parameter :: epsmac = 1d-18 ! machine precision number
!
!   Local variables
!
    integer, save :: ient = 0 ! number of entries in this subroutine
    integer       :: j        ! iteration counter
    integer       :: m        ! loop counter
    integer       :: maxit    ! maximum number of iterations
    integer       :: md       ! index of point m-1
    integer       :: mend     ! end index of loop over wl-points in x-direction
    integer       :: mu       ! index of point m+1
    integer       :: n        ! loop counter
    integer       :: nd       ! index of point n-1
    integer       :: ndm      ! pointer to m,n-1
    integer       :: nend     ! end index of loop over wl-points in y-direction
    integer       :: nm       ! pointer to m,n
    integer       :: nmd      ! pointer to m-1,n
    integer       :: nmu      ! pointer to m+1,n
    integer       :: nu       ! index of point n+1
    integer       :: num      ! pointer to m,n+1
    !
    real          :: a        ! inner product of conjugate vector and matrix times conjugate vector
    real          :: alpha    ! coefficients in expansion of the solution formed by conjugate vectors
    real          :: amod     ! parameter used in the average of ILUD and MILUD preconditioners
    real          :: beta     ! coefficients in Gram-Schmidt orthogonalization
    real          :: epslin   ! required accuracy in the linear solver
    real          :: reps     ! accuracy with respect to the initial residual used in the following termination criterion:
                              !
                              !  ||b-Ax || < reps*||b-Ax ||
                              !        j                0
                              !
    real          :: rho      ! inner product of the residual vector
    real          :: rhold    ! inner product of the residual vector from previous iteration
    real          :: rnorm    ! 2-norm of residual vector
    real          :: rnrm0    ! 2-norm of initial residual vector
    real          :: ueps     ! minimal accuracy based on machine precision
    !
    logical       :: STPNOW   ! indicates that program must stop
!
!   Structure
!
!   Description of the pseudo code
!
!   Source text
!
    if (ltrace) call strace (ient,'pcg')
    !
    reps  = pnums(21)
    maxit = nint(pnums(24))
    amod  = pnums(26)
    !
    ! not computed for points ml and nl at subdomain interfaces since, these end points are owned by the neighbouring subdomains
    !
    mend = ml - 1
    if ( LMXL ) mend = ml
    !
    nend = nl - 1
    if ( LMYL ) nend = nl
    !
    ! compute diagonal elements for scaling the system
    !
    do n = nfu, nl
       do m = mfu, ml
          !
          md = m - 1
          nd = n - 1
          !
          nm  = kgrpnt(m ,n )
          nmd = kgrpnt(md,n )
          ndm = kgrpnt(m ,nd)
          !
          diag(nm) = amat(nm,1) - amat(nm,2)*((1.-amod)*amat(nmd,3) + amod*(amat(nmd,3)+amat(nmd,5))) / diag(nmd)  &
                                - amat(nm,4)*((1.-amod)*amat(ndm,5) + amod*(amat(ndm,3)+amat(ndm,5))) / diag(ndm)
          !
       enddo
    enddo
    !
    ! scale the system with the inverse of the diagonal matrix
    !
    do n = nfu, nl
       do m = mfu, ml
          !
          nm = kgrpnt(m,n)
          !
          amat(nm,:) = amat(nm,:) / diag(nm)
          rhs (nm  ) = rhs (nm  ) / diag(nm)
          !
       enddo
    enddo
    !
    ! compute initial residual
    !
    do n = nfu, nend
       do m = mfu, mend
          !
          md = m - 1
          mu = m + 1
          nd = n - 1
          nu = n + 1
          !
          nm  = kgrpnt(m ,n )
          nmd = kgrpnt(md,n )
          nmu = kgrpnt(mu,n )
          ndm = kgrpnt(m ,nd)
          num = kgrpnt(m ,nu)
          !
          resd(nm) = amat(nm,1)*x(nm) + amat(nm,2)*x(nmd) + amat(nm,3)*x(nmu) + amat(nm,4)*x(ndm) + amat(nm,5)*x(num)
          !
       enddo
    enddo
    !
    do n = nfu, nend
       do m = mfu, mend
          !
          nm = kgrpnt(m,n)
          !
          resd(nm) = rhs (nm) - resd(nm)
          z   (nm) = resd(nm)
          !
       enddo
    enddo
    !
    call SWEXCHG ( z, kgrpnt, 1, 1 )
    if (STPNOW()) return
    !
    ! multiply vector by the inverse of a lower triangular matrix
    !
    do n = nfu, nend
       do m = mfu, mend
          !
          md = m - 1
          nd = n - 1
          !
          nm  = kgrpnt(m ,n )
          nmd = kgrpnt(md,n )
          ndm = kgrpnt(m ,nd)
          !
          z(nm) = z(nm) - amat(nm,2)*z(nmd) - amat(nm,4)*z(ndm)
          !
       enddo
    enddo
    !
    call SWEXCHG ( z, kgrpnt, 1, 1 )
    if (STPNOW()) return
    !
    ! multiply vector by the inverse of an upper triangular matrix
    !
    do n = nend, nfu, -1
       do m = mend, mfu, -1
          !
          mu = m + 1
          nu = n + 1
          !
          nm  = kgrpnt(m ,n )
          nmu = kgrpnt(mu,n )
          num = kgrpnt(m ,nu)
          !
          z(nm) = z(nm) - amat(nm,3)*z(nmu) - amat(nm,5)*z(num)
          !
       enddo
    enddo
    !
    do n = nfu, nend
       do m = mfu, mend
          !
          nm = kgrpnt(m,n)
          !
          s(nm) = z(nm)
          !
       enddo
    enddo
    !
    rho = 0.
    do n = nfu, nend
       do m = mfu, mend
          !
          nm = kgrpnt(m,n)
          !
          rho = rho + resd(nm)*z(nm)
          !
       enddo
    enddo
    !
    call SWREDUCE ( rho, 1, SWREAL, SWSUM )
    rnrm0 = sqrt(rho)
    !
    epslin = reps*rnrm0
    ueps   = 1000.*real(epsmac)*rnrm0
    !
    if ( epslin < ueps .and. rnrm0 > 0. ) then
       !
       if ( iamout > 0 .and. INODE == MASTER ) then
          write (PRINTF, '(a)') ' ++ pcg: the required accuracy is too small'
          write (PRINTF, '(a,e12.6)') '         required accuracy    = ',epslin
          write (PRINTF, '(a,e12.6)') '         appropriate accuracy = ',ueps
       endif
       !
       epslin = ueps
       !
    endif
    !
    do j = 1, maxit
       !
       call SWEXCHG ( s, kgrpnt, 1, 1 )
       if (STPNOW()) return
       !
       do n = nfu, nl
          do m = mfu, ml
             !
             md = m - 1
             mu = m + 1
             nd = n - 1
             nu = n + 1
             !
             nm  = kgrpnt(m ,n )
             nmd = kgrpnt(md,n )
             nmu = kgrpnt(mu,n )
             ndm = kgrpnt(m ,nd)
             num = kgrpnt(m ,nu)
             !
             r(nm) = amat(nm,1)*s(nm) + amat(nm,2)*s(nmd) + amat(nm,3)*s(nmu) + amat(nm,4)*s(ndm) + amat(nm,5)*s(num)
             !
          enddo
       enddo
       !
       a = 0.
       do n = nfu, nend
          do m = mfu, mend
             !
             nm = kgrpnt(m,n)
             !
             a = a + s(nm)*r(nm)
             !
          enddo
       enddo
       call SWREDUCE ( a, 1, SWREAL, SWSUM )
       !
       alpha = rho / a
       !
       do n = nfu, nend
          do m = mfu, mend
             !
             nm = kgrpnt(m,n)
             !
             x   (nm) = x   (nm) + alpha*s(nm)
             resd(nm) = resd(nm) - alpha*r(nm)
             z   (nm) = resd(nm)
             !
          enddo
       enddo
       !
       call SWEXCHG ( z, kgrpnt, 1, 1 )
       if (STPNOW()) return
       !
       ! multiply vector by the inverse of a lower triangular matrix
       !
       do n = nfu, nend
          do m = mfu, mend
             !
             md = m - 1
             nd = n - 1
             !
             nm  = kgrpnt(m ,n )
             nmd = kgrpnt(md,n )
             ndm = kgrpnt(m ,nd)
             !
             z(nm) = z(nm) - amat(nm,2)*z(nmd) - amat(nm,4)*z(ndm)
             !
          enddo
       enddo
       !
       call SWEXCHG ( z, kgrpnt, 1, 1 )
       if (STPNOW()) return
       !
       ! multiply vector by the inverse of an upper triangular matrix
       !
       do n = nend, nfu, -1
          do m = mend, mfu, -1
             !
             mu = m + 1
             nu = n + 1
             !
             nm  = kgrpnt(m ,n )
             nmu = kgrpnt(mu,n )
             num = kgrpnt(m ,nu)
             !
             z(nm) = z(nm) - amat(nm,3)*z(nmu) - amat(nm,5)*z(num)
             !
          enddo
       enddo
       !
       rhold = rho
       !
       rho = 0.
       do n = nfu, nend
          do m = mfu, mend
             !
             nm = kgrpnt(m,n)
             !
             rho = rho + resd(nm)*z(nm)
             !
          enddo
       enddo
       call SWREDUCE ( rho, 1, SWREAL, SWSUM )
       rnorm = sqrt(rho)
       !
       if ( iamout == 2 .and. INODE == MASTER ) then
          write (PRINTF, '(a,i3,a,e12.6)') ' ++ pcg: iter = ',j,'    res = ',rnorm
       endif
       !
       if ( .not. rnorm > epslin ) exit
       !
       beta = rho / rhold
       !
       do n = nfu, nend
          do m = mfu, mend
             !
             nm = kgrpnt(m,n)
             !
             s(nm) = z(nm) + beta*s(nm)
             !
          enddo
       enddo
       !
    enddo
    !
    ! investigate the reason for stopping
    !
    if ( rnorm > epslin .and. iamout > 0 .and. INODE == MASTER ) then
       !
       write (PRINTF, '(a)') ' ++ pcg: no convergence for solving system of equations'
       write (PRINTF, '(a,i3)'   ) '         total number of iterations     = ',j
       write (PRINTF, '(a,e12.6)') '         2-norm of the residual         = ',rnorm
       write (PRINTF, '(a,e12.6)') '         required accuracy              = ',epslin
       !
    else if ( iamout == 3 .and. INODE == MASTER ) then
       !
       write (PRINTF, '(a,e12.6)') ' ++ pcg: 2-norm of the initial residual = ',rnrm0
       write (PRINTF, '(a,i3)'   ) '         total number of iterations     = ',j
       write (PRINTF, '(a,e12.6)') '         2-norm of the residual         = ',rnorm
       !
    endif
    !
end subroutine pcg
!
subroutine pcg2 ( amat, rhs, x )
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
!    1.00, April 2010: New subroutine
!
!   Purpose
!
!   Solves system of equations by means of preconditioned conjugate gradient algorithm
!
!   Method
!
!   System of equations in double precision
!
!   The system is supposed to be symmetric and positive definite
!
!   A diagonal matrix D is computed in order to scale the system with the inverse of D
!
!   Based on the scaled system an incomplete decomposition is obtained, using the following rules:
!
!   A = LU - R
!
!   (a)  diag(L) = diag(U) = I,
!   (b)  the off-diagonal parts of L and U are equal to the corresponding parts of A
!   (c)  diag(LU) = diag(A)
!
!   This incomplete decomposition is used as the ILUD split preconditioner. If the last rule (c) is
!   replaced by
!
!   rowsum(LU) = rowsum(A)
!
!   the MILUD preconditioner is obtained. We also used a RILUD preconditioner, which is the
!   average of ILUD and MILUD. The contributions on the diagonal are weighted by the scalar amod,
!   amod = 0 corresponds with ILUD, whereas amod = 1 corresponds with MILUD.
!
!   Using this incomplete LU decomposition, the system Ax = b is split preconditioned as follows:
!
!    -1  -1     -1
!   L  AU  y = L  b    with  y = Ux
!
!   The well-known CG method is applied for the solution of preconditioned Ax = b
!
!   Modules used
!
    use ocpcomm4
    use SwashCommdata3, only: mcgrd, pnums
    use SwashFlowdata, only: mfu, ml, nfu, nl
    use m_genarr, only: kgrpnt
    use m_parall
    use SwashSolvedata
!
    implicit none
!
!   Argument variables
!
    real*8, dimension(mcgrd,5), intent(inout) :: amat ! coefficient matrix of the system of equations
    real*8, dimension(mcgrd  ), intent(inout) :: rhs  ! right-hand side of the system of equations
    real*8, dimension(mcgrd  ), intent(inout) :: x    ! solution of the system of equations
!
!   Parameter variables
!
    double precision, parameter :: epsmac = 1d-18 ! machine precision number
!
!   Local variables
!
    integer, save :: ient = 0 ! number of entries in this subroutine
    integer       :: j        ! iteration counter
    integer       :: m        ! loop counter
    integer       :: maxit    ! maximum number of iterations
    integer       :: md       ! index of point m-1
    integer       :: mend     ! end index of loop over wl-points in x-direction
    integer       :: mu       ! index of point m+1
    integer       :: n        ! loop counter
    integer       :: nd       ! index of point n-1
    integer       :: ndm      ! pointer to m,n-1
    integer       :: nend     ! end index of loop over wl-points in y-direction
    integer       :: nm       ! pointer to m,n
    integer       :: nmd      ! pointer to m-1,n
    integer       :: nmu      ! pointer to m+1,n
    integer       :: nu       ! index of point n+1
    integer       :: num      ! pointer to m,n+1
    !
    real*8        :: a        ! inner product of conjugate vector and matrix times conjugate vector
    real*8        :: alpha    ! coefficients in expansion of the solution formed by conjugate vectors
    real*8        :: amod     ! parameter used in the average of ILUD and MILUD preconditioners
    real*8        :: beta     ! coefficients in Gram-Schmidt orthogonalization
    real*8        :: epslin   ! required accuracy in the linear solver
    real*8        :: reps     ! accuracy with respect to the initial residual used in the following termination criterion:
                              !
                              !  ||b-Ax || < reps*||b-Ax ||
                              !        j                0
                              !
    real*8        :: rho      ! inner product of the residual vector
    real*8        :: rhold    ! inner product of the residual vector from previous iteration
    real*8        :: rnorm    ! 2-norm of residual vector
    real*8        :: rnrm0    ! 2-norm of initial residual vector
    real*8        :: ueps     ! minimal accuracy based on machine precision
    !
    logical       :: STPNOW   ! indicates that program must stop
!
!   Structure
!
!   Description of the pseudo code
!
!   Source text
!
    if (ltrace) call strace (ient,'pcg2')
    !
    reps  = dble(pnums(21))
    maxit = nint(pnums(24))
    amod  = dble(pnums(26))
    !
    ! not computed for points ml and nl at subdomain interfaces since, these end points are owned by the neighbouring subdomains
    !
    mend = ml - 1
    if ( LMXL ) mend = ml
    !
    nend = nl - 1
    if ( LMYL ) nend = nl
    !
    ! compute diagonal elements for scaling the system
    !
    do n = nfu, nl
       do m = mfu, ml
          !
          md = m - 1
          nd = n - 1
          !
          nm  = kgrpnt(m ,n )
          nmd = kgrpnt(md,n )
          ndm = kgrpnt(m ,nd)
          !
          diag2(nm) = amat(nm,1) - amat(nm,2)*((1d0-amod)*amat(nmd,3) + amod*(amat(nmd,3)+amat(nmd,5))) / diag2(nmd)  &
                                 - amat(nm,4)*((1d0-amod)*amat(ndm,5) + amod*(amat(ndm,3)+amat(ndm,5))) / diag2(ndm)
          !
       enddo
    enddo
    !
    ! scale the system with the inverse of the diagonal matrix
    !
    do n = nfu, nl
       do m = mfu, ml
          !
          nm = kgrpnt(m,n)
          !
          amat(nm,:) = amat(nm,:) / diag2(nm)
          rhs (nm  ) = rhs (nm  ) / diag2(nm)
          !
       enddo
    enddo
    !
    ! compute initial residual
    !
    do n = nfu, nend
       do m = mfu, mend
          !
          md = m - 1
          mu = m + 1
          nd = n - 1
          nu = n + 1
          !
          nm  = kgrpnt(m ,n )
          nmd = kgrpnt(md,n )
          nmu = kgrpnt(mu,n )
          ndm = kgrpnt(m ,nd)
          num = kgrpnt(m ,nu)
          !
          resd2(nm) = amat(nm,1)*x(nm) + amat(nm,2)*x(nmd) + amat(nm,3)*x(nmu) + amat(nm,4)*x(ndm) + amat(nm,5)*x(num)
          !
       enddo
    enddo
    !
    do n = nfu, nend
       do m = mfu, mend
          !
          nm = kgrpnt(m,n)
          !
          resd2(nm) = rhs  (nm) - resd2(nm)
          z3   (nm) = resd2(nm)
          !
       enddo
    enddo
    !
    call SWEXCHGD ( z3, kgrpnt, 1, 1 )
    if (STPNOW()) return
    !
    ! multiply vector by the inverse of a lower triangular matrix
    !
    do n = nfu, nend
       do m = mfu, mend
          !
          md = m - 1
          nd = n - 1
          !
          nm  = kgrpnt(m ,n )
          nmd = kgrpnt(md,n )
          ndm = kgrpnt(m ,nd)
          !
          z3(nm) = z3(nm) - amat(nm,2)*z3(nmd) - amat(nm,4)*z3(ndm)
          !
       enddo
    enddo
    !
    call SWEXCHGD ( z3, kgrpnt, 1, 1 )
    if (STPNOW()) return
    !
    ! multiply vector by the inverse of an upper triangular matrix
    !
    do n = nend, nfu, -1
       do m = mend, mfu, -1
          !
          mu = m + 1
          nu = n + 1
          !
          nm  = kgrpnt(m ,n )
          nmu = kgrpnt(mu,n )
          num = kgrpnt(m ,nu)
          !
          z3(nm) = z3(nm) - amat(nm,3)*z3(nmu) - amat(nm,5)*z3(num)
          !
       enddo
    enddo
    !
    do n = nfu, nend
       do m = mfu, mend
          !
          nm = kgrpnt(m,n)
          !
          s2(nm) = z3(nm)
          !
       enddo
    enddo
    !
    rho = 0d0
    do n = nfu, nend
       do m = mfu, mend
          !
          nm = kgrpnt(m,n)
          !
          rho = rho + resd2(nm)*z3(nm)
          !
       enddo
    enddo
    !
    call SWREDUCD ( rho, 1, SWSUM )
    rnrm0 = sqrt(rho)
    !
    epslin = reps*rnrm0
    ueps   = 1d3*epsmac*rnrm0
    !
    if ( epslin < ueps .and. rnrm0 > 0d0 ) then
       !
       if ( iamout > 0 .and. INODE == MASTER ) then
          write (PRINTF, '(a)') ' ++ pcg: the required accuracy is too small'
          write (PRINTF, '(a,e12.6)') '         required accuracy    = ',epslin
          write (PRINTF, '(a,e12.6)') '         appropriate accuracy = ',ueps
       endif
       !
       epslin = ueps
       !
    endif
    !
    do j = 1, maxit
       !
       call SWEXCHGD ( s2, kgrpnt, 1, 1 )
       if (STPNOW()) return
       !
       do n = nfu, nl
          do m = mfu, ml
             !
             md = m - 1
             mu = m + 1
             nd = n - 1
             nu = n + 1
             !
             nm  = kgrpnt(m ,n )
             nmd = kgrpnt(md,n )
             nmu = kgrpnt(mu,n )
             ndm = kgrpnt(m ,nd)
             num = kgrpnt(m ,nu)
             !
             r2(nm) = amat(nm,1)*s2(nm) + amat(nm,2)*s2(nmd) + amat(nm,3)*s2(nmu) + amat(nm,4)*s2(ndm) + amat(nm,5)*s2(num)
             !
          enddo
       enddo
       !
       a = 0d0
       do n = nfu, nend
          do m = mfu, mend
             !
             nm = kgrpnt(m,n)
             !
             a = a + s2(nm)*r2(nm)
             !
          enddo
       enddo
       call SWREDUCD ( a, 1, SWSUM )
       !
       alpha = rho / a
       !
       do n = nfu, nend
          do m = mfu, mend
             !
             nm = kgrpnt(m,n)
             !
             x    (nm) = x    (nm) + alpha*s2(nm)
             resd2(nm) = resd2(nm) - alpha*r2(nm)
             z3   (nm) = resd2(nm)
             !
          enddo
       enddo
       !
       call SWEXCHGD ( z3, kgrpnt, 1, 1 )
       if (STPNOW()) return
       !
       ! multiply vector by the inverse of a lower triangular matrix
       !
       do n = nfu, nend
          do m = mfu, mend
             !
             md = m - 1
             nd = n - 1
             !
             nm  = kgrpnt(m ,n )
             nmd = kgrpnt(md,n )
             ndm = kgrpnt(m ,nd)
             !
             z3(nm) = z3(nm) - amat(nm,2)*z3(nmd) - amat(nm,4)*z3(ndm)
             !
          enddo
       enddo
       !
       call SWEXCHGD ( z3, kgrpnt, 1, 1 )
       if (STPNOW()) return
       !
       ! multiply vector by the inverse of an upper triangular matrix
       !
       do n = nend, nfu, -1
          do m = mend, mfu, -1
             !
             mu = m + 1
             nu = n + 1
             !
             nm  = kgrpnt(m ,n )
             nmu = kgrpnt(mu,n )
             num = kgrpnt(m ,nu)
             !
             z3(nm) = z3(nm) - amat(nm,3)*z3(nmu) - amat(nm,5)*z3(num)
             !
          enddo
       enddo
       !
       rhold = rho
       !
       rho = 0d0
       do n = nfu, nend
          do m = mfu, mend
             !
             nm = kgrpnt(m,n)
             !
             rho = rho + resd2(nm)*z3(nm)
             !
          enddo
       enddo
       call SWREDUCD ( rho, 1, SWSUM )
       rnorm = sqrt(rho)
       !
       if ( iamout == 2 .and. INODE == MASTER ) then
          write (PRINTF, '(a,i3,a,e12.6)') ' ++ pcg: iter = ',j,'    res = ',rnorm
       endif
       !
       if ( .not. rnorm > epslin ) exit
       !
       beta = rho / rhold
       !
       do n = nfu, nend
          do m = mfu, mend
             !
             nm = kgrpnt(m,n)
             !
             s2(nm) = z3(nm) + beta*s2(nm)
             !
          enddo
       enddo
       !
    enddo
    !
    ! investigate the reason for stopping
    !
    if ( rnorm > epslin .and. iamout > 0 .and. INODE == MASTER ) then
       !
       write (PRINTF, '(a)') ' ++ pcg: no convergence for solving system of equations'
       write (PRINTF, '(a,i3)'   ) '         total number of iterations     = ',j
       write (PRINTF, '(a,e12.6)') '         2-norm of the residual         = ',rnorm
       write (PRINTF, '(a,e12.6)') '         required accuracy              = ',epslin
       !
    else if ( iamout == 3 .and. INODE == MASTER ) then
       !
       write (PRINTF, '(a,e12.6)') ' ++ pcg: 2-norm of the initial residual = ',rnrm0
       write (PRINTF, '(a,i3)'   ) '         total number of iterations     = ',j
       write (PRINTF, '(a,e12.6)') '         2-norm of the residual         = ',rnorm
       !
    endif
    !
end subroutine pcg2
!
subroutine sip ( amat, rhs, x )
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
!    1.00, April 2010: New subroutine
!
!   Purpose
!
!   Solves system of equations by means of Stone's SIP solver
!
!   Method
!
!   The system of equations is solved using an incomplete factorization technique called Strongly Implicit Procedure (SIP) as described in
!
!   H.L. Stone
!   Iterative solution of implicit approximations of multidimensional partial differential equations
!   SIAM J. of Numer. Anal., vol. 5, 530-558, 1968
!
!   This method constructs an incomplete lower-upper factorization that has the same sparsity as the original matrix.
!   Hereby, a parameter 0 <= alfa <= 1 is used, which should be around 0.91 (when alfa > 0.95, the method may diverge).
!   Furthermore, alfa = 0 means standard ILU decomposition.
!
!   Afterward, the resulting system is solved in an iterative manner by forward and backward substitutions.
!
!   Modules used
!
    use ocpcomm4
    use SwashCommdata3, only: mcgrd, pnums
    use SwashFlowdata, only: mfu, ml, nfu, nl
    use m_genarr, only: kgrpnt
    use m_parall
    use SwashSolvedata, only: iamout, cmat, resd
!
    implicit none
!
!   Argument variables
!
    real, dimension(mcgrd,5), intent(inout) :: amat ! coefficient matrix of the system of equations
    real, dimension(mcgrd  ), intent(inout) :: rhs  ! right-hand side of the system of equations
    real, dimension(mcgrd  ), intent(inout) :: x    ! solution of the system of equations
!
!   Parameter variables
!
    double precision, parameter :: epsmac = 1d-18   ! machine precision number
    real            , parameter :: small  = 1.e-15  ! a small number
!
!   Local variables
!
    integer, save :: ient = 0 ! number of entries in this subroutine
    integer       :: j        ! iteration counter
    integer       :: m        ! loop counter
    integer       :: maxit    ! maximum number of iterations
    integer       :: md       ! index of point m-1
    integer       :: mend     ! end index of loop over wl-points in x-direction
    integer       :: mm       ! counter
    integer       :: msta     ! start index of loop over wl-points in x-direction
    integer       :: mu       ! index of point m+1
    integer       :: n        ! loop counter
    integer       :: nd       ! index of point n-1
    integer       :: ndm      ! pointer to m,n-1
    integer       :: ndum     ! dummy loop counter
    integer       :: nend     ! end index of loop over wl-points in y-direction
    integer       :: nm       ! pointer to m,n
    integer       :: nmd      ! pointer to m-1,n
    integer       :: nmu      ! pointer to m+1,n
    integer       :: nn       ! counter
    integer       :: nsta     ! start index of loop over wl-points in y-direction
    integer       :: nu       ! index of point n+1
    integer       :: num      ! pointer to m,n+1
    !
    real          :: alfa     ! relaxation parameter
    real          :: bnorm    ! 2-norm of right-hand side vector
    real          :: epslin   ! required accuracy in the linear solver
    real          :: p1       ! auxiliary factor
    real          :: p2       ! auxiliary factor
    real          :: p3       ! auxiliary factor
    real          :: reps     ! accuracy with respect to the right-hand side used in the following termination criterion:
                              !
                              !  ||b-Ax || < reps*||b||
                              !        j
                              !
    real          :: rnorm    ! 2-norm of residual vector
    real          :: rnrm0    ! 2-norm of initial residual vector
    real          :: ueps     ! minimal accuracy based on machine precision
    !
    logical       :: STPNOW   ! indicates that program must stop
!
!   Structure
!
!   Description of the pseudo code
!
!   Source text
!
    if (ltrace) call strace (ient,'sip')
    !
    reps  = pnums(22)
    maxit = nint(pnums(25))
    alfa  = pnums(27)
    !
    ! construct L and U matrices (stored in cmat)
    !
    ! Note: execute loop over rows/columns in a pipelined parallel manner
    !
    if ( kpart == 4 ) then
       !
       msta = mfu
       mend = ml
       if ( .not.LMXL ) mend = mend - 1
       nsta = nfu
       nend = nl
       if ( .not.LMYL ) nend = nend - 1
       !
    else
       !
       msta = nfu
       mend = nl
       if ( .not.LMYL ) mend = mend - 1
       nsta = mfu
       nend = ml
       if ( .not.LMXL ) nend = nend - 1
       !
    endif
    !
    do ndum = nsta, nend+NPROC-1
       !
       nn = ndum - INODE +1
       !
       if ( INODE >= ndum-nend+1 .and. INODE <= ndum-nsta+1 ) then
          !
!TIMG          call SWTSTA(213)
          call SWRECVFLD( cmat, 2, msta, nn, 1, kgrpnt )
!TIMG          call SWTSTO(213)
          if (STPNOW()) return
          !
          do mm = msta, mend
             !
             if ( kpart == 4 ) then
                !
                m = mm
                n = nn
                !
             else
                !
                m = nn
                n = mm
                !
             endif
             !
             md = m - 1
             nd = n - 1
             !
             nm   = kgrpnt(m ,n )
             nmd  = kgrpnt(md,n )
             ndm  = kgrpnt(m ,nd)
             !
             p1 = alfa*cmat(nmd,2)
             p2 = alfa*cmat(ndm,1)
             !
             cmat(nm,4) = amat(nm,2)/(1.+p1)
             cmat(nm,5) = amat(nm,4)/(1.+p2)
             !
             p1 = p1*cmat(nm,4)
             p2 = p2*cmat(nm,5)
             p3 = amat(nm,1) + p1 + p2 - cmat(nm,4)*cmat(nmd,1) - cmat(nm,5)*cmat(ndm,2) + small
             !
             cmat(nm,3) = 1./p3
             cmat(nm,1) = (amat(nm,3)-p2)*cmat(nm,3)
             cmat(nm,2) = (amat(nm,5)-p1)*cmat(nm,3)
             !
          enddo
          !
!TIMG          call SWTSTA(213)
          call SWSENDFLD( cmat, 2, mend, nn, 1, kgrpnt )
!TIMG          call SWTSTO(213)
          if (STPNOW()) return
          !
       endif
       !
    enddo
    !
    ! compute 2-norm of right-hand side vector
    !
    bnorm = 0.
    !
    ! not computed for points ml and nl at subdomain interfaces since, these end points are owned by the neighbouring subdomains
    !
    mend = ml - 1
    if ( LMXL ) mend = ml
    !
    nend = nl - 1
    if ( LMYL ) nend = nl
    !
    do n = nfu, nend
       do m = mfu, mend
          !
          nm = kgrpnt(m,n)
          !
          bnorm = bnorm + rhs(nm)*rhs(nm)
          !
       enddo
    enddo
    !
    call SWREDUCE ( bnorm, 1, SWREAL, SWSUM )
    bnorm = sqrt(bnorm)
    !
    epslin = reps*bnorm
    ueps   = 1000.*real(epsmac)*bnorm
    !
    if ( epslin < ueps .and. bnorm > 0. ) then
       !
       if ( iamout > 0 .and. INODE == MASTER ) then
          write (PRINTF, '(a)') ' ++ sip: the required accuracy is too small'
          write (PRINTF, '(a,e12.6)') '         required accuracy    = ',epslin
          write (PRINTF, '(a,e12.6)') '         appropriate accuracy = ',ueps
       endif
       !
       epslin = ueps
       !
    endif
    !
    ! solve the system by forward and backward substitutions in an iterative manner
    !
    do j = 1, maxit
       !
       ! compute residual and its norm
       !
       rnorm = 0.
       !
       ! not computed for points ml and nl at subdomain interfaces since, these end points are owned by the neighbouring subdomains
       !
       mend = ml - 1
       if ( LMXL ) mend = ml
       !
       nend = nl - 1
       if ( LMYL ) nend = nl
       !
       do n = nfu, nend
          do m = mfu, mend
             !
             md = m - 1
             mu = m + 1
             nd = n - 1
             nu = n + 1
             !
             nm  = kgrpnt(m ,n )
             nmu = kgrpnt(mu,n )
             nmd = kgrpnt(md,n )
             num = kgrpnt(m ,nu)
             ndm = kgrpnt(m ,nd)
             !
             resd(nm) = rhs(nm) - amat(nm,1)*x(nm) - amat(nm,2)*x(nmd) - amat(nm,3)*x(nmu) - amat(nm,4)*x(ndm) - amat(nm,5)*x(num)
             rnorm = rnorm + resd(nm)*resd(nm)
             !
          enddo
       enddo
       !
       call SWEXCHG ( resd, kgrpnt, 1, 1 )
       if (STPNOW()) return
       !
       ! forward substitution
       !
       ! Note: execute loop over rows/columns in a pipelined parallel manner
       !
       if ( kpart == 4 ) then
          !
          msta = mfu
          mend = ml
          if ( .not.LMXL ) mend = mend - 1
          nsta = nfu
          nend = nl
          if ( .not.LMYL ) nend = nend - 1
          !
       else
          !
          msta = nfu
          mend = nl
          if ( .not.LMYL ) mend = mend - 1
          nsta = mfu
          nend = ml
          if ( .not.LMXL ) nend = nend - 1
          !
       endif
       !
       do ndum = nsta, nend+NPROC-1
          !
          nn = ndum - INODE +1
          !
          if ( INODE >= ndum-nend+1 .and. INODE <= ndum-nsta+1 ) then
             !
!TIMG             call SWTSTA(213)
             call SWRECVFLD( resd, 1, msta, nn, 1, kgrpnt )
!TIMG             call SWTSTO(213)
             if (STPNOW()) return
             !
             do mm = msta, mend
                !
                if ( kpart == 4 ) then
                   !
                   m = mm
                   n = nn
                   !
                else
                   !
                   m = nn
                   n = mm
                   !
                endif
                !
                md = m - 1
                nd = n - 1
                !
                nm   = kgrpnt(m ,n )
                nmd  = kgrpnt(md,n )
                ndm  = kgrpnt(m ,nd)
                !
                resd(nm) = cmat(nm,3) * (resd(nm) - cmat(nm,4)*resd(nmd) - cmat(nm,5)*resd(ndm))
                !
             enddo
             !
!TIMG             call SWTSTA(213)
             call SWSENDFLD( resd, 1, mend, nn, 1, kgrpnt )
!TIMG             call SWTSTO(213)
             if (STPNOW()) return
             !
          endif
          !
       enddo
       !
       call SWEXCHG ( resd, kgrpnt, 1, 1 )
       if (STPNOW()) return
       !
       ! backward substitution and update solution
       !
       ! Note: execute loop over rows/columns in a pipelined parallel manner
       !
       do ndum = nend+NPROC-1, nsta, -1
          !
          nn = ndum - NPROC + INODE
          !
          if ( NPROC+1-INODE >= ndum-nend+1 .and. NPROC+1-INODE <= ndum-nsta+1 ) then
             !
!TIMG             call SWTSTA(213)
             call SWRECVFLD( resd, 1, mend, nn, 2, kgrpnt )
!TIMG             call SWTSTO(213)
             if (STPNOW()) return
             !
             do mm = mend, msta, -1
                !
                if ( kpart == 4 ) then
                   !
                   m = mm
                   n = nn
                   !
                else
                   !
                   m = nn
                   n = mm
                   !
                endif
                !
                mu = m + 1
                nu = n + 1
                !
                nm   = kgrpnt(m ,n )
                nmu  = kgrpnt(mu,n )
                num  = kgrpnt(m ,nu)
                !
                resd(nm) = resd(nm) - cmat(nm,1)*resd(nmu) - cmat(nm,2)*resd(num)
                x   (nm) = x   (nm) + resd(nm)
                !
             enddo
             !
!TIMG             call SWTSTA(213)
             call SWSENDFLD( resd, 1, msta, nn, 2, kgrpnt )
!TIMG             call SWTSTO(213)
             if (STPNOW()) return
             !
          endif
          !
       enddo
       !
       call SWREDUCE ( rnorm, 1, SWREAL, SWSUM )
       !
       rnorm = sqrt(rnorm)
       if ( iamout == 3 .and. j == 1 ) rnrm0 = rnorm
       !
       if ( iamout == 2 .and. INODE == MASTER ) then
          write (PRINTF, '(a,i3,a,e12.6)') ' ++ sip: iter = ',j,'    res = ',rnorm
       endif
       !
       if ( .not. rnorm > epslin ) exit
       !
       call SWEXCHG ( x, kgrpnt, 1, 1 )
       if (STPNOW()) return
       !
    enddo
    !
    ! investigate the reason for stopping
    !
    if ( rnorm > epslin .and. iamout > 0 .and. INODE == MASTER ) then
       !
       write (PRINTF, '(a)') ' ++ sip: no convergence for solving system of equations'
       write (PRINTF, '(a,i3)'   ) '         total number of iterations     = ',j
       write (PRINTF, '(a,e12.6)') '         2-norm of the residual         = ',rnorm
       write (PRINTF, '(a,e12.6)') '         required accuracy              = ',epslin
       !
    else if ( iamout == 3 .and. INODE == MASTER ) then
       !
       write (PRINTF, '(a,e12.6)') ' ++ sip: 2-norm of the initial residual = ',rnrm0
       write (PRINTF, '(a,i3)'   ) '         total number of iterations     = ',j
       write (PRINTF, '(a,e12.6)') '         2-norm of the residual         = ',rnorm
       !
    endif
    !
end subroutine sip
!
subroutine ilu ( amat )
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
!    1.00, July 2010: New subroutine
!
!   Purpose
!
!   Computes classical incomplete factorization (RILU)
!
!   Method
!
!   An upper triangular matrix U and a lower triangular matrix L are computed according to the following rules:
!
!   A = LU - R
!
!   (a)  diag(L) = I,
!   (b)  the non-zero pattern of L and U are equal to the non-zero pattern of A
!   (c)  if A(i,j) <> 0 then L*U(i,j) = A(i,j)
!
!   This incomplete decomposition is used as the ILU right preconditioner. If the last rule (c) is replaced by
!
!   rowsum(LU) = rowsum(A)
!
!   the MILU preconditioner is obtained. We also used a RILU preconditioner, which is the average of ILU and MILU.
!   The contributions on the diagonal are weighted by the scalar amod: amod = 0 corresponds with ILU, whereas
!   amod = 1 corresponds with MILU.
!
!   Modules used
!
    use ocpcomm4
    use SwashCommdata3, only: mcgrd, qmax, oned, pnums
    use SwashFlowdata, only: mfu, ml, nfu, nl, nconct, ishif
    use m_genarr, only: kgrpnt
    use m_parall
    use SwashSolvedata
!
    implicit none
!
!   Argument variables
!
    real, dimension(mcgrd,qmax,nconct), intent(in) :: amat ! coefficient matrix of the system of equations
!
!   Local variables
!
    integer, save            :: ient = 0 ! number of entries in this subroutine
    integer                  :: j        ! loop counter
    integer                  :: jmax     ! =kmax-1 if Keller-box scheme is applied
                                         ! =2      if standard pressure gradient approximation is employed
    integer                  :: k        ! loop counter over vertical layers
    integer                  :: kd       ! index of layer k-1
    integer                  :: kdd      ! index of layer k-2
    integer                  :: kddd     ! index of layer k-3
    integer                  :: kdddd    ! index of layer k-4
    integer                  :: ku       ! index of layer k+1
    integer                  :: kuu      ! index of layer k+2
    integer                  :: m        ! loop counter
    integer                  :: md       ! index of point m-1
    integer                  :: mend     ! end index of loop over wl-points in x-direction
    integer                  :: mu       ! index of point m+1
    integer                  :: n        ! loop counter
    integer                  :: nd       ! index of point n-1
    integer                  :: nend     ! end index of loop over wl-points in y-direction
    integer                  :: ndm      ! pointer to m,n-1
    integer                  :: nm       ! pointer to m,n
    integer                  :: nmd      ! pointer to m-1,n
    integer                  :: nmu      ! pointer to m+1,n
    integer                  :: nu       ! index of point n+1
    integer                  :: num      ! pointer to m,n+1
    !
    real                     :: amod     ! parameter used in the average of ILU and MILU preconditioners
    real                     :: fac      ! factor meant for MILU approach
    real                     :: jl       ! =1. if index j   exists otherwise 0.
    real                     :: jld      ! =1. if index j-1 exists otherwise 0.
    real                     :: jldd     ! =1. if index j-2 exists otherwise 0.
    real                     :: jlu      ! =1. if index j+1 exists otherwise 0.
    real                     :: jlup     ! =0. if index jmax otherwise 1.
    real                     :: jluu     ! =1. if index j+2 exists otherwise 0.
    real                     :: jluup    ! =0. if index jmax-1 or jmax otherwise 1.
    real                     :: kld      ! =1. if index k-1 exists otherwise 0.
    real                     :: kldd     ! =1. if index k-2 exists otherwise 0.
    real                     :: klddd    ! =1. if index k-3 exists otherwise 0.
    real                     :: kldddd   ! =1. if index k-4 exists otherwise 0.
    real                     :: klu      ! =1. if index k+1 exists otherwise 0.
    real                     :: kluu     ! =1. if index k+2 exists otherwise 0.
    !
    integer, dimension(-2:2) :: mshif    ! local shifts of the layer index for a given point right from point of consideration
    integer, dimension(-2:2) :: nshif    ! local shifts of the layer index for a given point up from point of consideration
    !
    data mshif(-2:2) /18, 8, 3, 12, 22/
    data nshif(-2:2) /19, 9, 7, 13, 23/
!
!   Structure
!
!   Description of the pseudo code
!
!   Source text
!
    if (ltrace) call strace (ient,'ilu')
    !
    amod = pnums(27)
    !
    ! note: nconct-23 is equivalent with kmax-1 and 2 for Keller-box and standard pressure gradient approximations, respectively
    !
    jmax = nconct - 23
    !
    ! initialize some arrays
    !
    t = 0.
    w = 0.
    !
    ! not computed for points ml and nl at subdomain interfaces since, these end points are owned by the neighbouring subdomains
    !
    mend = ml - 1
    if ( LMXL ) mend = ml
    !
    nend = nl - 1
    if ( LMYL ) nend = nl
    !
    if ( oned ) then
       !
       do m = mfu, mend
          !
          nm = kgrpnt(m,1)
          !
          do k = 1, qmax
             !
             do j = 1, nconct
                !
                t(nm,k) = t(nm,k) + amat(nm,k,j)
                !
             enddo
             !
          enddo
          !
       enddo
       !
       do m = mfu, mend
          !
          md = m - 1
          mu = m + 1
          !
          nm  = kgrpnt(m ,1)
          nmd = kgrpnt(md,1)
          nmu = kgrpnt(mu,1)
          !
          do k = 1, qmax
             !
             kd    = max(k-1,1   )
             kdd   = max(k-2,1   )
             kddd  = max(k-3,-1  )
             kdddd = max(k-4,-1  )
             ku    = min(k+1,qmax)
             kuu   = min(k+2,qmax)
             !
             ! compute matrix U
             !
             do j = 0, jmax
                !
                jl = 1.
                if ( j > 2 ) jl = 0.
                !
                jlu = 1.
                if ( j > 1 ) jlu = 0.
                jluu = 1.
                if ( j > 0 ) jluu = 0.
                !
                jld = 1.
                if ( j > 3 ) jld = 0.
                jldd = 1.
                if ( j > 4 ) jldd = 0.
                !
                jlup = 1.
                if ( j == jmax ) jlup = 0.
                jluup = 1.
                if ( j == jmax-1 .or. j == jmax ) jluup = 0.
                !
                prec(nm,k,ishif(j)) = amat(nm,k,ishif(j)) -   jl*prec(nm,k, 2)*prec(nmd,k  ,mshif(min(j  ,2   ))) - jlup*prec(nm,k, 4)*prec(nm ,k-1,ishif(min(j+1,jmax)))  &
                                                          -  jlu*prec(nm,k,10)*prec(nmd,k-1,mshif(min(j+1,2   ))) -  jld*prec(nm,k,14)*prec(nmd,k+1,mshif(min(j-1,2   )))  &
                                                          -jluup*prec(nm,k,16)*prec(nm ,k-2,ishif(min(j+2,jmax))) - jluu*prec(nm,k,20)*prec(nmd,k-2,mshif(min(j+2,2   )))  &
                                                          - jldd*prec(nm,k,24)*prec(nmd,k+2,mshif(min(j-2,2   )))
                !
             enddo
             !
             prec(nm,k, 3) = amat(nm,k, 3) - prec(nm,k,4)*prec(nm,k-1,12) - prec(nm,k,16)*prec(nm,k-2,22)
             prec(nm,k, 8) = amat(nm,k, 8) - prec(nm,k,4)*prec(nm,k-1, 3) - prec(nm,k,16)*prec(nm,k-2,12)
             prec(nm,k,12) = amat(nm,k,12) - prec(nm,k,4)*prec(nm,k-1,22)
             prec(nm,k,18) = amat(nm,k,18) - prec(nm,k,4)*prec(nm,k-1, 8) - prec(nm,k,16)*prec(nm,k-2, 3)
             prec(nm,k,22) = amat(nm,k,22)
             !
             ! compute the row sum of U in case of MILU
             !
             w(nm,k) = prec(nm,k,3) + prec(nm,k,8) + prec(nm,k,12) + prec(nm,k,18) + prec(nm,k,22)
             !
             do j = 1, jmax
                !
                w(nm,k) = w(nm,k) + prec(nm,k,ishif(j))
                !
             enddo
             !
             fac = t(nm,k) - w(nm,k)
             fac = fac - prec(nm,k, 2)*w(nmd,k  ) - prec(nm,k, 4)*w(nm ,k-1) - prec(nm,k,10)*w(nmd,k-1) - prec(nm,k,14)*w(nmd,k+1)  &
                       - prec(nm,k,16)*w(nm ,k-2) - prec(nm,k,20)*w(nmd,k-2) - prec(nm,k,24)*w(nmd,k+2)
             !
             prec(nm,k,1) = (1.-amod) * prec(nm,k,1) + amod * fac
             w(nm,k)      = w(nm,k) + prec(nm,k,1)
             !
             if ( prec(nm,k,1) == 0. ) then
                prec(nm,k,1) = 1.
             else
                prec(nm,k,1) = 1./prec(nm,k,1)
             endif
             !
             ! compute matrix L
             !
             kldddd = 1.
             if ( k < 3 ) kldddd = 0.
             klddd = 1.
             if ( k < 2 ) klddd = 0.
             kldd = 1.
             if ( k == 1 .or. k == 2 ) kldd = 0.
             kld  = 1.
             if ( k == 1 ) kld = 0.
             klu  = 1.
             if ( k == qmax ) klu = 0.
             kluu = 1.
             if ( k == qmax-1 .or. k == qmax ) kluu = 0.
             !
             prec(nmu,k  , 2) = (      amat(nmu,k  , 2) - prec(nmu,k  ,10)*prec(nm ,k-1, 5) - prec(nmu,k  ,20)*prec(nm ,k-2,17) ) * prec(nm,k,1)
             prec(nm ,k+1, 4) = (  klu*amat(nm ,ku , 4) - prec(nm ,k+1, 2)*prec(nmd,k+1, 8) - prec(nm ,k+1,10)*prec(nmd,k  , 3)  &
                                                        - prec(nm ,k+1,14)*prec(nmd,k+2,18) - prec(nm ,k+1,16)*prec(nm ,k-1, 5)  &
                                                        - prec(nm ,k+1,20)*prec(nmd,k-1,12) )*prec(nm,k,1)
             prec(nmu,k-1,14) = (  kld*amat(nmu,kd ,14) - prec(nmu,k-1, 2)*prec(nm ,k-1, 5) - prec(nmu,k-1,10)*prec(nm ,k-2,17) ) * prec(nm,k,1)
             if ( jmax > 2 ) prec(nmu,k-1,14) = prec(nmu,k-1,14) - klddd*prec(nmu,k-1,20)*prec(nm,kddd,ishif(3)) * prec(nm,k,1)
             prec(nmu,k+1,10) = (  klu*amat(nmu,ku ,10) - prec(nmu,k+1,20)*prec(nm ,k-1, 5) )*prec(nm,k,1)
             prec(nm ,k+2,16) = ( kluu*amat(nm ,kuu,16) - prec(nm ,k+2, 2)*prec(nmd,k+2,18) - prec(nm ,k+2,10)*prec(nmd,k+1, 8)  &
                                                        - prec(nm ,k+2,20)*prec(nmd,k  , 3) )*prec(nm,k,1)
             prec(nmu,k-2,24) = ( kldd*amat(nmu,kdd,24) - prec(nmu,k-2,14)*prec(nm ,k-1, 5) - prec(nmu,k-2, 2)*prec(nm ,k-2,17) ) * prec(nm,k,1)
             if ( jmax > 2 ) prec(nmu,k-2,24) = prec(nmu,k-2,24) -  klddd*prec(nmu,k-2,10)*prec(nm,kddd ,ishif(3)) * prec(nm,k,1)
             if ( jmax > 3 ) prec(nmu,k-2,24) = prec(nmu,k-2,24) - kldddd*prec(nmu,k-2,20)*prec(nm,kdddd,ishif(4)) * prec(nm,k,1)
             prec(nmu,k+2,20) =   kluu*amat(nmu,kuu,20) * prec(nm,k,1)
             !
          enddo
          !
       enddo
       !
    else
       !
       do n = nfu, nend
          do m = mfu, mend
             !
             nm = kgrpnt(m,n)
             !
             do k = 1, qmax
                !
                do j = 1, nconct
                   !
                   t(nm,k) = t(nm,k) + amat(nm,k,j)
                   !
                enddo
                !
             enddo
             !
          enddo
       enddo
       !
       do n = nfu, nend
          do m = mfu, mend
             !
             md = m - 1
             mu = m + 1
             nd = n - 1
             nu = n + 1
             !
             nm  = kgrpnt(m ,n )
             nmd = kgrpnt(md,n )
             nmu = kgrpnt(mu,n )
             ndm = kgrpnt(m ,nd)
             num = kgrpnt(m ,nu)
             !
             do k = 1, qmax
                !
                kd    = max(k-1,1   )
                kdd   = max(k-2,1   )
                kddd  = max(k-3,-1  )
                kdddd = max(k-4,-1  )
                ku    = min(k+1,qmax)
                kuu   = min(k+2,qmax)
                !
                ! compute matrix U
                !
                do j = 0, jmax
                   !
                   jl = 1.
                   if ( j > 2 ) jl = 0.
                   !
                   jlu = 1.
                   if ( j > 1 ) jlu = 0.
                   jluu = 1.
                   if ( j > 0 ) jluu = 0.
                   !
                   jld = 1.
                   if ( j > 3 ) jld = 0.
                   jldd = 1.
                   if ( j > 4 ) jldd = 0.
                   !
                   jlup = 1.
                   if ( j == jmax ) jlup = 0.
                   jluup = 1.
                   if ( j == jmax-1 .or. j == jmax ) jluup = 0.
                   !
                   prec(nm,k,ishif(j)) = amat(nm,k,ishif(j)) -   jl*prec(nm,k, 2)*prec(nmd,k  ,mshif(min(j  ,2   ))) - jlup*prec(nm,k, 4)*prec(nm ,k-1,ishif(min(j+1,jmax)))  &
                                                             -  jlu*prec(nm,k,10)*prec(nmd,k-1,mshif(min(j+1,2   ))) -  jld*prec(nm,k,14)*prec(nmd,k+1,mshif(min(j-1,2   )))  &
                                                             -jluup*prec(nm,k,16)*prec(nm ,k-2,ishif(min(j+2,jmax))) - jluu*prec(nm,k,20)*prec(nmd,k-2,mshif(min(j+2,2   )))  &
                                                             - jldd*prec(nm,k,24)*prec(nmd,k+2,mshif(min(j-2,2   ))) -   jl*prec(nm,k, 6)*prec(ndm,k  ,nshif(min(j  ,2   )))  &
                                                             -  jlu*prec(nm,k,11)*prec(ndm,k-1,nshif(min(j+1,2   ))) -  jld*prec(nm,k,15)*prec(ndm,k+1,nshif(min(j-1,2   )))  &
                                                             - jluu*prec(nm,k,21)*prec(ndm,k-2,nshif(min(j+2,2   ))) - jldd*prec(nm,k,25)*prec(ndm,k+2,nshif(min(j-2,2   )))
                   !
                enddo
                !
                prec(nm,k, 3) = amat(nm,k, 3) - prec(nm,k,4)*prec(nm ,k-1,12) - prec(nm,k,16)*prec(nm,k-2,22)
                prec(nm,k, 7) = amat(nm,k, 7) - prec(nm,k,4)*prec(nm ,k-1,13) - prec(nm,k,16)*prec(nm,k-2,23)
                prec(nm,k, 8) = amat(nm,k, 8) - prec(nm,k,4)*prec(nm ,k-1, 3) - prec(nm,k,16)*prec(nm,k-2,12)
                prec(nm,k, 9) = amat(nm,k, 9) - prec(nm,k,4)*prec(nm ,k-1, 7) - prec(nm,k,16)*prec(nm,k-2,13)
                prec(nm,k,12) = amat(nm,k,12) - prec(nm,k,4)*prec(nm ,k-1,22)
                prec(nm,k,13) = amat(nm,k,13) - prec(nm,k,4)*prec(nm ,k-1,23)
                prec(nm,k,18) = amat(nm,k,18) - prec(nm,k,4)*prec(nm ,k-1, 8) - prec(nm,k,16)*prec(nm,k-2, 3)
                prec(nm,k,19) = amat(nm,k,19) - prec(nm,k,4)*prec(nm ,k-1, 9) - prec(nm,k,16)*prec(nm,k-2, 7)
                prec(nm,k,22) = amat(nm,k,22)
                prec(nm,k,23) = amat(nm,k,23)
                !
                ! compute the row sum of U in case of MILU
                !
                w(nm,k) = prec(nm,k, 3) + prec(nm,k, 7) + prec(nm,k, 8) + prec(nm,k, 9) + prec(nm,k,12) +  &
                          prec(nm,k,13) + prec(nm,k,18) + prec(nm,k,19) + prec(nm,k,22) + prec(nm,k,23)
                !
                do j = 1, jmax
                   !
                   w(nm,k) = w(nm,k) + prec(nm,k,ishif(j))
                   !
                enddo
                !
                fac = t(nm,k) - w(nm,k)
                fac = fac - prec(nm,k, 2)*w(nmd,k  ) - prec(nm,k, 4)*w(nm ,k-1) - prec(nm,k, 6)*w(ndm,k  ) - prec(nm,k,10)*w(nmd,k-1)  &
                          - prec(nm,k,11)*w(ndm,k-1) - prec(nm,k,14)*w(nmd,k+1) - prec(nm,k,15)*w(ndm,k+1) - prec(nm,k,16)*w(nm ,k-2)  &
                          - prec(nm,k,20)*w(nmd,k-2) - prec(nm,k,21)*w(ndm,k-2) - prec(nm,k,24)*w(nmd,k+2) - prec(nm,k,25)*w(ndm,k+2)
                !
                prec(nm,k,1) = (1.-amod) * prec(nm,k,1) + amod * fac
                w(nm,k)      = w(nm,k) + prec(nm,k,1)
                !
                if ( prec(nm,k,1) == 0. ) then
                   prec(nm,k,1) = 1.
                else
                   prec(nm,k,1) = 1./prec(nm,k,1)
                endif
                !
                ! compute matrix L
                !
                kldddd = 1.
                if ( k < 3 ) kldddd = 0.
                klddd = 1.
                if ( k < 2 ) klddd = 0.
                kldd = 1.
                if ( k == 1 .or. k == 2 ) kldd = 0.
                kld  = 1.
                if ( k == 1 ) kld = 0.
                klu  = 1.
                if ( k == qmax ) klu = 0.
                kluu = 1.
                if ( k == qmax-1 .or. k == qmax ) kluu = 0.
                !
                prec(nmu,k  , 2) = (      amat(nmu,k  , 2) - prec(nmu,k  ,10)*prec(nm ,k-1, 5) - prec(nmu,k  ,20)*prec(nm ,k-2,17) ) * prec(nm,k,1)
                prec(nm ,k+1, 4) = (  klu*amat(nm ,ku , 4) - prec(nm ,k+1, 2)*prec(nmd,k+1, 8) - prec(nm ,k+1, 6)*prec(ndm,k+1, 9)  &
                                                           - prec(nm ,k+1,10)*prec(nmd,k  , 3) - prec(nm ,k+1,11)*prec(ndm,k  , 7)  &
                                                           - prec(nm ,k+1,14)*prec(nmd,k+2,18) - prec(nm ,k+1,15)*prec(ndm,k+2,19)  &
                                                           - prec(nm ,k+1,16)*prec(nm ,k-1, 5) - prec(nm ,k+1,20)*prec(nmd,k-1,12)  &
                                                           - prec(nm ,k+1,21)*prec(ndm,k-1,13) )*prec(nm,k,1)
                prec(num,k  , 6) = (      amat(num,k  , 6) - prec(num,k  ,11)*prec(nm ,k-1, 5) - prec(num,k  ,21)*prec(nm ,k-2,17) ) * prec(nm,k,1)
                prec(nmu,k-1,14) = (  kld*amat(nmu,kd ,14) - prec(nmu,k-1, 2)*prec(nm ,k-1, 5) - prec(nmu,k-1,10)*prec(nm ,k-2,17) ) * prec(nm,k,1)
                if ( jmax > 2 ) prec(nmu,k-1,14) = prec(nmu,k-1,14) - klddd*prec(nmu,k-1,20)*prec(nm,kddd,ishif(3)) * prec(nm,k,1)
                prec(num,k-1,15) = (  kld*amat(num,kd ,15) - prec(num,k-1, 6)*prec(nm ,k-1, 5) - prec(num,k-1,11)*prec(nm ,k-2,17) ) * prec(nm,k,1)
                if ( jmax > 2 ) prec(num,k-1,15) = prec(num,k-1,15) - klddd*prec(num,k-1,21)*prec(nm,kddd,ishif(3)) * prec(nm,k,1)
                prec(nmu,k+1,10) = (  klu*amat(nmu,ku ,10) - prec(nmu,k+1,20)*prec(nm ,k-1, 5) )*prec(nm,k,1)
                prec(num,k+1,11) = (  klu*amat(num,ku ,11) - prec(num,k+1,21)*prec(nm ,k-1, 5) )*prec(nm,k,1)
                prec(nm ,k+2,16) = ( kluu*amat(nm ,kuu,16) - prec(nm ,k+2, 2)*prec(nmd,k+2,18) - prec(nm ,k+2, 6)*prec(ndm,k+2,19)  &
                                                           - prec(nm ,k+2,10)*prec(nmd,k+1, 8) - prec(nm ,k+2,11)*prec(ndm,k+1, 9)  &
                                                           - prec(nm ,k+2,20)*prec(nmd,k  , 3) - prec(nm ,k+2,21)*prec(ndm,k  , 7) ) * prec(nm,k,1)
                prec(nmu,k-2,24) = ( kldd*amat(nmu,kdd,24) - prec(nmu,k-2,14)*prec(nm ,k-1, 5) - prec(nmu,k-2, 2)*prec(nm ,k-2,17) ) * prec(nm,k,1)
                if ( jmax > 2 ) prec(nmu,k-2,24) = prec(nmu,k-2,24) -  klddd*prec(nmu,k-2,10)*prec(nm,kddd ,ishif(3)) * prec(nm,k,1)
                if ( jmax > 3 ) prec(nmu,k-2,24) = prec(nmu,k-2,24) - kldddd*prec(nmu,k-2,20)*prec(nm,kdddd,ishif(4)) * prec(nm,k,1)
                prec(num,k-2,25) = ( kldd*amat(num,kdd,25) - prec(num,k-2,15)*prec(nm ,k-1, 5) - prec(num,k-2, 6)*prec(nm ,k-2,17) ) * prec(nm,k,1)
                if ( jmax > 2 ) prec(num,k-2,25) = prec(num,k-2,25) -  klddd*prec(num,k-2,11)*prec(nm,kddd ,ishif(3)) * prec(nm,k,1)
                if ( jmax > 3 ) prec(num,k-2,25) = prec(num,k-2,25) - kldddd*prec(num,k-2,21)*prec(nm,kdddd,ishif(4)) * prec(nm,k,1)
                prec(nmu,k+2,20) =   kluu*amat(nmu,kuu,20) * prec(nm,k,1)
                prec(num,k+2,21) =   kluu*amat(num,kuu,21) * prec(nm,k,1)
                !
             enddo
             !
          enddo
       enddo
       !
    endif
    !
end subroutine ilu
!
subroutine iludr ( amat )
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
!    1.00, February 2012: New subroutine
!
!   Purpose
!
!   Computes incomplete factorization restricted to diagonal (RILUD)
!
!   Method
!
!   Based on the scaled system an incomplete decomposition is obtained, using the following rules:
!
!   A = LU - R
!
!   (a)  diag(L) = diag(U) = I,
!   (b)  the off-diagonal parts of L and U are equal to the corresponding parts of A
!   (c)  diag(LU) = diag(A)
!
!   This incomplete decomposition is used as the ILUD right preconditioner. If the last rule (c) is replaced by
!
!        rowsum(LU) = rowsum(A)
!
!   the MILUD preconditioner is obtained. We also used a RILUD preconditioner, which is the average of
!   ILUD and MILUD. The contributions on the diagonal are weighted by the scalar amod: amod = 0 corresponds
!   with ILUD, whereas amod = 1 corresponds with MILUD.
!
!   Modules used
!
    use ocpcomm4
    use SwashCommdata3, only: mcgrd, qmax, oned, pnums
    use SwashFlowdata, only: mfu, ml, nfu, nl, nconct, ishif
    use m_genarr, only: kgrpnt
    use SwashSolvedata
!
    implicit none
!
!   Argument variables
!
    real, dimension(mcgrd,qmax,nconct), intent(in) :: amat ! coefficient matrix of the system of equations
!
!   Local variables
!
    integer, save :: ient = 0 ! number of entries in this subroutine
    integer       :: j        ! loop counter
    integer       :: k        ! loop counter over vertical layers
    integer       :: kd       ! index of layer k-1
    integer       :: kdd      ! index of layer k-2
    integer       :: ku       ! index of layer k+1
    integer       :: kuu      ! index of layer k+2
    integer       :: m        ! loop counter
    integer       :: md       ! index of point m-1
    integer       :: n        ! loop counter
    integer       :: nd       ! index of point n-1
    integer       :: ndm      ! pointer to m,n-1
    integer       :: nm       ! pointer to m,n
    integer       :: nmd      ! pointer to m-1,n
    !
    real          :: amod     ! parameter used in the average of ILUD and MILUD preconditioners
!
!   Structure
!
!   Description of the pseudo code
!
!   Source text
!
    if (ltrace) call strace (ient,'iludr')
    !
    amod = pnums(27)
    !
    ! initialize an auxiliary array
    !
    w = 0.
    !
    if ( oned ) then
       !
       ! copy amat to prec and modify it such that it is an M matrix
       !
       do m = mfu, ml
          !
          nm = kgrpnt(m,1)
          !
          do k = 1, qmax
             !
             prec(nm,k,1) = amat(nm,k,1)
             !
             do j = 2, nconct
                !
                if ( amat(nm,k,j) > 0. ) then
                   prec(nm,k,j) = amat(nm,k,j)
                else
                   prec(nm,k,j) = 0.
                   prec(nm,k,1) = prec(nm,k,1) + amat(nm,k,j)
                endif
                !
             enddo
             !
          enddo
          !
       enddo
       !
       ! compute the row sum of U in case of MILUD
       !
       if ( amod /= 0. ) then
          !
          do m = mfu, ml
             !
             nm = kgrpnt(m,1)
             !
             do k = 1, qmax
                !
                w(nm,k) = prec(nm,k,3) + prec(nm,k,8) + prec(nm,k,12) + prec(nm,k,18) + prec(nm,k,22)
                !
                do j = 1, nconct-23
                   !
                   w(nm,k) = w(nm,k) + prec(nm,k,ishif(j))
                   !
                enddo
                !
             enddo
             !
          enddo
          !
       endif
       !
       do m = mfu, ml
          !
          md = m - 1
          !
          nm  = kgrpnt(m ,1)
          nmd = kgrpnt(md,1)
          !
          do k = 1, qmax
             !
             kd  = max(k-1,1   )
             kdd = max(k-2,1   )
             ku  = min(k+1,qmax)
             kuu = min(k+2,qmax)
             !
             diagb(nm,k) = prec(nm,k,1) - prec(nm,k, 2) * ( (1.-amod) * prec(nmd,k  , 3) + amod * w(nmd,k  ) ) / diagb(nmd,k  )  &
                                        - prec(nm,k, 4) * ( (1.-amod) * prec(nm ,kd , 5) + amod * w(nm ,k-1) ) / diagb(nm ,k-1)  &
                                        - prec(nm,k,10) * ( (1.-amod) * prec(nmd,kd ,12) + amod * w(nmd,k-1) ) / diagb(nmd,k-1)  &
                                        - prec(nm,k,14) * ( (1.-amod) * prec(nmd,ku , 8) + amod * w(nmd,k+1) ) / diagb(nmd,k+1)  &
                                        - prec(nm,k,16) * ( (1.-amod) * prec(nm ,kdd,17) + amod * w(nm ,k-2) ) / diagb(nm ,k-2)  &
                                        - prec(nm,k,20) * ( (1.-amod) * prec(nmd,kdd,22) + amod * w(nmd,k-2) ) / diagb(nmd,k-2)  &
                                        - prec(nm,k,24) * ( (1.-amod) * prec(nmd,kuu,18) + amod * w(nmd,k+2) ) / diagb(nmd,k+2)
             !
          enddo
          !
       enddo
       !
       do m = mfu, ml
          !
          nm = kgrpnt(m,1)
          !
          do k = 1, qmax
             !
             if ( diagb(nm,k) /= 0. ) then
                diagb(nm,k) = 1./diagb(nm,k)
             else
                diagb(nm,k) = 1.
             endif
             !
          enddo
          !
       enddo
       !
       ! scale the ILUD right preconditioner
       !
       do m = mfu, ml
          !
          nm = kgrpnt(m,1)
          !
          do k = 1, qmax
             !
             prec(nm,k,1) = diagb(nm,k)
             !
             do j = 2, nconct
                !
                prec(nm,k,j) = prec(nm,k,j) * diagb(nm,k)
                !
             enddo
             !
          enddo
          !
       enddo
       !
    else
       !
       ! copy amat to prec and modify it such that it is an M matrix
       !
       do n = nfu, nl
          do m = mfu, ml
             !
             nm = kgrpnt(m,n)
             !
             do k = 1, qmax
                !
                prec(nm,k,1) = amat(nm,k,1)
                !
                do j = 2, nconct
                   !
                   if ( amat(nm,k,j) > 0. ) then
                      prec(nm,k,j) = amat(nm,k,j)
                   else
                      prec(nm,k,j) = 0.
                      prec(nm,k,1) = prec(nm,k,1) + amat(nm,k,j)
                   endif
                   !
                enddo
                !
             enddo
          !
          enddo
       enddo
       !
       ! compute the row sum of U in case of MILUD
       !
       if ( amod /= 0. ) then
          !
          do n = nfu, nl
             do m = mfu, ml
                !
                nm = kgrpnt(m,n)
                !
                do k = 1, qmax
                   !
                   w(nm,k) = prec(nm,k, 3) + prec(nm,k, 7) + prec(nm,k, 8) + prec(nm,k, 9) + prec(nm,k,12) +  &
                             prec(nm,k,13) + prec(nm,k,18) + prec(nm,k,19) + prec(nm,k,22) + prec(nm,k,23)
                   !
                   do j = 1, nconct-23
                      !
                      w(nm,k) = w(nm,k) + prec(nm,k,ishif(j))
                      !
                   enddo
                   !
                enddo
                !
             enddo
          enddo
          !
       endif
       !
       ! compute the diagonal elements
       !
       do n = nfu, nl
          do m = mfu, ml
             !
             md = m - 1
             nd = n - 1
             !
             nm  = kgrpnt(m ,n )
             nmd = kgrpnt(md,n )
             ndm = kgrpnt(m ,nd)
             !
             do k = 1, qmax
                !
                kd  = max(k-1,1   )
                kdd = max(k-2,1   )
                ku  = min(k+1,qmax)
                kuu = min(k+2,qmax)
                !
                diagb(nm,k) = prec(nm,k,1) - prec(nm,k, 2) * ( (1.-amod) * prec(nmd,k  , 3) + amod * w(nmd,k  ) ) / diagb(nmd,k  )  &
                                           - prec(nm,k, 4) * ( (1.-amod) * prec(nm ,kd , 5) + amod * w(nm ,k-1) ) / diagb(nm ,k-1)  &
                                           - prec(nm,k, 6) * ( (1.-amod) * prec(ndm,k  , 7) + amod * w(ndm,k  ) ) / diagb(ndm,k  )  &
                                           - prec(nm,k,10) * ( (1.-amod) * prec(nmd,kd ,12) + amod * w(nmd,k-1) ) / diagb(nmd,k-1)  &
                                           - prec(nm,k,11) * ( (1.-amod) * prec(ndm,kd ,13) + amod * w(ndm,k-1) ) / diagb(ndm,k-1)  &
                                           - prec(nm,k,14) * ( (1.-amod) * prec(nmd,ku , 8) + amod * w(nmd,k+1) ) / diagb(nmd,k+1)  &
                                           - prec(nm,k,15) * ( (1.-amod) * prec(ndm,ku , 9) + amod * w(ndm,k+1) ) / diagb(ndm,k+1)  &
                                           - prec(nm,k,16) * ( (1.-amod) * prec(nm ,kdd,17) + amod * w(nm ,k-2) ) / diagb(nm ,k-2)  &
                                           - prec(nm,k,20) * ( (1.-amod) * prec(nmd,kdd,22) + amod * w(nmd,k-2) ) / diagb(nmd,k-2)  &
                                           - prec(nm,k,21) * ( (1.-amod) * prec(ndm,kdd,23) + amod * w(ndm,k-2) ) / diagb(ndm,k-2)  &
                                           - prec(nm,k,24) * ( (1.-amod) * prec(nmd,kuu,18) + amod * w(nmd,k+2) ) / diagb(nmd,k+2)  &
                                           - prec(nm,k,25) * ( (1.-amod) * prec(ndm,kuu,19) + amod * w(ndm,k+2) ) / diagb(ndm,k+2)
                !
             enddo
             !
          enddo
       enddo
       !
       do n = nfu, nl
          do m = mfu, ml
             !
             nm = kgrpnt(m,n)
             !
             do k = 1, qmax
                !
                if ( diagb(nm,k) /= 0. ) then
                   diagb(nm,k) = 1./diagb(nm,k)
                else
                   diagb(nm,k) = 1.
                endif
                !
             enddo
             !
          enddo
       enddo
       !
       ! scale the ILUD right preconditioner
       !
       do n = nfu, nl
          do m = mfu, ml
             !
             nm = kgrpnt(m,n)
             !
             do k = 1, qmax
                !
                prec(nm,k,1) = diagb(nm,k)
                !
                do j = 2, nconct
                   !
                   prec(nm,k,j) = prec(nm,k,j) * diagb(nm,k)
                   !
                enddo
                !
             enddo
             !
          enddo
       enddo
       !
    endif
    !
end subroutine iludr
!
subroutine iluds ( amat )
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
!    1.00, February 2012: New subroutine
!
!   Purpose
!
!   Computes incomplete factorization restricted to diagonal (RILUD)
!
!   Method
!
!   Based on the scaled system an incomplete decomposition is obtained, using the following rules:
!
!   A = LU - R
!
!   (a)  diag(L) = diag(U) = I,
!   (b)  the off-diagonal parts of L and U are equal to the corresponding parts of A
!   (c)  diag(LU) = diag(A)
!
!   This incomplete decomposition is used as the ILUD split preconditioner. If the last rule (c) is replaced by
!
!        rowsum(LU) = rowsum(A)
!
!   the MILUD preconditioner is obtained. We also used a RILUD preconditioner, which is the average of
!   ILUD and MILUD. The contributions on the diagonal are weighted by the scalar amod: amod = 0 corresponds
!   with ILUD, whereas amod = 1 corresponds with MILUD.
!
!   Modules used
!
    use ocpcomm4
    use SwashCommdata3, only: mcgrd, qmax, oned, pnums
    use SwashFlowdata, only: mfu, ml, nfu, nl, nconct, ishif
    use m_genarr, only: kgrpnt
    use SwashSolvedata
!
    implicit none
!
!   Argument variables
!
    real, dimension(mcgrd,qmax,nconct), intent(in) :: amat ! coefficient matrix of the system of equations
!
!   Local variables
!
    integer, save :: ient = 0 ! number of entries in this subroutine
    integer       :: j        ! loop counter
    integer       :: k        ! loop counter over vertical layers
    integer       :: kd       ! index of layer k-1
    integer       :: kdd      ! index of layer k-2
    integer       :: ku       ! index of layer k+1
    integer       :: kuu      ! index of layer k+2
    integer       :: m        ! loop counter
    integer       :: md       ! index of point m-1
    integer       :: n        ! loop counter
    integer       :: nd       ! index of point n-1
    integer       :: ndm      ! pointer to m,n-1
    integer       :: nm       ! pointer to m,n
    integer       :: nmd      ! pointer to m-1,n
    !
    real          :: amod     ! parameter used in the average of ILUD and MILUD preconditioners
!
!   Structure
!
!   Description of the pseudo code
!
!   Source text
!
    if (ltrace) call strace (ient,'iluds')
    !
    amod = pnums(27)
    !
    ! initialize an auxiliary array
    !
    w = 0.
    !
    if ( oned ) then
       !
       ! compute the row sum of U in case of MILUD
       !
       if ( amod /= 0. ) then
          !
          do m = mfu, ml
             !
             nm = kgrpnt(m,1)
             !
             do k = 1, qmax
                !
                w(nm,k) = amat(nm,k,3) + amat(nm,k,8) + amat(nm,k,12) + amat(nm,k,18) + amat(nm,k,22)
                !
                do j = 1, nconct-23
                   !
                   w(nm,k) = w(nm,k) + amat(nm,k,ishif(j))
                   !
                enddo
                !
             enddo
             !
          enddo
          !
       endif
       !
       do m = mfu, ml
          !
          md = m - 1
          !
          nm  = kgrpnt(m ,1)
          nmd = kgrpnt(md,1)
          !
          do k = 1, qmax
             !
             kd  = max(k-1,1   )
             kdd = max(k-2,1   )
             ku  = min(k+1,qmax)
             kuu = min(k+2,qmax)
             !
             diagb(nm,k) = amat(nm,k,1) - amat(nm,k, 2) * ( (1.-amod) * amat(nmd,k  , 3) + amod * w(nmd,k  ) ) / diagb(nmd,k  )  &
                                        - amat(nm,k, 4) * ( (1.-amod) * amat(nm ,kd , 5) + amod * w(nm ,k-1) ) / diagb(nm ,k-1)  &
                                        - amat(nm,k,10) * ( (1.-amod) * amat(nmd,kd ,12) + amod * w(nmd,k-1) ) / diagb(nmd,k-1)  &
                                        - amat(nm,k,14) * ( (1.-amod) * amat(nmd,ku , 8) + amod * w(nmd,k+1) ) / diagb(nmd,k+1)  &
                                        - amat(nm,k,16) * ( (1.-amod) * amat(nm ,kdd,17) + amod * w(nm ,k-2) ) / diagb(nm ,k-2)  &
                                        - amat(nm,k,20) * ( (1.-amod) * amat(nmd,kdd,22) + amod * w(nmd,k-2) ) / diagb(nmd,k-2)  &
                                        - amat(nm,k,24) * ( (1.-amod) * amat(nmd,kuu,18) + amod * w(nmd,k+2) ) / diagb(nmd,k+2)
             !
          enddo
          !
       enddo
       !
       do m = mfu, ml
          !
          nm = kgrpnt(m,1)
          !
          do k = 1, qmax
             !
             if ( diagb(nm,k) == 0. ) then
                diagb(nm,k) = 1.
             else
                diagb(nm,k) = 1./diagb(nm,k)
             endif
             !
          enddo
          !
       enddo
       !
    else
       !
       ! compute the row sum of U in case of MILUD
       !
       if ( amod /= 0. ) then
          !
          do n = nfu, nl
             do m = mfu, ml
                !
                nm = kgrpnt(m,n)
                !
                do k = 1, qmax
                   !
                   w(nm,k) = amat(nm,k, 3) + amat(nm,k, 7) + amat(nm,k, 8) + amat(nm,k, 9) + amat(nm,k,12) +  &
                             amat(nm,k,13) + amat(nm,k,18) + amat(nm,k,19) + amat(nm,k,22) + amat(nm,k,23)
                   !
                   do j = 1, nconct-23
                      !
                      w(nm,k) = w(nm,k) + amat(nm,k,ishif(j))
                      !
                   enddo
                   !
                enddo
                !
             enddo
          enddo
          !
       endif
       !
       ! compute the diagonal elements
       !
       do n = nfu, nl
          do m = mfu, ml
             !
             md = m - 1
             nd = n - 1
             !
             nm  = kgrpnt(m ,n )
             nmd = kgrpnt(md,n )
             ndm = kgrpnt(m ,nd)
             !
             do k = 1, qmax
                !
                kd  = max(k-1,1   )
                kdd = max(k-2,1   )
                ku  = min(k+1,qmax)
                kuu = min(k+2,qmax)
                !
                diagb(nm,k) = amat(nm,k,1) - amat(nm,k, 2) * ( (1.-amod) * amat(nmd,k  , 3) + amod * w(nmd,k  ) ) / diagb(nmd,k  )  &
                                           - amat(nm,k, 4) * ( (1.-amod) * amat(nm ,kd , 5) + amod * w(nm ,k-1) ) / diagb(nm ,k-1)  &
                                           - amat(nm,k, 6) * ( (1.-amod) * amat(ndm,k  , 7) + amod * w(ndm,k  ) ) / diagb(ndm,k  )  &
                                           - amat(nm,k,10) * ( (1.-amod) * amat(nmd,kd ,12) + amod * w(nmd,k-1) ) / diagb(nmd,k-1)  &
                                           - amat(nm,k,11) * ( (1.-amod) * amat(ndm,kd ,13) + amod * w(ndm,k-1) ) / diagb(ndm,k-1)  &
                                           - amat(nm,k,14) * ( (1.-amod) * amat(nmd,ku , 8) + amod * w(nmd,k+1) ) / diagb(nmd,k+1)  &
                                           - amat(nm,k,15) * ( (1.-amod) * amat(ndm,ku , 9) + amod * w(ndm,k+1) ) / diagb(ndm,k+1)  &
                                           - amat(nm,k,16) * ( (1.-amod) * amat(nm ,kdd,17) + amod * w(nm ,k-2) ) / diagb(nm ,k-2)  &
                                           - amat(nm,k,20) * ( (1.-amod) * amat(nmd,kdd,22) + amod * w(nmd,k-2) ) / diagb(nmd,k-2)  &
                                           - amat(nm,k,21) * ( (1.-amod) * amat(ndm,kdd,23) + amod * w(ndm,k-2) ) / diagb(ndm,k-2)  &
                                           - amat(nm,k,24) * ( (1.-amod) * amat(nmd,kuu,18) + amod * w(nmd,k+2) ) / diagb(nmd,k+2)  &
                                           - amat(nm,k,25) * ( (1.-amod) * amat(ndm,kuu,19) + amod * w(ndm,k+2) ) / diagb(ndm,k+2)
                !
             enddo
             !
          enddo
       enddo
       !
       do n = nfu, nl
          do m = mfu, ml
             !
             nm = kgrpnt(m,n)
             !
             do k = 1, qmax
                !
                if ( diagb(nm,k) == 0. ) then
                   diagb(nm,k) = 1.
                else
                   diagb(nm,k) = 1./diagb(nm,k)
                endif
                !
             enddo
             !
          enddo
       enddo
       !
    endif
    !
end subroutine iluds
!
subroutine bicgstab ( amat, rhs, x )
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
!    1.00, July 2010: New subroutine
!
!   Purpose
!
!   Solves preconditioned system of equations by means of BiCGSTAB method
!
!   Method
!
!   Using the incomplete LU decomposition, the system Ax = b is split preconditioned as follows:
!
!    -1  -1     -1
!   L  AU  y = L  b    with  y = Ux
!
!   so that the Eisenstat trick is possible for the cheap matrix-vector multiplication
!
!   or the system Ax = b is right preconditioned as follows:
!
!     -1 -1
!   AU  L  y = b    with  y = LUx
!
!   The well-known BiCGSTAB method is applied for the solution of preconditioned Ax = b
!
!   H.A. van der Vorst
!   Bi-CGSTAB: a fast and smoothly converging variant of Bi-CG for the solution of nonsymmetric linear systems
!   SIAM J. Sci. Stat. Comput., vol. 13, 631-644, 1992
!
!   Modules used
!
    use ocpcomm4
    use SwashCommdata3, only: mcgrd, qmax, oned, pnums
    use SwashFlowdata, only: mfu, ml, nfu, nl, nconct, ishif
    use m_genarr, only: kgrpnt
    use m_parall
    use SwashSolvedata
!
    implicit none
!
!   Argument variables
!
    real, dimension(mcgrd,qmax,nconct), intent(inout) :: amat ! coefficient matrix of the system of equations
    real, dimension(mcgrd,qmax       ), intent(inout) :: rhs  ! right-hand side of the system of equations
    real, dimension(mcgrd,qmax       ), intent(inout) :: x    ! solution of the system of equations
!
!   Parameter variables
!
    double precision, parameter :: epsmac = 1d-18 ! machine precision number
!
!   Local variables
!
    integer            :: i        ! loop counter
    integer, save      :: ient = 0 ! number of entries in this subroutine
    integer            :: j        ! iteration counter
    integer            :: k        ! loop counter over vertical layers
    integer            :: kd       ! index of layer k-1
    integer            :: kdd      ! index of layer k-2
    integer            :: ku       ! index of layer k+1
    integer            :: kuu      ! index of layer k+2
    integer            :: m        ! loop counter
    integer            :: maxit    ! maximum number of iterations
    integer            :: md       ! index of point m-1
    integer            :: mend     ! end index of loop over wl-points in x-direction
    integer            :: mu       ! index of point m+1
    integer            :: n        ! loop counter
    integer            :: nd       ! index of point n-1
    integer            :: ndm      ! pointer to m,n-1
    integer            :: nend     ! end index of loop over wl-points in y-direction
    integer            :: nm       ! pointer to m,n
    integer            :: nmd      ! pointer to m-1,n
    integer            :: nmu      ! pointer to m+1,n
    integer            :: nu       ! index of point n+1
    integer            :: num      ! pointer to m,n+1
    !
    real               :: alpha    ! factor in the BiCGSTAB algorithm, i.e. rho/sigma or (res0,r)/(res0,Ap)
    real               :: beta     ! factor in the BiCGSTAB algorithm, i.e. alpha/gamma
    real               :: bnorm    ! 2-norm of right-hand side vector
    real               :: crel     ! relaxation parameter
    real               :: epslin   ! required accuracy in the linear solver
    real               :: gamma    ! factor in the BiCGSTAB algorithm
    real, dimension(2) :: reps     ! accuracies with respect to the right-hand side and initial residual used in the following termination criterion:
                                   !
                                   !  ||b-Ax || < reps(1)*||b|| + reps(2)*||b-Ax ||
                                   !        j                                   0
                                   !
    real               :: rho      ! inner product of quasi-residual vector and residual vector
    real               :: rnorm    ! 2-norm of residual vector
    real               :: rnrm0    ! 2-norm of initial residual vector
    real, dimension(2) :: rtmp     ! temporary array for communication
    real               :: sigma    ! inner product of quasi-residual vector and Ap
    real               :: tnorm    ! 2-norm of auxiliary vector t
    real               :: unorm    ! 2-norm of auxiliary vector u
    real               :: ueps     ! minimal accuracy based on machine precision
    !
    logical            :: STPNOW   ! indicates that program must stop
!
!   Structure
!
!   Description of the pseudo code
!
!   Source text
!
    if (ltrace) call strace (ient,'bicgstab')
    !
    reps(1) = pnums(22)
    reps(2) = pnums(23)
    maxit   = nint(pnums(25))
    !
    ! not computed for points ml and nl at subdomain interfaces since, these end points are owned by the neighbouring subdomains
    !
    mend = ml - 1
    if ( LMXL ) mend = ml
    !
    nend = nl - 1
    if ( LMYL ) nend = nl
    !
    if ( icond == 1 .or. icond == 2 ) then
       !
       ! scale the system with the inverse of the diagonal in case of ILUD
       !
       if ( oned ) then
          !
          do m = mfu, ml
             !
             nm = kgrpnt(m,1)
             !
             do k = 1, qmax
                !
                amat(nm,k,:) = amat(nm,k,:) * diagb(nm,k)
                rhs (nm,k  ) = rhs (nm,k  ) * diagb(nm,k)
                !
             enddo
             !
          enddo
          !
       else
          !
          do n = nfu, nl
             do m = mfu, ml
                !
                nm = kgrpnt(m,n)
                !
                do k = 1, qmax
                   !
                   amat(nm,k,:) = amat(nm,k,:) * diagb(nm,k)
                   rhs (nm,k  ) = rhs (nm,k  ) * diagb(nm,k)
                   !
                enddo
                !
             enddo
          enddo
          !
       endif
       !
    endif
    !
    if ( icond == 1 ) then
       !
       call SWEXCHG ( rhs, kgrpnt, 1, qmax )
       if (STPNOW()) return
       !
       ! left preconditioning of the right-hand side
       !
       if ( oned ) then
          !
          do m = mfu, mend
             !
             md = m - 1
             !
             nm  = kgrpnt(m ,1)
             nmd = kgrpnt(md,1)
             !
             do k = 1, qmax
                !
                kd  = max(k-1,1   )
                kdd = max(k-2,1   )
                ku  = min(k+1,qmax)
                kuu = min(k+2,qmax)
                !
                rhs(nm,k) = rhs(nm,k) - amat(nm,k, 2)*rhs(nmd,k  ) - amat(nm,k, 4)*rhs(nm ,kd ) - amat(nm,k,10)*rhs(nmd,kd ) - amat(nm,k,14)*rhs(nmd,ku)  &
                                      - amat(nm,k,16)*rhs(nm ,kdd) - amat(nm,k,20)*rhs(nmd,kdd) - amat(nm,k,24)*rhs(nmd,kuu)
                !
             enddo
             !
          enddo
          !
       else
          !
          do n = nfu, nend
             do m = mfu, mend
                !
                md = m - 1
                nd = n - 1
                !
                nm  = kgrpnt(m ,n )
                nmd = kgrpnt(md,n )
                ndm = kgrpnt(m ,nd)
                !
                do k = 1, qmax
                   !
                   kd  = max(k-1,1   )
                   kdd = max(k-2,1   )
                   ku  = min(k+1,qmax)
                   kuu = min(k+2,qmax)
                   !
                   rhs(nm,k) = rhs(nm,k) - amat(nm,k, 2)*rhs(nmd,k  ) - amat(nm,k, 4)*rhs(nm ,kd ) - amat(nm,k, 6)*rhs(ndm,k  ) - amat(nm,k,10)*rhs(nmd,kd )  &
                                         - amat(nm,k,11)*rhs(ndm,kd ) - amat(nm,k,14)*rhs(nmd,ku ) - amat(nm,k,15)*rhs(ndm,ku ) - amat(nm,k,16)*rhs(nm ,kdd)  &
                                         - amat(nm,k,20)*rhs(nmd,kdd) - amat(nm,k,21)*rhs(ndm,kdd) - amat(nm,k,24)*rhs(nmd,kuu) - amat(nm,k,25)*rhs(ndm,kuu)
                   !
                enddo
                !
             enddo
          enddo
          !
       endif
       !
       ! right preconditioning of the solution
       !
       call SWEXCHG ( x, kgrpnt, 1, qmax )
       if (STPNOW()) return
       !
       if ( oned ) then
          !
          do m = mfu, mend
             !
             mu = m + 1
             !
             nm  = kgrpnt(m ,1)
             nmu = kgrpnt(mu,1)
             !
             do k = 1, qmax
                !
                kd  = max(k-1,1   )
                kdd = max(k-2,1   )
                ku  = min(k+1,qmax)
                kuu = min(k+2,qmax)
                !
                x(nm,k) = x(nm,k) + amat(nm,k,3)*x(nmu,k) + amat(nm,k,8)*x(nmu,kd) + amat(nm,k,12)*x(nmu,ku) + amat(nm,k,18)*x(nmu,kdd) + amat(nm,k,22)*x(nmu,kuu)
                !
                ! note: nconct-23 is equivalent with kmax-1 and 2 for Keller-box and standard pressure gradient approximations, respectively
                !
                do i = 1, nconct-23
                   !
                   x(nm,k) = x(nm,k) + amat(nm,k,ishif(i))*x(nm,min(k+i,qmax))
                   !
                enddo
                !
             enddo
             !
          enddo
          !
       else
          !
          do n = nfu, nend
             do m = mfu, mend
                !
                mu = m + 1
                nu = n + 1
                !
                nm  = kgrpnt(m ,n )
                nmu = kgrpnt(mu,n )
                num = kgrpnt(m ,nu)
                !
                do k = 1, qmax
                   !
                   kd  = max(k-1,1   )
                   kdd = max(k-2,1   )
                   ku  = min(k+1,qmax)
                   kuu = min(k+2,qmax)
                   !
                   x(nm,k) = x(nm,k) + amat(nm,k, 3)*x(nmu,k  ) + amat(nm,k, 7)*x(num,k  ) + amat(nm,k, 8)*x(nmu,kd )  &
                                     + amat(nm,k, 9)*x(num,kd ) + amat(nm,k,12)*x(nmu,ku ) + amat(nm,k,13)*x(num,ku )  &
                                     + amat(nm,k,18)*x(nmu,kdd) + amat(nm,k,19)*x(num,kdd) + amat(nm,k,22)*x(nmu,kuu)  &
                                     + amat(nm,k,23)*x(num,kuu)
                   !
                   ! note: nconct-23 is equivalent with kmax-1 and 2 for Keller-box and standard pressure gradient approximations, respectively
                   !
                   do i = 1, nconct-23
                      !
                      x(nm,k) = x(nm,k) + amat(nm,k,ishif(i))*x(nm,min(k+i,qmax))
                      !
                   enddo
                   !
                enddo
                !
             enddo
          enddo
          !
       endif
       !
    endif
    !
    j     = 0
    rho   = 1.
    alpha = 1.
    gamma = 1.
    crel  = 0.7
    !
    ! initialize some arrays
    !
    p = 0.
    v = 0.
    w = 0.
    !
    if ( oned ) then
       !
       ! compute initial residual vector and its 2-norm
       !
       if ( icond == 1 ) then
          call ave ( amat, x, res, vt2 )
       else if ( icond == 2 .or. icond == 3 ) then
          call avm ( amat, x, res )
       endif
       !
       rnorm = 0.
       bnorm = 0.
       !
       do m = mfu, mend
          !
          nm = kgrpnt(m,1)
          !
          do k = 1, qmax
             !
             res (nm,k) = rhs(nm,k) - res(nm,k)
             res0(nm,k) = res(nm,k)
             rnorm = rnorm + res(nm,k)*res(nm,k)
             bnorm = bnorm + rhs(nm,k)*rhs(nm,k)
             !
          enddo
          !
       enddo
       !
       !
       rtmp(1) = rnorm
       rtmp(2) = bnorm
       call SWREDUCE ( rtmp, 2, SWREAL, SWSUM )
       rnorm = sqrt(rtmp(1))
       bnorm = sqrt(rtmp(2))
       rnrm0 = rnorm
       !
       if ( iamout == 2 .and. INODE == MASTER ) then
          write (PRINTF, '(a,i3,a,e12.6)') ' ++ bicgstab: iter = ',j,'    res = ',rnorm
       endif
       !
       epslin = reps(1)*bnorm + reps(2)*rnrm0
       ueps   = 1000.*real(epsmac)*bnorm
       !
       if ( epslin < ueps .and. rnorm > 0. ) then
          !
          if ( iamout > 0 .and. INODE == MASTER ) then
             write (PRINTF, '(a)') ' ++ bicgstab: the required accuracy is too small'
             write (PRINTF, '(a,e12.6)') '              required accuracy    = ',epslin
             write (PRINTF, '(a,e12.6)') '              appropriate accuracy = ',ueps
          endif
          !
          epslin = ueps
          !
       endif
       !
 10    if ( rnorm > epslin .and. j < maxit ) then
          !
          j = j + 1
          !
          beta = alpha / (rho*gamma)
          !
          rho = 0.
          !
          do m = mfu, mend
             !
             nm = kgrpnt(m,1)
             !
             do k = 1, qmax
                !
                rho = rho + res0(nm,k)*res(nm,k)
                !
             enddo
             !
          enddo
          !
          call SWREDUCE ( rho, 1, SWREAL, SWSUM )
          !
          beta = beta * rho
          !
          do m = mfu, mend
             !
             nm = kgrpnt(m,1)
             !
             do k = 1, qmax
                !
                p(nm,k) = res(nm,k) + beta * ( p(nm,k) - gamma*v(nm,k) )
                w(nm,k) = p(nm,k)
                !
             enddo
             !
          enddo
          !
          if ( icond == 1 ) then
             !
             call SWEXCHG ( p, kgrpnt, 1, qmax )
             if (STPNOW()) return
             !
             call ave ( amat, p, v, vt2 )
             !
          else if ( icond == 2 .or. icond == 3 ) then
             !
             call SWEXCHG ( w, kgrpnt, -1, qmax+nconct-23 )
             if (STPNOW()) return
             !
             call ivl ( prec, w )
             !
             call SWEXCHG ( w, kgrpnt, -1, qmax+nconct-23 )
             if (STPNOW()) return
             !
             call ivu ( prec, w )
             !
             call SWEXCHG ( w, kgrpnt, -1, qmax+nconct-23 )
             if (STPNOW()) return
             !
             call avm ( amat, w(1,1), v )
             !
          endif
          !
          sigma = 0.
          !
          do m = mfu, mend
             !
             nm = kgrpnt(m,1)
             !
             do k = 1, qmax
                !
                sigma = sigma + res0(nm,k)*v(nm,k)
                !
             enddo
             !
          enddo
          !
          call SWREDUCE ( sigma, 1, SWREAL, SWSUM )
          !
          if ( .not. sigma /= 0. ) then
             !
             if ( iamout > 0 .and. INODE == MASTER ) then
                write (PRINTF,'(a)') ' ++ bicgstab: process is halted due to a division by zero'
                write (PRINTF,'(a)') '              cause: parameter sigma = 0'
             endif
             !
             goto 20
             !
          endif
          !
          alpha = rho / sigma
          !
          unorm = 0.
          !
          do m = mfu, mend
             !
             nm = kgrpnt(m,1)
             !
             do k = 1, qmax
                !
                u(nm,k) = res(nm,k) - alpha * v(nm,k)
                x(nm,k) = x  (nm,k) + alpha * w(nm,k)
                w(nm,k) = u  (nm,k)
                unorm = unorm + u(nm,k)*u(nm,k)
                !
             enddo
             !
          enddo
          !
          call SWREDUCE ( unorm, 1, SWREAL, SWSUM )
          unorm = sqrt(unorm)
          if ( unorm < epslin ) goto 20
          !
          if ( icond == 1 ) then
             !
             call SWEXCHG ( u, kgrpnt, 1, qmax )
             if (STPNOW()) return
             !
             call ave ( amat, u, t, vt2 )
             !
          else if ( icond == 2 .or. icond == 3 ) then
             !
             call SWEXCHG ( w, kgrpnt, -1, qmax+nconct-23 )
             if (STPNOW()) return
             !
             call ivl ( prec, w )
             !
             call SWEXCHG ( w, kgrpnt, -1, qmax+nconct-23 )
             if (STPNOW()) return
             !
             call ivu ( prec, w )
             !
             call SWEXCHG ( w, kgrpnt, -1, qmax+nconct-23 )
             if (STPNOW()) return
             !
             call avm ( amat, w(1,1), t )
             !
          endif
          !
          if ( unorm /= 0. ) then
             !
             gamma = 0.
             tnorm = 0.
             !
             do m = mfu, mend
                !
                nm = kgrpnt(m,1)
                !
                do k = 1, qmax
                   !
                   gamma = gamma + t(nm,k)*u(nm,k)
                   tnorm = tnorm + t(nm,k)*t(nm,k)
                   !
                enddo
                !
             enddo
             !
             rtmp(1) = gamma
             rtmp(2) = tnorm
             call SWREDUCE ( rtmp, 2, SWREAL, SWSUM )
             gamma = rtmp(1)
             tnorm = rtmp(2)
             !
             tnorm = sqrt(tnorm)
             gamma = gamma / (unorm*tnorm)
             gamma = sign(1.,gamma) * max(abs(gamma),crel) * unorm / tnorm
             !
          else
             !
             gamma = 1.
             !
          endif
          !
          if ( .not. gamma /= 0. ) then
             !
             if ( iamout > 0 .and. INODE == MASTER ) then
                write (PRINTF,'(a)') ' ++ bicgstab: process is halted due to a division by zero'
                write (PRINTF,'(a)') '              cause: parameter gamma = 0'
             endif
             !
             goto 20
             !
          endif
          !
          rnorm = 0.
          !
          do m = mfu, mend
             !
             nm = kgrpnt(m,1)
             !
             do k = 1, qmax
                !
                x  (nm,k) = x(nm,k) + gamma * w(nm,k)
                res(nm,k) = u(nm,k) - gamma * t(nm,k)
                rnorm = rnorm + res(nm,k)*res(nm,k)
                !
             enddo
             !
          enddo
          !
          call SWREDUCE ( rnorm, 1, SWREAL, SWSUM )
          rnorm = sqrt(rnorm)
          !
          if ( iamout == 2 .and. INODE == MASTER ) then
             write (PRINTF,'(a,i3,a,e12.6)') ' ++ bicgstab: iter = ',j,'    res = ',rnorm
          endif
          !
          goto 10
          !
       endif
       !
 20    continue
       !
    else
       !
       ! compute initial residual vector and its 2-norm
       !
       if ( icond == 1 ) then
          call ave ( amat, x, res, vt2 )
       else if ( icond == 2 .or. icond == 3 ) then
          call avm ( amat, x, res )
       endif
       !
       rnorm = 0.
       bnorm = 0.
       !
       do n = nfu, nend
          do m = mfu, mend
             !
             nm = kgrpnt(m,n)
             !
             do k = 1, qmax
                !
                res (nm,k) = rhs(nm,k) - res(nm,k)
                res0(nm,k) = res(nm,k)
                rnorm = rnorm + res(nm,k)*res(nm,k)
                bnorm = bnorm + rhs(nm,k)*rhs(nm,k)
                !
             enddo
             !
          enddo
       enddo
       !
       rtmp(1) = rnorm
       rtmp(2) = bnorm
       call SWREDUCE ( rtmp, 2, SWREAL, SWSUM )
       rnorm = sqrt(rtmp(1))
       bnorm = sqrt(rtmp(2))
       rnrm0 = rnorm
       !
       if ( iamout == 2 .and. INODE == MASTER ) then
          write (PRINTF, '(a,i3,a,e12.6)') ' ++ bicgstab: iter = ',j,'    res = ',rnorm
       endif
       !
       epslin = reps(1)*bnorm + reps(2)*rnrm0
       ueps   = 1000.*real(epsmac)*bnorm
       !
       if ( epslin < ueps .and. rnorm > 0. ) then
          !
          if ( iamout > 0 .and. INODE == MASTER ) then
             write (PRINTF, '(a)') ' ++ bicgstab: the required accuracy is too small'
             write (PRINTF, '(a,e12.6)') '              required accuracy    = ',epslin
             write (PRINTF, '(a,e12.6)') '              appropriate accuracy = ',ueps
          endif
          !
          epslin = ueps
          !
       endif
       !
 30    if ( rnorm > epslin .and. j < maxit ) then
          !
          j = j + 1
          !
          beta = alpha / (rho*gamma)
          !
          rho = 0.
          !
          do n = nfu, nend
             do m = mfu, mend
                !
                nm = kgrpnt(m,n)
                !
                do k = 1, qmax
                   !
                   rho = rho + res0(nm,k)*res(nm,k)
                   !
                enddo
                !
             enddo
          enddo
          !
          call SWREDUCE ( rho, 1, SWREAL, SWSUM )
          !
          beta = beta * rho
          !
          do n = nfu, nend
             do m = mfu, mend
                !
                nm = kgrpnt(m,n)
                !
                do k = 1, qmax
                   !
                   p(nm,k) = res(nm,k) + beta * ( p(nm,k) - gamma*v(nm,k) )
                   w(nm,k) = p(nm,k)
                   !
                enddo
                !
             enddo
          enddo
          !
          if ( icond == 1 ) then
             !
             call SWEXCHG ( p, kgrpnt, 1, qmax )
             if (STPNOW()) return
             !
             call ave ( amat, p, v, vt2 )
             !
          else if ( icond == 2 .or. icond == 3 ) then
             !
             call SWEXCHG ( w, kgrpnt, -1, qmax+nconct-23 )
             if (STPNOW()) return
             !
             call ivl ( prec, w )
             !
             call SWEXCHG ( w, kgrpnt, -1, qmax+nconct-23 )
             if (STPNOW()) return
             !
             call ivu ( prec, w )
             !
             call SWEXCHG ( w, kgrpnt, -1, qmax+nconct-23 )
             if (STPNOW()) return
             !
             call avm ( amat, w(1,1), v )
             !
          endif
          !
          sigma = 0.
          !
          do n = nfu, nend
             do m = mfu, mend
                !
                nm = kgrpnt(m,n)
                !
                do k = 1, qmax
                   !
                   sigma = sigma + res0(nm,k)*v(nm,k)
                   !
                enddo
                !
             enddo
          enddo
          !
          call SWREDUCE ( sigma, 1, SWREAL, SWSUM )
          !
          if ( .not. sigma /= 0. ) then
             !
             if ( iamout > 0 .and. INODE == MASTER ) then
                write (PRINTF,'(a)') ' ++ bicgstab: process is halted due to a division by zero'
                write (PRINTF,'(a)') '              cause: parameter sigma = 0'
             endif
             !
             goto 40
             !
          endif
          !
          alpha = rho / sigma
          !
          unorm = 0.
          !
          do n = nfu, nend
             do m = mfu, mend
                !
                nm = kgrpnt(m,n)
                !
                do k = 1, qmax
                   !
                   u(nm,k) = res(nm,k) - alpha * v(nm,k)
                   x(nm,k) = x  (nm,k) + alpha * w(nm,k)
                   w(nm,k) = u  (nm,k)
                   unorm = unorm + u(nm,k)*u(nm,k)
                   !
                enddo
                !
             enddo
          enddo
          !
          call SWREDUCE ( unorm, 1, SWREAL, SWSUM )
          unorm = sqrt(unorm)
          if ( unorm < epslin ) goto 40
          !
          if ( icond == 1 ) then
             !
             call SWEXCHG ( u, kgrpnt, 1, qmax )
             if (STPNOW()) return
             !
             call ave ( amat, u, t, vt2 )
             !
          else if ( icond == 2 .or. icond == 3 ) then
             !
             call SWEXCHG ( w, kgrpnt, -1, qmax+nconct-23 )
             if (STPNOW()) return
             !
             call ivl ( prec, w )
             !
             call SWEXCHG ( w, kgrpnt, -1, qmax+nconct-23 )
             if (STPNOW()) return
             !
             call ivu ( prec, w )
             !
             call SWEXCHG ( w, kgrpnt, -1, qmax+nconct-23 )
             if (STPNOW()) return
             !
             call avm ( amat, w(1,1), t )
             !
          endif
          !
          if ( unorm /= 0. ) then
             !
             gamma = 0.
             tnorm = 0.
             !
             do n = nfu, nend
                do m = mfu, mend
                   !
                   nm = kgrpnt(m,n)
                   !
                   do k = 1, qmax
                      !
                      gamma = gamma + t(nm,k)*u(nm,k)
                      tnorm = tnorm + t(nm,k)*t(nm,k)
                      !
                   enddo
                   !
                enddo
             enddo
             !
             rtmp(1) = gamma
             rtmp(2) = tnorm
             call SWREDUCE ( rtmp, 2, SWREAL, SWSUM )
             gamma = rtmp(1)
             tnorm = rtmp(2)
             !
             tnorm = sqrt(tnorm)
             gamma = gamma / (unorm*tnorm)
             gamma = sign(1.,gamma) * max(abs(gamma),crel) * unorm / tnorm
             !
          else
             !
             gamma = 1.
             !
          endif
          !
          if ( .not. gamma /= 0. ) then
             !
             if ( iamout > 0 .and. INODE == MASTER ) then
                write (PRINTF,'(a)') ' ++ bicgstab: process is halted due to a division by zero'
                write (PRINTF,'(a)') '              cause: parameter gamma = 0'
             endif
             !
             goto 40
             !
          endif
          !
          rnorm = 0.
          !
          do n = nfu, nend
             do m = mfu, mend
                !
                nm = kgrpnt(m,n)
                !
                do k = 1, qmax
                   !
                   x  (nm,k) = x(nm,k) + gamma * w(nm,k)
                   res(nm,k) = u(nm,k) - gamma * t(nm,k)
                   rnorm = rnorm + res(nm,k)*res(nm,k)
                   !
                enddo
                !
             enddo
          enddo
          !
          call SWREDUCE ( rnorm, 1, SWREAL, SWSUM )
          rnorm = sqrt(rnorm)
          !
          if ( iamout == 2 .and. INODE == MASTER ) then
             write (PRINTF,'(a,i3,a,e12.6)') ' ++ bicgstab: iter = ',j,'    res = ',rnorm
          endif
          !
          goto 30
          !
       endif
       !
 40    continue
       !
    endif
    !
    ! investigate the reason for stopping
    !
    if ( j > 0 ) rnorm = min(rnorm,unorm)
    !
    if ( rnorm > epslin .and. iamout > 0 .and. INODE == MASTER ) then
       !
       write (PRINTF, '(a)') ' ++ bicgstab: no convergence for solving system of equations'
       write (PRINTF, '(a,i3)'   ) '              total number of iterations     = ',j
       write (PRINTF, '(a,e12.6)') '              2-norm of the residual         = ',rnorm
       write (PRINTF, '(a,e12.6)') '              required accuracy              = ',epslin
       !
    else if ( iamout == 3 .and. INODE == MASTER ) then
       !
       write (PRINTF, '(a,e12.6)') ' ++ bicgstab: 2-norm of the initial residual = ',rnrm0
       write (PRINTF, '(a,i3)'   ) '              total number of iterations     = ',j
       write (PRINTF, '(a,e12.6)') '              2-norm of the residual         = ',rnorm
       !
    endif
    !
    if ( icond == 1 ) then
       !
       ! rescale the solution
       !
       call SWEXCHG ( x, kgrpnt, 1, qmax )
       if (STPNOW()) return
       !
       if ( oned ) then
          !
          do m = ml, mfu, -1
             !
             mu = m + 1
             !
             nm  = kgrpnt(m ,1)
             nmu = kgrpnt(mu,1)
             !
             do k = qmax, 1, -1
                !
                kd  = max(k-1,1   )
                kdd = max(k-2,1   )
                ku  = min(k+1,qmax)
                kuu = min(k+2,qmax)
                !
                x(nm,k) = x(nm,k) - amat(nm,k,3)*x(nmu,k) - amat(nm,k,8)*x(nmu,kd) - amat(nm,k,12)*x(nmu,ku) - amat(nm,k,18)*x(nmu,kdd) - amat(nm,k,22)*x(nmu,kuu)
                !
                ! note: nconct-23 is equivalent with kmax-1 and 2 for Keller-box and standard pressure gradient approximations, respectively
                !
                do i = 1, nconct-23
                   !
                   x(nm,k) = x(nm,k) - amat(nm,k,ishif(i))*x(nm,min(k+i,qmax))
                   !
                enddo
                !
             enddo
             !
          enddo
          !
       else
          !
          do n = nl, nfu, -1
             do m = ml, mfu, -1
                !
                mu = m + 1
                nu = n + 1
                !
                nm  = kgrpnt(m ,n )
                nmu = kgrpnt(mu,n )
                num = kgrpnt(m ,nu)
                !
                do k = qmax, 1, -1
                   !
                   kd  = max(k-1,1   )
                   kdd = max(k-2,1   )
                   ku  = min(k+1,qmax)
                   kuu = min(k+2,qmax)
                   !
                   x(nm,k) = x(nm,k) - amat(nm,k, 3)*x(nmu,k  ) - amat(nm,k, 7)*x(num,k  ) - amat(nm,k, 8)*x(nmu,kd )  &
                                     - amat(nm,k, 9)*x(num,kd ) - amat(nm,k,12)*x(nmu,ku ) - amat(nm,k,13)*x(num,ku )  &
                                     - amat(nm,k,18)*x(nmu,kdd) - amat(nm,k,19)*x(num,kdd) - amat(nm,k,22)*x(nmu,kuu)  &
                                     - amat(nm,k,23)*x(num,kuu)
                   !
                   do i = 1, nconct-23
                      !
                      x(nm,k) = x(nm,k) - amat(nm,k,ishif(i))*x(nm,min(k+i,qmax))
                      !
                   enddo
                   !
                enddo
                !
             enddo
          enddo
          !
       endif
       !
    endif
    !
end subroutine bicgstab
!
subroutine tridiag ( a, b, c, d, x, kgrpnt )
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
!    1.00, Augustus 2014: New subroutine
!
!   Purpose
!
!   Solves tri-diagonal system of equations
!
!   Method
!
!   System of equations in single precision
!
!   Modules used
!
    use ocpcomm4
    use SwashCommdata3, only: mcgrd, mxc
    use SwashFlowdata, only: mf, mfu, ml
    use m_parall
!
    implicit none
!
!   Argument variables
!
    integer, dimension(mxc), intent(in)    :: kgrpnt ! index table containing the address of each (active) grid point
                                                     ! =1: not active grid point
                                                     ! >1: active grid point
    !
    real, dimension(mcgrd), intent(in   )  :: a      ! lower diagonal of tri-diagonal matrix
    real, dimension(mcgrd), intent(inout)  :: b      ! main diagonal of tri-diagonal matrix
    real, dimension(mcgrd), intent(in   )  :: c      ! upper diagonal of tri-diagonal matrix
    real, dimension(mcgrd), intent(inout)  :: d      ! right-hand side of tri-diagonal system of equations
    real, dimension(mcgrd), intent(out  )  :: x      ! solution of the system of equations
!
!   Local variables
!
    integer, save :: ient = 0 ! number of entries in this subroutine
    integer       :: m        ! loop counter
    integer       :: md       ! index of point m-1
    integer       :: mend     ! end index of loop over wl-points in x-direction
    integer       :: mu       ! index of point m+1
    integer       :: nm       ! pointer to m
    integer       :: nmd      ! pointer to m-1
    integer       :: nmf      ! pointer to mf
    integer       :: nmfu     ! pointer to mfu
    integer       :: nml      ! pointer to ml
    integer       :: nmld     ! pointer to mld
    integer       :: nmu      ! pointer to m+1
    !
    real          :: fac      ! a factor
    !
    logical       :: STPNOW   ! indicates that program must stop
!
!   Structure
!
!   Description of the pseudo code
!
!   Source text
!
    if (ltrace) call strace (ient,'tridiag')
    !
    nmf  = kgrpnt(mf   )
    nmfu = kgrpnt(mfu  )
    nml  = kgrpnt(ml   )
    nmld = kgrpnt(ml -1)
    !
    if ( NPROC == 1 .or. NPROC == 2 ) then
       !
       ! The next piece of code works for both serial and parallel runs:
       !
       !    in case of serial run double sweep is carried out
       !    in case of parallel run twisted factorization technique is employed which is perfectly parallelizable for two processors only
       !
       if ( INODE == 1 ) then
          !
          ! forward sweep (elimination)
          !
          ! not computed for point ml at subdomain interface since, this end point is owned by the neighbouring subdomain
          !
          mend = ml - 1
          if ( LMXL ) mend = ml
          !
          do m = mfu+1, mend
             !
             md = m - 1
             !
             nm  = kgrpnt(m )
             nmd = kgrpnt(md)
             !
             fac   = a(nm) / b(nmd)
             b(nm) = b(nm) - fac*c(nmd)
             d(nm) = d(nm) - fac*d(nmd)
             !
          enddo
          !
       else
          !
          ! backward sweep (elimination)
          !
          do m = ml-1, mfu, -1
             !
             mu = m + 1
             !
             nm  = kgrpnt(m )
             nmu = kgrpnt(mu)
             !
             fac   = c(nm) / b(nmu)
             b(nm) = b(nm) - fac*a(nmu)
             d(nm) = d(nm) - fac*d(nmu)
             !
          enddo
          !
       endif
       !
       ! exchange coefficients a, b, c and d with neighbouring subdomains
       !
       call SWEXCHG ( a, kgrpnt, 1, 1 )
       if (STPNOW()) return
       call SWEXCHG ( b, kgrpnt, 1, 1 )
       if (STPNOW()) return
       call SWEXCHG ( c, kgrpnt, 1, 1 )
       if (STPNOW()) return
       call SWEXCHG ( d, kgrpnt, 1, 1 )
       if (STPNOW()) return
       !
       if ( INODE == 1 ) then
          !
          ! forward sweep in point ml (elimination)
          !
          if ( .not.LMXL ) then
             !
             fac    = a(nml) / b(nmld)
             b(nml) = b(nml) - fac*c(nmld)
             d(nml) = d(nml) - fac*d(nmld)
             !
          endif
          !
          ! backward sweep (substitution)
          !
          x(nml) = d(nml) / b(nml)
          !
          do m = ml-1, mfu, -1
             !
             mu = m + 1
             !
             nm  = kgrpnt(m )
             nmu = kgrpnt(mu)
             !
             x(nm) = ( d(nm)-c(nm)*x(nmu) ) / b(nm)
             !
          enddo
          !
       else
          !
          ! backward sweep in point mf (elimination)
          !
          fac    = c(nmf) / b(nmfu)
          b(nmf) = b(nmf) - fac*a(nmfu)
          d(nmf) = d(nmf) - fac*d(nmfu)
          !
          ! forward sweep (substitution)
          !
          x(nmf) = d(nmf) / b(nmf)
          !
          do m = mfu, ml
             !
             md = m - 1
             !
             nm  = kgrpnt(m )
             nmd = kgrpnt(md)
             !
             x(nm) = ( d(nm)-a(nm)*x(nmd) ) / b(nm)
             !
          enddo
          !
       endif
       !
    else
       !
       ! Bondeli's divide and conquer algorithm
       !
       call dac ( a, b, c, d, x, kgrpnt )
       if (STPNOW()) return
       !
    endif
    !
end subroutine tridiag
!
subroutine tridiag2 ( a, b, c, d, x, kgrpnt )
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
!    1.00, Augustus 2014: New subroutine
!
!   Purpose
!
!   Solves tri-diagonal system of equations
!
!   Method
!
!   System of equations in double precision
!
!   Modules used
!
    use ocpcomm4
    use SwashCommdata3, only: mcgrd, mxc
    use SwashFlowdata, only: mf, mfu, ml
    use m_parall
!
    implicit none
!
!   Argument variables
!
    integer, dimension(mxc), intent(in)     :: kgrpnt ! index table containing the address of each (active) grid point
                                                      ! =1: not active grid point
                                                      ! >1: active grid point
    !
    real*8, dimension(mcgrd), intent(in   ) :: a      ! lower diagonal of tri-diagonal matrix
    real*8, dimension(mcgrd), intent(inout) :: b      ! main diagonal of tri-diagonal matrix
    real*8, dimension(mcgrd), intent(in   ) :: c      ! upper diagonal of tri-diagonal matrix
    real*8, dimension(mcgrd), intent(inout) :: d      ! right-hand side of tri-diagonal system of equations
    real*8, dimension(mcgrd), intent(out  ) :: x      ! solution of the system of equations
!
!   Local variables
!
    integer, save :: ient = 0 ! number of entries in this subroutine
    integer       :: m        ! loop counter
    integer       :: md       ! index of point m-1
    integer       :: mend     ! end index of loop over wl-points in x-direction
    integer       :: mu       ! index of point m+1
    integer       :: nm       ! pointer to m
    integer       :: nmd      ! pointer to m-1
    integer       :: nmf      ! pointer to mf
    integer       :: nmfu     ! pointer to mfu
    integer       :: nml      ! pointer to ml
    integer       :: nmld     ! pointer to mld
    integer       :: nmu      ! pointer to m+1
    !
    real*8        :: fac      ! a factor
    !
    logical       :: STPNOW   ! indicates that program must stop
!
!   Structure
!
!   Description of the pseudo code
!
!   Source text
!
    if (ltrace) call strace (ient,'tridiag2')
    !
    nmf  = kgrpnt(mf   )
    nmfu = kgrpnt(mfu  )
    nml  = kgrpnt(ml   )
    nmld = kgrpnt(ml -1)
    !
    if ( NPROC == 1 .or. NPROC == 2 ) then
       !
       ! The next piece of code works for both serial and parallel runs:
       !
       !    in case of serial run double sweep is carried out
       !    in case of parallel run twisted factorization technique is employed which is perfectly parallelizable for two processors only
       !
       if ( INODE == 1 ) then
          !
          ! forward sweep (elimination)
          !
          ! not computed for point ml at subdomain interface since, this end point is owned by the neighbouring subdomain
          !
          mend = ml - 1
          if ( LMXL ) mend = ml
          !
          do m = mfu+1, mend
             !
             md = m - 1
             !
             nm  = kgrpnt(m )
             nmd = kgrpnt(md)
             !
             fac   = a(nm) / b(nmd)
             b(nm) = b(nm) - fac*c(nmd)
             d(nm) = d(nm) - fac*d(nmd)
             !
          enddo
          !
       else
          !
          ! backward sweep (elimination)
          !
          do m = ml-1, mfu, -1
             !
             mu = m + 1
             !
             nm  = kgrpnt(m )
             nmu = kgrpnt(mu)
             !
             fac   = c(nm) / b(nmu)
             b(nm) = b(nm) - fac*a(nmu)
             d(nm) = d(nm) - fac*d(nmu)
             !
          enddo
          !
       endif
       !
       ! exchange coefficients a, b, c and d with neighbouring subdomains
       !
       call SWEXCHGD ( a, kgrpnt, 1, 1 )
       if (STPNOW()) return
       call SWEXCHGD ( b, kgrpnt, 1, 1 )
       if (STPNOW()) return
       call SWEXCHGD ( c, kgrpnt, 1, 1 )
       if (STPNOW()) return
       call SWEXCHGD ( d, kgrpnt, 1, 1 )
       if (STPNOW()) return
       !
       if ( INODE == 1 ) then
          !
          ! forward sweep in point ml (elimination)
          !
          if ( .not.LMXL ) then
             !
             fac    = a(nml) / b(nmld)
             b(nml) = b(nml) - fac*c(nmld)
             d(nml) = d(nml) - fac*d(nmld)
             !
          endif
          !
          ! backward sweep (substitution)
          !
          x(nml) = d(nml) / b(nml)
          !
          do m = ml-1, mfu, -1
             !
             mu = m + 1
             !
             nm  = kgrpnt(m )
             nmu = kgrpnt(mu)
             !
             x(nm) = ( d(nm)-c(nm)*x(nmu) ) / b(nm)
             !
          enddo
          !
       else
          !
          ! backward sweep in point mf (elimination)
          !
          fac    = c(nmf) / b(nmfu)
          b(nmf) = b(nmf) - fac*a(nmfu)
          d(nmf) = d(nmf) - fac*d(nmfu)
          !
          ! forward sweep (substitution)
          !
          x(nmf) = d(nmf) / b(nmf)
          !
          do m = mfu, ml
             !
             md = m - 1
             !
             nm  = kgrpnt(m )
             nmd = kgrpnt(md)
             !
             x(nm) = ( d(nm)-a(nm)*x(nmd) ) / b(nm)
             !
          enddo
          !
       endif
       !
    else
       !
       ! Bondeli's divide and conquer algorithm
       !
       call dac2 ( a, b, c, d, x, kgrpnt )
       if (STPNOW()) return
       !
    endif
    !
end subroutine tridiag2
!
subroutine dac ( a, b, c, d, x, kgrpnt )
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
!    1.00, September 2011: New subroutine
!
!   Purpose
!
!   Solves tri-diagonal system of equations by means of Bondeli's DAC solver
!
!   Method
!
!   System of equations in single precision
!
!   The tri-diagonal system of equations is solved on a parallel computer using a divide and conquer method as described in
!
!   S. Bondeli
!   Divide and conquer: a parallel algorithm for the solution of a tri-diagonal linear system of equations
!   Parallel Computing, vol. 17, 419-434, 1991
!
!   This method is highly parallelizable and requires almost twice the number of operations as double sweep algorithm
!
!   Modules used
!
    use ocpcomm4
    use SwashCommdata3, only: mcgrd, mxc
    use SwashFlowdata, only: mfu, ml
    use m_parall
    use SwashSolvedata, only: z1, z2
!
    implicit none
!
!   Argument variables
!
    integer, dimension(mxc), intent(in)    :: kgrpnt ! index table containing the address of each (active) grid point
                                                     ! =1: not active grid point
                                                     ! >1: active grid point
    !
    real, dimension(mcgrd), intent(inout)  :: a      ! lower diagonal of tri-diagonal matrix
    real, dimension(mcgrd), intent(inout)  :: b      ! main diagonal of tri-diagonal matrix
    real, dimension(mcgrd), intent(inout)  :: c      ! upper diagonal of tri-diagonal matrix
    real, dimension(mcgrd), intent(inout)  :: d      ! right-hand side of tri-diagonal system of equations
    real, dimension(mcgrd), intent(inout)  :: x      ! solution of the system of equations
!
!   Local variables
!
    integer                    :: i        ! loop counter
    integer, save              :: ient = 0 ! number of entries in this subroutine
    integer                    :: ip       ! loop counter
    integer                    :: j        ! counter
    integer                    :: m        ! loop counter
    integer                    :: md       ! index of point m-1
    integer                    :: mend     ! end index of loop over wl-points in x-direction
    integer                    :: mu       ! index of point m+1
    integer                    :: nm       ! pointer to m
    integer                    :: nmd      ! pointer to m-1
    integer                    :: nmend    ! pointer to mend
    integer                    :: nmfu     ! pointer to mfu
    integer                    :: nmu      ! pointer to m+1
    !
    real, dimension(2*NPROC-2) :: alpha    ! solution of small system of equations
    real, dimension(8,NPROC)   :: arrc     ! auxiliary array for collecting data
    real, dimension(8)         :: arrl     ! auxiliary array containing some data of each subdomain
    real, dimension(2*NPROC-2) :: beta     ! right-hand side of small system of equations
    real                       :: fac      ! a factor
    real, dimension(2*NPROC-2) :: r        ! lower diagonal of small system of equations
    real, dimension(2*NPROC-2) :: s        ! main diagonal of small system of equations
    real, dimension(2*NPROC-2) :: t        ! upper diagonal of small system of equations
    !
    logical                    :: STPNOW   ! indicates that program must stop
!
!   Structure
!
!   Description of the pseudo code
!
!   Source text
!
    if (ltrace) call strace (ient,'dac')
    !
    nmfu = kgrpnt(mfu)
    !
    ! initialize the small system of equations
    !
    alpha = 0.
    beta  = 0.
    r     = 0.
    s     = 1.
    t     = 0.
    !
    ! not computed for point ml at subdomain interface since, this end point is owned by the neighbouring subdomain
    !
    mend = ml - 1
    if ( LMXL ) mend = ml
    nmend = kgrpnt(mend)
    !
    ! initialize z2 and z1 as first and last unit vector, respectively
    !
    z1        = 0.
    z2        = 0.
    z1(nmend) = 1.
    z2(nmfu ) = 1.
    !
    ! solve three local systems of equations by means of double sweep
    !
    do m = mfu+1, mend
       !
       md = m - 1
       !
       nm  = kgrpnt(m )
       nmd = kgrpnt(md)
       !
       fac    = a(nm) / b(nmd)
       b (nm) = b (nm) - fac*c (nmd)
       d (nm) = d (nm) - fac*d (nmd)
       z1(nm) = z1(nm) - fac*z1(nmd)
       z2(nm) = z2(nm) - fac*z2(nmd)
       !
    enddo
    !
    x (nmend) = d (nmend) / b(nmend)
    z1(nmend) = z1(nmend) / b(nmend)
    z2(nmend) = z2(nmend) / b(nmend)
    !
    do m = mend-1, mfu, -1
       !
       mu = m + 1
       !
       nm  = kgrpnt(m )
       nmu = kgrpnt(mu)
       !
       x (nm) = ( d (nm)-c(nm)*x (nmu) ) / b(nm)
       z1(nm) = ( z1(nm)-c(nm)*z1(nmu) ) / b(nm)
       z2(nm) = ( z2(nm)-c(nm)*z2(nmu) ) / b(nm)
       !
    enddo
    !
    ! gather first and last elements of matrix and solutions from all subdomains to the master
    !
    arrl(1) = a (nmfu ) + 1.e-8
    arrl(2) = c (nmend) + 1.e-8
    arrl(3) = x (nmend)
    arrl(4) = x (nmfu )
    arrl(5) = z1(nmend)
    arrl(6) = z1(nmfu )
    arrl(7) = z2(nmend)
    arrl(8) = z2(nmfu )
    !
    call SWGATHER ( arrc, 8*NPROC, arrl, 8, SWREAL )
    if (STPNOW()) return
    !
    ! find solution of the small system of 2*(NPROC-1) equations by the master
    !
    if ( INODE == MASTER ) then
       !
       ! build the small system of equations
       !
       fac     =  arrc(1,2)
       s   (1) =  arrc(5,1)*fac
       t   (1) =  1.
       beta(1) = -arrc(3,1)*fac
       !
       fac             =  arrc(2,NPROC-1)
       s   (2*NPROC-2) =  arrc(8,NPROC)*fac
       r   (2*NPROC-2) =  1.
       beta(2*NPROC-2) = -arrc(4,NPROC)*fac
       !
       do ip = 2, NPROC-1
          !
          ! odd index
          !
          j = 2*ip-1
          !
          fac     =  arrc(1,ip+1)
          s   (j) =  arrc(5,ip)*fac
          r   (j) =  arrc(7,ip)*fac
          t   (j) =  1.
          beta(j) = -arrc(3,ip)*fac
          !
          ! even index
          !
          j = 2*ip-2
          !
          fac     =  arrc(2,ip-1)
          s   (j) =  arrc(8,ip)*fac
          r   (j) =  1.
          t   (j) =  arrc(6,ip)*fac
          beta(j) = -arrc(4,ip)*fac
          !
       enddo
       !
       ! solve the small system of equations
       !
       do j = 2, 2*NPROC-2
          !
          fac     = r(j) / s(j-1)
          s   (j) = s   (j) - fac*t   (j-1)
          beta(j) = beta(j) - fac*beta(j-1)
          !
       enddo
       !
       alpha(2*NPROC-2) = beta(2*NPROC-2) / s(2*NPROC-2)
       !
       do j = 2*NPROC-3, 1, -1
          !
          alpha(j) = ( beta(j)-t(j)*alpha(j+1) ) / s(j)
          !
       enddo
       !
    endif
    !
    ! scatter alpha to all nodes
    !
    call SWBROADC ( alpha, 2*(NPROC-1), SWREAL )
    if (STPNOW()) return
    !
    ! compute the solution
    !
    if ( INODE == 1 ) then
       !
       do m = mfu, mend
          nm = kgrpnt(m)
          x(nm) = x(nm) + alpha(1)*z1(nm)
       enddo
       !
    else if ( INODE == NPROC ) then
       !
       do m = mfu, mend
          nm = kgrpnt(m)
          x(nm) = x(nm) + alpha(2*NPROC-2)*z2(nm)
       enddo
       !
    else
       !
       do m = mfu, mend
          nm = kgrpnt(m)
          x(nm) = x(nm) + alpha(2*INODE-2)*z2(nm) + alpha(2*INODE-1)*z1(nm)
       enddo
       !
    endif
    !
end subroutine dac
!
subroutine dac2 ( a, b, c, d, x, kgrpnt )
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
!    1.00, September 2011: New subroutine
!
!   Purpose
!
!   Solves tri-diagonal system of equations by means of Bondeli's DAC solver
!
!   Method
!
!   System of equations in double precision
!
!   The tri-diagonal system of equations is solved on a parallel computer using a divide and conquer method as described in
!
!   S. Bondeli
!   Divide and conquer: a parallel algorithm for the solution of a tri-diagonal linear system of equations
!   Parallel Computing, vol. 17, 419-434, 1991
!
!   This method is highly parallelizable and requires almost twice the number of operations as double sweep algorithm
!
!   Modules used
!
    use ocpcomm4
    use SwashCommdata3, only: mcgrd, mxc
    use SwashFlowdata, only: mfu, ml
    use m_parall
    use SwashSolvedata, only: z12, z22
!
    implicit none
!
!   Argument variables
!
    integer, dimension(mxc), intent(in)     :: kgrpnt ! index table containing the address of each (active) grid point
                                                      ! =1: not active grid point
                                                      ! >1: active grid point
    !
    real*8, dimension(mcgrd), intent(inout) :: a      ! lower diagonal of tri-diagonal matrix
    real*8, dimension(mcgrd), intent(inout) :: b      ! main diagonal of tri-diagonal matrix
    real*8, dimension(mcgrd), intent(inout) :: c      ! upper diagonal of tri-diagonal matrix
    real*8, dimension(mcgrd), intent(inout) :: d      ! right-hand side of tri-diagonal system of equations
    real*8, dimension(mcgrd), intent(inout) :: x      ! solution of the system of equations
!
!   Local variables
!
    integer                      :: i        ! loop counter
    integer, save                :: ient = 0 ! number of entries in this subroutine
    integer                      :: ip       ! loop counter
    integer                      :: j        ! counter
    integer                      :: m        ! loop counter
    integer                      :: md       ! index of point m-1
    integer                      :: mend     ! end index of loop over wl-points in x-direction
    integer                      :: mu       ! index of point m+1
    integer                      :: nm       ! pointer to m
    integer                      :: nmd      ! pointer to m-1
    integer                      :: nmend    ! pointer to mend
    integer                      :: nmfu     ! pointer to mfu
    integer                      :: nmu      ! pointer to m+1
    !
    real*8, dimension(2*NPROC-2) :: alpha    ! solution of small system of equations
    real*8, dimension(8,NPROC)   :: arrc     ! auxiliary array for collecting data
    real*8, dimension(8)         :: arrl     ! auxiliary array containing some data of each subdomain
    real*8, dimension(2*NPROC-2) :: beta     ! right-hand side of small system of equations
    real*8                       :: fac      ! a factor
    real*8, dimension(2*NPROC-2) :: r        ! lower diagonal of small system of equations
    real*8, dimension(2*NPROC-2) :: s        ! main diagonal of small system of equations
    real*8, dimension(2*NPROC-2) :: t        ! upper diagonal of small system of equations
    !
    logical                      :: STPNOW   ! indicates that program must stop
!
!   Structure
!
!   Description of the pseudo code
!
!   Source text
!
    if (ltrace) call strace (ient,'dac2')
    !
    nmfu = kgrpnt(mfu)
    !
    ! initialize the small system of equations
    !
    alpha = 0d0
    beta  = 0d0
    r     = 0d0
    s     = 1d0
    t     = 0d0
    !
    ! not computed for point ml at subdomain interface since, this end point is owned by the neighbouring subdomain
    !
    mend = ml - 1
    if ( LMXL ) mend = ml
    nmend = kgrpnt(mend)
    !
    ! initialize z22 and z12 as first and last unit vector, respectively
    !
    z12        = 0d0
    z22        = 0d0
    z12(nmend) = 1d0
    z22(nmfu ) = 1d0
    !
    ! solve three local systems of equations by means of double sweep
    !
    do m = mfu+1, mend
       !
       md = m - 1
       !
       nm  = kgrpnt(m )
       nmd = kgrpnt(md)
       !
       fac    = a(nm) / b(nmd)
       b  (nm) = b  (nm) - fac*c  (nmd)
       d  (nm) = d  (nm) - fac*d  (nmd)
       z12(nm) = z12(nm) - fac*z12(nmd)
       z22(nm) = z22(nm) - fac*z22(nmd)
       !
    enddo
    !
    x  (nmend) = d  (nmend) / b(nmend)
    z12(nmend) = z12(nmend) / b(nmend)
    z22(nmend) = z22(nmend) / b(nmend)
    !
    do m = mend-1, mfu, -1
       !
       mu = m + 1
       !
       nm  = kgrpnt(m )
       nmu = kgrpnt(mu)
       !
       x  (nm) = ( d  (nm)-c(nm)*x  (nmu) ) / b(nm)
       z12(nm) = ( z12(nm)-c(nm)*z12(nmu) ) / b(nm)
       z22(nm) = ( z22(nm)-c(nm)*z22(nmu) ) / b(nm)
       !
    enddo
    !
    ! gather first and last elements of matrix and solutions from all subdomains to the master
    !
    arrl(1) = a  (nmfu ) + 1.d-8
    arrl(2) = c  (nmend) + 1.d-8
    arrl(3) = x  (nmend)
    arrl(4) = x  (nmfu )
    arrl(5) = z12(nmend)
    arrl(6) = z12(nmfu )
    arrl(7) = z22(nmend)
    arrl(8) = z22(nmfu )
    !
    call SWGATHERD ( arrc, 8*NPROC, arrl, 8 )
    if (STPNOW()) return
    !
    ! find solution of the small system of 2*(NPROC-1) equations by the master
    !
    if ( INODE == MASTER ) then
       !
       ! build the small system of equations
       !
       fac     =  arrc(1,2)
       s   (1) =  arrc(5,1)*fac
       t   (1) =  1d0
       beta(1) = -arrc(3,1)*fac
       !
       fac             =  arrc(2,NPROC-1)
       s   (2*NPROC-2) =  arrc(8,NPROC)*fac
       r   (2*NPROC-2) =  1d0
       beta(2*NPROC-2) = -arrc(4,NPROC)*fac
       !
       do ip = 2, NPROC-1
          !
          ! odd index
          !
          j = 2*ip-1
          !
          fac     =  arrc(1,ip+1)
          s   (j) =  arrc(5,ip)*fac
          r   (j) =  arrc(7,ip)*fac
          t   (j) =  1d0
          beta(j) = -arrc(3,ip)*fac
          !
          ! even index
          !
          j = 2*ip-2
          !
          fac     =  arrc(2,ip-1)
          s   (j) =  arrc(8,ip)*fac
          r   (j) =  1d0
          t   (j) =  arrc(6,ip)*fac
          beta(j) = -arrc(4,ip)*fac
          !
       enddo
       !
       ! solve the small system of equations
       !
       do j = 2, 2*NPROC-2
          !
          fac     = r(j) / s(j-1)
          s   (j) = s   (j) - fac*t   (j-1)
          beta(j) = beta(j) - fac*beta(j-1)
          !
       enddo
       !
       alpha(2*NPROC-2) = beta(2*NPROC-2) / s(2*NPROC-2)
       !
       do j = 2*NPROC-3, 1, -1
          !
          alpha(j) = ( beta(j)-t(j)*alpha(j+1) ) / s(j)
          !
       enddo
       !
    endif
    !
    ! scatter alpha to all nodes
    !
    call SWBROADCD ( alpha, 2*(NPROC-1) )
    if (STPNOW()) return
    !
    ! compute the solution
    !
    if ( INODE == 1 ) then
       !
       do m = mfu, mend
          nm = kgrpnt(m)
          x(nm) = x(nm) + alpha(1)*z12(nm)
       enddo
       !
    else if ( INODE == NPROC ) then
       !
       do m = mfu, mend
          nm = kgrpnt(m)
          x(nm) = x(nm) + alpha(2*NPROC-2)*z22(nm)
       enddo
       !
    else
       !
       do m = mfu, mend
          nm = kgrpnt(m)
          x(nm) = x(nm) + alpha(2*INODE-2)*z22(nm) + alpha(2*INODE-1)*z12(nm)
       enddo
       !
    endif
    !
end subroutine dac2
!
subroutine newton1D ( a, b, c, d, lo, up, x, kgrpnt )
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
!    1.00: Dirk Rijnsdorp
!
!   Updates
!
!    1.00, August 2014: New subroutine
!
!   Purpose
!
!   Solves a piecewise linear tri-diagonal system by means of the nested Newton iteration method
!
!   Method
!
!   The piecewise linear system is of the following form
!
!   max[l,min(u,x)] + Ax = b
!
!   with A a symmetric and positive definite matrix and l, u, and b are known vectors
!
!   L. Brugnano and V. Casulli
!   Iterative solution of piecewise linear system and applications to flows in porous media
!   SIAM J. Sci. Comput., vol. 31, 1858-1873, 2009
!
!   Modules used
!
    use ocpcomm4
    use SwashCommdata3, only: mcgrd, mxc
    use SwashFlowdata, only: mfu, ml
    use m_parall
    use SwashSolvedata, only: p0, p1, q0, q1, ba, da, sola, vt
!
    implicit none
!
!   Argument variables
!
    integer, dimension(mxc) , intent(in ) :: kgrpnt ! index table containing the address of each (active) grid point
                                                    ! =1: not active grid point
                                                    ! >1: active grid point
    !
    real*8, dimension(mcgrd), intent(in ) :: a      ! lower diagonal of tri-diagonal matrix
    real*8, dimension(mcgrd), intent(in ) :: b      ! main diagonal of tri-diagonal matrix
    real*8, dimension(mcgrd), intent(in ) :: c      ! upper diagonal of tri-diagonal matrix
    real*8, dimension(mcgrd), intent(in ) :: d      ! right-hand side of tri-diagonal system of equations
    real*8, dimension(mcgrd), intent(in ) :: lo     ! lower bound of solution
    real*8, dimension(mcgrd), intent(in ) :: up     ! upper bound of solution
    real  , dimension(mcgrd), intent(out) :: x      ! solution of the system of equations
!
!   Parameter variables
!
    integer, parameter :: maxit = 100    ! maximum number of iterations
    real*8 , parameter :: reps  = 1.d-14 ! convergence criterion
!
!   Local variables
!
    integer, save :: ient = 0 ! number of entries in this subroutine
    integer       :: j        ! iteration counter of outer iteration
    integer       :: k        ! iteration counter of inner iteration
    integer       :: m        ! loop counter
    integer       :: mend     ! end index of loop over wl-points in x-direction
    integer       :: nm       ! pointer to m
    !
    real*8        :: rnorm    ! 2-norm of vector
!
!   Structure
!
!   Description of the pseudo code
!
!   Source text
!
    if (ltrace) call strace (ient,'newton1D')
    !
    ! not computed for point ml at subdomain interface since, this end point is owned by the neighbouring subdomain
    !
    mend = ml - 1
    if ( LMXL ) mend = ml
    !
    ! initialize solution
    !
    sola = 0d0
    !
    ! set diagonal matrix q to zero
    !
    q0 = 0d0
    !
    ! nested Newton iterations
    !
    do j = 1, maxit
       !
       ! set diagonal matrix p to identity
       !
       p0 = 1d0
       !
       do k = 1, maxit
          !
          ! set main diagonal and right-hand side
          !
          ba = b + p0 - q0
          da = d - lo + p0*lo - q0*up
          !
          ! solve tri-diagonal system of equations by means of double sweep
          !
          call tridiag2 ( a, ba, c, da, sola, kgrpnt )
          !
          ! update diagonal matrix p
          !
          where ( .not. sola < lo )
             p1 = 1d0
          elsewhere
             p1 = 0d0
          end where
          !
          ! do a check for convergence
          !
          if ( any( p1 < q0 ) ) call msgerr (1, 'no convergence in nested Newton iteration')
          !
          ! check convergence
          !
          vt = ( p1 - p0 ) * ( sola - lo )
          !
          rnorm = 0d0
          !
          do m = mfu, mend
             !
             nm = kgrpnt(m)
             !
             rnorm = rnorm + vt(nm)*vt(nm)
             !
          enddo
          !
          call SWREDUCD ( rnorm, 1, SWSUM )
          rnorm = sqrt(rnorm)
          !
          if ( .not. rnorm > reps ) exit
          !
          p0 = p1
          !
       enddo
       !
       ! update diagonal matrix q
       !
       where ( sola > up )
          q1 = 1d0
       elsewhere
          q1 = 0d0
       end where
       !
       ! check convergence
       !
       vt = ( q1 - q0 ) * ( sola - up )
       !
       rnorm = 0d0
       !
       do m = mfu, mend
          !
          nm = kgrpnt(m)
          !
          rnorm = rnorm + vt(nm)*vt(nm)
          !
       enddo
       !
       call SWREDUCD ( rnorm, 1, SWSUM )
       rnorm = sqrt(rnorm)
       !
       if ( .not. rnorm > reps ) exit
       !
       q0 = q1
       !
    enddo
    !
    x = real(sola)
    !
    if ( ITEST >= 60 .and. INODE == MASTER ) write (PRINTF, 101) k, j
    !
 101 format (' ++ nested Newton: number of inner iterations = ',i2,/,   &
             '                   number of outer iterations = ',i2)
    !
end subroutine newton1D
!
subroutine newton2D ( amat, rhs, lo, up, x )
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
!    1.00: Dirk Rijnsdorp
!
!   Updates
!
!    1.00, August 2014: New subroutine
!
!   Purpose
!
!   Solves a piecewise linear penta-diagonal system by means of the nested Newton iteration method
!
!   Method
!
!   The piecewise linear system is of the following form
!
!   max[l,min(u,x)] + Ax = b
!
!   with A a symmetric and positive definite matrix and l, u, and b are known vectors
!
!   L. Brugnano and V. Casulli
!   Iterative solution of piecewise linear system and applications to flows in porous media
!   SIAM J. Sci. Comput., vol. 31, 1858-1873, 2009
!
!   Modules used
!
    use ocpcomm4
    use SwashCommdata3, only: mcgrd, pnums
    use SwashFlowdata, only: mfu, ml, nfu, nl
    use m_genarr, only: kgrpnt
    use m_parall
    use SwashSolvedata, only: p0, p1, q0, q1, amata, rhsa, sola, vt
!
    implicit none
!
!   Argument variables
!
    real*8, dimension(mcgrd,5), intent(in ) :: amat ! coefficient matrix of the system of equations
    real*8, dimension(mcgrd  ), intent(in ) :: lo   ! lower bound of solution
    real*8, dimension(mcgrd  ), intent(in ) :: rhs  ! right-hand side of the system of equations
    real*8, dimension(mcgrd  ), intent(in ) :: up   ! upper bound of solution
    real  , dimension(mcgrd  ), intent(out) :: x    ! solution of the system of equations
!
!   Local variables
!
    integer, save :: ient = 0 ! number of entries in this subroutine
    integer       :: j        ! iteration counter of outer iteration
    integer       :: k        ! iteration counter of inner iteration
    integer       :: m        ! loop counter
    integer       :: maxit    ! maximum number of iterations
    integer       :: mend     ! end index of loop over wl-points in x-direction
    integer       :: n        ! loop counter
    integer       :: nend     ! end index of loop over wl-points in y-direction
    integer       :: nm       ! pointer to m,n
    !
    real*8        :: reps     ! convergence criterion
    real*8        :: rnorm    ! 2-norm of vector
!
!   Structure
!
!   Description of the pseudo code
!
!   Source text
!
    if (ltrace) call strace (ient,'newton2D')
    !
    reps  = dble(pnums(21))
    maxit = nint(pnums(24))
    !
    ! not computed for points ml and nl at subdomain interfaces since, these end points are owned by the neighbouring subdomains
    !
    mend = ml - 1
    if ( LMXL ) mend = ml
    !
    nend = nl - 1
    if ( LMYL ) nend = nl
    !
    ! initialize solution
    !
    sola = dble(x)
    !
    ! set diagonal matrix q to zero
    !
    q0 = 0d0
    !
    ! nested Newton iterations
    !
    do j = 1, maxit
       !
       ! set diagonal matrix p to identity
       !
       p0 = 1d0
       !
       do k = 1, maxit
          !
          ! copy original coefficient matrix
          !
          amata = amat
          !
          ! set main diagonal and right-hand side
          !
          amata(:,1) = amat(:,1) + p0 - q0
          rhsa       = rhs       - lo + p0*lo - q0*up
          !
          ! solve penta-diagonal system of equations by means of preconditioned CG method
          !
          call pcg2 ( amata, rhsa, sola )
          !
          ! update diagonal matrix p
          !
          where ( .not. sola < lo )
             p1 = 1d0
          elsewhere
             p1 = 0d0
          end where
          !
          ! do a check for convergence
          !
          if ( any( p1 < q0 ) ) call msgerr (1, 'no convergence in nested Newton iteration')
          !
          ! check convergence
          !
          vt = ( p1 - p0 ) * ( sola - lo )
          !
          rnorm = 0d0
          !
          do n = nfu, nend
             do m = mfu, mend
                !
                nm = kgrpnt(m,n)
                !
                rnorm = rnorm + vt(nm)*vt(nm)
                !
             enddo
          enddo
          !
          call SWREDUCD ( rnorm, 1, SWSUM )
          rnorm = sqrt(rnorm)
          !
          if ( .not. rnorm > reps ) exit
          !
          p0 = p1
          !
       enddo
       !
       ! update diagonal matrix q
       !
       where ( sola > up )
          q1 = 1d0
       elsewhere
          q1 = 0d0
       end where
       !
       ! check convergence
       !
       vt = ( q1 - q0 ) * ( sola - up )
       !
       rnorm = 0d0
       !
       do n = nfu, nend
          do m = mfu, mend
             !
             nm = kgrpnt(m,n)
             !
             rnorm = rnorm + vt(nm)*vt(nm)
             !
          enddo
       enddo
       !
       call SWREDUCD ( rnorm, 1, SWSUM )
       rnorm = sqrt(rnorm)
       !
       if ( .not. rnorm > reps ) exit
       !
       q0 = q1
       !
    enddo
    !
    x = real(sola)
    !
    if ( ITEST >= 60 .and. INODE == MASTER ) write (PRINTF, 101) k, j
    !
 101 format (' ++ nested Newton: number of inner iterations = ',i2,/,   &
             '                   number of outer iterations = ',i2)
    !
end subroutine newton2D
