!
!  SWASH service routines for flexible mesh
!
!  Contents of this file
!
!     perot
!     chkdiv
!     pcgu  (single precision)
!     pcgu2 (double precision)
!     iluu
!     iludu
!     bicgstabu
!     bicgstab3
!     newtonU
!     csrf
!
subroutine perot ( u, ks, ke )
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
!    1.00, January 2017: New subroutine
!
!   Purpose
!
!   Computes the cell velocity vector in terms of face normal components using Perot's formula
!
!   Modules used
!
    use ocpcomm4
    use SwanGriddata
    use SwanGridobjects
    use SwashFlowdata, only: uvc
!
    implicit none
!
!   Argument variables
!
    integer,                       intent(in) :: ke       ! last index in vertical direction
    integer,                       intent(in) :: ks       ! first index in vertical direction
    real, dimension(nfaces,ks:ke), intent(in) :: u        ! flow velocity at faces
!
!   Local variables
!
    integer                                   :: icell    ! loop counter over cells
    integer                                   :: icelll   ! left cell of present face
    integer                                   :: icellr   ! right cell of present face
    integer, save                             :: ient = 0 ! number of entries in this subroutine
    integer                                   :: iface    ! face index
    integer                                   :: jf       ! loop counter
    integer                                   :: k        ! loop counter in vertical direction
    !
    real                                      :: area     ! area of present cell
    real                                      :: lf       ! length of face
    real                                      :: rsgn     ! sign for indicating face orientation
    real                                      :: xc       ! x-coordinate of cell circumcenter
    real                                      :: xf       ! x-coordinate of face center
    real                                      :: yc       ! y-coordinate of cell circumcenter
    real                                      :: yf       ! y-coordinate of face center
    !
    type(celltype), dimension(:), pointer     :: cell     ! datastructure for cells with their attributes
    type(facetype), dimension(:), pointer     :: face     ! datastructure for faces with their attributes
!
!   Structure
!
!   Description of the pseudo code
!
!   Source text
!
    if (ltrace) call strace (ient,'perot')
    !
    ! point to cell and face objects
    !
    cell => gridobject%cell_grid
    face => gridobject%face_grid
    !
    uvc = 0.
    !
    do icell = 1, ncells
       !
       ! get coordinates of the cell circumcenter
       !
       xc = cell(icell)%attr(CELLCCX)
       yc = cell(icell)%attr(CELLCCY)
       !
       do k = ks, ke
          !
          do jf = 1, cell(icell)%nof
             !
             ! face identifier
             !
             iface = cell(icell)%face(jf)%atti(FACEID)
             !
             ! get length of the face
             !
             lf = face(iface)%attr(FACELEN)
             !
             ! get coordinates of the midface
             !
             xf = face(iface)%attr(FACEMX)
             yf = face(iface)%attr(FACEMY)
             !
             ! consider left and right cells of the face
             !
             icelll = face(iface)%atti(FACECL)
             icellr = face(iface)%atti(FACECR)
             !
             ! take into account orientation of the face
             !
             if ( icell == icelll ) then
                rsgn =  1.
             else if ( icell == icellr ) then
                rsgn = -1.
             endif
             !
             ! compute the cell velocity vector using Perot's formula
             !
             uvc(icell,k,1) = uvc(icell,k,1) + rsgn * lf * ( xf - xc ) * u(iface,k)
             uvc(icell,k,2) = uvc(icell,k,2) + rsgn * lf * ( yf - yc ) * u(iface,k)
             !
          enddo
          !
       enddo
       !
       area = cell(icell)%attr(CELLAREA)
       !
       uvc(icell,:,:) = uvc(icell,:,:) / area
       !
    enddo
    !
end subroutine perot
!
subroutine chkdiv
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
!    1.00, January 2017: New subroutine
!
!   Purpose
!
!   Computes horizontal divergence operators for checking purposes
!
!   Modules used
!
    use ocpcomm4
    use SwashCommdata3
    use SwanGriddata
    use SwanGridobjects
    use SwashFlowdata
!
    implicit none
!
!   Local variables
!
    integer                                   :: icell    ! loop counter over cells
    integer                                   :: icelll   ! left cell of present face
    integer                                   :: icellr   ! right cell of present face
    integer, save                             :: ient = 0 ! number of entries in this subroutine
    integer                                   :: iface    ! face index
    integer                                   :: jf       ! loop counter
    !
    real                                      :: area     ! area of present cell
    real                                      :: lf       ! length of face
    real                                      :: rsgn     ! sign for indicating face orientation
    !
    type(celltype), dimension(:), pointer     :: cell     ! datastructure for cells with their attributes
    type(facetype), dimension(:), pointer     :: face     ! datastructure for faces with their attributes
!
!   Structure
!
!   Description of the pseudo code
!
!   Source text
!
    if (ltrace) call strace (ient,'chkdiv')
    !
    ! point to cell and face objects
    !
    cell => gridobject%cell_grid
    face => gridobject%face_grid
    !
    ! compute local divergence of flow velocity and mass flux
    !
    do icell = 1, ncells
       !
       divu(icell,:) = 0.
       divq(icell,:) = 0.
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
          divu(icell,:) = divu(icell,:) + rsgn * lf * u1(iface,:)
          !
          if ( kmax == 1 ) then
             divq(icell,1) = divq(icell,1) + rsgn * lf * hu (iface  ) * u1(iface,1)
          else
             divq(icell,:) = divq(icell,:) + rsgn * lf * hku(iface,:) * u1(iface,:)
          endif
          !
       enddo
       !
       area = cell(icell)%attr(CELLAREA)
       !
       divu(icell,:) = divu(icell,:) / area
       !
       if ( kmax == 1 ) then
          divq(icell,1) = divq(icell,1) / area / hs (icell  )
       else
          divq(icell,:) = divq(icell,:) / area / hks(icell,:)
       endif
       !
    enddo
    !
end subroutine chkdiv
!
subroutine pcgu ( amat, rhs, x )
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
!    1.00, February 2020: New subroutine
!
!   Purpose
!
!   Solves system of equations by means of preconditioned conjugate gradient algorithm
!
!   Method
!
!   System of equations assembled on triangular mesh
!
!   System of equations in single precision and in CSR format
!
!   The system is supposed to be symmetric and positive definite
!
!   An upper triangular matrix U and a lower triangular matrix L are computed according to the following rules:
!
!   A = LU - R
!
!   (a)  diag(L) = I,
!   (b)  the non-zero pattern of L and U are equal to the non-zero pattern of A
!   (c)  if A(i,j) <> 0 then L*U(i,j) = A(i,j)
!
!   This incomplete decomposition is used as the ILU preconditioner. If the last rule (c) is replaced by
!
!   rowsum(LU) = rowsum(A)
!
!   the MILU preconditioner is obtained. We also used a RILU preconditioner, which is the average of ILU and MILU.
!   The contributions on the diagonal are weighted by the scalar amod: amod = 0 corresponds with ILU, whereas
!   amod = 1 corresponds with MILU.
!
!   Auke van der Ploeg (1994)
!   Preconditioning for sparse matrices with applications
!   PhD thesis, Delft University of Technology, Delft, the Netherlands
!
!   The well-known CG method is applied for the solution of preconditioned Ax = b
!
!   Modules used
!
    use ocpcomm4
    use SwashCommdata3, only: pnums
    use SwashSolvedata
    use SwanGriddata
!
    implicit none
!
!   Argument variables
!
    real, dimension(ncells,0:3), intent(in   ) :: amat ! coefficient matrix of the system of equations
    real, dimension(ncells    ), intent(in   ) :: rhs  ! right-hand side of the system of equations
    real, dimension(ncells    ), intent(inout) :: x    ! solution of the system of equations
!
!   Parameter variables
!
    double precision, parameter :: epsmac = 1d-18 ! machine precision number
!
!   Local variables
!
    integer       :: icell    ! loop counter over cells
    integer, save :: ient = 0 ! number of entries in this subroutine
    integer       :: ip       ! row pointer
    integer       :: j        ! iteration counter
    integer       :: jc       ! loop counter
    integer       :: k        ! counter
    integer       :: m        ! counter
    integer       :: maxit    ! maximum number of iterations
    integer       :: n        ! loop counter
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
!   Structure
!
!   Description of the pseudo code
!
!   Source text
!
    if (ltrace) call strace (ient,'pcgu')
    !
    reps  = pnums(21)
    maxit = nint(pnums(24))
    amod  = pnums(26)
    !
    ! restore matrix using the CSR format
    !
    do icell = 1, ncells
       !
       ip = ias(icell)
       !
       do jc = 1, 4
          !
          k = irws(icell,jc)
          !
          if ( k /= 0 ) axs2(ip+jc-1) = amat(icell,k-1)
          !
       enddo
       !
    enddo
    !
    ! compute incomplete LU factorization
    !
    prec2(1:) = axs2
    iplu      = 0
    !
    do icell = 2, ncells
       !
       prec2(0) = 0.
       !
       do jc = ias(icell), ias(icell+1)-1
          iplu(jas(jc)) = jc
       enddo
       !
       do jc = ias(icell), dis(icell)-1
          !
          m = jas(jc)
          !
          prec2(jc) = prec2(jc) / prec2(dis(m))
          !
          do n = dis(m)+1, ias(m+1)-1
             k = iplu(jas(n))
             prec2(k) = prec2(k) - prec2(jc) * prec2(n)
          enddo
          !
       enddo
       !
       prec2(dis(icell)) = prec2(dis(icell)) + amod * prec2(0)
       !
       do jc = ias(icell), ias(icell+1)-1
          iplu(jas(jc)) = 0
       enddo
       !
    enddo
    !
    ! compute initial residual
    !
    do icell = 1, ncells
       !
       resd(icell) = 0.
       !
       do jc = ias(icell), ias(icell+1)-1
          resd(icell) = resd(icell) + axs2(jc) * x(jas(jc))
       enddo
       !
    enddo
    !
    resd = rhs - resd
    !
    ! solve LU z = r
    !
    do icell = 1, ncells
       !
       z(icell) = resd(icell)
       !
       do jc = ias(icell), dis(icell)-1
          !
          z(icell) = z(icell) - prec2(jc) * z(jas(jc))
          !
       enddo
       !
    enddo
    !
    do icell = ncells, 1, -1
       !
       do jc = dis(icell)+1, ias(icell+1)-1
          !
          z(icell) = z(icell) - prec2(jc) * z(jas(jc))
          !
       enddo
       !
       z(icell) = z(icell) / prec2(dis(icell))
       !
    enddo
    !
    s = z
    !
    rho = 0.
    do icell = 1, ncells
       !
       rho = rho + resd(icell)*z(icell)
       !
    enddo
    !
    rnrm0 = sqrt(rho)
    !
    epslin = reps*rnrm0
    ueps   = 1000.*real(epsmac)*rnrm0
    !
    if ( epslin < ueps .and. rnrm0 > 0. ) then
       !
       if ( iamout > 0 ) then
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
       do icell = 1, ncells
          !
          r(icell) = 0.
          !
          do jc = ias(icell), ias(icell+1)-1
             r(icell) = r(icell) + axs2(jc) * s(jas(jc))
          enddo
          !
       enddo
       !
       a = 0.
       do icell = 1, ncells
          !
          a = a + s(icell)*r(icell)
          !
       enddo
       !
       alpha = rho / a
       !
       x    = x    + alpha*s
       resd = resd - alpha*r
       !
       ! solve LU z = r
       !
       do icell = 1, ncells
          !
          z(icell) = resd(icell)
          !
          do jc = ias(icell), dis(icell)-1
             !
             z(icell) = z(icell) - prec2(jc) * z(jas(jc))
             !
          enddo
          !
       enddo
       !
       do icell = ncells, 1, -1
          !
          do jc = dis(icell)+1, ias(icell+1)-1
             !
             z(icell) = z(icell) - prec2(jc) * z(jas(jc))
             !
          enddo
          !
          z(icell) = z(icell) / prec2(dis(icell))
          !
       enddo
       !
       rhold = rho
       !
       rho = 0.
       do icell = 1, ncells
          !
          rho = rho + resd(icell)*z(icell)
          !
       enddo
       rnorm = sqrt(rho)
       !
       if ( iamout == 2 ) then
          write (PRINTF, '(a,i3,a,e12.6)') ' ++ pcg: iter = ',j,'    res = ',rnorm
       endif
       !
       if ( .not. rnorm > epslin ) exit
       !
       beta = rho / rhold
       !
       s = z + beta*s
       !
    enddo
    !
    ! investigate the reason for stopping
    !
    if ( rnorm > epslin .and. iamout > 0 ) then
       !
       write (PRINTF, '(a)') ' ++ pcg: no convergence for solving system of equations'
       write (PRINTF, '(a,i3)'   ) '         total number of iterations     = ',j
       write (PRINTF, '(a,e12.6)') '         2-norm of the residual         = ',rnorm
       write (PRINTF, '(a,e12.6)') '         required accuracy              = ',epslin
       !
    else if ( iamout == 3 ) then
       !
       write (PRINTF, '(a,e12.6)') ' ++ pcg: 2-norm of the initial residual = ',rnrm0
       write (PRINTF, '(a,i3)'   ) '         total number of iterations     = ',j
       write (PRINTF, '(a,e12.6)') '         2-norm of the residual         = ',rnorm
       !
    endif
    !
end subroutine pcgu
!
subroutine pcgu2 ( amat, rhs, x )
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
!    1.00, February 2020: New subroutine
!
!   Purpose
!
!   Solves system of equations by means of preconditioned conjugate gradient algorithm
!
!   Method
!
!   System of equations assembled on triangular mesh
!
!   System of equations in double precision and in CSR format
!
!   The system is supposed to be symmetric and positive definite
!
!   An upper triangular matrix U and a lower triangular matrix L are computed according to the following rules:
!
!   A = LU - R
!
!   (a)  diag(L) = I,
!   (b)  the non-zero pattern of L and U are equal to the non-zero pattern of A
!   (c)  if A(i,j) <> 0 then L*U(i,j) = A(i,j)
!
!   This incomplete decomposition is used as the ILU preconditioner. If the last rule (c) is replaced by
!
!   rowsum(LU) = rowsum(A)
!
!   the MILU preconditioner is obtained. We also used a RILU preconditioner, which is the average of ILU and MILU.
!   The contributions on the diagonal are weighted by the scalar amod: amod = 0 corresponds with ILU, whereas
!   amod = 1 corresponds with MILU.
!
!   Auke van der Ploeg (1994)
!   Preconditioning for sparse matrices with applications
!   PhD thesis, Delft University of Technology, Delft, the Netherlands
!
!   The well-known CG method is applied for the solution of preconditioned Ax = b
!
!   Modules used
!
    use ocpcomm4
    use SwashCommdata3, only: pnums
    use SwashSolvedata
    use SwanGriddata
!
    implicit none
!
!   Argument variables
!
    real*8, dimension(ncells,0:3), intent(in   ) :: amat ! coefficient matrix of the system of equations
    real*8, dimension(ncells    ), intent(in   ) :: rhs  ! right-hand side of the system of equations
    real*8, dimension(ncells    ), intent(inout) :: x    ! solution of the system of equations
!
!   Parameter variables
!
    double precision, parameter :: epsmac = 1d-18 ! machine precision number
!
!   Local variables
!
    integer       :: icell    ! loop counter over cells
    integer, save :: ient = 0 ! number of entries in this subroutine
    integer       :: ip       ! row pointer
    integer       :: j        ! iteration counter
    integer       :: jc       ! loop counter
    integer       :: k        ! counter
    integer       :: m        ! counter
    integer       :: maxit    ! maximum number of iterations
    integer       :: n        ! loop counter
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
!   Structure
!
!   Description of the pseudo code
!
!   Source text
!
    if (ltrace) call strace (ient,'pcgu2')
    !
    reps  = dble(pnums(21))
    maxit = nint(pnums(24))
    amod  = dble(pnums(26))
    !
    ! restore matrix using the CSR format
    !
    do icell = 1, ncells
       !
       ip = ias(icell)
       !
       do jc = 1, 4
          !
          k = irws(icell,jc)
          !
          if ( k /= 0 ) axs(ip+jc-1) = amat(icell,k-1)
          !
       enddo
       !
    enddo
    !
    ! compute incomplete LU factorization
    !
    prec3(1:) = axs
    iplu      = 0
    !
    do icell = 2, ncells
       !
       prec3(0) = 0.
       !
       do jc = ias(icell), ias(icell+1)-1
          iplu(jas(jc)) = jc
       enddo
       !
       do jc = ias(icell), dis(icell)-1
          !
          m = jas(jc)
          !
          prec3(jc) = prec3(jc) / prec3(dis(m))
          !
          do n = dis(m)+1, ias(m+1)-1
             k = iplu(jas(n))
             prec3(k) = prec3(k) - prec3(jc) * prec3(n)
          enddo
          !
       enddo
       !
       prec3(dis(icell)) = prec3(dis(icell)) + amod * prec3(0)
       !
       do jc = ias(icell), ias(icell+1)-1
          iplu(jas(jc)) = 0
       enddo
       !
    enddo
    !
    ! compute initial residual
    !
    do icell = 1, ncells
       !
       resd2(icell) = 0d0
       !
       do jc = ias(icell), ias(icell+1)-1
          resd2(icell) = resd2(icell) + axs(jc) * x(jas(jc))
       enddo
       !
    enddo
    !
    resd2 = rhs - resd2
    !
    ! solve LU z = r
    !
    do icell = 1, ncells
       !
       z3(icell) = resd2(icell)
       !
       do jc = ias(icell), dis(icell)-1
          !
          z3(icell) = z3(icell) - prec3(jc) * z3(jas(jc))
          !
       enddo
       !
    enddo
    !
    do icell = ncells, 1, -1
       !
       do jc = dis(icell)+1, ias(icell+1)-1
          !
          z3(icell) = z3(icell) - prec3(jc) * z3(jas(jc))
          !
       enddo
       !
       z3(icell) = z3(icell) / prec3(dis(icell))
       !
    enddo
    !
    s2 = z3
    !
    rho = 0d0
    do icell = 1, ncells
       !
       rho = rho + resd2(icell)*z3(icell)
       !
    enddo
    !
    rnrm0 = sqrt(rho)
    !
    epslin = reps*rnrm0
    ueps   = 1d3*epsmac*rnrm0
    !
    if ( epslin < ueps .and. rnrm0 > 0d0 ) then
       !
       if ( iamout > 0 ) then
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
       do icell = 1, ncells
          !
          r2(icell) = 0d0
          !
          do jc = ias(icell), ias(icell+1)-1
             r2(icell) = r2(icell) + axs(jc) * s2(jas(jc))
          enddo
          !
       enddo
       !
       a = 0d0
       do icell = 1, ncells
          !
          a = a + s2(icell)*r2(icell)
          !
       enddo
       !
       alpha = rho / a
       !
       x     = x     + alpha*s2
       resd2 = resd2 - alpha*r2
       !
       ! solve LU z = r
       !
       do icell = 1, ncells
          !
          z3(icell) = resd2(icell)
          !
          do jc = ias(icell), dis(icell)-1
             !
             z3(icell) = z3(icell) - prec3(jc) * z3(jas(jc))
             !
          enddo
          !
       enddo
       !
       do icell = ncells, 1, -1
          !
          do jc = dis(icell)+1, ias(icell+1)-1
             !
             z3(icell) = z3(icell) - prec3(jc) * z3(jas(jc))
             !
          enddo
          !
          z3(icell) = z3(icell) / prec3(dis(icell))
          !
       enddo
       !
       rhold = rho
       !
       rho = 0d0
       do icell = 1, ncells
          !
          rho = rho + resd2(icell)*z3(icell)
          !
       enddo
       rnorm = sqrt(rho)
       !
       if ( iamout == 2 ) then
          write (PRINTF, '(a,i3,a,e12.6)') ' ++ pcg: iter = ',j,'    res = ',rnorm
       endif
       !
       if ( .not. rnorm > epslin ) exit
       !
       beta = rho / rhold
       !
       s2 = z3 + beta*s2
       !
    enddo
    !
    ! investigate the reason for stopping
    !
    if ( rnorm > epslin .and. iamout > 0 ) then
       !
       write (PRINTF, '(a)') ' ++ pcg: no convergence for solving system of equations'
       write (PRINTF, '(a,i3)'   ) '         total number of iterations     = ',j
       write (PRINTF, '(a,e12.6)') '         2-norm of the residual         = ',rnorm
       write (PRINTF, '(a,e12.6)') '         required accuracy              = ',epslin
       !
    else if ( iamout == 3 ) then
       !
       write (PRINTF, '(a,e12.6)') ' ++ pcg: 2-norm of the initial residual = ',rnrm0
       write (PRINTF, '(a,i3)'   ) '         total number of iterations     = ',j
       write (PRINTF, '(a,e12.6)') '         2-norm of the residual         = ',rnorm
       !
    endif
    !
end subroutine pcgu2
!
subroutine iluu
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
!    1.00, February 2023: New subroutine
!
!   Purpose
!
!   Computes classical incomplete factorization
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
!   This incomplete decomposition is used as the ILU preconditioner. If the last rule (c) is replaced by
!
!   rowsum(LU) = rowsum(A)
!
!   the MILU preconditioner is obtained. We also used a RILU preconditioner, which is the average of ILU and MILU.
!   The contributions on the diagonal are weighted by the scalar amod: amod = 0 corresponds with ILU, whereas
!   amod = 1 corresponds with MILU.
!
!   Auke van der Ploeg (1994)
!   Preconditioning for sparse matrices with applications
!   PhD thesis, Delft University of Technology, Delft, the Netherlands
!
!   Modules used
!
    use ocpcomm4
    use SwashCommdata3, only: pnums
    use SwashSolvedata
    use SwanGriddata, only: ncells
!
    implicit none
!
!   Local variables
!
    integer       :: icell    ! loop counter over cells
    integer, save :: ient = 0 ! number of entries in this subroutine
    integer       :: j        ! counter
    integer       :: k        ! counter
    integer       :: m        ! loop counter
    integer       :: n        ! loop counter
    !
    real          :: amod     ! parameter used in the average of ILU and MILU preconditioners
!
!   Structure
!
!   Description of the pseudo code
!
!   Source text
!
    if (ltrace) call strace (ient,'iluu')
    !
    amod = pnums(27)
    !
    prec2(1:) = ax2
    iplu      = 0
    !
    do icell = 2, ncells
       !
       prec2(0) = 0.
       !
       do m = ia(icell), ia(icell+1)-1
          iplu(ja(m)) = m
       enddo
       !
       do m = ia(icell), di(icell)-1
          !
          j = ja(m)
          !
          prec2(m) = prec2(m) / prec2(di(j))
          !
          do n = di(j)+1, ia(j+1)-1
             k = iplu(ja(n))
             prec2(k) = prec2(k) - prec2(m) * prec2(n)
          enddo
          !
       enddo
       !
       prec2(di(icell)) = prec2(di(icell)) + amod * prec2(0)
       !
       do m = ia(icell), ia(icell+1)-1
          iplu(ja(m)) = 0
       enddo
       !
    enddo
    !
end subroutine iluu
!
subroutine iludu
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
!    1.00, February 2023: New subroutine
!
!   Purpose
!
!   Computes incomplete factorization restricted to diagonal
!
!   Method
!
!   An incomplete decomposition is obtained using the following rules (Pommerell, 1992):
!
!   A = LU - R
!
!   (a)  diag(L) = diag(U) = I,
!   (b)  the off-diagonal parts of L and U are equal to the corresponding parts of A
!   (c)  diag(LU) = diag(A)
!
!   Modules used
!
    use ocpcomm4
    use SwashSolvedata
    use SwanGriddata, only: ncells
!
    implicit none
!
!   Local variables
!
    integer       :: icell    ! loop counter over cells
    integer, save :: ient = 0 ! number of entries in this subroutine
    integer       :: j        ! loop counter
    integer       :: k        ! loop counter
    !
    real          :: rval     ! auxiliary real
    !
    logical       :: found    ! a nonzero matrix element found
!
!   Structure
!
!   Description of the pseudo code
!
!   Source text
!
    if (ltrace) call strace (ient,'iludu')
    !
    do icell = 1, ncells
       !
       diag(icell) = ax2(di(icell))
       !
    enddo
    !
    do icell = 1, ncells
       !
       diag(icell) = 1. / diag(icell)
       !
       do j = di(icell)+1, ia(icell+1)-1
          !
          found = .false.
          !
          do k = ia(ja(j)), di(ja(j))-1
             !
             if ( ja(k) == icell ) then
                !
                found = .true.
                rval  = ax2(k)
                !
             endif
             !
          enddo
          !
          if ( found ) diag(ja(j)) = diag(ja(j)) - rval * diag(icell) * ax2(j)
          !
       enddo
       !
    enddo
    !
end subroutine iludu
!
subroutine bicgstabu ( amat, rhs, x )
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
!    1.00, February 2020: New subroutine
!
!   Purpose
!
!   Solves preconditioned system of equations by means of BiCGSTAB method
!
!   Method
!
!   System of equations assembled on triangular mesh
!
!   System of equations in single precision
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
    use SwashCommdata3, only: pnums, lprecon
    use SwashSolvedata
    use SwanGriddata
!
    implicit none
!
!   Argument variables
!
    real, dimension(ncells,0:3), intent(in   ) :: amat ! coefficient matrix of the system of equations
    real, dimension(ncells    ), intent(in   ) :: rhs  ! right-hand side of the system of equations
    real, dimension(ncells    ), intent(inout) :: x    ! solution of the system of equations
!
!   Parameter variables
!
    double precision, parameter :: epsmac = 1d-18 ! machine precision number
!
!   Local variables
!
    integer            :: icell    ! loop counter over cells
    integer, save      :: ient = 0 ! number of entries in this subroutine
    integer            :: ip       ! row pointer
    integer            :: j        ! iteration counter
    integer            :: jc       ! loop counter
    integer            :: k        ! counter
    integer            :: maxit    ! maximum number of iterations
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
    real               :: rval     ! auxiliary real
    real               :: sigma    ! inner product of quasi-residual vector and Ap
    real               :: tnorm    ! 2-norm of auxiliary vector t
    real               :: unorm    ! 2-norm of auxiliary vector u
    real               :: ueps     ! minimal accuracy based on machine precision
!
!   Structure
!
!   Description of the pseudo code
!
!   Source text
!
    if (ltrace) call strace (ient,'bicgstabu')
    !
    reps(1) = pnums(22)
    reps(2) = pnums(23)
    maxit   = nint(pnums(25))
    !
    ! restore matrix using the CSR format
    !
    do icell = 1, ncells
       !
       ip = ia(icell)
       !
       do jc = 1, 4
          !
          k = irw(icell,jc)
          !
          if ( k /= 0 ) ax2(ip+jc-1) = amat(icell,k-1)
          !
       enddo
       !
    enddo
    !
    ! compute incomplete LU factorization
    !
    if ( lprecon ) then
       !
       if ( icond == 2 ) then
          !
          ! ILU restricted to diagonal (Pommerell, 1992)
          !
          call iludu
          !
       else if ( icond == 3 ) then
          !
          ! ILU(0)
          !
          call iluu
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
    p  = 0.
    v  = 0.
    w2 = 0.
    !
    ! compute initial residual vector and its 2-norm
    !
    do icell = 1, ncells
       !
       res(icell,1) = 0.
       !
       do jc = ia(icell), ia(icell+1)-1
          res(icell,1) = res(icell,1) + ax2(jc) * x(ja(jc))
       enddo
       !
    enddo
    !
    rnorm = 0.
    bnorm = 0.
    !
    do icell = 1, ncells
       !
       res (icell,1) = rhs(icell  ) - res(icell,1)
       res0(icell,1) = res(icell,1)
       rnorm = rnorm + res(icell,1)*res(icell,1)
       bnorm = bnorm + rhs(icell  )*rhs(icell  )
       !
    enddo
    !
    rnorm = sqrt(rnorm)
    bnorm = sqrt(bnorm)
    rnrm0 = rnorm
    !
    if ( iamout == 2 ) then
       write (PRINTF, '(a,i3,a,e12.6)') ' ++ bicgstab: iter = ',j,'    res = ',rnorm
    endif
    !
    epslin = reps(1)*bnorm + reps(2)*rnrm0
    ueps   = 1000.*real(epsmac)*bnorm
    !
    if ( epslin < ueps .and. rnorm > 0. ) then
       !
       if ( iamout > 0 ) then
          write (PRINTF, '(a)') ' ++ bicgstab: the required accuracy is too small'
          write (PRINTF, '(a,e12.6)') '              required accuracy    = ',epslin
          write (PRINTF, '(a,e12.6)') '              appropriate accuracy = ',ueps
       endif
       !
       epslin = ueps
       !
    endif
    !
 10 if ( rnorm > epslin .and. j < maxit ) then
       !
       j = j + 1
       !
       beta = alpha / (rho*gamma)
       !
       rho = 0.
       !
       do icell = 1, ncells
          !
          rho = rho + res0(icell,1)*res(icell,1)
          !
       enddo
       !
       beta = beta * rho
       !
       do icell = 1, ncells
          !
          p(icell,1) = res(icell,1) + beta * ( p(icell,1) - gamma*v(icell,1) )
          !
       enddo
       !
       ! solve LU w = p
       !
       if ( icond == 2 ) then
          !
          ! elimination (lower triangular matrix)
          !
          do icell = 1, ncells
             !
             rval = 0.
             !
             do jc = ia(icell), di(icell)-1
                !
                rval = rval + ax2(jc) * w2(ja(jc))
                !
             enddo
             !
             w2(icell) = diag(icell) * ( p(icell,1) - rval )
             !
          enddo
          !
          ! substitution (upper triangular matrix)
          !
          do icell = ncells, 1, -1
             !
             rval = 0.
             !
             do jc = di(icell)+1, ia(icell+1)-1
                !
                rval = rval + ax2(jc) * w2(ja(jc))
                !
             enddo
             !
             w2(icell) = w2(icell) - diag(icell) * rval
             !
          enddo
          !
       else if ( icond == 3 ) then
          !
          ! elimination (lower triangular matrix)
          !
          do icell = 1, ncells
             !
             w2(icell) = p(icell,1)
             !
             do jc = ia(icell), di(icell)-1
                !
                w2(icell) = w2(icell) - prec2(jc) * w2(ja(jc))
                !
             enddo
             !
          enddo
          !
          ! substitution (upper triangular matrix)
          !
          do icell = ncells, 1, -1
             !
             do jc = di(icell)+1, ia(icell+1)-1
                !
                w2(icell) = w2(icell) - prec2(jc) * w2(ja(jc))
                !
             enddo
             !
             w2(icell) = w2(icell) / prec2(di(icell))
             !
          enddo
          !
       endif
       !
       do icell = 1, ncells
          !
          v(icell,1) = 0.
          !
          do jc = ia(icell), ia(icell+1)-1
             v(icell,1) = v(icell,1) + ax2(jc) * w2(ja(jc))
          enddo
          !
       enddo
       !
       sigma = 0.
       !
       do icell = 1, ncells
          !
          sigma = sigma + res0(icell,1)*v(icell,1)
          !
       enddo
       !
       if ( .not. sigma /= 0. ) then
          !
          if ( iamout > 0 ) then
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
       do icell = 1, ncells
          !
          u(icell,1) = res(icell,1) - alpha * v (icell,1)
          x(icell  ) = x  (icell  ) + alpha * w2(icell)
          unorm = unorm + u(icell,1)*u(icell,1)
          !
       enddo
       !
       unorm = sqrt(unorm)
       if ( unorm < epslin ) goto 20
       !
       ! solve LU w = u
       !
       if ( icond == 2 ) then
          !
          ! elimination (lower triangular matrix)
          !
          do icell = 1, ncells
             !
             rval = 0.
             !
             do jc = ia(icell), di(icell)-1
                !
                rval = rval + ax2(jc) * w2(ja(jc))
                !
             enddo
             !
             w2(icell) = diag(icell) * ( u(icell,1) - rval )
             !
          enddo
          !
          ! substitution (upper triangular matrix)
          !
          do icell = ncells, 1, -1
             !
             rval = 0.
             !
             do jc = di(icell)+1, ia(icell+1)-1
                !
                rval = rval + ax2(jc) * w2(ja(jc))
                !
             enddo
             !
             w2(icell) = w2(icell) - diag(icell) * rval
             !
          enddo
          !
       else if ( icond == 3 ) then
          !
          ! elimination (lower triangular matrix)
          !
          do icell = 1, ncells
             !
             w2(icell) = u(icell,1)
             !
             do jc = ia(icell), di(icell)-1
                !
                w2(icell) = w2(icell) - prec2(jc) * w2(ja(jc))
                !
             enddo
             !
          enddo
          !
          ! substitution (upper triangular matrix)
          !
          do icell = ncells, 1, -1
             !
             do jc = di(icell)+1, ia(icell+1)-1
                !
                w2(icell) = w2(icell) - prec2(jc) * w2(ja(jc))
                !
             enddo
             !
             w2(icell) = w2(icell) / prec2(di(icell))
             !
          enddo
          !
       endif
       !
       do icell = 1, ncells
          !
          t(icell,1) = 0.
          !
          do jc = ia(icell), ia(icell+1)-1
             t(icell,1) = t(icell,1) + ax2(jc) * w2(ja(jc))
          enddo
          !
       enddo
       !
       if ( unorm /= 0. ) then
          !
          gamma = 0.
          tnorm = 0.
          !
          do icell = 1, ncells
             !
             gamma = gamma + t(icell,1)*u(icell,1)
             tnorm = tnorm + t(icell,1)*t(icell,1)
             !
          enddo
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
          if ( iamout > 0 ) then
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
       do icell = 1, ncells
          !
          x  (icell  ) = x(icell  ) + gamma * w2(icell)
          res(icell,1) = u(icell,1) - gamma * t(icell,1)
          rnorm = rnorm + res(icell,1)*res(icell,1)
          !
       enddo
       rnorm = sqrt(rnorm)
       !
       if ( iamout == 2 ) then
          write (PRINTF,'(a,i3,a,e12.6)') ' ++ bicgstab: iter = ',j,'    res = ',rnorm
       endif
       !
       goto 10
       !
    endif
    !
 20 continue
    !
    ! investigate the reason for stopping
    !
    if ( j > 0 ) rnorm = min(rnorm,unorm)
    !
    if ( rnorm > epslin .and. iamout > 0 ) then
       !
       write (PRINTF, '(a)') ' ++ bicgstab: no convergence for solving system of equations'
       write (PRINTF, '(a,i3)'   ) '              total number of iterations     = ',j
       write (PRINTF, '(a,e12.6)') '              2-norm of the residual         = ',rnorm
       write (PRINTF, '(a,e12.6)') '              required accuracy              = ',epslin
       !
    else if ( iamout == 3 ) then
       !
       write (PRINTF, '(a,e12.6)') ' ++ bicgstab: 2-norm of the initial residual = ',rnrm0
       write (PRINTF, '(a,i3)'   ) '              total number of iterations     = ',j
       write (PRINTF, '(a,e12.6)') '              2-norm of the residual         = ',rnorm
       !
    endif
    !
end subroutine bicgstabu
!
subroutine bicgstab3 ( amat, rhs, x )
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
!    1.00, February 2020: New subroutine
!
!   Purpose
!
!   Solves multiple layer system of equations by means of preconditioned BiCGSTAB method
!
!   Method
!
!   System of equations assembled on triangular mesh
!
!   System of equations in single precision and in CSR format
!
!   An upper triangular matrix U and a lower triangular matrix L are computed according to the following rules:
!
!   A = LU - R
!
!   (a)  diag(L) = I,
!   (b)  the non-zero pattern of L and U are equal to the non-zero pattern of A
!   (c)  if A(i,j) <> 0 then L*U(i,j) = A(i,j)
!
!   This incomplete decomposition is used as the ILU preconditioner. If the last rule (c) is replaced by
!
!   rowsum(LU) = rowsum(A)
!
!   the MILU preconditioner is obtained. We also used a RILU preconditioner, which is the average of ILU and MILU.
!   The contributions on the diagonal are weighted by the scalar amod: amod = 0 corresponds with ILU, whereas
!   amod = 1 corresponds with MILU.
!
!   Auke van der Ploeg (1994)
!   Preconditioning for sparse matrices with applications
!   PhD thesis, Delft University of Technology, Delft, the Netherlands
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
    use SwashCommdata3, only: qmax, pnums, lprecon
    use SwashSolvedata
    use SwashFlowdata, only: nconct, ishif
    use SwanGriddata
!
    implicit none
!
!   Argument variables
!
    real, dimension(ncells,qmax,0:nconct), intent(in   ) :: amat ! coefficient matrix of the system of equations
    real, dimension(ncells,qmax         ), intent(in   ) :: rhs  ! right-hand side of the system of equations
    real, dimension(ncells,qmax         ), intent(inout) :: x    ! solution of the system of equations
!
!   Parameter variables
!
    double precision, parameter :: epsmac = 1d-18 ! machine precision number
!
!   Local variables
!
    integer                 :: icell    ! loop counter over cells
    integer, save           :: ient = 0 ! number of entries in this subroutine
    integer, dimension(0:3) :: irow     ! order of matrix columns in each matrix row
    integer                 :: j        ! iteration counter
    integer                 :: jc       ! counter
    integer                 :: jf       ! loop counter
    integer                 :: k        ! loop counter over vertical layers
    integer                 :: m        ! counter
    integer                 :: n        ! counter
    integer                 :: maxit    ! maximum number of iterations
    !
    real                    :: alpha    ! factor in the BiCGSTAB algorithm, i.e. rho/sigma or (res0,r)/(res0,Ap)
    real                    :: amod     ! parameter used in the average of ILU and MILU preconditioners
    real                    :: beta     ! factor in the BiCGSTAB algorithm, i.e. alpha/gamma
    real                    :: bnorm    ! 2-norm of right-hand side vector
    real                    :: crel     ! relaxation parameter
    real                    :: epslin   ! required accuracy in the linear solver
    real                    :: gamma    ! factor in the BiCGSTAB algorithm
    real, dimension(2)      :: reps     ! accuracies with respect to the right-hand side and initial residual used in the following termination criterion:
                                        !
                                        !  ||b-Ax || < reps(1)*||b|| + reps(2)*||b-Ax ||
                                        !        j                                   0
                                        !
    real                    :: rho      ! inner product of quasi-residual vector and residual vector
    real                    :: rnorm    ! 2-norm of residual vector
    real                    :: rnrm0    ! 2-norm of initial residual vector
    real                    :: rval     ! auxiliary real
    real                    :: sigma    ! inner product of quasi-residual vector and Ap
    real                    :: tnorm    ! 2-norm of auxiliary vector t
    real                    :: unorm    ! 2-norm of auxiliary vector u
    real                    :: ueps     ! minimal accuracy based on machine precision
!
!   Structure
!
!   Description of the pseudo code
!
!   Source text
!
    if (ltrace) call strace (ient,'bicgstab3')
    !
    ! set the parameters for the iterative solver
    !
    reps(1) = pnums(22)
    reps(2) = pnums(23)
    maxit   = nint(pnums(25))
    amod    = pnums(27)
    !
    ! restore matrix using the CSR format
    !
    do k = 1, qmax
       do icell = 1, ncells
          !
          m = ia(icell+(k-1)*ncells)
          !
          n = 0
          !
          do jf = 1, 4
             !
             jc = irw(icell,jf)
             !
             if ( jc /= 0 ) then
                irow(n) = jc - 1
                n = n + 1
             endif
             !
          enddo
          !
          if ( k > 2 ) then
             !
             do jf = 0, n-1
                !
                jc = irow(jf)
                !
                ax2(m+jf) = amat(icell,k,4+jc)
                !
             enddo
             !
             m = m + n
             !
          endif
          !
          if ( k > 1 ) then
             !
             do jf = 0, n-1
                !
                jc = irow(jf)
                !
                ax2(m+jf) = amat(icell,k,8+jc)
                !
             enddo
             !
             m = m + n
             !
          endif
          !
          do jf = 0, n-1
             !
             jc = irow(jf)
             !
             ax2(m+jf) = amat(icell,k,jc)
             !
          enddo
          !
          m = m + n
          !
          if ( k < qmax ) then
             !
             do jf = 0, n-1
                !
                jc = irow(jf)
                !
                ax2(m+jf) = amat(icell,k,12+jc)
                !
             enddo
             !
             m = m + n
             !
          endif
          !
          if ( k < qmax-1 ) then
             !
             do jf = 0, n-1
                !
                jc = irow(jf)
                !
                ax2(m+jf) = amat(icell,k,16+jc)
                !
             enddo
             !
             m = m + n
             !
          endif
          !
          do j = 3, nconct-17
             !
             if ( k+j < qmax+1 ) then
                !
                ax2(m) = amat(icell,k,ishif(j))
                !
                m = m + 1
                !
             endif
             !
          enddo
          !
       enddo
    enddo
    !
    ! compute incomplete LU factorization
    !
    if ( lprecon ) then
       !
       prec4(1:) = ax2
       iplu3     = 0
       !
       do k = 1, qmax
          do icell = 1, ncells
             !
             m = icell + (k-1)*ncells
             !
             if ( m == 1 ) cycle
             !
             prec4(0) = 0.
             !
             do jc = ia(m), ia(m+1)-1
                iplu3(ja(jc)) = jc
             enddo
             !
             do jc = ia(m), di(m)-1
                !
                j = ja(jc)
                !
                prec4(jc) = prec4(jc) / prec4(di(j))
                !
                do jf = di(j)+1, ia(j+1)-1
                   n = iplu3(ja(jf))
                   prec4(n) = prec4(n) - prec4(jc) * prec4(jf)
                enddo
                !
             enddo
             !
             prec4(di(m)) = prec4(di(m)) + amod * prec4(0)
             !
             do jc = ia(m), ia(m+1)-1
                iplu3(ja(jc)) = 0
             enddo
             !
          enddo
       enddo
       !
    endif
    !
    do k = 1, qmax
       do icell = 1, ncells
          m = icell + (k-1)*ncells
          sol3(m) = x(icell,k)
       enddo
    enddo
    !
    j     = 0
    rho   = 1.
    alpha = 1.
    gamma = 1.
    crel  = 0.7
    !
    ! initialize some arrays
    !
    p  = 0.
    v  = 0.
    w2 = 0.
    !
    ! compute initial residual vector and its 2-norm
    !
    do k = 1, qmax
       do icell = 1, ncells
          !
          m = icell + (k-1)*ncells
          !
          res(icell,k) = 0.
          !
          do jc = ia(m), ia(m+1)-1
             res(icell,k) = res(icell,k) + ax2(jc) * sol3(ja(jc))
          enddo
          !
       enddo
    enddo
    !
    rnorm = 0.
    bnorm = 0.
    !
    do k = 1, qmax
       do icell = 1, ncells
          !
          res (icell,k) = rhs(icell,k) - res(icell,k)
          res0(icell,k) = res(icell,k)
          rnorm = rnorm + res(icell,k)*res(icell,k)
          bnorm = bnorm + rhs(icell,k)*rhs(icell,k)
          !
       enddo
    enddo
    !
    rnorm = sqrt(rnorm)
    bnorm = sqrt(bnorm)
    rnrm0 = rnorm
    !
    if ( iamout == 2 ) then
       write (PRINTF, '(a,i3,a,e12.6)') ' ++ bicgstab: iter = ',j,'    res = ',rnorm
    endif
    !
    epslin = reps(1)*bnorm + reps(2)*rnrm0
    ueps   = 1000.*real(epsmac)*bnorm
    !
    if ( epslin < ueps .and. rnorm > 0. ) then
       !
       if ( iamout > 0 ) then
          write (PRINTF, '(a)') ' ++ bicgstab: the required accuracy is too small'
          write (PRINTF, '(a,e12.6)') '              required accuracy    = ',epslin
          write (PRINTF, '(a,e12.6)') '              appropriate accuracy = ',ueps
       endif
       !
       epslin = ueps
       !
    endif
    !
 10 if ( rnorm > epslin .and. j < maxit ) then
       !
       j = j + 1
       !
       beta = alpha / (rho*gamma)
       !
       rho = 0.
       !
       do k = 1, qmax
          do icell = 1, ncells
             !
             rho = rho + res0(icell,k)*res(icell,k)
             !
          enddo
       enddo
       !
       beta = beta * rho
       !
       do k = 1, qmax
          do icell = 1, ncells
             !
             p(icell,k) = res(icell,k) + beta * ( p(icell,k) - gamma*v(icell,k) )
             !
          enddo
       enddo
       !
       ! solve LU w = p
       !
       ! elimination (lower triangular matrix)
       !
       do k = 1, qmax
          do icell = 1, ncells
             !
             m = icell + (k-1)*ncells
             !
             w2(m) = p(icell,k)
             !
             do jc = ia(m), di(m)-1
                !
                w2(m) = w2(m) - prec4(jc) * w2(ja(jc))
                !
             enddo
             !
          enddo
       enddo
       !
       ! substitution (upper triangular matrix)
       !
       do k = qmax, 1, -1
          do icell = ncells, 1, -1
             !
             m = icell + (k-1)*ncells
             !
             do jc = di(m)+1, ia(m+1)-1
                !
                w2(m) = w2(m) - prec4(jc) * w2(ja(jc))
                !
             enddo
             !
             w2(m) = w2(m) / prec4(di(m))
             !
          enddo
       enddo
       !
       do k = 1, qmax
          do icell = 1, ncells
             !
             m = icell + (k-1)*ncells
             !
             v(icell,k) = 0.
             !
             do jc = ia(m), ia(m+1)-1
                v(icell,k) = v(icell,k) + ax2(jc) * w2(ja(jc))
             enddo
             !
          enddo
       enddo
       !
       sigma = 0.
       !
       do k = 1, qmax
          do icell = 1, ncells
             !
             sigma = sigma + res0(icell,k)*v(icell,k)
             !
          enddo
       enddo
       !
       if ( .not. sigma /= 0. ) then
          !
          if ( iamout > 0 ) then
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
       do k = 1, qmax
          do icell = 1, ncells
             !
             m = icell + (k-1)*ncells
             !
             u(icell,k) = res(icell,k) - alpha * v(icell,k)
             sol3(m)    = sol3(m)      + alpha * w2(m)
             unorm = unorm + u(icell,k)*u(icell,k)
             !
          enddo
       enddo
       !
       unorm = sqrt(unorm)
       if ( unorm < epslin ) goto 20
       !
       ! solve LU w = u
       !
       ! elimination (lower triangular matrix)
       !
       do k = 1, qmax
          do icell = 1, ncells
             !
             m = icell + (k-1)*ncells
             !
             w2(m) = u(icell,k)
             !
             do jc = ia(m), di(m)-1
                !
                w2(m) = w2(m) - prec4(jc) * w2(ja(jc))
                !
             enddo
             !
          enddo
       enddo
       !
       ! substitution (upper triangular matrix)
       !
       do k = qmax, 1, -1
          do icell = ncells, 1, -1
             !
             m = icell + (k-1)*ncells
             !
             do jc = di(m)+1, ia(m+1)-1
                !
                w2(m) = w2(m) - prec4(jc) * w2(ja(jc))
                !
             enddo
             !
             w2(m) = w2(m) / prec4(di(m))
             !
          enddo
       enddo
       !
       do k = 1, qmax
          do icell = 1, ncells
             !
             m = icell + (k-1)*ncells
             !
             t(icell,k) = 0.
             !
             do jc = ia(m), ia(m+1)-1
                t(icell,k) = t(icell,k) + ax2(jc) * w2(ja(jc))
             enddo
             !
          enddo
       enddo
       !
       if ( unorm /= 0. ) then
          !
          gamma = 0.
          tnorm = 0.
          !
          do k = 1, qmax
             do icell = 1, ncells
                !
                gamma = gamma + t(icell,k)*u(icell,k)
                tnorm = tnorm + t(icell,k)*t(icell,k)
                !
             enddo
          enddo
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
          if ( iamout > 0 ) then
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
       do k = 1, qmax
          do icell = 1, ncells
             !
             m = icell + (k-1)*ncells
             !
             sol3(m)      = sol3(m)    + gamma * w2(m)
             res(icell,k) = u(icell,k) - gamma * t(icell,k)
             rnorm = rnorm + res(icell,k)*res(icell,k)
             !
          enddo
       enddo
       rnorm = sqrt(rnorm)
       !
       if ( iamout == 2 ) then
          write (PRINTF,'(a,i3,a,e12.6)') ' ++ bicgstab: iter = ',j,'    res = ',rnorm
       endif
       !
       goto 10
       !
    endif
    !
 20 continue
    !
    do k = 1, qmax
       do icell = 1, ncells
          m = icell + (k-1)*ncells
          x(icell,k) = sol3(m)
       enddo
    enddo
    !
    ! investigate the reason for stopping
    !
    if ( j > 0 ) rnorm = min(rnorm,unorm)
    !
    if ( rnorm > epslin .and. iamout > 0 ) then
       !
       write (PRINTF, '(a)') ' ++ bicgstab: no convergence for solving system of equations'
       write (PRINTF, '(a,i3)'   ) '              total number of iterations     = ',j
       write (PRINTF, '(a,e12.6)') '              2-norm of the residual         = ',rnorm
       write (PRINTF, '(a,e12.6)') '              required accuracy              = ',epslin
       !
    else if ( iamout == 3 ) then
       !
       write (PRINTF, '(a,e12.6)') ' ++ bicgstab: 2-norm of the initial residual = ',rnrm0
       write (PRINTF, '(a,i3)'   ) '              total number of iterations     = ',j
       write (PRINTF, '(a,e12.6)') '              2-norm of the residual         = ',rnorm
       !
    endif
    !
end subroutine bicgstab3
!
subroutine newtonU ( amat, rhs, x )
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
!    1.00, April 2023: New subroutine
!
!   Purpose
!
!   Solves a piecewise linear system by means of the Newton-type iteration method
!
!   Method
!
!   The piecewise linear system is of the following form
!
!   p(x)x + Ax = b
!
!   with p a Heaviside function, as follows
!
!           | 1 if x >= 0
!   p(x) = <
!           | 0 if x  < 0
!
!   and A a symmetric and positive semidefinite matrix
!
!   V. Casulli
!   A high-resolution wetting and drying algorithm for free-surface hydrodynamics
!   IJNMF, vol. 60, 391-408, 2009
!
!   Modules used
!
    use ocpcomm4
    use SwashCommdata3, only: epsdry, pnums
    use SwashFlowdata, only: hs
    use SwashSolvedata, only: iamout, amata, rhsa, sol0, sol1, vt, resd2
    use SwanGriddata
    use SwanGridobjects
!
    implicit none
!
!   Argument variables
!
    real*8, dimension(ncells,0:3), intent(in ) :: amat ! coefficient matrix of the system of equations
    real*8, dimension(ncells    ), intent(in ) :: rhs  ! right-hand side of the system of equations
    real  , dimension(ncells    ), intent(out) :: x    ! solution of the system of equations
!
!   Parameter variables
!
    real*8, parameter :: epsmac = 1d-18 ! machine precision number
!
!   Local variables
!
    integer       :: icell    ! loop counter over cells
    integer       :: icelll   ! left cell of present face
    integer       :: icellr   ! right cell of present face
    integer, save :: ient = 0 ! number of entries in this subroutine
    integer       :: iface    ! face index
    integer       :: j        ! iteration counter
    integer       :: jf       ! loop counter
    integer       :: maxit    ! maximum number of iterations
    !
    real          :: area     ! area of present cell
    !
    real*8        :: a        ! cell area
    real*8        :: epslin   ! required accuracy for Newton-type iteration
    real*8        :: h        ! water depth
    real*8        :: p        ! derivative of v (=area)
    real*8        :: reps     ! convergence criterion for matrix solver
    real*8        :: rnorm    ! 2-norm of residual vector
    real*8        :: rnrm0    ! 2-norm of initial residual vector
    real*8        :: s        ! solution (= difference in water level over time step)
    real*8        :: ueps     ! minimal accuracy based on machine precision
    real*8        :: v        ! cell volume
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
    if (ltrace) call strace (ient,'newtonU')
    !
    ! point to cell and face objects
    !
    cell => gridobject%cell_grid
    face => gridobject%face_grid
    !
    reps  = dble(pnums(21))
    maxit = nint(pnums(24))
    !
    ! copy original coefficient matrix
    !
    amata = amat
    !
    ! initialize solution
    !
    sol1 = 0d0
    !
    ! compute initial residual
    !
    rnrm0 = 0d0
    do icell = 1, ncells
       rnrm0 = rnrm0 + rhs(icell)*rhs(icell)
    enddo
    rnrm0 = sqrt(rnrm0)
    !
    epslin = reps*rnrm0
    ueps   = 1d3*epsmac*rnrm0
    !
    if ( epslin < ueps .and. rnrm0 > 0d0 ) then
       !
       if ( iamout > 0 ) then
          write (PRINTF, '(a)') ' ++ Newton: the required accuracy is too small'
          write (PRINTF, '(a,e12.6)') '            required accuracy    = ',epslin
          write (PRINTF, '(a,e12.6)') '            appropriate accuracy = ',ueps
       endif
       !
       epslin = ueps
       !
    endif
    !
    iterloop: do j = 1, maxit
       !
       sol0 = sol1
       !
       ! set main diagonal and right-hand side due to Newton
       !
       do icell = 1, ncells
          !
          area = cell(icell)%attr(CELLAREA)
          a    = dble(area)
          !
          s = sol0(icell)
          h = dble(hs(icell))
          !
          if ( s < -h + epsdry ) then
             !
             ! water level drops below bottom level
             !
             v = a * epsdry
             p = 0d0
             !
          else
             !
             ! volume cell is non-negative
             !
             v = a * h
             p = a
             !
          endif
          !
          amata(icell,0) = amat(icell,0) + p
          rhsa (icell  ) = rhs (icell  ) - v
          !
       enddo
       !
       ! solve system of equations by means of preconditioned CG method
       !
       call pcgu2 ( amata, rhsa, sol1 )
       !
       ! compute residual
       !
       do icell = 1, ncells
          !
          resd2(icell) = amata(icell,0) * sol1(icell)
          !
          ! loop over faces of the cell
          !
          do jf = 1, cell(icell)%nof
             !
             ! face identifier
             !
             iface = cell(icell)%face(jf)%atti(FACEID)
             !
             ! consider left and right cells of current face
             !
             icelll = face(iface)%atti(FACECL)
             icellr = face(iface)%atti(FACECR)
             !
             if ( icell == icelll .and. icellr /= 0 ) then
                resd2(icell) = resd2(icell) + amata(icell,jf) * sol1(icellr)
             else if ( icell == icellr .and. icelll /= 0 ) then
                resd2(icell) = resd2(icell) + amata(icell,jf) * sol1(icelll)
             endif
             !
          enddo
          !
       enddo
       !
       resd2 = rhsa - resd2
       !
       rnorm = 0d0
       do icell = 1, ncells
          rnorm = rnorm + resd2(icell)*resd2(icell)
       enddo
       rnorm = sqrt(rnorm)
       !
       if ( iamout == 2 ) then
          write (PRINTF, '(a,i3,a,e12.6)') ' ++ Newton: iter = ',j,'    res = ',rnorm
       endif
       !
       ! check convergence
       !
       if ( .not. rnorm > epslin ) exit iterloop
       !
       !
    enddo iterloop
    !
    x = real(sol1)
    !
    ! investigate the reason for iteration stop
    !
    if ( rnorm > epslin .and. iamout > 0 ) then
       !
       write (PRINTF, '(a)') ' ++ Newton: no convergence for solving system of equations'
       write (PRINTF, '(a,i3)'   ) '            total number of iterations     = ',j
       write (PRINTF, '(a,e12.6)') '            2-norm of the residual         = ',rnorm
       write (PRINTF, '(a,e12.6)') '            required accuracy              = ',epslin
       !
    else if ( iamout == 3 ) then
       !
       write (PRINTF, '(a,e12.6)') ' ++ Newton: 2-norm of the initial residual = ',rnrm0
       write (PRINTF, '(a,i3)'   ) '            total number of iterations     = ',j
       write (PRINTF, '(a,e12.6)') '            2-norm of the residual         = ',rnorm
       !
    endif
    !
end subroutine newtonU
!
subroutine csrf ( ia, ja, irw, di, kmax, nod )
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
!    1.00, February 2023: New subroutine
!
!   Purpose
!
!   Sets up the Compressed Sparse Row format
!
!   Modules used
!
    use ocpcomm4
    use SwanGriddata
    use SwanGridobjects
!
    implicit none
!
!   Argument variables
!
    integer,                                 intent(in ) :: kmax ! actual number of vertical layers
    integer,                                 intent(in ) :: nod  ! actual number of off-diagonals in matrix row
!
    integer, dimension(ncells*kmax        ), intent(out) :: di  ! array containing the positions of the main diagonal elements per row
    integer, dimension(ncells*kmax+1      ), intent(out) :: ia  ! array containing the positions that start a row
    integer, dimension(ncells,4           ), intent(out) :: irw ! array containing the order of matrix columns per row
    integer, dimension((nod+1)*ncells*kmax), intent(out) :: ja  ! array containing the column indices of matrix elements
!
!   Local variables
!
    integer               :: icell    ! loop counter over cells
    integer               :: icelll   ! left cell of present face
    integer               :: icellr   ! right cell of present face
    integer, dimension(4) :: icol     ! column numbers in each matrix row
    integer, save         :: ient = 0 ! number of entries in this subroutine
    integer               :: iface    ! face index
    integer, dimension(4) :: irow     ! order of matrix columns in each matrix row
    integer               :: itmp     ! temporary stored integer for swapping
    integer               :: j        ! loop counter
    integer               :: jf       ! loop counter over faces of cell
    integer               :: k        ! loop counter over vertical layers
    integer               :: l        ! counter
    integer, dimension(1) :: lm       ! location of minimum value in array icol
    integer               :: m        ! counter
    integer               :: msiz     ! matrix size
    integer               :: n        ! counter
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
    if (ltrace) call strace (ient,'csrf')
    !
    ! point to cell and face objects
    !
    cell => gridobject%cell_grid
    face => gridobject%face_grid
    !
    msiz = ncells * kmax
    !
    ! initialization
    !
    ia  = 0
    ja  = 0
    irw = 0
    di  = 0
    !
    n = 0
    !
    do k = 1, kmax
       !
       do icell = 1, ncells
          !
          icol = ncells + 1
          !
          icol(1) = icell
          !
          m = 1
          !
          ! loop over faces of the cell
          !
          do jf = 1, cell(icell)%nof
             !
             ! face identifier
             !
             iface = cell(icell)%face(jf)%atti(FACEID)
             !
             ! consider left and right cells of current face
             !
             icelll = face(iface)%atti(FACECL)
             icellr = face(iface)%atti(FACECR)
             !
             if ( face(iface)%atti(FMARKER) == 0 ) then   ! internal face
                !
                if ( icell == icelll ) then
                   icol(1+jf) = icellr
                else if ( icell == icellr ) then
                   icol(1+jf) = icelll
                endif
                !
                m = m + 1
                !
             endif
             !
          enddo
          !
          ! sort column numbers
          !
          do j = 1, 4
             irow(j) = j
          enddo
          !
          if ( m > 1 ) then
             !
             do j = 1, 3
                !
                lm = minloc(icol(j:4))
                l  = lm(1) + j-1
                !
                if ( l /= j ) then
                   !
                   itmp    = icol(j)
                   icol(j) = icol(l)
                   icol(l) = itmp
                   !
                   itmp    = irow(j)
                   irow(j) = irow(l)
                   irow(l) = itmp
                   !
                endif
                !
             enddo
             !
          endif
          !
          ! fill CSR arrays
          !
          irw(icell,1:m) = irow(1:m)
          !
          ia(icell+(k-1)*ncells) = n + 1
          !
          if ( k > 2 ) then
             !
             ja(n+1:n+m) = icol(1:m) + (k-3)*ncells
             !
             n = n + m
             !
          endif
          !
          if ( k > 1 ) then
             !
             ja(n+1:n+m) = icol(1:m) + (k-2)*ncells
             !
             n = n + m
             !
          endif
          !
          ja(n+1:n+m) = icol(1:m) + (k-1)*ncells
          !
          do j = 1, m
             !
             if ( icol(j) == icell ) di(icell+(k-1)*ncells) = n + j
             !
          enddo
          !
          n = n + m
          !
          if ( k < kmax ) then
             !
             ja(n+1:n+m) = icol(1:m) + k*ncells
             !
             n = n + m
             !
          endif
          !
          if ( k < kmax-1 ) then
             !
             ja(n+1:n+m) = icol(1:m) + (k+1)*ncells
             !
             n = n + m
             !
          endif
          !
          do j = 3, nod-17
             !
             if ( k+j < kmax+1 ) then
                !
                n = n + 1
                !
                ja(n) = icell + (k+j-1)*ncells
                !
             endif
             !
          enddo
          !
       enddo
       !
    enddo
    !
    ia(msiz+1) = ia(1) + n
    !
end subroutine csrf
