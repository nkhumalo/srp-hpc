      submodule (wfn1) wfn1_4index
      implicit none
      contains
!>
!>----------------------------------------------------------------------
!>
!>    \brief 4-Index transform same spin
!>
!>    Do a 4-index transformation where all indeces are transformed
!>    with the natural orbitals of the same spin.
!>
!>    The algorithm is simply:
!>
!>    - g(a,c,b,d)*n(a,i) = g2(i,c)
!>    - g2(i,c)*n(c,j)    = gt(i,j,b,d)
!>
!>    - gt(i,j,b,d)*n(b,k) = g2(k,d)
!>    - g2(k,d)*n(d,l)     = go(i,j,k,l)
!>
!>    where each multiplication is a matrix-matrix multiplication.
!>    To make this algorithm work we have to choose:
!>
!>    - g\_4mo = go
!>    - g\_4t1 = gt
!>
      subroutine wfn1_4index_aa(wfn1_param,g_na,g_4ao,g_4mo,g_4t1,
     &                          g_bfbf)
      implicit none
!>    Constant parameters for this calculation
      type(wfn1_prm), intent(in) :: wfn1_param
!>    GA handle for the natural orbitals
      integer,        intent(in) :: g_na
!>    GA handle for the 2-electron integrals in AO basis
      integer,        intent(in) :: g_4ao
!>    GA handle for the 2-electron integrals in MO basis
      integer,        intent(in) :: g_4mo
!>    GA handle for a 4-dimensional temporary array
      integer,        intent(in) :: g_4t1
!>    GA handle for a 2-dimensional temporary array
      integer,        intent(in) :: g_bfbf
!>
!>    local variables
!>
      integer :: nbf ! number of basis functions
      integer :: nmo ! number of molecular orbitals
      integer :: ii  ! counter
      integer :: ij  ! counter
      integer :: ik  ! counter
      integer :: il  ! counter
      integer :: ia  ! counter
      integer :: ib  ! counter
      integer :: ic  ! counter
      integer :: id  ! counter
!>
      integer :: alo(4), ahi(4)
      integer :: blo(4), bhi(4)
      integer :: clo(4), chi(4)
!>
!>    Code
!>
      nbf = wfn1_param%nbf
      nmo = wfn1_param%nmo
      do id = 1, nbf
        do ib = 1, nbf
!>
          alo(1) = 1
          alo(2) = 1
          alo(3) = -1
          alo(4) = -1
          ahi(1) = nmo
          ahi(2) = nbf
          ahi(3) = -2
          ahi(4) = -2
!>
          blo(1) = 1
          blo(2) = 1
          blo(3) = ib
          blo(4) = id
          bhi(1) = nbf
          bhi(2) = nbf
          bhi(3) = ib
          bhi(4) = id
!>
          clo(1) = 1
          clo(2) = 1
          clo(3) = -1
          clo(4) = -1
          chi(1) = nmo
          chi(2) = nbf
          chi(3) = -2
          chi(4) = -2
!>
          call nga_matmul_patch('t','n',1.0d0,0.0d0,
     &                          g_na,  alo,ahi,
     &                          g_4ao, blo,bhi,
     &                          g_bfbf,clo,chi)
!>
          alo(1) = 1
          alo(2) = 1
          alo(3) = -1
          alo(4) = -1
          ahi(1) = nmo
          ahi(2) = nbf
          ahi(3) = -2
          ahi(4) = -2
!>
          blo(1) = 1
          blo(2) = 1
          blo(3) = -1
          blo(4) = -1
          bhi(1) = nbf
          bhi(2) = nmo
          bhi(3) = -2
          bhi(4) = -2
!>
          clo(1) = 1
          clo(2) = 1
          clo(3) = ib
          clo(4) = id
          chi(1) = nmo
          chi(2) = nmo
          chi(3) = ib
          chi(4) = id
!>
          call nga_matmul_patch('n','n',1.0d0,0.0d0,
     &                          g_bfbf,alo,ahi,
     &                          g_na,  blo,bhi,
     &                          g_4t1, clo,chi)
!>
        enddo
      enddo
!>
      do ij = 1, nmo
        do ii = 1, nmo
!>
          alo(1) = 1
          alo(2) = 1
          alo(3) = -1
          alo(4) = -1
          ahi(1) = nmo
          ahi(2) = nbf
          ahi(3) = -2
          ahi(4) = -2
!>
          blo(1) = ii
          blo(2) = ij
          blo(3) = 1
          blo(4) = 1
          bhi(1) = ii
          bhi(2) = ij
          bhi(3) = nbf
          bhi(4) = nbf
!>
          clo(1) = 1
          clo(2) = 1
          clo(3) = -1
          clo(4) = -1
          chi(1) = nmo
          chi(2) = nbf
          chi(3) = -2
          chi(4) = -2
!>
          call nga_matmul_patch('t','n',1.0d0,0.0d0,
     &                          g_na,  alo,ahi,
     &                          g_4t1, blo,bhi,
     &                          g_bfbf,clo,chi)
!>
          alo(1) = 1
          alo(2) = 1
          alo(3) = -1
          alo(4) = -1
          ahi(1) = nmo
          ahi(2) = nbf
          ahi(3) = -2
          ahi(4) = -2
!>
          blo(1) = 1
          blo(2) = 1
          blo(3) = -1
          blo(4) = -1
          bhi(1) = nbf
          bhi(2) = nmo
          bhi(3) = -2
          bhi(4) = -2
!>
          clo(1) = ii
          clo(2) = ij
          clo(3) = 1
          clo(4) = 1
          chi(1) = ii
          chi(2) = ij
          chi(3) = nmo
          chi(4) = nmo
!>
          call nga_matmul_patch('n','n',1.0d0,0.0d0,
     &                          g_bfbf,alo,ahi,
     &                          g_na,  blo,bhi,
     &                          g_4mo, clo,chi)
!>
        enddo
      enddo
!>
      end subroutine wfn1_4index_aa
!>
!>----------------------------------------------------------------------
!>
!>    \brief 4-Index transform opposite spin
!>
!>    Do a 4-index transformation where two indeces are transformed
!>    with the natural orbitals of the one spin, and the other two
!>    are transformed with the natural orbitals of the opposite spin.
!>
!>    The algorithm is simply:
!>
!>    - g(a,c,b,d)*na(a,i)  = g2(i,c)
!>    - g2(i,c)*na(c,j)     = gt(i,j,b,d)
!>
!>    - gt(i,j,b,d)*nb(b,k) = g2(k,d)
!>    - g2(k,d)*nb(d,l)     = go(i,j,k,l)
!>
!>    where each multiplication is a matrix-matrix multiplication.
!>    To make this algorithm work we have to choose:
!>
!>    - g\_4mo = go
!>    - g\_4t1 = gt
!>
      subroutine wfn1_4index_ab(wfn1_param,g_na,g_nb,g_4ao,g_4mo,g_4t1,
     &                          g_bfbf)
      implicit none
!>    Constant parameters for this calculation
      type(wfn1_prm), intent(in) :: wfn1_param
!>    GA handle for the natural orbitals of sigma spin
      integer,        intent(in) :: g_na
!>    GA handle for the natural orbitals of sigma' spin
      integer,        intent(in) :: g_nb
!>    GA handle for the 2-electron integrals in AO basis
      integer,        intent(in) :: g_4ao
!>    GA handle for the 2-electron integrals in MO basis
      integer,        intent(in) :: g_4mo
!>    GA handle for a 4-dimensional temporary array
      integer,        intent(in) :: g_4t1
!>    GA handle for a 2-dimensional temporary array
      integer,        intent(in) :: g_bfbf
!>
!>    local variables
!>
      integer :: nbf ! number of basis functions
      integer :: nmo ! number of molecular orbitals
      integer :: ii  ! counter
      integer :: ij  ! counter
      integer :: ik  ! counter
      integer :: il  ! counter
      integer :: ia  ! counter
      integer :: ib  ! counter
      integer :: ic  ! counter
      integer :: id  ! counter
!>
      integer :: alo(4), ahi(4)
      integer :: blo(4), bhi(4)
      integer :: clo(4), chi(4)
!>
!>    Code
!>
      nbf = wfn1_param%nbf
      nmo = wfn1_param%nmo
      do id = 1, nbf
        do ib = 1, nbf
!>
          alo(1) = 1
          alo(2) = 1
          alo(3) = -1
          alo(4) = -1
          ahi(1) = nmo
          ahi(2) = nbf
          ahi(3) = -2
          ahi(4) = -2
!>
          blo(1) = 1
          blo(2) = 1
          blo(3) = ib
          blo(4) = id
          bhi(1) = nbf
          bhi(2) = nbf
          bhi(3) = ib
          bhi(4) = id
!>
          clo(1) = 1
          clo(2) = 1
          clo(3) = -1
          clo(4) = -1
          chi(1) = nmo
          chi(2) = nbf
          chi(3) = -2
          chi(4) = -2
!>
          call nga_matmul_patch('t','n',1.0d0,0.0d0,
     &                          g_na,  alo,ahi,
     &                          g_4ao, blo,bhi,
     &                          g_bfbf,clo,chi)
!>
          alo(1) = 1
          alo(2) = 1
          alo(3) = -1
          alo(4) = -1
          ahi(1) = nmo
          ahi(2) = nbf
          ahi(3) = -2
          ahi(4) = -2
!>
          blo(1) = 1
          blo(2) = 1
          blo(3) = -1
          blo(4) = -1
          bhi(1) = nbf
          bhi(2) = nmo
          bhi(3) = -2
          bhi(4) = -2
!>
          clo(1) = 1
          clo(2) = 1
          clo(3) = ib
          clo(4) = id
          chi(1) = nmo
          chi(2) = nmo
          chi(3) = ib
          chi(4) = id
!>
          call nga_matmul_patch('n','n',1.0d0,0.0d0,
     &                          g_bfbf,alo,ahi,
     &                          g_na,  blo,bhi,
     &                          g_4t1, clo,chi)
!>
        enddo
      enddo
!>
      do ij = 1, nmo
        do ii = 1, nmo
!>
          alo(1) = 1
          alo(2) = 1
          alo(3) = -1
          alo(4) = -1
          ahi(1) = nmo
          ahi(2) = nbf
          ahi(3) = -2
          ahi(4) = -2
!>
          blo(1) = ii
          blo(2) = ij
          blo(3) = 1
          blo(4) = 1
          bhi(1) = ii
          bhi(2) = ij
          bhi(3) = nbf
          bhi(4) = nbf
!>
          clo(1) = 1
          clo(2) = 1
          clo(3) = -1
          clo(4) = -1
          chi(1) = nmo
          chi(2) = nbf
          chi(3) = -2
          chi(4) = -2
!>
          call nga_matmul_patch('t','n',1.0d0,0.0d0,
     &                          g_nb,  alo,ahi,
     &                          g_4t1, blo,bhi,
     &                          g_bfbf,clo,chi)
!>
          alo(1) = 1
          alo(2) = 1
          alo(3) = -1
          alo(4) = -1
          ahi(1) = nmo
          ahi(2) = nbf
          ahi(3) = -2
          ahi(4) = -2
!>
          blo(1) = 1
          blo(2) = 1
          blo(3) = -1
          blo(4) = -1
          bhi(1) = nbf
          bhi(2) = nmo
          bhi(3) = -2
          bhi(4) = -2
!>
          clo(1) = ii
          clo(2) = ij
          clo(3) = 1
          clo(4) = 1
          chi(1) = ii
          chi(2) = ij
          chi(3) = nmo
          chi(4) = nmo
!>
          call nga_matmul_patch('n','n',1.0d0,0.0d0,
     &                          g_bfbf,alo,ahi,
     &                          g_nb,  blo,bhi,
     &                          g_4mo, clo,chi)
!>
        enddo
      enddo
!>
      end subroutine wfn1_4index_ab
!>
!>----------------------------------------------------------------------
!>
      end submodule wfn1_4index
