      submodule (wfn1) wfn1_clc_ortho
      implicit none
      contains
c
c-----------------------------------------------------------------------
c
C>
C>    \brief Orthonormalize the vectors 
C>
C>    This subroutine ortho-normalizes all vectors making up the 
C>    wavefunction, i.e. the natural orbitals as well as the correlation
C>    functions. 
C>
      module subroutine wfn1_ortho(wfn1_mat,wfn1_wave)
      implicit none
!include "wfn1_matcnst.fh"
!include "wfn1_wfn.fh"
      type(wfn1_matcnst), intent(in) :: wfn1_mat ! matrix constants
      type(wfn1_wfn),  intent(inout) :: wfn1_wave ! the wave function
c
      integer :: idum  ! dummy variable
c
      idum = 0
      call wfn1_orthog(wfn1_wave%g_ca,idum,.false.)
      call wfn1_orthog(wfn1_wave%g_cb,idum,.false.)
      call wfn1_orthog(wfn1_wave%g_na,wfn1_mat%g_s,.true.)
      call wfn1_orthog(wfn1_wave%g_nb,wfn1_mat%g_s,.true.)
c
      end subroutine wfn1_ortho
c
c-----------------------------------------------------------------------
c
!> \brief Orthogonalize a set of complex or double vectors
!>
!> Orthogonalize a set of vectors in either an orthogonal basis or 
!> a non-orthogonal basis:
!> - ometric = .false. the vectors are given in an orthogonal basis
!>   and `g_over` is not used
!> - ometric = .true.  the vectors are given in a non-orthogonal bais
!>   in which `g_over` contains the overlap matrix
!> The vector may be real or complex:
!> - if `g_vecs` is real then g_over is also real
!> - if `g_vecs` is complex then g_over is also complex and Hermitian
!>
      module subroutine wfn1_orthog(g_vecs, g_over, ometric)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "msgids.fh"
      integer, intent(in) :: g_vecs, g_over
      logical, intent(in) :: ometric
c
      integer ulo, uhi, itypv, ityps, ni, nu, uchunk, ld
      integer mu, mv
      integer l_tmp, l_s, l_over, l_w
      integer k_tmp, k_s, k_over, k_w
c
      real(kind=8),    allocatable :: dvec(:,:)
      real(kind=8),    allocatable :: dovl(:,:)
      real(kind=8),    allocatable :: dsij(:)
      real(kind=8),    allocatable :: dw(:)
c
      complex(kind=8), allocatable :: zvec(:,:)
      complex(kind=8), allocatable :: zovl(:,:)
      complex(kind=8), allocatable :: zsij(:)
      complex(kind=8), allocatable :: zw(:)
c     
c     Redistribute the input matrix ... block the leading
c     dimension, leave second dimension undistributed
c     ... each process has a(ulo:uhi,1:ni).  If
c
c     nu = no. of basis functions
c     ni = no. of vectors
c
c     overlap should be (nu,nu)
c
      call ga_inquire(g_vecs, itypv, nu, ni)
      if (ometric) then
        call ga_inquire(g_over, ityps, mu, mv)
        if (itypv.ne.ityps) then
          call errquit('wfn1_orthog: vectors and overlap not of same '//
     &                 'type',UERR)
        endif
      endif
      uchunk = max(8, (nu-1)/ga_nnodes()+1)
      ulo  = ga_nodeid()*uchunk + 1
      uhi  = min(ulo + uchunk - 1, nu)
      ld   = uhi - ulo + 1
      if (ulo .gt. uhi) then
         ulo = 0
         uhi = -1
         ld  = 1
      end if
*      write(6,*) ' node ulo uhi ', ga_nodeid(), ulo, uhi
*      call util_flush(6)
*      call ga_sync
*      call ga_print(g_vecs)
*      if (ometric) call ga_print(g_over)
c
c     Allocate local workspace
c
      if      (itypv.eq.MT_DBL) then
        if (ulo.gt.0) then
          allocate(dvec(ulo:uhi,1:ni))
          if (ometric) then
            allocate(dovl(1:nu,ulo:uhi))
          endif
          allocate(dsij(1:ni))
          allocate(dw(1:nu))
        endif
      else if (itypv.eq.MT_DCPL) then
        if (ulo.gt.0) then
          allocate(zvec(ulo:uhi,1:ni))
          if (ometric) then
            allocate(zovl(1:nu,ulo:uhi))
          endif
          allocate(zsij(1:ni))
          allocate(zw(1:nu))
        endif
      else
        call errquit("wfn1_orthog: unsupported type",UERR)
      endif
c
      call ga_sync()
      if (uhi.ge.ulo) then
        if      (itypv.eq.MT_DBL) then
          call ga_get(g_vecs, ulo, uhi, 1, ni, dvec, ld)
          if (ometric) call ga_get(g_over, 1, nu, ulo, uhi, dovl, nu)
        else if (itypv.eq.MT_DCPL) then
          call ga_get(g_vecs, ulo, uhi, 1, ni, zvec, ld)
          if (ometric) call ga_get(g_over, 1, nu, ulo, uhi, zovl, nu)
        else
          call errquit("wfn1_orthog: unsupported type",UERR)
        endif
      endif
c
c     Do the deed
c
      if (itypv.eq.MT_DBL) then
        call wfn1_orthog_dvecs(dvec, ld, ni, ulo, uhi, dsij,
     &       dovl, nu, dw, ometric)
      else if (itypv.eq.MT_DCPL) then
        call wfn1_orthog_zvecs(zvec, ld, ni, ulo, uhi, zsij,
     &       zovl, nu, zw, ometric)
      else
        call errquit("wfn1_orthog: unsupported type",UERR)
      endif
c
c     Put results back
c
      if (itypv.eq.MT_DBL) then
        if (uhi .ge. ulo) then
          call ga_put(g_vecs, ulo, uhi, 1, ni, dvec, ld)
        endif
      else if (itypv.eq.MT_DCPL) then
        if (uhi .ge. ulo) then
          call ga_put(g_vecs, ulo, uhi, 1, ni, zvec, ld)
        endif
      else
        call errquit("wfn1_orthog: unsupported type",UERR)
      endif
c
c     Tidy up memory
c
      if (itypv.eq.MT_DBL) then
        if (ulo.gt.0) then
          deallocate(dvec)
          if (ometric) then
            deallocate(dovl)
          endif
          deallocate(dsij)
          deallocate(dw)
        endif
      else if (itypv.eq.MT_DCPL) then
        if (ulo.gt.0) then
          deallocate(zvec)
          if (ometric) then
            deallocate(zovl)
          endif
          deallocate(zsij)
          deallocate(zw)
        endif
      else
        call errquit("wfn1_orthog: unsupported type",UERR)
      endif
c
      call ga_sync
c
      end subroutine wfn1_orthog
c
c-----------------------------------------------------------------------
c
      module subroutine wfn1_orthog_dvecs(vecs, ld, ni, ulo, uhi, s,
     &                                    o, nu, w, ometric)
      implicit none
#include "errquit.fh"
#include "msgids.fh"
c     
      integer, intent(in) :: ld, ni, ulo, uhi, nu
      double precision, intent(inout) :: 
     &     vecs(ulo:(ulo+ld-1),1:ni), ! [input/output] Local patch of vectors
     &     s(ni),               ! [scratch] workspace for making overlap
     &     w(nu)                ! [scratch] workspace for making overlap
      double precision, intent(in) :: 
     &     o(1:nu,ulo:(ulo+ld-1))  ! [input] metric matrix if (ometric)
      logical, intent(in) :: ometric
c     
      integer i, j, u, npass
      double precision si, scale
c     
c     Orthogonalize columns of a matrix distributed so that
c     each process has vecs(ulo:uhi,1:ni) ... uses global sums only.
c     
      do i = 1, ni
        npass = 0
   10   npass = npass + 1
c     
c       If have a metric then first form overlap*vec(i)
c     
        if (ometric) then
          call dfill(nu, 0.0d0, w, 1)
          if ((uhi-ulo+1) .gt. 0) 
     &      call dgemv('n', nu, (uhi-ulo+1), 1.0d0, o, nu,
     &           vecs(ulo,i), 1, 0.0d0, w, 1)
          call ga_dgop(msg_orthog_1, w, nu, '+') ! Actually only need ulo:hi local
        else
          do u = ulo, uhi
            w(u) = vecs(u,i)
          end do
        end if
c     
c       Now form overlap between vector i and vectors 1...i
c     
        call dfill(i, 0.0d0, s, 1)
        if ((uhi-ulo+1) .gt. 0) then
          call dgemv('t', (uhi-ulo+1), i, 1.0d0,
     &         vecs(ulo,1), ld, w(ulo), 1, 0.0d0, s(1), 1)
        endif
        call ga_dgop(msg_orthog_2, s, i, '+')
c
c       Apply the rotation
c     
        if ((uhi-ulo+1).gt.0 .and. i.gt.1) then
          call dgemv('n', (uhi-ulo+1), i-1, -1.0d0,
     &         vecs(ulo,1), ld, s(1), 1, 1.0d0, vecs(ulo,i), 1)
        end if
c     
c       Renormalize vector i
c     
        si = s(i)
        do j = 1, i-1
          si = si - s(j)*s(j)
        end do
c     
c       If the vector norm changed a lot then repeat
c     
        if (i .gt. 1) then
          scale = si/s(i)
          if (scale .lt. 0.9d0) then
            if (npass .lt. 3) then
              goto 10
            else
              call errquit('wfn1_orthog_dvecs: failed to orthog vector',
     &                     i, GA_ERR)
            end if
          end if
        end if
c
        if (si .eq. 0.0d0) then
          call errquit('wfn1_orthog_dvecs: hard zero', i, GA_ERR)
        endif
        if (si .lt. 0.0d0) then
          write(6,*) 'Negative value orthog: ',si
          call errquit('wfn1_orthog_dvecs: negative', i, 0)
        endif
        scale = 1.0d0/sqrt(si)
        do u = ulo, uhi
          vecs(u,i) = vecs(u,i) * scale
        end do
c     
      end do
c
      end subroutine wfn1_orthog_dvecs
c
c-----------------------------------------------------------------------
c
      module subroutine wfn1_orthog_zvecs(vecs, ld, ni, ulo, uhi, s,
     &                                    o, nu, w, ometric)
      implicit none
#include "errquit.fh"
#include "msgids.fh"
c     
      integer, intent(in) :: ld, ni, ulo, uhi, nu
      complex(kind=8), intent(inout) :: 
     &     vecs(ulo:(ulo+ld-1),1:ni), ! [input/output] Local patch of vectors
     &     s(ni),               ! [scratch] workspace for making overlap
     &     w(nu)                ! [scratch] workspace for making overlap
      complex(kind=8), intent(in) :: 
     &     o(1:nu,ulo:(ulo+ld-1)) ! [input] metric matrix if (ometric)
      logical, intent(in) :: ometric
c     
      integer i, j, u, npass
      complex(kind=8) :: si, scale
c     
c     Orthogonalize columns of a matrix distributed so that
c     each process has vecs(ulo:uhi,1:ni) ... uses global sums only.
c     
      do i = 1, ni
        npass = 0
   10   npass = npass + 1
c     
c       If have a metric then first form overlap*vec(i)
c     
        if (ometric) then
          w = complex(0.0d0,0.0d0)
          if ((uhi-ulo+1) .gt. 0) 
     &      call zgemv('n', nu, (uhi-ulo+1), 1.0d0, o, nu,
     &           vecs(ulo,i), 1, 0.0d0, w, 1)
          call ga_zgop(msg_orthog_1, w, nu, '+') ! Actually only need ulo:hi local
          do j = 1, nu
            w(j) = conjg(w(j))
          enddo
        else
          do u = ulo, uhi
            w(u) = conjg(vecs(u,i))
          end do
        end if
c     
c       Now form overlap between vector i and vectors 1...i
c     
        s = complex(0.0d0,0.0d0)
        if ((uhi-ulo+1) .gt. 0) then
          call zgemv('t', (uhi-ulo+1), i, 1.0d0,
     &         vecs(ulo,1), ld, w(ulo), 1, 0.0d0, s(1), 1)
        endif
        call ga_zgop(msg_orthog_2, s, i, '+')
c
c       Apply the rotation
c     
        if ((uhi-ulo+1).gt.0 .and. i.gt.1) then
          call zgemv('n', (uhi-ulo+1), i-1, -1.0d0,
     &         vecs(ulo,1), ld, s(1), 1, 1.0d0, vecs(ulo,i), 1)
        end if
c     
c       Renormalize vector i
c     
        si = s(i)
        do j = 1, i-1
          si = si - s(j)*conjg(s(j))
        end do
c     
c       If the vector norm changed a lot then repeat
c     
        if (i .gt. 1) then
          scale = si/s(i)
          if (abs(scale) .lt. 0.9d0) then
            if (npass .lt. 3) then
              goto 10
            else
              call errquit('wfn1_orthog_zvecs: failed to orthog vector',
     &                     i, GA_ERR)
            end if
          end if
        end if
c
        if (abs(si) .eq. 0.0d0) then
          call errquit('wfn1_orthog_zvecs: hard zero', i, GA_ERR)
        endif
        scale = complex(1.0d0,0.0d0)/sqrt(si)
        do u = ulo, uhi
          vecs(u,i) = vecs(u,i) * scale
        end do
c     
      end do
c
      end subroutine wfn1_orthog_zvecs
c
c-----------------------------------------------------------------------
c
C>
C>    \brief Test routine for wfn1_ortho
C>
C>    This routine simply calls wfn1_ortho but then checks the results
C>    to make sure that the orbitals returned are actually orthonormal.
C>
      module subroutine wfn1_tst_ortho(wfn1_mat,wfn1_wave,g_bfmo,
     &                                 gd_momo,gz_momo)
      implicit none
#include "stdio.fh"
#include "global.fh"
!include "wfn1_matcnst.fh"
!include "wfn1_wfn.fh"
      type(wfn1_matcnst), intent(in) :: wfn1_mat ! matrix constants
      type(wfn1_wfn),  intent(inout) :: wfn1_wave ! the wave function
c
      integer, intent(inout) :: g_bfmo  ! temporary matrix nbf x nmo
      integer, intent(inout) :: gd_momo ! temporary MT_DBL matrix nmo x nmo
      integer, intent(inout) :: gz_momo ! temporary MT_DCPL matrix nmo x nmo
c
      integer :: nbf   ! the number of basis functions
      integer :: nmo   ! the number of MOs
      integer :: itype ! the GA type
      integer :: iproc ! the current processor rank
      integer :: ii    ! counter
c
      integer :: ilo, ihi ! local matrix block limits
      integer :: jlo, jhi ! local matrix block limits
c
      double precision :: dnormT  ! norm of temporary matrix
      complex(kind=8)  :: znormT  ! norm of temporary matrix
c
      call wfn1_ortho(wfn1_mat,wfn1_wave)
c
      iproc = ga_nodeid()
      call ga_inquire(wfn1_wave%g_na,itype,nbf,nmo)
c
c     Check alpha natural orbitals Na^T S Na = I
c
      call ga_dgemm('n','n',nbf,nmo,nbf,1.0d0,wfn1_mat%g_s,
     &              wfn1_wave%g_na,0.0d0,g_bfmo)
      call ga_dgemm('t','n',nmo,nmo,nbf,1.0d0,wfn1_wave%g_na,
     &              g_bfmo,0.0d0,gd_momo)
      call ga_distribution(gd_momo,iproc,ilo,ihi,jlo,jhi)
      do ii = max(ilo,jlo),min(ihi,jhi)
        call ga_acc(gd_momo,ii,ii,ii,ii,-1.0d0,1,1.0d0)
      enddo
      call ga_norm_infinity(gd_momo,dnormT)
      if (dnormT.gt.1.0d-10) then
        if (iproc.eq.0) then
          write(LuOut,'(" Error in wfn1_ortho: Na: ",f24.15)')
     &          dnormT
        endif
c       call ga_print(wfn1_wave%g_na)
c       call ga_print(gd_momo)
      endif
c
c     Check beta natural orbitals Nb^T S Nb = I
c
      call ga_dgemm('n','n',nbf,nmo,nbf,1.0d0,wfn1_mat%g_s,
     &              wfn1_wave%g_nb,0.0d0,g_bfmo)
      call ga_dgemm('t','n',nmo,nmo,nbf,1.0d0,wfn1_wave%g_nb,
     &              g_bfmo,0.0d0,gd_momo)
      call ga_distribution(gd_momo,iproc,ilo,ihi,jlo,jhi)
      do ii = max(ilo,jlo),min(ihi,jhi)
        call ga_acc(gd_momo,ii,ii,ii,ii,-1.0d0,1,1.0d0)
      enddo
      call ga_norm_infinity(gd_momo,dnormT)
      if (dnormT.gt.1.0d-10) then
        if (iproc.eq.0) then
          write(LuOut,'(" Error in wfn1_ortho: Nb: ",f24.15)')
     &          dnormT
        endif
c       call ga_print(wfn1_wave%g_nb)
c       call ga_print(gd_momo)
      endif
c
c     Check alpha correlation functions Ca^T Ca = I
c
      call ga_zgemm('t','n',nmo,nmo,nmo,complex(1.0d0,0.0d0),
     &              wfn1_wave%g_ca,wfn1_wave%g_ca,complex(0.0d0,0.0d0),
     &              gz_momo)
      call ga_distribution(gz_momo,iproc,ilo,ihi,jlo,jhi)
      do ii = max(ilo,jlo),min(ihi,jhi)
        call ga_acc(gz_momo,ii,ii,ii,ii,complex(-1.0d0,0.0d0),
     &              1,complex(1.0d0,0.0d0))
      enddo
      call ga_norm_infinity(gz_momo,znormT)
      if (abs(znormT).gt.1.0d-10) then
        if (iproc.eq.0) then
          write(LuOut,'(" Error in wfn1_ortho: Ca: ",f24.15)')
     &          znormT
        endif
c       call ga_print(wfn1_wave%g_ca)
c       call ga_print(gz_momo)
      endif
c
c     Check beta correlation functions Cb^T Cb = I
c
      call ga_zgemm('t','n',nmo,nmo,nmo,complex(1.0d0,0.0d0),
     &              wfn1_wave%g_cb,wfn1_wave%g_cb,complex(0.0d0,0.0d0),
     &              gz_momo)
      call ga_distribution(gz_momo,iproc,ilo,ihi,jlo,jhi)
      do ii = max(ilo,jlo),min(ihi,jhi)
        call ga_acc(gz_momo,ii,ii,ii,ii,complex(-1.0d0,0.0d0),
     &              1,complex(1.0d0,0.0d0))
      enddo
      call ga_norm_infinity(gz_momo,znormT)
      if (abs(znormT).gt.1.0d-10) then
        if (iproc.eq.0) then
          write(LuOut,'(" Error in wfn1_ortho: Cb: ",f24.15)')
     &          znormT
        endif
c       call ga_print(wfn1_wave%g_cb)
c       call ga_print(gz_momo)
      endif
      end subroutine wfn1_tst_ortho
c
c-----------------------------------------------------------------------
c
      end submodule wfn1_clc_ortho
