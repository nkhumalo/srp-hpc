c-----------------------------------------------------------------------
c
C> \brief Make derivative matrices symmetric
C>
      subroutine wfn1_symm(wfn1_param,wfn1_deriv,
     +                     v_dna,v_dnb,v_dca,v_dcb)
      implicit none
c
#include "mafdecls.fh"
#include "errquit.fh"
#include "global.fh"
!include "wfn1_param.fh"
!include "wfn1_drv.fh"
c
!>    The molecular system parameters
      type(wfn1_prm), intent(in)    :: wfn1_param
!>    The energy derivative matrices
      type(wfn1_drv), intent(inout) :: wfn1_deriv
c
!>    The \$\alpha\f$ natural orbital energies
      double precision, intent(out) :: v_dna(1:wfn1_param%nmo)
!>    The \$\beta\f$ natural orbital energies
      double precision, intent(out) :: v_dnb(1:wfn1_param%nmo)
!>    The \$\alpha\f$ correlation function energies
      double precision, intent(out) :: v_dca(1:wfn1_param%nmo)
!>    The \$\beta\f$ correlation function energies
      double precision, intent(out) :: v_dcb(1:wfn1_param%nmo)
c
!>    My processor rank
      integer :: iproc
!>    Lower and upper block limits
      integer :: ilo, ihi, jlo, jhi
!>    Counters
      integer :: ii, jj
c
!>    Buffers
      double precision, allocatable :: buf(:,:), bufji(:,:)
      complex(kind=8),  allocatable :: zbuf(:,:), zbufji(:,:)
c
      iproc = ga_nodeid()
c
      v_dna = 0.0d0
      v_dnb = 0.0d0
c
      call ga_distribution(wfn1_deriv%g_dna,iproc,ilo,ihi,jlo,jhi)
      allocate(buf(ilo:ihi,jlo:jhi))
      allocate(bufji(jlo:jhi,ilo:ihi))
      call ga_get(wfn1_deriv%g_dna,ilo,ihi,jlo,jhi,buf,ihi-ilo+1)
      call ga_get(wfn1_deriv%g_dna,jlo,jhi,ilo,ihi,bufji,jhi-jlo+1)
      do jj = jlo, jhi
        do ii = ilo, ihi
          if      (ii.lt.jj) then
            buf(ii,jj) = 0.5d0*(buf(ii,jj)+bufji(jj,ii))
          else if (ii.eq.jj) then
            v_dna(ii)  = buf(ii,jj)
          else if (ii.gt.jj) then
            buf(ii,jj) = 0.5d0*(buf(ii,jj)+bufji(jj,ii))
          endif
        enddo
      enddo
      call ga_sync()
      call ga_put(wfn1_deriv%g_dna,ilo,ihi,jlo,jhi,buf,ihi-ilo+1)
      deallocate(bufji)
      deallocate(buf)
c
      call ga_distribution(wfn1_deriv%g_dnb,iproc,ilo,ihi,jlo,jhi)
      allocate(buf(ilo:ihi,jlo:jhi))
      allocate(bufji(jlo:jhi,ilo:ihi))
      call ga_get(wfn1_deriv%g_dnb,ilo,ihi,jlo,jhi,buf,ihi-ilo+1)
      call ga_get(wfn1_deriv%g_dnb,jlo,jhi,ilo,ihi,bufji,jhi-jlo+1)
      do jj = jlo, jhi
        do ii = ilo, ihi
          if      (ii.lt.jj) then
            buf(ii,jj) = 0.5d0*(buf(ii,jj)+bufji(jj,ii))
          else if (ii.eq.jj) then
            v_dnb(ii)  = buf(ii,jj)
          else if (ii.gt.jj) then
            buf(ii,jj) = 0.5d0*(buf(ii,jj)+bufji(jj,ii))
          endif
        enddo
      enddo
      call ga_sync()
      call ga_put(wfn1_deriv%g_dnb,ilo,ihi,jlo,jhi,buf,ihi-ilo+1)
      deallocate(bufji)
      deallocate(buf)
c
      v_dca = 0.0d0
      v_dcb = 0.0d0
c
      call ga_distribution(wfn1_deriv%g_dca,iproc,ilo,ihi,jlo,jhi)
      allocate(zbuf(ilo:ihi,jlo:jhi))
      allocate(zbufji(jlo:jhi,ilo:ihi))
      call ga_get(wfn1_deriv%g_dca,ilo,ihi,jlo,jhi,zbuf,ihi-ilo+1)
      call ga_get(wfn1_deriv%g_dca,jlo,jhi,ilo,ihi,zbufji,jhi-jlo+1)
      do jj = jlo, jhi
        do ii = ilo, ihi
          if      (ii.lt.jj) then
            zbuf(ii,jj) = 0.5d0*(zbuf(ii,jj)+conjg(zbufji(jj,ii)))
          else if (ii.eq.jj) then
            v_dca(ii)  = zbuf(ii,jj)
          else if (ii.gt.jj) then
            zbuf(ii,jj) = 0.5d0*(zbuf(ii,jj)+conjg(zbufji(jj,ii)))
          endif
        enddo
      enddo
      call ga_sync()
      call ga_put(wfn1_deriv%g_dca,ilo,ihi,jlo,jhi,zbuf,ihi-ilo+1)
      deallocate(zbufji)
      deallocate(zbuf)
c
      call ga_distribution(wfn1_deriv%g_dcb,iproc,ilo,ihi,jlo,jhi)
      allocate(zbuf(ilo:ihi,jlo:jhi))
      allocate(zbufji(jlo:jhi,ilo:ihi))
      call ga_get(wfn1_deriv%g_dcb,ilo,ihi,jlo,jhi,zbuf,ihi-ilo+1)
      call ga_get(wfn1_deriv%g_dcb,jlo,jhi,ilo,ihi,zbufji,jhi-jlo+1)
      do jj = jlo, jhi
        do ii = ilo, ihi
          if      (ii.lt.jj) then
            zbuf(ii,jj) = 0.5d0*(zbuf(ii,jj)+conjg(zbufji(jj,ii)))
          else if (ii.eq.jj) then
            v_dcb(ii)  = zbuf(ii,jj)
          else if (ii.gt.jj) then
            zbuf(ii,jj) = 0.5d0*(zbuf(ii,jj)+conjg(zbufji(jj,ii)))
          endif
        enddo
      enddo
      call ga_sync()
      call ga_put(wfn1_deriv%g_dcb,ilo,ihi,jlo,jhi,zbuf,ihi-ilo+1)
      deallocate(zbufji)
      deallocate(zbuf)
c
      call ga_dgop(1402,v_dna,wfn1_param%nmo,"+")
      call ga_dgop(1403,v_dnb,wfn1_param%nmo,"+")
      call ga_dgop(1404,v_dca,wfn1_param%nmo,"+")
      call ga_dgop(1405,v_dcb,wfn1_param%nmo,"+")
c
      end subroutine wfn1_symm
