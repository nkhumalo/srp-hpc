!     submodule (wfn1) wfn1_gss_crr
!     implicit none
!     contains
c-----------------------------------------------------------------------
c
C>
C> \brief Guess the correlation functions 
C>
C> Currently it is believed that it is important that the correlation
C> function coefficients are non-zero everywhere so that all functions
C> will mix if fractional occupation numbers are beneficial. At the
C> same time it is important that the correlation functions 
C> coefficients decay as a function of energy so that the calculation
C> starts reasonably close to the expected answer. 
C>
C> To guess the correlation function coefficients the eigenvalues of
C> the natural orbitals are used. For each orbital it is assumed that
C> the coefficients decay exponentially as a function of the energy
C> difference with the current function. In the case guess functions
C> are degenerate some additional tricks are needed to ensure 
C> orthogonality within the degenerate set. Basically we assume that the
C> weights are multiplied with a plane wave style phase factor.
C>
C> Considering the coefficients as a function of the natural orbital
C> energies we assume they are proportional to 
C> \f{eqnarray}{
C>   C^s_i \propto \exp(|\epsilon_s-\epsilon_i|/2T)
C> \f}
C> In addition a degenerate set is defined as the set of functions such
C> all functions within the set are connected to all other functions
C> through a path
C> \f$i\rightarrow k = i\rightarrow j \land j\rightarrow k\f$ such that
C> the probability \f$\exp(|\epsilon_j-\epsilon_i|/T)\f$ associated with
C> \f$i\rightarrow j\f$ exceeds \f$1/2\f$. This means there should be
C> at least one other vector in the degenerate set such that
C> \f{eqnarray}{
C>   -|\epsilon_j-\epsilon_i| \geq T \ln(1/2)
C> \f}
C> Within a degenerate set the coefficients are multiplied with a 
C> phase factor that depends on the distance to the lowest energy
C> natural orbital within the degenerate set.
C>
C> In addition it is assumed that the natural orbitals that constitute
C> a degenerate set appear in a consecutive sequence of orbitals.
C>
C> See also wfn1_nxt_init_occ
C>
      subroutine wfn1_guess_corr(wfn1_param,wfn1_wave,eig,neig)
      implicit none
#include "mafdecls.fh"
#include "global.fh"
!include "wfn1_param.fh"
!include "wfn1_wfn.fh"
C>    The calculation parameters
      type(wfn1_prm),   intent(in)    :: wfn1_param
C>    The wave function
      type(wfn1_wfn),   intent(inout) :: wfn1_wave
C>    The number of eigenvalues
      integer,          intent(in)    :: neig
C>    The eigenvalues of the natural orbitals
      double precision, intent(in)    :: eig(neig)
c
      integer nproc   ! the number of ranks
      integer iproc   ! the number of the current rank
      integer ii      ! counter
      integer jj      ! counter
      integer kk      ! counter
      integer ll      ! counter
      integer ilo,ihi ! loop limits
      integer jlo,jhi ! loop limits
c
      integer nea     ! local number of alpha electrons
      integer neb     ! local number of alpha electrons
      integer nmo     ! local number of molecular orbitals
      integer ncorea  ! local number of alpha core orbitals
      integer ncoreb  ! local number of beta core orbitals
      integer nvirta  ! local number of alpha virtual orbitals
      integer nvirtb  ! local number of beta virtual orbitals
c
      double precision :: Tinitial
c
!     complex(kind=8), allocatable :: buf(:,:)  ! buffer 
!     complex(kind=8), allocatable :: rndm(:)   ! random numbers of norm 1.0
      real(kind=8),    allocatable :: rndm(:)   ! random numbers of norm 1.0
!     complex(kind=8), allocatable :: dnrm(:)   ! vector norms
!     double precision :: pi
!     double precision :: delta
!     double precision :: x
c
      integer :: ideglo ! the bottom end of a degenerate set
      integer :: ideghi ! the top end of a degenerate set
      integer :: ndeg   ! the number of degenerate set functions
c
      nproc = ga_nnodes()
      iproc = ga_nodeid()
c
!     pi = acos(-1.0d0)
      nea = wfn1_param%nea
      neb = wfn1_param%neb
      nmo = wfn1_param%nmo
      ncorea   = wfn1_param%ncorea
      ncoreb   = wfn1_param%ncoreb
      nvirta   = wfn1_param%nvirta
      nvirtb   = wfn1_param%nvirtb
      Tinitial = wfn1_param%Tinitial
!     delta    = 0.0d0
c
c     Set up random numbers of norm 1.0
c
      allocate(rndm(1:nmo))
      if (ga_nodeid().eq.0) then
c       allocate(rrndm(1:2,1:nmo))
        call random_number(rndm)
        do ii = 1, nmo
c         rndm(ii) = complex(2.0d0*rrndm(1,ii)-1.0d0,
c    &                       2.0d0*rrndm(2,ii)-1.0d0)
c         rndm(ii) = rndm(ii)/abs(rndm(ii))
          rndm(ii) = 2.0d0*rndm(ii)-1.0d0
          rndm(ii) = rndm(ii)/abs(rndm(ii))
        enddo
c       deallocate(rrndm)
      endif
      call ga_brdcst(100,rndm,ma_sizeof(MT_DBL,nmo,MT_BYTE),0)
      call ga_sync()
c
c     Setup the alpha correlation functions
c
      call wfn1_guess_corr_spin(nmo,ncorea,nvirta,wfn1_wave%g_ca,
     &                          eig,neig,rndm,Tinitial)
c
c     Setup the beta correlation functions
c
      call wfn1_guess_corr_spin(nmo,ncoreb,nvirtb,wfn1_wave%g_cb,
     &                          eig,neig,rndm,Tinitial)
c
      call ga_sync()
      deallocate(rndm)
c
c     Orthonormalize the vectors
c
      call wfn1_orthog(wfn1_wave%g_ca,ii,.false.)
      call wfn1_orthog(wfn1_wave%g_cb,ii,.false.)
c
      end subroutine
c
c-----------------------------------------------------------------------
c
      subroutine wfn1_guess_corr_spin(nmo,ncore,nvirt,g_c,
     &                                       eig,neig,rndm,Tt)
      implicit none
#include "global.fh"
      integer,          intent(in) :: nmo
      integer,          intent(in) :: ncore
      integer,          intent(in) :: nvirt
      integer,          intent(in) :: g_c
      integer,          intent(in) :: neig
      double precision, intent(in) :: eig(neig)
!     complex(kind=8),  intent(in) :: rndm(:)
      real   (kind=8),  intent(in) :: rndm(:)
      double precision, intent(in) :: Tt      ! temperature
c
c     Local variables
c
      integer :: ilo, ihi, jlo, jhi ! delimiters of local patch of g_c
      integer :: ii, jj             ! counters
      integer :: iproc              ! this processor's rank
c
      integer :: ideglo ! the bottom end of a degenerate set
      integer :: ideghi ! the top end of a degenerate set
      integer :: ndeg   ! the number of degenerate set functions
c
      double precision :: delta
      double precision :: pi
      double precision :: x
c
c     complex(kind=8), allocatable :: buf(:,:)  ! buffer 
      real   (kind=8), allocatable :: buf(:,:)  ! buffer 
      real(kind=8),    allocatable :: dnrm(:)   ! vector norms
c
c     Setup the alpha correlation functions
c
      iproc = ga_nodeid()
      call ga_distribution(g_c,iproc,ilo,ihi,jlo,jhi)
      allocate(buf(ilo:ihi,jlo:jhi))
      allocate(dnrm(1:nmo))
c
c     Put initial values in (just weights, no phases)
c
      delta = 0.0d0
      pi = acos(-1.0d0)
c     buf = complex(0.0d0,0.0d0)
      buf = 0.0d0
      do jj = jlo, jhi
        if (jj.le.ncore) then
          if (jj.ge.ilo.and.jj.le.ihi) then
c           buf(jj,jj) = complex(1.0d0,0.0d0)
            buf(jj,jj) = 1.0d0
          endif
        else if (jj.gt.nmo-nvirt) then
          if (jj.ge.ilo.and.jj.le.ihi) then
c           buf(jj,jj) = complex(1.0d0,0.0d0)
            buf(jj,jj) = 1.0d0
          endif
        else
          if (Tt.lt.1.0d-10) then
            do ii = max(ncore+1,ilo), min(ihi,nmo-nvirt)
              if (ii.eq.jj) then
c               buf(ii,jj) = complex(1.0d0,0.0d0)
                buf(ii,jj) = 1.0d0
              else
c               buf(ii,jj) = complex(0.0d0,0.0d0)
                buf(ii,jj) = 1.0d0
              endif
            enddo
          else
            do ii = max(ncore+1,ilo), min(ihi,nmo-nvirt)
c             buf(ii,jj) = complex(exp(-abs(eig(ii)-eig(jj))/
c    &                             (2*Tt)),0.0d0)
              buf(ii,jj) = exp(-abs(eig(ii)-eig(jj))/(2*Tt))
            enddo
          endif
        endif
      enddo
c
c     Next adjust the phases to avoid all correlation functions within
c     a degenerate set becoming linearly dependent
c
      ideglo = 0
      ideghi = 0
      do while (ideglo.lt.nmo)
        ideglo = ideglo + 1
        ideghi = ideglo
        do while(ideghi.gt.ncore.and.ideghi.le.nmo-nvirt.and.
     &           (-abs(eig(ideghi+1)-eig(ideglo)).gt.Tt*log(0.5d0)
     &            .or.
     &            -abs(eig(ideghi+1)-eig(ideghi)).gt.Tt*log(0.5d0)
     &          ))
          ideghi = ideghi + 1
        enddo ! ideghi
        ndeg = ideghi-ideglo+1
        if (ndeg.eq.1) then
          delta = 0.0d0
        else
          delta = pi/ndeg
        endif
        do jj = max(ideglo,jlo), min(ideghi,jhi)
          do ii = ilo, ihi
            x = (ii-ideglo+0.5d0)*(jj-ideglo)*delta
            if (abs(cos(x)).lt.0.1d0) then
              buf(ii,jj) = 0.0d0
            else if (cos(x).lt.0.0d0) then
              buf(ii,jj) = -buf(ii,jj)
            endif
          enddo ! ii
        enddo ! jj
        ideglo = ideghi
      enddo ! ideglo
c
c     Normalize vectors
c
      dnrm = 0.0d0
      do jj = jlo, jhi
        do ii = ilo, ihi
c         dnrm(jj) = dnrm(jj) + buf(ii,jj)*conjg(buf(ii,jj))
          dnrm(jj) = dnrm(jj) + buf(ii,jj)*buf(ii,jj)
        enddo
      enddo
      call ga_dgop(1600,dnrm,nmo,"+")
      do jj = jlo, jhi
        dnrm(jj) = 1.0d0/dsqrt(dnrm(jj))
      enddo
      do jj = jlo, jhi
        do ii = ilo, ihi
          buf(ii,jj) = dnrm(jj) * buf(ii,jj) * rndm(ii)
        enddo
      enddo
      if (ihi.ge.ilo.and.jhi.ge.jlo) then
        call ga_put(g_c,ilo,ihi,jlo,jhi,buf,ihi-ilo+1)
      else
c       This process hold no part of g_c
      endif
      deallocate(dnrm)
      deallocate(buf)
c
      end subroutine
c
c-----------------------------------------------------------------------
c
!     end submodule wfn1_gss_crr
