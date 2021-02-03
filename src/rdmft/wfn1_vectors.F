!> \brief Copy of ddscf/vectors.F but extended for complex SCF vectors
!>
!> These subroutines check the datatype of the g\_vecs global arrays
!> and adapt their function accordingly. The basic approach is to
!> consider every complex number as a pair of doubles. Obviously for 
!> compatibility with the rest of NWChem this code needs to be kept in
!> sync with the code in ddscf/vectors.F.
!>
      submodule (wfn1) wfn1_io_movecs
      implicit none
      contains
c
c-----------------------------------------------------------------------
c
      module function wfn1_movecs_read_header(filename, title,
     &                basis_name, scftype, nbf, nsets, nmo, ldnmo)
     &   result(wfn1_return)
      implicit none
#include "errquit.fh"
#include "global.fh"
#include "tcgmsg.fh"
#include "msgtypesf.h"
#include "msgids.fh"
#include "cscfps.fh"
#include "inp.fh"
#include "stdio.fh"
#include "mafdecls.fh"
c
c     Temporary routine
c
      character(len=*), intent(in)  :: filename    ! File to read header from
      character(len=*), intent(out) :: title       ! Returns title of job that created vectors
      character(len=*), intent(out) :: basis_name  ! Returns name of basis set
      character(len=*), intent(out) :: scftype     ! Returns the SCF type of the vectors
      integer, intent(out) :: nbf               ! Returns no. of functions in basis
      integer, intent(out) :: nsets             ! Returns no. of functions in each set
      integer, intent(in)  :: ldnmo             ! Inputs size of nmo
      integer, intent(out) :: nmo(ldnmo)        ! Returns no. of vectors in each set
c
      logical :: wfn1_return
c
      integer unitno            ! Unit no. for reading
      parameter (unitno = 67)   ! These need to be managed !!!

      integer lentit
      integer lenbas
      integer ok, i
      integer ioserr,inntsize
      character*20 scftype20
c
      wfn1_return = .false.
      if (oscfps) call pstat_on(ps_vecio)
      ok = 0
      inntsize=MA_sizeof(MT_INT,1,MT_BYTE)
      if (ga_nodeid() .eq. 0) then
         open(unitno, status='old', form='unformatted', file=filename,
     $        iostat=ioserr,err=1000)
         read(unitno, err=1001, end=2001) ! SKIP convergence info
         read(unitno, err=1001, end=2001) scftype20
         scftype = scftype20
         read(unitno, err=1001, end=2001) lentit
         if (len(title) .lt. lentit) call errquit
     $        ('movecs_read_header: title too short', lentit, INPUT_ERR)
         title = ' '
         read(unitno, err=1001, end=2001) title(1:lentit)
         read(unitno, err=1001, end=2001) lenbas
         if (len(basis_name) .lt. lenbas) call errquit
     $        ('movecs_read_header: basis_name too short', lenbas,
     &       INPUT_ERR)
         basis_name = ' '
         read(unitno, err=1001, end=2001) basis_name(1:lenbas)
         read(unitno, err=1001, end=2001) nsets
         read(unitno, err=1001, end=2001) nbf
         if (nsets .gt. ldnmo) then
            write(6,*) ' movecs_read_header: ldnmo too small ',
     $           nsets, ldnmo
            close(unitno, err=1002)
            goto 10
         endif
         read(unitno, err=1001, end=2001) (nmo(i),i=1,nsets)
         close(unitno, err=1002)
         ok = 1
      endif
c
 10   call ga_brdcst(Msg_Vec_Stat+MSGINT, ok, inntsize, 0) ! Propagate status
      wfn1_return = ok .eq. 1
      if (.not. wfn1_return) return
c
      call util_char_ga_brdcst(Msg_Vec_Ttl, title, 0)
      call util_char_ga_brdcst(Msg_Vec_Nam, basis_name, 0)
      call util_char_ga_brdcst(Msg_Vec_Type, scftype, 0)
      call ga_brdcst(Msg_Vec_NBF+MSGINT, nbf, inntsize, 0)
      call ga_brdcst(Msg_Vec_Set+MSGINT, nsets, inntsize, 0)
      call ga_brdcst(Msg_Vec_NMO+MSGINT, nmo, inntsize*nsets, 0)
c
      if (oscfps) call pstat_off(ps_vecio)
c
      return
c
 1000 write(6,*) ' wfn1_movecs_read_header: failed to open ',
     $     filename(1:inp_strlen(filename)),
     A ' IERR = ', ioserr
      call util_flush(luout)
      ok = 0
      goto 10
c
 1001 write(6,*) ' wfn1_movecs_read_header: failing reading from ',
     $     filename(1:inp_strlen(filename))
      call util_flush(luout)
      ok = 0
      close(unitno, err=1002)
      goto 10
c
 2001 write(6,*) ' wfn1_movecs_read_header:eof: failing reading from '
     &    , filename(1:inp_strlen(filename))
      call util_flush(luout)
      ok = 0
      close(unitno, err=1002)
      goto 10
c
 1002 write(6,*) ' wfn1_movecs_read_header: failed to close',
     $     filename(1:inp_strlen(filename))
      call util_flush(luout)
      ok = 0
      goto 10
c
      end function wfn1_movecs_read_header
c
c-----------------------------------------------------------------------
c
      module function wfn1_movecs_read(filename, iset, occ, evals,
     &                g_vecs)
     &  result(wfn1_return)
      implicit none
#include "errquit.fh"
#include "global.fh"
#include "tcgmsg.fh"
#include "msgtypesf.h"
#include "mafdecls.fh"
#include "msgids.fh"
#include "cscfps.fh"
#include "inp.fh"
#include "util.fh"
#include "stdio.fh"
c
      character(len=*), intent(in) :: filename
      integer, intent(in) :: iset           ! No. (1,2) of set of vectors to read
      real(kind=8), intent(out) :: occ(*)   ! Must be at least nbf long (not nmo)
      real(kind=8), intent(out) :: evals(*) ! Must be at least nbf long (not nmo)
      integer, intent(in) :: g_vecs
c
      logical :: wfn1_return
c
      integer :: nsets             ! No. of sets of vectors
      integer :: nbf               ! No. of functions in basis
      integer :: nmo(2)            ! No. of vectors in each set
      integer :: ok, jset, i, j
      integer :: l_vecs, k_vecs
      integer :: itype
      integer :: unitno
      parameter (unitno = 67)
      integer :: inntsize,ddblsize,dcplsize
c     
      l_vecs = -1               ! An invalid MA handle
      wfn1_return = .false.
c
      if (oscfps) call pstat_on(ps_vecio)
      inntsize=MA_sizeof(MT_INT,1,MT_BYTE)
      ddblsize=MA_sizeof(MT_DBL,1,MT_BYTE)
      dcplsize=MA_sizeof(MT_DCPL,1,MT_BYTE)
      call ga_sync()
      call ga_enquire(g_vecs,itype,i,j)
      ok = 0
      if (ga_nodeid() .eq. 0) then
         open(unitno, status='old', form='unformatted', file=filename,
     $        err=1000)
c
c        Skip over uninteresting bits of the header
c
         read(unitno, err=1001, end=1001) ! convergence info
         read(unitno, err=1001, end=1001) ! scftype
         read(unitno, err=1001, end=1001) ! lentit
         read(unitno, err=1001, end=1001) ! title
         read(unitno, err=1001, end=1001) ! lenbas
         read(unitno, err=1001, end=1001) ! basis_name
         read(unitno, err=1001, end=1001) nsets
         read(unitno, err=1001, end=1001) nbf
         read(unitno, err=1001, end=1001) (nmo(i),i=1,nsets)
c
         if (.not. ma_push_get(itype,nbf,'wfn1_movecs_read',
     &                         l_vecs,k_vecs))
     &        call errquit('wfn1_movecs_read: ma failed', nbf, MA_ERR)
c
c        Skip over unwanted sets
c
         do jset = 1, iset-1
            read(unitno, err=1001, end=1001)
            read(unitno, err=1001, end=1001)
            do i = 1, nmo(jset)
               read(unitno, err=1001, end=1001)
            enddo
         enddo
         read(unitno, err=1001, end=1001) (occ(j),j=1,nbf)
         read(unitno, err=1001, end=1001) (evals(j),j=1,nbf)
         if (itype.eq.MT_DBL) then
           do i = 1, nmo(iset)
*            read(unitno, err=1001, end=1001)
*     $           (dbl_mb(k_vecs+j), j=0,nbf-1)
             call sread(unitno, dbl_mb(k_vecs), nbf)
             call ga_put(g_vecs, 1, nbf, i, i, dbl_mb(k_vecs), nbf)
           enddo
         else if (itype.eq.MT_DCPL) then
           do i = 1, nmo(iset)
*            read(unitno, err=1001, end=1001)
*     $           (dbl_mb(k_vecs+j), j=0,nbf-1)
             call sread(unitno, dcpl_mb(k_vecs), 2*nbf)
             call ga_put(g_vecs, 1, nbf, i, i, dcpl_mb(k_vecs), nbf)
           enddo
         else
           call errquit('wfn1_movecs read: invalid type',itype,GA_ERR)
         endif
 9       close(unitno,err=1002)
         ok = 1
      endif
c
 10   continue
      if (l_vecs .ne. -1) then
        if (.not. ma_pop_stack(l_vecs)) call errquit
     $      ('wfn1_movecs_read: pop failed', l_vecs, MA_ERR)
      endif
      call ga_sync()
      call ga_brdcst(Msg_Vec_Stat+MSGINT, ok, inntsize, 0) ! Propagate status
      if (ok .eq. 1) then
        call ga_brdcst(Msg_Vec_NBF+MSGINT, nbf, inntsize, 0)
        call ga_brdcst(Msg_Vec_EVal+MSGDBL, evals, ddblsize*nbf, 0)
        call ga_brdcst(Msg_Vec_Occ+MSGDBL, occ, ddblsize*nbf, 0)
      endif
c
      wfn1_return = ok .eq. 1
      if (ga_nodeid() .eq. 0 .and. wfn1_return .and.
     $     util_print('vectors i/o', print_high)) then
         write(6,22) filename(1:inp_strlen(filename))
 22      format(/' Read molecular orbitals from ',a/)
         call util_flush(luout)
      endif
      if (oscfps) call pstat_off(ps_vecio)
      return
c
 1000 write(6,*) ' wfn1_movecs_read: failed to open ',
     $     filename(1:inp_strlen(filename))
      call util_flush(luout)
      ok = 0
      goto 10
c
 1001 write(6,*) ' wfn1_movecs_read: failing reading from ',
     $     filename(1:inp_strlen(filename))
      call util_flush(luout)
      ok = 0
      close(unitno,err=1002)
      goto 10
c
 1002 write(6,*) ' wfn1_movecs_read: failed to close',
     $     filename(1:inp_strlen(filename))
      call util_flush(luout)
      ok = 0
      goto 10
c
      end function wfn1_movecs_read
c
c-----------------------------------------------------------------------
c
      module function wfn1_movecs_gather(filename, iset, nwant, want,
     &   occ, evals, g_vecs) result(wfn1_return)
      implicit none
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
c
      character(len=*), intent(in) :: filename ! name of movecs file
      integer, intent(in) :: iset               ! no. (1,2) of set of vectors to read
      integer, intent(in) :: nwant              ! length of gather vector
      integer, intent(in) :: want(nwant)        ! list of mos to be gathered
      real(kind=8), intent(out) :: occ(nwant)   ! occupation numbers
      real(kind=8), intent(out) :: evals(nwant) ! eigenvalues
      integer, intent(in) :: g_vecs             ! handle for ga to hold results
c
      logical :: wfn1_return
c
      integer nsets                 ! no. of sets of vectors
      integer nbf                   ! no. of functions in basis
      integer ldnmo                 ! size of nmo
      parameter (ldnmo = 2)
      integer nmo(ldnmo)            ! no. of vectors in each set
      integer i, l_fullocc, k_fullocc, l_fulleval, k_fulleval, g_tmp
      character*255 title, basis
      integer minchunk              ! for data distirbution in ga
      parameter (minchunk = -1)     ! let ga determine even distrib.
      character*20 scftype
      integer :: itype, ii, jj
c
!     logical movecs_read_header, movecs_read
!     external movecs_read_header, movecs_read
c
      wfn1_return = .false.
c     
      if ( .not. wfn1_movecs_read_header(filename, title, basis,
     $     scftype, nbf, nsets, nmo, ldnmo) ) call errquit(
     $   'wfn1_movecs_gather: failed to read movecs file', 1, INPUT_ERR)
c
c     allocate memory to be used while reading info
C     movecs_read expects fullocc and fulleval arrays to be nbf,
C     not nmo(iset)!
c
      if (.not. ma_push_get(mt_dbl,nbf,'movecs_gather',
     $   l_fullocc, k_fullocc)) call errquit(
     $   'wfn1_movecs_gather: ma failed', nmo(iset), MA_ERR)
c
      if (.not. ma_push_get(mt_dbl,nbf,'movecs_gather',
     $   l_fulleval, k_fulleval)) call errquit(
     $   'wfn1_movecs_gather: ma failed', nmo(iset), MA_ERR)
c
      call ga_enquire(g_vecs,itype,ii,jj)
      if ( .not. ga_create( itype, nbf, nmo(iset), 'full vectors',
     $   minchunk, minchunk, g_tmp) ) call errquit(
     $   'wfn1_movecs_gather: ga_create failed', nbf*nmo(iset) , GA_ERR)
c
c     read the full data
c
      if ( .not. wfn1_movecs_read(filename, iset, dbl_mb(k_fullocc),
     $   dbl_mb(k_fulleval), g_tmp) ) call errquit(
     $   'wfn1_movecs_gather: failed reading movecs file', 2, DISK_ERR)
c
c     gather the data into the requested subset
c
      do i = 1, nwant
         occ(i)   = dbl_mb( k_fullocc  + want(i) - 1)
         evals(i) = dbl_mb( k_fulleval + want(i) - 1)
         call ga_copy_patch('n', g_tmp, 1, nbf, want(i), want(i),
     $      g_vecs, 1, nbf, i, i)
      enddo
c
c     free temporaries
c
      if ( .not. ga_destroy(g_tmp)) call errquit(
     $   'wfn1_movecs_gather: ga_destroy failed', g_tmp, GA_ERR)
      if (.not. ma_pop_stack(l_fulleval)) call errquit(
     $   'wfn1_movecs_gather: pop failed', l_fulleval, MA_ERR)
      if (.not. ma_pop_stack(l_fullocc)) call errquit(
     $   'wfn1_movecs_gather: pop failed', l_fullocc,
     &       MA_ERR)
c
      wfn1_return = .true.
c
      end function wfn1_movecs_gather
c
c-----------------------------------------------------------------------
c
      module function wfn1_movecs_write(rtdb, basis, filename, scftype,
     &       title, nbf, nsets, nmo, occ, ldocc, evals, ldevals, g_vecs)
     &  result(wfn1_return)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "tcgmsg.fh"
#include "msgtypesf.h"
#include "inp.fh"
#include "msgids.fh"
#include "cscfps.fh"
#include "util.fh"
#include "bas.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "stdio.fh"
#include "bq.fh"
c
c     Temporary routine
c
      integer, intent(in) :: rtdb              ! RTDB handle (-1 if not accessible)
      integer, intent(in) :: basis             ! Basis handle(-1 if not accessible)
      character(len=*), intent(in) :: filename ! File to write to
      character(len=*), intent(in) :: scftype  ! Type of SCF (dft, scf, mcscf)
      character(len=*), intent(in) :: title    ! Title of job that created vectors
      integer, intent(in) :: nbf               ! No. of functions in basis
      integer, intent(in) :: nsets             ! No. of sets of vectors
      integer, intent(in) :: nmo(nsets)        ! No. of vectors in each set
      integer, intent(in) :: ldocc
      integer, intent(in) :: ldevals
      real(kind=8), intent(in) :: occ(ldocc, nsets)     ! Occupation numbers 
      real(kind=8), intent(in) :: evals(ldevals, nsets) ! Eigenvalues
      integer, intent(in) :: g_vecs(nsets)              ! Global array with eigen-vectors
      logical wfn1_return
c
      integer unitno            ! Unit no. for writing
      parameter (unitno = 67)   ! These need to be managed !!!
      integer bq_handle         ! Handle of the BQ fragment
      integer lentit
      integer lenbas
      integer l_vecs, k_vecs
      integer ok, iset, i, j
      integer geom, ma_type, nelem
      integer nmem     ! the number of real*8 words available
      integer nvec     ! the number of vectors that fit in local memory
      integer nmo_top  ! the top MO currently considered
      character*26 date
      character*32 geomsum, basissum, bqsum, key
      character*20 scftype20    ! DECL MUST match movecs_converged routine
      character*128 basis_name, trans_name
      double precision energy, enrep
      integer inntsize,kcount
      integer :: itype, ii, jj
c
      wfn1_return = .false.
      if (oscfps) call pstat_on(ps_vecio)
      inntsize=MA_sizeof(MT_INT,1,MT_BYTE)
      call ga_sync()
c
c     Generate info to go into convergence and restart line
c
      if (basis.ne.-1 .and. rtdb.ne.-1) then
         if (.not. bas_checksum(basis, basissum))
     $      call errquit('wfn1_movecs_write: bad basis',0, BASIS_ERR)
         if (.not. bas_geom(basis, geom))
     $      call errquit('wfn1_movecs_write: bad basis',0, BASIS_ERR)
         if (.not. geom_checksum(geom, geomsum))
     $      call errquit('wfn1_movecs_write: bad geometry',0, GEOM_ERR)
         if (.not. bas_name(basis, basis_name, trans_name))
     $      call errquit('wfn1_movecs_write: bad basis',0, BASIS_ERR)
         if (trans_name .ne. ' ') basis_name = trans_name
         if (.not. bq_get_active(bq_handle)) then
              bqsum = ' '
         else
           if (.not. bq_checksum(bq_handle,bqsum)) then
              bqsum = ' '
           endif
         endif
c         
         key = ' '
         write(key,'(a,'':converged'')') scftype(1:inp_strlen(scftype))
         if (.not. rtdb_get_info(rtdb, key, ma_type, nelem, date))
     $        date = ' '
      else
         basissum = ' '
         geomsum  = ' '
         bqsum    = ' '
         date     = ' '
      endif
      scftype20 = scftype
c
c
      ok = 0
c
c Write scf energy at the end of the movecs file. If energy
c is not in rtdb, it is given a value of zero.
c
      call ga_enquire(g_vecs,itype,ii,jj)
      write(key,'(a,'':energy'')') scftype(1:inp_strlen(scftype))
      if (.not. rtdb_get(rtdb, key, mt_dbl, 1, energy)) then
        energy=0.0d0
      endif
c
c get nuclear repulsion energy
      if (.not. geom_nuc_rep_energy(geom, enrep))
     $  call errquit('unable to get nuclear rep energy',geom, GEOM_ERR)
c
      if (ga_nodeid() .eq. 0) then
         if (nsets .gt. 2) then
            write(6,*) ' wfn1_movecs_write: nsets > 2 ', nsets
            goto 10
         endif
         open(unitno, status='unknown', form='unformatted',
     $        file=filename, err=1000)
c
c        Information about convergence
c
         write(unitno, err=1001) basissum, geomsum, bqsum, scftype20,
     $                           date
c
c        Check that read routines are both consistent with this
c
         write(unitno, err=1001) scftype20
         lentit = max(1,inp_strlen(title)) ! 0 length record confuses f2c
         write(unitno, err=1001) lentit
         write(unitno, err=1001) title(1:lentit)
         lenbas = max(1,inp_strlen(basis_name))
         write(unitno, err=1001) lenbas
         write(unitno, err=1001) basis_name(1:lenbas)
         write(unitno, err=1001) nsets
         write(unitno, err=1001) nbf
         write(unitno, err=1001) (nmo(i),i=1,nsets)
         nmem = ma_inquire_avail(itype)
         nvec = max(1,min(nbf,nmem/nbf))
         if (.not. ma_push_get(itype,nbf*nvec,'movecs_write',l_vecs,
     K        k_vecs))
     $        call errquit('wfn1_movecs_write: ma failed', nbf, MA_ERR)
         do iset = 1, nsets
           write(unitno, err=1001) (occ(j,iset),j=1,nbf)
           write(unitno, err=1001) (evals(j,iset),j=1,nbf)
           if (itype.eq.MT_DBL) then
             do i = 1, nmo(iset), nvec
               nmo_top = min(nmo(iset),i+nvec-1)
               call ga_get(g_vecs(iset), 1, nbf, i, nmo_top,
     D              dbl_mb(k_vecs),nbf)
               kcount = k_vecs
               do j = i, nmo_top
                 call swrite(unitno, dbl_mb(kcount), nbf)
                 kcount = kcount + nbf
               enddo
             enddo
           else if (itype.eq.MT_DCPL) then
             do i = 1, nmo(iset), nvec
               nmo_top = min(nmo(iset),i+nvec-1)
               call ga_get(g_vecs(iset), 1, nbf, i, nmo_top,
     D              dcpl_mb(k_vecs),nbf)
               kcount = k_vecs
               do j = i, nmo_top
                 call swrite(unitno, dcpl_mb(kcount), 2*nbf)
                 kcount = kcount + nbf
               enddo
             enddo
           else
             call errquit('wfn1_movecs_write: invalid type',itype,UERR)
           endif

         enddo
         if (.not. ma_pop_stack(l_vecs))
     $        call errquit('wfn1_movecs_write: ma pop failed', l_vecs,
     &       MA_ERR)
c
c        Write scf energy at the end of the movecs file. If energy
c        is not in rtdb, it is given a value of zero.
c
         write(unitno, err=1001) energy, enrep
c         
         close(unitno,err=1002)
         ok = 1
      endif
c
 10   call ga_brdcst(Msg_Vec_Stat+MSGINT, ok, inntsize, 0) ! Propagate status
c
      wfn1_return = ok .eq. 1
      if (ga_nodeid() .eq. 0 .and.
     $     util_print('vectors i/o', print_high)) then
         write(6,22) filename(1:inp_strlen(filename))
 22      format(/' Wrote molecular orbitals to ',a/)
         call util_flush(luout)
      endif
      call ga_sync()
      if (oscfps) call pstat_off(ps_vecio)
      return
c
 1000 write(6,*) ' wfn1_movecs_write: failed to open ',
     $     filename(1:inp_strlen(filename))
      call util_flush(luout)
      ok = 0
      goto 10
c
 1001 write(6,*) ' wfn1_movecs_write: failing writing to ',
     $     filename(1:inp_strlen(filename))
      call util_flush(luout)
      ok = 0
      close(unitno,err=1002)
      goto 10
c
 1002 write(6,*) ' wfn1_movecs_write: failed to close',
     $     filename(1:inp_strlen(filename))
      call util_flush(luout)
      ok = 0
      goto 10
c
      end function wfn1_movecs_write
c
c-----------------------------------------------------------------------
c
      module subroutine wfn1_movecs_swap(rtdb, module, scftype,g_vecs,
     &       occ,eval)
      implicit none
#include "errquit.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "util.fh"
#include "inp.fh"
#include "stdio.fh"
c
      character(len=*), intent(in) :: scftype
      integer, intent(in) :: rtdb, g_vecs(*)
      real(kind=8), intent(inout) :: occ(*), eval(*)
      character(len=*), intent(in) :: module
c
c     In the database may be an entry '<module>:swap' (with
c     alpha/beta appended for UHF) which is a list of pairs
c     of vectors to swap.  Apply these to the input vectors, 
c     occupation numbers and eigenvalues.  In the case of
c     UHF g_vecs is assumed to be an array with 2 handles.
c     When finished delete the database entries so that
c     a restart does not re-swap the vectors.
c
      double precision tmp
      integer itype, nelem, handle, index
      integer l_veci, k_veci, l_vecj, k_vecj, dim1, dim2
      integer pair, i, j, ioff, joff, off, set, nset, mlen
      logical oprint,dorot
      double precision mixang,cosa,sina,pi
      integer k_tj,l_tj,k_ti,l_ti
      complex(kind=8), parameter :: zzero = complex(0.0d0,0.0d0)
      complex(kind=8) :: zcosa, zsina
      character*40 name
      character*8 text
c
      oprint = util_print('mo guess', print_default)
      mlen = inp_strlen(module)
c
      if (scftype .eq. 'UHF') then
         nset = 2
      else
         nset = 1
      endif
c     check if we want the mix rotation
      if (rtdb_get(rtdb, 'mixang', mt_dbl, 1, mixang)) then
         dorot=.true.
         pi=acos(-1d0)
         cosa=cos(mixang*pi/180d0)
         sina=sin(mixang*pi/180d0)
         zcosa=complex(cosa,0.0d0)
         zsina=complex(sina,0.0d0)
c        write(6,*) ' nseeeet ',nset
      else
         dorot=.false.
      endif
c
      do set = 1, nset
         name = ' '
         if (scftype .eq. 'UHF') then
            if (set.eq.1) then
               write(name,'(a,a)') module(1:mlen),':swap alpha'
               text = 'alpha'
            else
               write(name,'(a,a)') module(1:mlen),':swap beta'
               text = 'beta'
            endif
         else
            write(name,'(a,a)') module(1:mlen),':swap'
            text = scftype
         endif
         if (.not. rtdb_ma_get(rtdb, name, itype, nelem, 
     $        handle)) goto 1000
         if (.not. ma_get_index(handle,index)) call errquit
     $        ('movecs_swap: get index failed', handle, MA_ERR)
         if (itype.ne.mt_int) call errquit
     $        ('movecs_swap: invalid list of pairs',nelem, MA_ERR)
c
         call ga_sync()
         call ga_inquire(g_vecs, itype, dim1, dim2)
         if (ga_nodeid() .eq. 0) then
            if (oprint) write(6,*)
            if (.not. ma_push_get(itype, dim1, 'swap', l_veci, k_veci))
     $           call errquit('movecs_swap: no scratch space', dim1,
     &       MA_ERR)
            if (.not. ma_push_get(itype, dim1, 'swap', l_vecj, k_vecj))
     $           call errquit('movecs_swap: no scratch space', dim1,
     &       MA_ERR)
            do pair = 1, nelem, 2
              i = int_mb(index+pair-1)
              j = int_mb(index+pair  )
              if (oprint) write(6,1) text(1:inp_strlen(text)), i, j
 1            format(' Swapping ',a,' orbitals ', 2i5)
              if (i.lt.0 .or. i.gt.dim2) call errquit
     $             ('movecs_swap: invalid vector ', i, INPUT_ERR)
              if (j.lt.0 .or. j.gt.dim2) call errquit
     $             ('movecs_swap: invalid vector ', j, INPUT_ERR)
c     
              if (itype.eq.MT_DBL) then
                call ga_get(g_vecs(set), 1, dim1, i, i, 
     D               dbl_mb(k_veci),dim1)
                call ga_get(g_vecs(set), 1, dim1, j, j, 
     D               dbl_mb(k_vecj),dim1)
c
                if(dorot) then 
                  if (.not. ma_push_get(mt_dbl, dim1, 'swap', 
     J                 l_ti, k_ti))
     $                 call errquit('moswap: pushget fail', dim1, 0)
                  if (.not. ma_push_get(mt_dbl, dim1, 'swap', 
     J                 l_tj, k_tj))
     $                 call errquit('moswap: pushget fail', dim1, 0)
                  if(set.eq.1) then
c                   alpha
c                   psi_1 =  cos*psi_1 + sin*psi_2
                    call dcopy(dim1, 0.d0, 0, dbl_mb(k_ti), 1)
                    call daxpy(dim1,cosa,dbl_mb(k_veci),1,
     &                   dbl_mb(k_ti),1)
                    call daxpy(dim1,sina,dbl_mb(k_vecj),1,
     &                   dbl_mb(k_ti),1)
c                   psi_2 = -sin*psi_1 + cos*psi_2
                    call dcopy(dim1, 0.d0, 0, dbl_mb(k_tj), 1)
                    call daxpy(dim1,-sina,dbl_mb(k_veci),1,
     &                   dbl_mb(k_tj),1)
                    call daxpy(dim1,cosa,dbl_mb(k_vecj),1,
     &                   dbl_mb(k_tj),1)
c
                    call ga_put(g_vecs(set),1,dim1,i,i,
     &                   dbl_mb(k_ti),1)
                    call ga_put(g_vecs(set),1,dim1,j,j,
     &                   dbl_mb(k_tj),1)
                  else
c
c                   for beta we swap 1 and 2
c
                    call ga_get(g_vecs(2), 1, dim1, i, i, 
     &                   dbl_mb(k_veci),1)
                    call ga_get(g_vecs(2), 1, dim1, j, j, 
     &                   dbl_mb(k_vecj),1)
c                   psi_1 =  -sin*psi_1 + cos*psi_2
                    call dcopy(dim1, 0.d0, 0, dbl_mb(k_ti), 1)
                    call daxpy(dim1,-sina,dbl_mb(k_veci),1,
     &                   dbl_mb(k_ti),1)
                    call daxpy(dim1,cosa,dbl_mb(k_vecj),1,
     &                   dbl_mb(k_ti),1)
c                   psi_2 = cos*psi_1 + sin*psi_2
                    call dcopy(dim1, 0.d0, 0, dbl_mb(k_tj), 1)
                    call daxpy(dim1,cosa,dbl_mb(k_veci),1,
     &                   dbl_mb(k_tj),1)
                    call daxpy(dim1,sina,dbl_mb(k_vecj),1,
     &                   dbl_mb(k_tj),1)
c
                    call ga_put(g_vecs(2), 1, dim1, i, i, 
     &                   dbl_mb(k_ti),1)
                    call ga_put(g_vecs(2), 1, dim1, j, j, 
     &                   dbl_mb(k_tj),1)
                  endif
                  if (.not. ma_chop_stack(l_ti)) call errquit(
     M                 'moswap: chop stack?', 0,MA_ERR)
                else
                  call ga_put(g_vecs(set), 1, dim1, i, i, 
     D                 dbl_mb(k_vecj),dim1)
                  call ga_put(g_vecs(set), 1, dim1, j, j, 
     D                 dbl_mb(k_veci),dim1)
                endif
              else if (itype.eq.MT_DCPL) then
                call ga_get(g_vecs(set), 1, dim1, i, i, 
     D               dcpl_mb(k_veci),dim1)
                call ga_get(g_vecs(set), 1, dim1, j, j, 
     D               dcpl_mb(k_vecj),dim1)
c
                if(dorot) then 
                  if (.not. ma_push_get(mt_dcpl, dim1, 'swap', 
     J                 l_ti, k_ti))
     $                 call errquit('moswap: pushget fail', dim1, 0)
                  if (.not. ma_push_get(mt_dcpl, dim1, 'swap', 
     J                 l_tj, k_tj))
     $                 call errquit('moswap: pushget fail', dim1, 0)
                  if(set.eq.1) then
c                   alpha
c                   psi_1 =  cos*psi_1 + sin*psi_2
                    call zcopy(dim1, zzero, 0, dcpl_mb(k_ti), 1)
                    call zaxpy(dim1,zcosa,dcpl_mb(k_veci),1,
     &                   dcpl_mb(k_ti),1)
                    call zaxpy(dim1,zsina,dcpl_mb(k_vecj),1,
     &                   dcpl_mb(k_ti),1)
c                   psi_2 = -sin*psi_1 + cos*psi_2
                    call zcopy(dim1, zzero, 0, dcpl_mb(k_tj), 1)
                    call zaxpy(dim1,-zsina,dcpl_mb(k_veci),1,
     &                   dcpl_mb(k_tj),1)
                    call zaxpy(dim1,zcosa,dcpl_mb(k_vecj),1,
     &                   dcpl_mb(k_tj),1)
c
                    call ga_put(g_vecs(set),1,dim1,i,i,
     &                   dcpl_mb(k_ti),1)
                    call ga_put(g_vecs(set),1,dim1,j,j,
     &                   dcpl_mb(k_tj),1)
                  else
c
c                   for beta we swap 1 and 2
c
                    call ga_get(g_vecs(2), 1, dim1, i, i, 
     &                   dcpl_mb(k_veci),1)
                    call ga_get(g_vecs(2), 1, dim1, j, j, 
     &                   dcpl_mb(k_vecj),1)
c                   psi_1 =  -sin*psi_1 + cos*psi_2
                    call zcopy(dim1, zzero, 0, dcpl_mb(k_ti), 1)
                    call zaxpy(dim1,-zsina,dcpl_mb(k_veci),1,
     &                   dcpl_mb(k_ti),1)
                    call zaxpy(dim1,zcosa,dcpl_mb(k_vecj),1,
     &                   dcpl_mb(k_ti),1)
c                   psi_2 = cos*psi_1 + sin*psi_2
                    call zcopy(dim1, zzero, 0, dcpl_mb(k_tj), 1)
                    call zaxpy(dim1,zcosa,dcpl_mb(k_veci),1,
     &                   dcpl_mb(k_tj),1)
                    call zaxpy(dim1,zsina,dcpl_mb(k_vecj),1,
     &                   dcpl_mb(k_tj),1)
c
                    call ga_put(g_vecs(2), 1, dim1, i, i, 
     &                   dcpl_mb(k_ti),1)
                    call ga_put(g_vecs(2), 1, dim1, j, j, 
     &                   dcpl_mb(k_tj),1)
                  endif
                  if (.not. ma_chop_stack(l_ti)) call errquit(
     M                 'moswap: chop stack?', 0,MA_ERR)
                else
                  call ga_put(g_vecs(set), 1, dim1, i, i, 
     D                 dcpl_mb(k_vecj),dim1)
                  call ga_put(g_vecs(set), 1, dim1, j, j, 
     D                 dcpl_mb(k_veci),dim1)
                endif
              else
                call errquit('wfn1_swap: invalid type',itype,UERR)
              endif ! itype
c
            enddo
            if (oprint) call util_flush(luout)
            if (.not. ma_pop_stack(l_vecj)) call errquit('ms:ma?', 0,
     &       MA_ERR)
            if (.not. ma_pop_stack(l_veci)) call errquit('ms:ma?', 0,
     &       MA_ERR)
         endif
c
         do pair = 1, nelem, 2
           i = int_mb(index+pair-1)
           j = int_mb(index+pair  )
           ioff = i + (set-1)*dim1
           joff = j + (set-1)*dim1
           tmp = occ(ioff)
           occ(ioff) = occ(joff)
           occ(joff) = tmp
           tmp = eval(ioff)
           eval(ioff) = eval(joff)
           eval(joff) = tmp
         enddo
c
         call ga_sync
c
         if (.not. rtdb_delete(rtdb, name)) call errquit
     $        ('movecs_swap: delete of swap entry failed', 0, RTDB_ERR)
c
         if (.not. ma_free_heap(handle)) call errquit('ms:ma?', 0,
     &       MA_ERR)
 1000    continue
      enddo
c
      end subroutine wfn1_movecs_swap
c
c-----------------------------------------------------------------------
c
      module function wfn1_movecs_converged(rtdb, basis, scftype,
     &        filename) result(wfn1_converged)
      implicit none
#include "errquit.fh"
#include "global.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "inp.fh"
#include "bas.fh"
#include "geom.fh"
#include "tcgmsg.fh"
#include "bq.fh"
      integer, intent(in) :: rtdb              ! [input]
      integer, intent(in) :: basis             ! [input]
      character(len=*), intent(in) :: scftype  ! Current SCF type
      character(len=*), intent(in) :: filename ! MOvecs filename
      logical :: wfn1_converged
c
c     Return true/false if the MO vectors in the file correspond
c     to converged vectors of the current SCF type with the
c     current set of options.
c
c     1) The movecs file contains basis checksum, geometry checksum,
c     .  scftype, date on convergence flag in the database at the
c     .  time the vectors were written.
c
c     2) The database contains a convergence flag that is deleted
c     .  whenever any input is changed or options automatically 
c     .  changed.
c
c     The calculation is converged if all of the following are true
c     
c     a) convergence flag is present in the database and the date 
c     .  matches the date in the movecs file
c     b) scftypes match
c     c) basis and geometry checksums match
c
      character*20 filescftype  ! MUST MATCH WRITE ROUTINE DECL
      character*32 filebasissum, filegeomsum, basissum, geomsum
      character*32 filebqsum, bqsum
      logical oconverged
      character*40 key
      character*26 date, filedate
      integer geom, ma_type, nelem
      integer bq_handle         ! Handle of the BQ fragment
      integer unitno            ! Unit no. for reading
      integer inntsize
      parameter (unitno = 67)   ! These need to be managed !!!
      logical status, oldmode
c
      status = .false.
      wfn1_converged = status
      inntsize=MA_sizeof(MT_INT,1,MT_BYTE)
c
      oldmode =  rtdb_parallel(.false.)
      if (ga_nodeid() .eq. 0) then
         open(unitno, status='old', form='unformatted',
     $        file=filename, err=1000)
         read(unitno, end=1001, err=2001)
     $        filebasissum, filegeomsum, filebqsum, filescftype,
     $        filedate
         goto 3001
c
c        If the previous read statement failed we may have a movecs
c        file from an old version of the code. So try reading the 
c        old format (but this means that the vectors cannot be
c        considered as converged). Non-convergence of the vectors
c        is detected from the filebqsum value which we set to something
c        unusual so it will not match anything obvious.
c
 2001    continue
         filebqsum = " "
         rewind(unitno, err=1001)
         read(unitno, end=1001, err=1001)
     $     filebasissum, filegeomsum, filescftype, filedate
 3001    continue
         close(unitno, err=1002)
c     
c     a)
c     
         key = ' '
         write(key,'(a,'':converged'')') scftype(1:inp_strlen(scftype))
         if (.not. rtdb_get(rtdb, key, mt_log, 1, oconverged)) goto 10
         if (.not. oconverged) goto 10
         if (.not. rtdb_get_info(rtdb, key, ma_type, nelem, date))
     $        call errquit('movecs_converged: rtdb corrupt?',0,
     &       RTDB_ERR)
         if (date .ne. filedate) goto 10
c     
c     b) 
c     
         if (.not. inp_compare(.false., scftype, filescftype)) goto 10
c     
c     c) 
c     
         if (.not. bas_geom(basis, geom)) 
     $        call errquit('movecs_converged: bad basis',0, BASIS_ERR)
         if (.not. bas_checksum(basis, basissum))
     $        call errquit('movecs_converged: bad basis',0, BASIS_ERR)
         if (.not. geom_checksum(geom, geomsum))
     $        call errquit('movecs_converged: bad geometry',0,
     &       BASIS_ERR)
         if (.not. bq_get_active(bq_handle)) then
              bqsum = ' '
         else
           if (.not. bq_checksum(bq_handle,bqsum)) then
              bqsum = ' '
           endif
         endif
         if (basissum.ne.filebasissum .or. geomsum.ne.filegeomsum .or.
     $       bqsum.ne.filebqsum) 
     $        goto 10
c     
c     Cool ... it is converged!
c
         status = .true.
      endif
c
 10   oldmode =  rtdb_parallel(oldmode)
      call ga_brdcst(1311, status, inntsize, 0)
      wfn1_converged = status
      return
c
 1000 call errquit('wfn1_movecs_converged: error opening file',0,
     &             DISK_ERR)
 1001 call errquit('wfn1_movecs_converged: error reading file',0,
     &             DISK_ERR)
 1002 call errquit('wfn1_movecs_converged: error closing file',0,
     &             DISK_ERR)
c
      end function wfn1_movecs_converged
c
c-----------------------------------------------------------------------
c
      end submodule wfn1_io_movecs
