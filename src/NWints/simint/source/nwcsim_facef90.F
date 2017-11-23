#ifdef USE_SIMINT
      subroutine nwcsim_init(rtdb,nbas,bases,num_der)
      use nwcsim_shell
      implicit none
#include "bas.fh"
#include "geom.fh"
#include "basP.fh"
#include "basdeclsP.fh"
#include "geomP.fh"
#include "geobasmapP.fh"
#include "apiP.fh"
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "stdio.fh"
c      use iso_c_binding
      integer rtdb,nbas,bases(3)
      integer :: num_der
c      type(c_simint_multi_shellpair), target :: msh
      double precision :: alpha(3), coef(3)
      integer :: ibasis,basis,bas
      logical cart_2_sphe
      integer :: geom, natoms
      integer :: max_der, max_ang
      integer :: mxmx
      integer :: iat,icset,ic1,ic2, nprim
      integer :: l, ncontr, isphere
      character*16 tag 
      integer l_coef,k_coef,l_exp,k_exp
      double precision coord(3)
      double precision q
      character*14 pname
      pname="nwcsim_f90init"
c     stick to nbas=1 for now
         if(nbas.gt.2) call errquit(
     C        ' simint interface not ready for no basis gt 2',0,0)
c     init
      max_der=0
#ifdef SIMINT_GRADIENT
      max_der=1
#endif
      if(num_der.gt.max_der) call errquit(
     C        ' simint interface not ready for derivative ',
     D num_der,0)
      if(.not.nwcsim_initialized) then
      call simint_init()
c     needed for 1-e integrals
      call igamma_init()

      nwcsim_initialized=.true.
c
c     create shells simint structure
c
      nwcsim_nbas=nbas
      do ibasis=1,nbas
         basis = bases(ibasis) 
         bas = basis + BASIS_HANDLE_OFFSET
         nwcsim_bas(ibasis)=bas
c     offset=-565
         cart_2_sphe=bas_spherical(bas)
         if (.not.bas_nprim_cn_max(basis, mxmx))
     C        call errquit(' bas_nprim_cn_max failed ',0,BASIS_ERR)
         if (.not. MA_Push_Get(MT_Dbl,mxmx,'simexp',l_exp,k_exp))
     C        call errquit(' simint wrk ',mxmx,MA_ERR)
         if (.not. MA_Push_Get(MT_Dbl,mxmx,'simcoef',l_coef,k_coef))
     C        call errquit(' simint wrk ',mxmx, MA_ERR)
         geom  = ibs_geom(bas)
c
         natoms   =  ncenter(geom)
         nwcsim_noshell(bas)=0
         max_ang=0
         do iat=1,natoms
            if (.not.bas_ce2cnr(basis,iat,ic1,ic2))
     &           call errquit(pname//'Exiting ',11, BASIS_ERR)
            do icset = ic1,ic2
               if (.not.bas_continfo(basis, icset,
     &              l, nprim, ncontr, isphere))
     &              call errquit(pname//'Exiting ',3, BASIS_ERR)
               max_ang=max(l,max_ang)
               if (.not.bas_get_exponent(basis, icset, dbl_mb(k_exp)))
     &              call errquit(pname//'Exiting ',7, BASIS_ERR)
               if(.not. bas_get_coeff(basis,icset,dbl_mb(k_coef)))
     &              call errquit(pname//'Exiting ',8, BASIS_ERR)
               if (.not. geom_cent_get(geom, iat, tag,
     &              coord, q))call errquit
               nwcsim_noshell(bas)=nwcsim_noshell(bas)+1
               call simint_initialize_shell(
     S              smnt_sh(nwcsim_noshell(bas),bas))
               
               call simint_create_shell(nprim, l , 
     C              coord(1), coord(2), coord(3),
     &              dbl_mb(k_exp), dbl_mb(k_coef), 
     L              smnt_sh(nwcsim_noshell(bas),bas))
c     dummy shell for 3c-
               dbl_mb(k_exp)=0d0
               dbl_mb(k_coef)=1d0
               call simint_create_shell(1, 0 , 
     C              coord(1), coord(2), coord(3),
     &              dbl_mb(k_exp), dbl_mb(k_coef), 
     L              zero_sh(nwcsim_noshell(bas),bas))
            enddo
         enddo
         if(.not.ma_chop_stack(l_exp)) call errquit(
     E        ' pop stack failed ',l_exp,0)
      enddo ! basis loop
      endif
c
c     normalization not needed since nwchem and simint use the same
cedo      call simint_normalize_shells(nwcsim_noshell,smnt_sh)
c
c     memory allocation
c

      mem_2e4c = max(mem_2e4c,
     S     simint_eri_worksize(num_der, max_ang))
      mem_2e3c = mem_2e4c
      if(num_der.eq.1) then
         memb_2e4c = mem_2e4c + mem_2e4c/5
      else
         memb_2e4c = mem_2e4c + mem_2e4c/10 ! +10% to be safe
      endif
      if(ga_nodeid().eq.0) then
         write(luout,*)
         write(luout,*) ' Using Simint Integral package'
         write(luout,*) ' simint: mem_2e4c   ',mem_2e4c
         write(luout,*) ' simint: memb_2e4c  ',memb_2e4c
         write(luout,*) ' 0 der',simint_eri_worksize(0, max_ang)
         write(luout,*) ' 1 der',simint_eri_worksize(1, max_ang)
         write(luout,*)
         call util_flush(luout)
      endif

      return
      end
      subroutine nwcsim_terminate()
      use nwcsim_shell
      implicit none
#include "errquit.fh"
      integer nbas     ! [in]
      integer bases(3) ! [in]
c
      integer icsh,ibasis,bas
      character*14 pname
      pname="nwcsim_termina"
c
      if(nwcsim_initialized) then
      do ibasis=1,nwcsim_nbas
         bas=nwcsim_bas(ibasis)
            do icsh=1,nwcsim_noshell(bas)
               call simint_free_shell(smnt_sh(icsh,bas))
               call simint_free_shell(zero_sh(icsh,bas))
            enddo
      enddo
      endif
      call simint_finalize()
      nwcsim_initialized=.false.
      return
      end
      logical function nwcsim_initdone()
      use nwcsim_shell
      implicit none
c
      nwcsim_initdone=nwcsim_initialized
      return
      end
#else 
c
c     stubs when simint is not used
c
      subroutine nwcsim_init(rtdb,nbas,bases)
      implicit none
#include "errquit.fh"
      integer rtdb,nbas,bases
      return
      end
      subroutine nwcsim_terminate(nbas,bases)
      implicit none
#include "errquit.fh"
      integer nbas,bases
      return
      end
      logical function  nwcsim_initdone()
      implicit none
c
      nwcsim_initdone=.false.
      return
      end
#endif
      