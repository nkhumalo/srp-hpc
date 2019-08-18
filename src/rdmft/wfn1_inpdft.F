c-----------------------------------------------------------------------
c
!> \brief load the DFT parameters
!>
!> This subroutine loads the DFT parameters from the runtime database
!> and stores them in the appropriate DFT common blocks. This step
!> enables using the xc_getv routine to compute the Kohn-Sham matrix
!> contributions.
!>
!> Note that this code has its own solver strategy and therefore not
!> all parameters of the DFT code are relevant here. What we do not 
!> use is simply not retrieved.
!>
!> In particular we currently always use 100% exact exchange, and so
!> exchange functionals are suppressed.
!>
      subroutine wfn1_inpdft(wfn1_param,wfn1_inst)
      implicit none
#include "errquit.fh"
#include "stdio.fh"
#include "rtdb.fh"
#include "geom.fh"
#include "bas.fh"
#include "mafdecls.fh"
#include "../nwdft/include/cdft.fh"
#include "../nwdft/grid/cgridfile.fh"
#include "wfn1_param.fh"
!> General calculation parameters
      type(wfn1_prm),     intent(in) :: wfn1_param
!> Calculation instances such as the runtime database, geometry, etc.
      type(wfn1_prminst), intent(in) :: wfn1_inst
c
c     Local
c
      integer      :: rtdb
      integer      :: icenter, jcenter
      integer      :: icharge
      integer      :: lcharge
      integer      :: icoord
      integer      :: lcoord
      integer      :: itags
      integer      :: ltags
      integer      :: mxnang
      integer      :: mxnrad
      logical      :: isbq
      logical      :: lnewtype
      logical      :: same_atom
      logical      :: same_bq
      character*40 :: blank_40c
      character*11 pname
c
      logical, external :: atom_tag_check
      logical, external :: xc_gotxc
c
      pname = 'wfn1_inpdft'
      rtdb  = wfn1_inst%rtdb
      geom  = wfn1_inst%geom
c
c     The WFN1 code only supports unrestricted wavefunctions
c
cDEBUG
      write(*,*)'HvD: wfn1_inpdft'
cDEBUG
      ipol         = 2
      mult         = wfn1_param%nea-wfn1_param%neb+1
      noc(1)       = wfn1_param%nea
      noc(2)       = wfn1_param%neb
      xfac         = 0.0d0
      cfac         = 0.0d0
      lcfac        = .false.
      nlcfac       = .false.
      lxfac        = .false.
      nlxfac       = .false.
      xccomb       = .false.
      xcfit        = .false.
      cdfit        = .false.
      leb          = .false.
      grid_written = .false.
c
c     if (.not. rtdb_get(rtdb, 'dft:no_prune', mt_log, 1,
c    &                   no_prune))no_prune = .false.
      if (.not. rtdb_get(rtdb, 'dft:nquad_task', mt_int, 1, nquad_task))
     &then
         nquad_task = 1
         if (.not. rtdb_put(rtdb, 'dft:nquad_task', mt_int, 1,
     &                      nquad_task))
     &     call errquit('wfn1_inpdft: dft:nquad_task put failed', 123,
     &                  RTDB_ERR)
      endif
      if (.not. rtdb_get(rtdb, 'dft:iangquad', mt_int, 1, iangquad))then
         iangquad=0
         if (.not. rtdb_put(rtdb, 'dft:iangquad', mt_int, 1, iangquad))
     &      call errquit('wfn1_inpdft: dft:iangquad put failed', 129,
     &                   RTDB_ERR)
      endif
c
c     Geometry
c
      if (.not. geom_ncent(geom, ncenters))
     &  call errquit(pname//' geom_ncent failed',73, GEOM_ERR)
c
c     allocate space for atomic coordinates and charges
c
      if (.not. Ma_Push_Get(MT_Dbl,ncenters*3,'coordinates',lcoord,
     &  icoord))call errquit(pname//'failed to alloc coordinates',0,
     &                       MA_ERR)
      if (.not. Ma_Push_Get(MT_Dbl,ncenters,'charges',lcharge,
     &  icharge))call errquit(pname//'failed to alloc charges',0,MA_ERR)
      if (.not. Ma_Push_Get(MT_Byte, ncenters*16, 'K center tags',
     &  ltags, itags))call errquit(
     &  pname//'failed to alloc center tags',0, MA_ERR)
c     
      if (.not. geom_cart_get(geom, ncenters, Byte_MB(itags),
     &                        Dbl_MB(icoord), Dbl_MB(icharge)))
     &  call errquit(pname//'geom_cart_get failed',74, GEOM_ERR)
c
c     set up number of atom types and atom type array iatype(icenter)
c
      ntypes = 0
      do icenter = 1, ncenters
c
c       is this a new type of atom?
c
        isbq = geom_isbq(geom,icenter)
        lnewtype = .true.
        do jcenter = 1, icenter - 1
          same_atom = Dbl_MB(icharge + icenter - 1) .eq.
     &                Dbl_MB(icharge + jcenter - 1)
          same_bq = geom_isbq(geom,jcenter) .and. isbq
          same_atom = same_atom .or. same_bq
          if (same_atom .and.
     &        atom_tag_check(Byte_MB(itags + (icenter - 1)*16),
     &                       Byte_MB(itags + (jcenter - 1)*16))
     &       )then   ! same atom type
               lnewtype = .false.
               iatype(icenter) = iatype(jcenter)
               goto 100
          endif
        enddo
  100   continue
        if (lnewtype)then
          ntypes = ntypes + 1
          iatype(icenter) = ntypes
        endif
      enddo
      if (ntypes.gt.dft_ntags_bsmx)then
        write(LuOut,*) 'RDINPUT:  Too many types of atoms.'
        call errquit(' rdinput: raise dft_ntags_bsmx',2, INPUT_ERR)
      end if
c
c     Basis set
c
      AO_bas_han = wfn1_inst%basis
      if (.not. bas_numbf(AO_bas_han, nbf))
     &   call errquit(pname//'basis set error:', 86, BASIS_ERR)
      nbf_ao = nbf
      if (.not. bas_nprim_cn_max(AO_bas_han,nbf_ao_mxprim))
     &   call errquit(pname//'basis set error:', 86, BASIS_ERR)
      if (.not. bas_high_angular(AO_bas_han,nbf_ao_mxang))
     &   call errquit(pname//'basis set error:', 86, BASIS_ERR)
      if (.not. bas_ncontr_cn_max(AO_bas_han,nbf_ao_mxcont))
     &   call errquit(pname//'basis set error:', 86, BASIS_ERR)
      if (.not. bas_nbf_cn_max(AO_bas_han,nbf_ao_mxnbf_cn))
     &   call errquit(pname//'basis set error:', 86, BASIS_ERR)
      if (.not. bas_nbf_ce_max(AO_bas_han,nbf_ao_mxnbf_ce))
     &   call errquit(pname//'basis set error:', 86, BASIS_ERR)
      if (.not. bas_numcont(AO_bas_han,nshells_ao))
     &   call errquit(pname//'basis set error:', 86, BASIS_ERR)
c
c     Quadrature Grids
c
c     - General grid definition
c
      if (.not. rtdb_cget(rtdb, 'dft:gridtype', 1, gridtype)) then
        gridtype = 'user'
        if (.not. rtdb_cput(rtdb, 'dft:gridtype', 1, gridtype)) then
          call errquit('wfn1_inpdft: put of gridtype failed',0,
     &                 RTDB_ERR)
        endif
      endif
c
c     - Angular grid (Lebedev by default)
c
      if (.not. rtdb_get(rtdb, 'dft:lebedev', mt_log, 1, leb)) then
        leb=.true.
        if (.not.rtdb_put(rtdb, 'dft:lebedev', mt_log, 1, leb))
     &    call errquit('wfn1_inpdft: put of leb failed',0,
     &                 RTDB_ERR)
      endif
c
c     - Radial grid
c
      if (.not. rtdb_cget(rtdb,'dft:radgrid', 1,radgrid)) then
        if(leb) then
          radgrid='murakn'
        else
          radgrid='eumac'
        endif
        if (.not. rtdb_cput(rtdb,'dft:radgrid', 1,radgrid))
     &    call errquit('wfn1_inpdft: put radgrid failed',7, RTDB_ERR)
      endif
      if (.not.rtdb_get(rtdb, 'dft:trunc_neigh', mt_log, 1,
     &  trunc_neigh)) trunc_neigh=.true.
c
c     - Weighting scheme (erf1 by default)
c
      if (.not.rtdb_get(rtdb, 'dft:ssw', mt_log, 1, lssw)) then
        lssw=.true.
        if (.not.rtdb_put(rtdb, 'dft:ssw', mt_log, 1, lssw))
     &    call errquit('wfn1_inpdft: rtdb_put failed', 1710, RTDB_ERR)
        whichssw='erf1'
        if (.not.rtdb_cput(rtdb, 'dft:ssweights',  1, whichssw))
     &    call errquit('wfn1_inpdft: rtdb_cput failed', 1720,
     &                 RTDB_ERR)
      else
        if(lssw) then
          if (.not.rtdb_cget(rtdb, 'dft:ssweights',  1, whichssw))
     &      call errquit('wfn1_inpdft: rtdb_cget failed',1700,
     &                   RTDB_ERR)
        endif
      endif
c     endif
c
c     Functionals
c
      if (.not. rtdb_get(rtdb, 'dft:cfac', mt_dbl, numfunc, cfac))then
        if (.not. rtdb_put(rtdb, 'dft:cfac', mt_dbl, numfunc, cfac))
     &    call errquit('wfn1_inpdft: dft:cfac put failed', 21, RTDB_ERR)
      endif
c
      if (.not. rtdb_get(rtdb, 'dft:lcfac', mt_log, numfunc, lcfac))then
        if (.not. rtdb_put(rtdb, 'dft:lcfac', mt_log, numfunc, lcfac))
     &    call errquit('wfn1_inpdft: dft:lcfac put failed', 9, RTDB_ERR)
      endif
c
      if (.not. rtdb_get(rtdb, 'dft:nlcfac', mt_log, numfunc, nlcfac))
     &then
        if (.not. rtdb_put(rtdb, 'dft:nlcfac', mt_log, numfunc, nlcfac))
     &      call errquit('wfn1_inpdft: dft:nlcfac put failed', 10,
     &                   RTDB_ERR)
      endif
c
      if (.not. rtdb_put(rtdb, 'dft:xfac', mt_dbl, numfunc, xfac))
     &  call errquit('wfn1_inpdft: dft:xfac put failed', 211, RTDB_ERR)
      if (.not. rtdb_put(rtdb, 'dft:lxfac', mt_log, numfunc, lxfac))
     &  call errquit('wfn1_inpdft: dft:lxfac put failed', 11, RTDB_ERR)
      if (.not. rtdb_put(rtdb, 'dft:nlxfac', mt_log, numfunc, nlxfac))
     &  call errquit('wfn1_inpdft: dft:nlxfac put failed', 12, RTDB_ERR)
      if (.not. rtdb_put(rtdb, 'dft:xccomb', mt_log, numfunc, xccomb))
     &  call errquit('wfn1_inpdft: dft:xccomb put failed', 13, RTDB_ERR)
c
      if (.not. rtdb_put(rtdb, 'dft:ipol', mt_int, 1, ipol))
     &  call errquit('wfn1_inpdft: dft:ipol put failed', 0, RTDB_ERR)
      if (.not. rtdb_put(rtdb, 'dft:noc', mt_int, 2, noc))
     &  call errquit('wfn1_inpdft: dft:noc put failed', 0, RTDB_ERR)
c
      blank_40c = '                                        '
      cname  = blank_40c
      xname  = blank_40c
      xcname = blank_40c
c
c
c     Correlation functional name defaults.
c
      cname(1)='            VWN V Correlation Functional'
      cname(2)='    Lee-Yang-Parr Correlation Functional'
      cname(3)='      Perdew 1981 Correlation Functional'
      cname(4)='      Perdew 1986 Correlation Functional'
      cname(5)='      Perdew 1991 Correlation Functional'
      cname(6)='  Perdew 1991 LDA Correlation Functional'
      cname(7)='        VWN I RPA Correlation Functional'
      cname(8)='            VWN I Correlation Functional'
      cname(9)='           VWN II Correlation Functional'
      cname(10)='          VWN III Correlation Functional'
      cname(11)='           VWN IV Correlation Functional'
      cname(12)=' PerdewBurkeErnz. Correlation Functional'
      cname(13)='            HCTH  Correlation Functional'
      cname(14)='        Becke 1997 Correlation Potential'
      cname(15)='      Becke 1997-1 Correlation Potential'
      cname(16)='         HCTH120  Correlation Functional'
      cname(17)='         HCTH147  Correlation Functional'
      cname(18)='       Becke 1998 Correlation Functional'
      cname(19)='   Becke 1997GGA-1 Correlation Potential'
      cname(20)='         HCTH407  Correlation Functional'
      cname(21)='         HCTHP1/4 Correlation Functional'
      cname(22)='            OPTC     Exchange Functional'
      cname(23)='      Becke 1997-2 Correlation Potential'
      cname(24)='  Filatov Thiel 1997 Method C Functional'
      cname(25)='   PKZB99 metaGGA Correlation Functional'
      cname(26)='         HCTH407+ Correlation Functional'
      cname(27)='   TPSS03 metaGGA Correlation Functional'
      cname(28)='      Becke 1997-3 Correlation Potential'
      cname(29)='             VS98 Correlation Functional'
      cname(30)='             M06-L Correlation Potential'
      cname(31)='               B95 Correlation Potential'
      cname(32)='            PW6B95 Correlation Potential'
      cname(33)='             PWB6K Correlation Potential'
      cname(34)='               M05 Correlation Potential'
      cname(35)='            M05-2X Correlation Potential'
      cname(36)='        OP_Becke88 Correlation Potential'
      cname(37)='          OP_PBE96 Correlation Potential'
      cname(38)='            M06-HF Correlation Potential'
      cname(39)='               M06 Correlation Potential'
      cname(40)='            M06-2X Correlation Potential'
      cname(41)='      Becke 1997-D Correlation Potential'
      cname(42)='              dlDF Correlation Potential'
      cname(46)='    simplified PBE Correlation Potential'
      cname(47)='        ------  NO Correlation ---------'
      cname(48)='            M08-HX Correlation Potential'
      cname(49)='            M08-SO Correlation Potential'
      cname(50)='               M11 Correlation Potential'
      cname(51)='             M11-L Correlation Potential'
      cname(52)='             SOGGA Correlation Potential'
      cname(53)='           SOGGA11 Correlation Potential'
      cname(54)='         SOGGA11-X Correlation Potential'
      cname(64)='revTPSSvariant PBE Correlation Potential'
c
c     Exchange functional name defaults.
c
      xname(1)='           Hartree-Fock (Exact) Exchange'
      xname(2)='              Slater Exchange Functional'
      xname(3)='          Becke 1988 Exchange Functional'
      xname(4)='               HCTH  Exchange Functional'
      xname(5)='          Becke 1997 Exchange Functional'
      xname(6)='        Becke 1997-1 Exchange Functional'
      xname(7)='PerdewBurkeErnzerhof Exchange Functional'
      xname(8)='         Gill 1996   Exchange Functional'
      xname(9)='       Perdew 1991   Exchange Functional'
      xname(10)='            HCTH120  Exchange Functional'
      xname(11)='            HCTH147  Exchange Functional'
      xname(12)='         Becke 1998  Exchange Functional'
      xname(13)='     Becke 1997GGA-1 Exchange Functional'
      xname(14)='            HCTH407  Exchange Functional'
      xname(15)='           HCTHP1/4  Exchange Functional'
      xname(16)='            OPTX     Exchange Functional'
      xname(17)='       Becke 1997-2  Exchange Functional'
      xname(18)='        PKZB metaGGA Exchange Functional'
      xname(19)='              mPW91  Exchange Functional'
      xname(20)='           HCTH407+  Exchange Functional'
      xname(21)='        TPSS metaGGA Exchange Functional'
      xname(22)='       Becke 1997-3  Exchange Functional'
      xname(25)='  Filatov Thiel 1997 Method X Functional'
      xname(26)='          PW6B95 GGA Exchange Functional'
      xname(27)='           PWB6K GGA Exchange Functional'
      xname(28)='         M05 metaGGA Exchange Functional'
      xname(29)='      M05-2X metaGGA Exchange Functional'
      xname(30)='                RPBE Exchange Functional'
      xname(31)='              RevPBE Exchange Functional'
      xname(32)='                dlDF Exchange Functional'
      xname(33)='                VS98 Exchange Functional'
      xname(34)='               M06-L Exchange Functional'
      xname(35)='              M06-HF Exchange Functional'
      xname(36)='                 M06 Exchange Functional'
      xname(37)='              M06-2X Exchange Functional'
      xname(38)='         LC-BNL 2007 Exchange Functional'
      xname(39)='        Becke 1997-D Exchange Functional'
      xname(40)='         CAM-Becke88 Exchange Functional'
      xname(41)='           CAM-PBE96 Exchange Functional'
      xname(42)='           CAM-LSD   Exchange Functional'
      xname(43)='           LC-wPBE(h)Exchange Functional'
      xname(46)='               SSB-1 Exchange Functional'
      xname(47)='                 KT1 gradient correction'
      xname(48)='              M08-HX gradient correction'
      xname(49)='              M08-SO gradient correction'
      xname(50)='                 M11 gradient correction'
      xname(51)='               M11-L gradient correction'
      xname(52)='               SOGGA gradient correction'
      xname(53)='             SOGGA11 gradient correction'
      xname(54)='           SOGGA11-X gradient correction'
      xname(55)='          Becke 1986 Exchange Functional'
      xname(56)='         Perdew 1986 Exchange Functional'
      xname(60)='      S12g(GGA-part) Exchange Functional'
      xname(61)='      S12h(GGA-part) Exchange Functional'
      xname(62)='  CAM-S12g(GGA-part) Exchange Functional'
      xname(63)='  CAM-S12h(GGA-part) Exchange Functional'
      xname(64)='         MVS metaGGA Exchange Functional'
c
c     Exchange-Correlation combination functional name defaults.
c
      xcname(1)='Adiabatic Connection Method XC Potential'
      xcname(2)='               B3LYP Method XC Potential'
      xcname(3)=' Becke half-and-half Method XC Potential'
      xcname(4)='          Becke 1997 Method XC Potential'
      xcname(5)='        Becke 1997-1 Method XC Potential'
      xcname(6)='              HCTH  Method XC Functional'
      xcname(7)='           HCTH120  Method XC Functional'
      xcname(8)='           HCTH147  Method XC Functional'
      xcname(9)='         Becke 1998 Method XC Functional'
      xcname(10)='               PBE0 Method XC Functional'
      xcname(11)='     Becke 1997GGA-1 Method XC Potential'
      xcname(12)='            HCTH407 Method XC Functional'
      xcname(13)='           HCTHP1/4 Method XC Functional'
      xcname(14)='         Becke 97-2 Method XC Functional'
      xcname(15)=' Filatov Thiel 1997 Method XC Functional'
      xcname(16)='              MPW1K Method XC Functional'
      xcname(17)='           HCTH407+ Method XC Functional'
      xcname(18)='        TPSS hybrid Method XC Functional'
      xcname(19)='         Becke 97-3 Method XC Functional'
      xcname(21)='            MPW1B95 Method XC Functional'
      xcname(22)='             MPWB1K Method XC Functional'
      xcname(23)='              B1B95 Method XC Functional'
      xcname(24)='               BB1K Method XC Functional'
      xcname(25)='             PW6B95 Method XC Functional'
      xcname(26)='              PWB6K Method XC Functional'
      xcname(27)='                BOP Method XC Functional'
      xcname(28)='              PBEOP Method XC Functional'
      xcname(29)='                M05 Method XC Functional'
      xcname(30)='             M05-2X Method XC Functional'
      xcname(31)='               dlDF Method XC Functional'
      xcname(32)='               VS98 Method XC Functional'
      xcname(33)='              M06-L Method XC Functional'
      xcname(34)='             M06-HF Method XC Functional'
      xcname(35)='                M06 Method XC Functional'
      xcname(36)='             M06-2X Method XC Functional'
      xcname(37)='         Becke 97-D Method XC Functional'
      xcname(38)='             LC-BNL Method XC Functional'
      xcname(46)='              SSB-D Method XC Functional'
      xcname(47)='                KT2 Method XC Functional'
      xcname(48)='             M08-HX Method XC Functional'
      xcname(49)='             M08-SO Method XC Functional'
      xcname(50)='                M11 Method XC Functional'
      xcname(51)='              M11-L Method XC Functional'
      xcname(52)='              SOGGA Method XC Functional'
      xcname(53)='            SOGGA11 Method XC Functional'
      xcname(54)='          SOGGA11-X Method XC Functional'
      xcname(55)='              BHLYP Method XC Functional'
      xcname(56)='              B3P86 Method XC Functional'
      xcname(57)='             B3PW91 Method XC Functional'
      xcname(58)='              PBE96 Method XC Functional'
      xcname(59)='              HSE03 Method XC Functional'
      xcname(60)='               S12g Method XC Functional'
      xcname(61)='               S12h Method XC Functional'
      xcname(62)='           CAM-S12g Method XC Functional'
      xcname(63)='           CAM-S12h Method XC Functional'
      xcname(64)=' MadeVerySimple MVS Method XC Functional'
      xcname(65)='              HLE16 Method XC Functional'
      if (.not. rtdb_cput(rtdb, 'dft:cname', numfunc, cname))
     &  call errquit('wfn1_inpdft: dft:cname put failed', 302, RTDB_ERR)
      if (.not. rtdb_cput(rtdb, 'dft:xname', numfunc, xname))
     &  call errquit('wfn1_inpdft: dft:xname put failed', 303, RTDB_ERR)
      if (.not. rtdb_cput(rtdb, 'dft:xcname', numfunc, xcname))
     &  call errquit('wfn1_inpdft: dft:xcname put failed',304, RTDB_ERR)
c
      if (.not. MA_Pop_Stack(ltags))
     &   call errquit('wfn1_inpdft: pop stack failed.',0, MA_ERR)
      if (.not. MA_Pop_Stack(lcharge))
     &   call errquit('wfn1_inpdft: pop stack failed.',0, MA_ERR)
      if (.not. MA_Pop_Stack(lcoord))
     &   call errquit('wfn1_inpdft: pop stack failed.',0, MA_ERR)
c
      mxnrad = 0
      mxnang = 0
      if (.not. rtdb_get(rtdb, 'dft:mxnrad', mt_int, 1, mxnrad)) then
         if (xc_gotxc())
     &       call errquit('gridv0: rtdb_get failed', 119, RTDB_ERR)
      endif
      if (.not. rtdb_get(rtdb, 'dft:mxnang', mt_int, 1, mxnang)) then
         if (xc_gotxc())
     &       call errquit('gridv0: rtdb_get failed', 119, RTDB_ERR)
      endif
cDEBUG
      write(*,*)'HvD: nquad_task = ',nquad_task
      write(*,*)'HvD: iangquad   = ',iangquad
      write(*,*)'HvD: gridtype   = ',gridtype
      write(*,*)'HvD: leb        = ',leb
      write(*,*)'HvD: radgrid    = ',radgrid
      write(*,*)'HvD: trunc_neigh= ',trunc_neigh
      write(*,*)'HvD: lssw       = ',lssw
      write(*,*)'HvD: whichssw   = ',whichssw
      write(*,*)'HvD: mxnrad     = ',mxnrad
      write(*,*)'HvD: mxnang     = ',mxnang
cDEBUG
c
      end subroutine wfn1_inpdft
c
c-----------------------------------------------------------------------
