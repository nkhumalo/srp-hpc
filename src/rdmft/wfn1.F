!     submodule (wfn1) wfn1_for_nw
!     implicit none
!     contains
      subroutine wfn1_input(rtdb)
      implicit none
!>
!> \brief Read the input for the RDMFT calculations
!>
      integer, intent(inout) :: rtdb !> [Input] The runtime database
      end subroutine wfn1_input
!>
!>----------------------------------------------------------------------
!>
!> \brief Minimize the free energy
!>
!> The wavefunction we are using is outlined in [1]. The energy
!> expression is quite different though. The reason is that attempts to
!> reconstruct a valid 2-electron density matrix from an arbitrary
!> 1-electron density matrix never work. Instead when starting from 
!> any wavefunction the 2-electron density matrix is N-representable
!> by construction (this includes the wavefunction from [1]). Hence
!> it would seem that the correlation energy can be more easily
!> formulated in other ways than in terms of 1-electron density
!> matrices. For this reason RDMFT was abandoned and we explore a
!> rather different approach here.
!>
!> \return `true`: if the energy minimization converged, and `false`
!> otherwise.
!>
!> ### References ###
!>
!> [1] Hubertus van Dam,
!>     "Generalization of the Kohn-Sham system enabling representing
!>     arbitary one electron density matrices",
!>     arXiv: <a href="http://arxiv.org/abs/1512.06120">
!>     http://arxiv.org/abs/1512.06120</a>
!>
      function wfn1_task_energy(rtdb) result(wfn1_success)
      implicit none
!include "wfn1_param.fh"
!include "wfn1_wfn.fh"
!include "wfn1_drv.fh"
!include "wfn1_matcnst.fh"
!include "wfn1_size.fh"
#include "wfn1_trace.fh"
!include "wfn1_en.fh"
#include "bas.fh"
#include "errquit.fh"
#include "geom.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "stdio.fh"
#include "util.fh"
! For filesize and memsize
#include "cfock.fh"
      integer, intent(inout) :: rtdb !> [Input] The runtime database
      logical :: wfn1_success
c
c     Local variables
c
      logical lstatus ! the status of the runtime database
      logical oprint  ! should something be printed
      logical oprerr  ! should error be printed
      integer iproc   ! the processor rank
      integer nproc   ! the number of processors
      integer nelec   ! the total number of electrons
      integer nopen   ! the number of unpaired electrons
      integer nclosed ! the number of closed shell orbitals
      integer maxit   ! the maximum number of iterations
      integer curit   ! the current iteration number
c
      logical oconverged ! is this calculation converged
c
      double precision charge     ! the molecular charge
      double precision charge_nuc ! the total nuclear charge
c
      character*12  pname ! the subroutine name
      character*256 title ! the job title
      character*30  tag   ! the RTDB key word
c
      integer g_bfbf1    ! nbf x nbf scratch global array
      integer g_bfbf2    ! nbf x nbf scratch global array
      integer g_bfbf3    ! nbf x nbf scratch global array
      integer g_bfbf4    ! nbf x nbf scratch global array
      integer g_bfmo     ! nbf x nmo scratch global array
      integer g_momo1    ! nmo x nmo scratch global array
      integer g_momo2    ! nmo x nmo scratch global array
      integer g_momo3    ! nmo x nmo scratch global array
      integer g_momo4    ! nmo x nmo scratch global array
c
      integer ndep_nmo          ! the number of linearly independent
                                ! vectors set as a user parameter
      double precision ndep_tol ! the tolerance on the eigenvalues of
                                ! the overlap matrix set as a user
                                ! parameter
      double precision ndep_tol_default ! default setting for ndep_tol
      parameter (ndep_tol_default = 1.0d-5)
c
      double precision conv_tol ! the convergence tolerance
c
      double precision, allocatable :: eig(:)
      double precision, allocatable :: v_dna(:)
      double precision, allocatable :: v_dnb(:)
      double precision, allocatable :: v_dca(:)
      double precision, allocatable :: v_dcb(:)
      double precision, allocatable :: v_dnna(:)
      double precision, allocatable :: v_dnnb(:)
      double precision, allocatable :: v_dcna(:)
      double precision, allocatable :: v_dcnb(:)
      type(wfn1_error)              :: error
c
      type(wfn1_prminst) :: wfn1_inst    ! paramater instances
      type(wfn1_prm)     :: wfn1_param   ! calculation parameters
      type(wfn1_matcnst) :: wfn1_mat     ! matrix constants
      type(wfn1_wfn)     :: wfn1_wave    ! the current wave function
      type(wfn1_drv)     :: wfn1_deriv   ! the derivative of the current
                                         ! wave function
      type(wfn1_en)      :: wfn1_energy  ! the energy terms
c
      type(wfn1_wfn)     :: wfn1_trial_wv  ! the trial wave function
      type(wfn1_drv)     :: wfn1_trial_drv ! the derivative wrt the trial
                                           ! wave function
      type(wfn1_drv)     :: wfn1_part_drv ! the non-correlation, non-
                                          ! thermal part of derivative
                                          ! wrt the wave function
      type(wfn1_en)      :: wfn1_trial_en  ! the trial energy
c
      type(wfn1_wfn) :: wfn1_t_wv  ! temporary variable
      type(wfn1_drv) :: wfn1_t_drv ! temporary variable
c
      logical, external :: int_normalize
!     logical, external :: wfn1_solve_line
!     logical, external :: wfn1_solve_diag
!     logical, external :: wfn1_solve_mc
c
      double precision aa_hl, ab_hl
      double precision fac_redist
      double precision fac_aa,         pow_x1_aa, pow_x2_aa,
     +                                 pow_y1_aa, pow_z1_aa,
     +                                 pow_y2_aa, pow_z2_aa,
     +                                 pow_y3_aa, pow_z3_aa,
     +                                 pow_y4_aa, pow_z4_aa,
     +                                 pow_y5_aa, pow_z5_aa
      double precision fac_ab, del_ab, pow_x1_ab, pow_x2_ab,
     +                                 pow_y1_ab, pow_z1_ab,
     +                                 pow_y2_ab, pow_z2_ab,
     +                                 pow_y3_ab, pow_z3_ab,
     +                                 pow_y4_ab, pow_z4_ab,
     +                                 pow_y5_ab, pow_z5_ab
      double precision                 powhf_x1_aa, powhf_x2_aa,
     +                                 powhf_y1_aa, powhf_z1_aa,
     +                                 powhf_y2_aa, powhf_z2_aa,
     +                                 powhf_y3_aa, powhf_z3_aa,
     +                                 powhf_y4_aa, powhf_z4_aa,
     +                                 powhf_y5_aa, powhf_z5_aa
      double precision                 powhf_x1_ab, powhf_x2_ab,
     +                                 powhf_y1_ab, powhf_z1_ab,
     +                                 powhf_y2_ab, powhf_z2_ab,
     +                                 powhf_y3_ab, powhf_z3_ab,
     +                                 powhf_y4_ab, powhf_z4_ab,
     +                                 powhf_y5_ab, powhf_z5_ab
c
      iproc = ga_nodeid()
      nproc = ga_nnodes()
c
      wfn1_inst%rtdb = rtdb
      pname = "wfn1_energy:"
      wfn1_success = .false.
      lstatus = rtdb_parallel(.true.)
      call util_print_push
      call util_print_rtdb_load(wfn1_inst%rtdb, 'wfn1') 
c
c     == Print Header ==
c
      oprint = util_print('information',print_low)
      if (oprint) then
        if (.not. rtdb_cget(wfn1_inst%rtdb, 'title', 1, title))
     &       title = ' '
        if(iproc.eq.0) then
          write(LuOut,*)
          call util_print_centered(
     &         LuOut, 'NWChem WFN1 Module', 40, .true.)
          write(LuOut,*)
          write(LuOut,*)
          if (title .ne. ' ') then
            call util_print_centered(LuOut, title, 40, .false.)
            write(LuOut,*)
            write(LuOut,*)
          endif
        endif
      endif
c
c     == Rig some settings that we need ==
c
      if (.not. rtdb_put(rtdb, 'bgj:scf_type', MT_INT, 1, 1))
     &  call errquit(pname//" put of bgj:scf_type failed",0, RTDB_ERR)
c     ---- Switch to DIRECT ----
      filesize = -1
      memsize  = -1
      if (.not. rtdb_put(rtdb,'int2e:filesize', MT_INT, 1, filesize))
     $  call errquit(pname//" failed to store filesize",0, RTDB_ERR)
      if (.not. rtdb_put(rtdb,'int2e:memsize', MT_INT, 1, memsize))
     $  call errquit(pname//" failed to store filesize",0, RTDB_ERR)
c
c     == Load geometry and basis set ==
c
      if (.not.geom_create(wfn1_inst%geom,'geometry'))
     &  call errquit(pname//' geom_create failed',0,GEOM_ERR)
      if (.not.geom_rtdb_load(wfn1_inst%rtdb,wfn1_inst%geom,'geometry'))
     &  call errquit(pname//' no geometry found',0,RTDB_ERR)
      if (.not.bas_create(wfn1_inst%basis,'ao basis'))
     &  call errquit(pname//' bas_create failed',0, BASIS_ERR)
      if (.not.bas_rtdb_load(wfn1_inst%rtdb,wfn1_inst%geom,
     &                       wfn1_inst%basis,'ao basis'))
     &  call errquit(pname//' no ao basis found',0, RTDB_ERR)
      if (.not.bas_numbf(wfn1_inst%basis,wfn1_param%nbf))
     &  call errquit(pname//" no basis functions",0,BASIS_ERR)
      if (.not.int_normalize(wfn1_inst%rtdb,wfn1_inst%basis))
     &  call errquit(pname//' int_normalize failed', 0, INT_ERR)
c
c     == Load linearly independence settings ==
c
      if (rtdb_get(rtdb, 'lindep:n_dep', MT_INT, 1, ndep_nmo)) then
        if (oprint) then
          write(6, 2) ndep_nmo
 2        format(/' Forcing the number of linearly dependent ',
     &            'vectors to be', i6)
          call util_flush(6)
        endif
      else
        ndep_nmo = -1
      endif
      if (rtdb_get(rtdb, 'lindep:tol', MT_DBL, 1, ndep_tol)) then
        if (oprint) then
          write(6, 3) ndep_tol
 3        format(/' Setting tolerance for linear dependence to',
     &          1p,d9.2)
          call util_flush(6)
        endif
      else
         ndep_tol = -1.0d0
      endif
c
      if ((ndep_tol.ne.-1.0d0) .and. (ndep_nmo.ne.-1)) call errquit
     &  (pname//' only one of lindep:n_dep and lindep:tol permitted',0,
     &  UNKNOWN_ERR)
c
      if (ndep_nmo .eq. -1) then
        if (ndep_tol .le. -1.0d0) ndep_tol = ndep_tol_default
      endif
c
c     == Print the geometry ==
c
      if(iproc.eq.0.and.util_print('geometry', print_high)) then
        if(.not.geom_print(wfn1_inst%geom)) call errquit(
     &       pname//'geom print failed',0, GEOM_ERR)
      endif
c
c     == Evaluate the nuclear repulsion energy ==
c
      lstatus = geom_nuc_rep_energy(wfn1_inst%geom, wfn1_param%Enuc)
      if (.not.lstatus)
     &  call errquit(pname//'Error in nuclear repulsion',1, GEOM_ERR)
c
c     == Work out numbers of electrons ==
c
      if (.not.geom_nuc_charge(wfn1_inst%geom,charge_nuc))
     &  call errquit(pname//" geometry nuclear charge failed",
     &               0,GEOM_ERR)
      if (.not.rtdb_get(wfn1_inst%rtdb,"charge",MT_DBL,1,charge))
     &  charge = 0.0d0
      if (.not. rtdb_get(wfn1_inst%rtdb,'wfn1:print_error',MT_LOG,1,
     &                   oprerr))
     &  oprerr = .false.
      if (.not. rtdb_get(wfn1_inst%rtdb,'wfn1:nopen',MT_INT,1,nopen))
     &  nopen = 0
      if (.not. rtdb_get(wfn1_inst%rtdb,'wfn1:ncorea',MT_INT,1,
     &                   wfn1_param%ncorea))
     &  wfn1_param%ncorea = 0
      if (.not. rtdb_get(wfn1_inst%rtdb,'wfn1:ncoreb',MT_INT,1,
     &                   wfn1_param%ncoreb))
     &  wfn1_param%ncoreb = 0
      if (.not. rtdb_get(wfn1_inst%rtdb,'wfn1:nvirta',MT_INT,1,
     &                   wfn1_param%nvirta))
     &  wfn1_param%nvirta = 0
      if (.not. rtdb_get(wfn1_inst%rtdb,'wfn1:nvirtb',MT_INT,1,
     &                   wfn1_param%nvirtb))
     &  wfn1_param%nvirtb = 0
      if (.not. rtdb_get(wfn1_inst%rtdb,'wfn1:fraca',MT_DBL,1,
     &                   wfn1_param%fraca))
     &  wfn1_param%fraca = 1.0d0
      if (.not. rtdb_get(wfn1_inst%rtdb,'wfn1:fracb',MT_DBL,1,
     &                   wfn1_param%fracb))
     &  wfn1_param%fracb = 1.0d0
      if (wfn1_param%fraca.lt.0.0d0) then
        call errquit(pname//" negative fractional alpha electron",
     &               0,INPUT_ERR)
      endif
      if (wfn1_param%fraca.gt.1.0d0) then
        call errquit(pname//" fractional alpha electron greater "//
     &               "than one",0,INPUT_ERR)
      endif
      if (wfn1_param%fracb.lt.0.0d0) then
        call errquit(pname//" negative fractional beta electron",
     &               0,INPUT_ERR)
      endif
      if (wfn1_param%fracb.gt.1.0d0) then
        call errquit(pname//" fractional beta electron greater "//
     &               "than one",0,INPUT_ERR)
      endif
      nelec = int(charge_nuc - charge)
      if (nelec.lt.0) then
        call errquit(pname//" negative number of electrons",
     &               nelec,INPUT_ERR)
      endif
      nclosed =  (nelec-nopen)/2
      wfn1_param%nea = nclosed + nopen
      wfn1_param%neb = nclosed
      if (nopen.eq.0.and.wfn1_param%nea+wfn1_param%neb.ne.nelec) then
        nopen = nelec-2*nclosed
        wfn1_param%nea = nclosed + nopen
      endif
c
c     == Get the correlation energy expression
c
      if (.not.rtdb_cget(wfn1_inst%rtdb,"wfn1:corr_expr",
     &                   len(wfn1_param%corr_expr),
     &                   wfn1_param%corr_expr)) then
         wfn1_param%corr_expr="wfn1_mx"
      endif
      if (.not. rtdb_get(wfn1_inst%rtdb,'wfn1:fac_a',MT_DBL,1,
     &                   wfn1_param%fac_a))
     &  wfn1_param%fac_a = 0.0d0
      if (wfn1_param%fac_a.lt.0.0d0) then
        call errquit(pname//" factor A less than 0", 0,INPUT_ERR)
      endif
      if (wfn1_param%fac_a.gt.1.0d0) then
        call errquit(pname//" factor A greater than 1", 0,INPUT_ERR)
      endif
      if (.not. rtdb_get(wfn1_inst%rtdb,'wfn1:aa',MT_DBL,1,
     &                   fac_aa))
     &  fac_aa = 0.0d0
      if (.not. rtdb_get(wfn1_inst%rtdb,'wfn1:aa_hl',MT_DBL,1,
     &                   aa_hl))
     &  aa_hl = 0.0d0
      if (.not. rtdb_get(wfn1_inst%rtdb,'wfn1:ab',MT_DBL,1,
     &                   fac_ab))
     &  fac_ab = 0.0d0
      if (.not. rtdb_get(wfn1_inst%rtdb,'wfn1:ab_hl',MT_DBL,1,
     &                   ab_hl))
     &  ab_hl = 0.0d0
      if (.not. rtdb_get(wfn1_inst%rtdb,'wfn1:fac_aa',MT_DBL,1,
     &                   fac_aa))
     &  fac_aa    = 0.0d0
      if (.not. rtdb_get(wfn1_inst%rtdb,'wfn1:fac_redist',MT_DBL,1,
     &                   fac_redist))
     &  fac_redist= 0.0d0
      if (.not. rtdb_get(wfn1_inst%rtdb,'wfn1:pow_x1_aa',MT_DBL,1,
     &                   pow_x1_aa))
     &  pow_x1_aa  = 0.0d0
      if (.not. rtdb_get(wfn1_inst%rtdb,'wfn1:pow_x2_aa',MT_DBL,1,
     &                   pow_x2_aa))
     &  pow_x2_aa  = 1.0d0
      if (.not. rtdb_get(wfn1_inst%rtdb,'wfn1:pow_y1_aa',MT_DBL,1,
     &                   pow_y1_aa))
     &  pow_y1_aa = 1.0d0
      if (.not. rtdb_get(wfn1_inst%rtdb,'wfn1:pow_y2_aa',MT_DBL,1,
     &                   pow_y2_aa))
     &  pow_y2_aa = 1.0d0
      if (.not. rtdb_get(wfn1_inst%rtdb,'wfn1:pow_y3_aa',MT_DBL,1,
     &                   pow_y3_aa))
     &  pow_y3_aa = 1.0d0
      if (.not. rtdb_get(wfn1_inst%rtdb,'wfn1:pow_y4_aa',MT_DBL,1,
     &                   pow_y4_aa))
     &  pow_y4_aa = 1.0d0
      if (.not. rtdb_get(wfn1_inst%rtdb,'wfn1:pow_y5_aa',MT_DBL,1,
     &                   pow_y5_aa))
     &  pow_y5_aa = 1.0d0
      if (.not. rtdb_get(wfn1_inst%rtdb,'wfn1:pow_z1_aa',MT_DBL,1,
     &                   pow_z1_aa))
     &  pow_z1_aa = 1.0d0
      if (.not. rtdb_get(wfn1_inst%rtdb,'wfn1:pow_z2_aa',MT_DBL,1,
     &                   pow_z2_aa))
     &  pow_z2_aa = 0.0d0
      if (.not. rtdb_get(wfn1_inst%rtdb,'wfn1:pow_z3_aa',MT_DBL,1,
     &                   pow_z3_aa))
     &  pow_z3_aa = 0.0d0
      if (.not. rtdb_get(wfn1_inst%rtdb,'wfn1:pow_z4_aa',MT_DBL,1,
     &                   pow_z4_aa))
     &  pow_z4_aa = 0.0d0
      if (.not. rtdb_get(wfn1_inst%rtdb,'wfn1:pow_z5_aa',MT_DBL,1,
     &                   pow_z5_aa))
     &  pow_z5_aa = 0.0d0
      if (.not. rtdb_get(wfn1_inst%rtdb,'wfn1:powhf_y1_aa',MT_DBL,1,
     &                   powhf_y1_aa))
     &  powhf_y1_aa = 1.0d0
      if (.not. rtdb_get(wfn1_inst%rtdb,'wfn1:powhf_y2_aa',MT_DBL,1,
     &                   powhf_y2_aa))
     &  powhf_y2_aa = 1.0d0
      if (.not. rtdb_get(wfn1_inst%rtdb,'wfn1:powhf_y3_aa',MT_DBL,1,
     &                   powhf_y3_aa))
     &  powhf_y3_aa = 1.0d0
      if (.not. rtdb_get(wfn1_inst%rtdb,'wfn1:powhf_y4_aa',MT_DBL,1,
     &                   powhf_y4_aa))
     &  powhf_y4_aa = 1.0d0
      if (.not. rtdb_get(wfn1_inst%rtdb,'wfn1:powhf_y5_aa',MT_DBL,1,
     &                   powhf_y5_aa))
     &  powhf_y5_aa = 1.0d0
      if (.not. rtdb_get(wfn1_inst%rtdb,'wfn1:powhf_z1_aa',MT_DBL,1,
     &                   powhf_z1_aa))
     &  powhf_z1_aa = 1.0d0
      if (.not. rtdb_get(wfn1_inst%rtdb,'wfn1:powhf_z2_aa',MT_DBL,1,
     &                   powhf_z2_aa))
     &  powhf_z2_aa = 0.0d0
      if (.not. rtdb_get(wfn1_inst%rtdb,'wfn1:powhf_z3_aa',MT_DBL,1,
     &                   powhf_z3_aa))
     &  powhf_z3_aa = 0.0d0
      if (.not. rtdb_get(wfn1_inst%rtdb,'wfn1:powhf_z4_aa',MT_DBL,1,
     &                   powhf_z4_aa))
     &  powhf_z4_aa = 0.0d0
      if (.not. rtdb_get(wfn1_inst%rtdb,'wfn1:powhf_z5_aa',MT_DBL,1,
     &                   powhf_z5_aa))
     &  powhf_z5_aa = 0.0d0
      if (.not. rtdb_get(wfn1_inst%rtdb,'wfn1:fac_ab',MT_DBL,1,
     &                   fac_ab))
     &  fac_ab    = 0.0d0
      if (.not. rtdb_get(wfn1_inst%rtdb,'wfn1:del_ab',MT_DBL,1,
     &                   del_ab))
     &  del_ab    = 0.0d0
      if (.not. rtdb_get(wfn1_inst%rtdb,'wfn1:pow_x1_ab',MT_DBL,1,
     &                   pow_x1_ab))
     &  pow_x1_ab  = 0.0d0
      if (.not. rtdb_get(wfn1_inst%rtdb,'wfn1:pow_x2_ab',MT_DBL,1,
     &                   pow_x2_ab))
     &  pow_x2_ab  = 1.0d0
      if (.not. rtdb_get(wfn1_inst%rtdb,'wfn1:pow_y1_ab',MT_DBL,1,
     &                   pow_y1_ab))
     &  pow_y1_ab = 1.0d0
      if (.not. rtdb_get(wfn1_inst%rtdb,'wfn1:pow_y2_ab',MT_DBL,1,
     &                   pow_y2_ab))
     &  pow_y2_ab = 1.0d0
      if (.not. rtdb_get(wfn1_inst%rtdb,'wfn1:pow_y3_ab',MT_DBL,1,
     &                   pow_y3_ab))
     &  pow_y3_ab = 1.0d0
      if (.not. rtdb_get(wfn1_inst%rtdb,'wfn1:pow_y4_ab',MT_DBL,1,
     &                   pow_y4_ab))
     &  pow_y4_ab = 1.0d0
      if (.not. rtdb_get(wfn1_inst%rtdb,'wfn1:pow_y5_ab',MT_DBL,1,
     &                   pow_y5_ab))
     &  pow_y5_ab = 1.0d0
      if (.not. rtdb_get(wfn1_inst%rtdb,'wfn1:pow_z1_ab',MT_DBL,1,
     &                   pow_z1_ab))
     &  pow_z1_ab = 0.0d0
      if (.not. rtdb_get(wfn1_inst%rtdb,'wfn1:pow_z2_ab',MT_DBL,1,
     &                   pow_z2_ab))
     &  pow_z2_ab = 0.0d0
      if (.not. rtdb_get(wfn1_inst%rtdb,'wfn1:pow_z3_ab',MT_DBL,1,
     &                   pow_z3_ab))
     &  pow_z3_ab = 0.0d0
      if (.not. rtdb_get(wfn1_inst%rtdb,'wfn1:pow_z4_ab',MT_DBL,1,
     &                   pow_z4_ab))
     &  pow_z4_ab = 0.0d0
      if (.not. rtdb_get(wfn1_inst%rtdb,'wfn1:pow_z5_ab',MT_DBL,1,
     &                   pow_z5_ab))
     &  pow_z5_ab = 0.0d0
      if (.not. rtdb_get(wfn1_inst%rtdb,'wfn1:powhf_y1_ab',MT_DBL,1,
     &                   powhf_y1_ab))
     &  powhf_y1_ab = 1.0d0
      if (.not. rtdb_get(wfn1_inst%rtdb,'wfn1:powhf_y2_ab',MT_DBL,1,
     &                   powhf_y2_ab))
     &  powhf_y2_ab = 1.0d0
      if (.not. rtdb_get(wfn1_inst%rtdb,'wfn1:powhf_y3_ab',MT_DBL,1,
     &                   powhf_y3_ab))
     &  powhf_y3_ab = 1.0d0
      if (.not. rtdb_get(wfn1_inst%rtdb,'wfn1:powhf_y4_ab',MT_DBL,1,
     &                   powhf_y4_ab))
     &  powhf_y4_ab = 1.0d0
      if (.not. rtdb_get(wfn1_inst%rtdb,'wfn1:powhf_y5_ab',MT_DBL,1,
     &                   powhf_y5_ab))
     &  powhf_y5_ab = 1.0d0
      if (.not. rtdb_get(wfn1_inst%rtdb,'wfn1:powhf_z1_ab',MT_DBL,1,
     &                   powhf_z1_ab))
     &  powhf_z1_ab = 0.0d0
      if (.not. rtdb_get(wfn1_inst%rtdb,'wfn1:powhf_z2_ab',MT_DBL,1,
     &                   powhf_z2_ab))
     &  powhf_z2_ab = 0.0d0
      if (.not. rtdb_get(wfn1_inst%rtdb,'wfn1:powhf_z3_ab',MT_DBL,1,
     &                   powhf_z3_ab))
     &  powhf_z3_ab = 0.0d0
      if (.not. rtdb_get(wfn1_inst%rtdb,'wfn1:powhf_z4_ab',MT_DBL,1,
     &                   powhf_z4_ab))
     &  powhf_z4_ab = 0.0d0
      if (.not. rtdb_get(wfn1_inst%rtdb,'wfn1:powhf_z5_ab',MT_DBL,1,
     &                   powhf_z5_ab))
     &  powhf_z5_ab = 0.0d0
      call wfn1_set_entropy_fac(fac_aa,fac_ab,fac_redist,del_ab,
     &                          pow_x1_aa,pow_x2_aa,pow_x1_ab,pow_x2_ab,
     &                          pow_y1_aa,pow_y1_ab,pow_y2_aa,pow_y2_ab,
     &                          pow_y3_aa,pow_y3_ab,pow_y4_aa,pow_y4_ab,
     &                          pow_y5_aa,pow_y5_ab,
     &                          pow_z1_aa,pow_z1_ab,pow_z2_aa,pow_z2_ab,
     &                          pow_z3_aa,pow_z3_ab,pow_z4_aa,pow_z4_ab,
     &                          pow_z5_aa,pow_z5_ab,
     &                          powhf_y1_aa,powhf_y1_ab,powhf_y2_aa,
     &                          powhf_y2_ab,powhf_y3_aa,powhf_y3_ab,
     &                          powhf_y4_aa,powhf_y4_ab,powhf_y5_aa,
     &                          powhf_y5_ab,
     &                          powhf_z1_aa,powhf_z1_ab,powhf_z2_aa,
     &                          powhf_z2_ab,powhf_z3_aa,powhf_z3_ab,
     &                          powhf_z4_aa,powhf_z4_ab,powhf_z5_aa,
     &                          powhf_z5_ab)
c
c     == Get the convergence criterion and maximum number of iterations
c
      if (.not.rtdb_cget(wfn1_inst%rtdb,"wfn1:solver",
     &                   len(wfn1_param%solver),
     &                   wfn1_param%solver)) then
         wfn1_param%solver="diagonalize"
      endif
      if (.not.rtdb_get(wfn1_inst%rtdb,"wfn1:conv_tol",MT_DBL,
     &                  1,conv_tol)) conv_tol = 1.0d-6
      if (.not.rtdb_get(wfn1_inst%rtdb,"wfn1:maxit",MT_INT,
     &                  1,maxit)) maxit = 50
      if (.not.rtdb_get(wfn1_inst%rtdb,"wfn1:t_bath",MT_DBL,
     &                  1,wfn1_param%Tbath)) wfn1_param%Tbath = 0.0d0
      if (wfn1_param%Tbath < 0.0d0) then
        call errquit(pname//" Bath temperature must be at least 0.0",
     &               0,INPUT_ERR)
      endif
      if (.not.rtdb_get(wfn1_inst%rtdb,"wfn1:t_initial",MT_DBL,
     &                  1,wfn1_param%Tinitial)) then
        wfn1_param%Tinitial = wfn1_param%Tbath
      endif
      if (wfn1_param%Tinitial < 0.0d0) then
        call errquit(pname//" Initial temperature must be at least 0.0",
     &               0,INPUT_ERR)
      endif
      if (.not.rtdb_get(wfn1_inst%rtdb,"wfn1:damp_n",MT_DBL,
     &                  1,wfn1_param%damp_n)) wfn1_param%damp_n = 1.0d-3
c    &                  1,wfn1_param%damp_n)) wfn1_param%damp_n = 1.0d-2
      if (.not.rtdb_get(wfn1_inst%rtdb,"wfn1:damp_c",MT_DBL,
     &                  1,wfn1_param%damp_c)) wfn1_param%damp_c = 1.0d-3
c    &                  1,wfn1_param%damp_c)) wfn1_param%damp_c = 2.0d-3
      if (.not.rtdb_get(wfn1_inst%rtdb,"wfn1:denom_min",MT_DBL,
     &                  1,wfn1_param%denom_min)) 
     &                                     wfn1_param%denom_min = 1.0d-2
c    &                                     wfn1_param%denom_min = 1.0d-8
c    &                                     wfn1_param%denom_min = 0.01d0
      if (.not.rtdb_get(wfn1_inst%rtdb,"wfn1:pfac",MT_DBL,
     &                  1,wfn1_param%pfac)) wfn1_param%pfac = 256.0d0
      if (.not.rtdb_get(wfn1_inst%rtdb,"wfn1:rlim",MT_DBL,
     &                  1,wfn1_param%rlim)) wfn1_param%rlim = 1.0d-2
      wfn1_param%convergence = conv_tol
      wfn1_param%maxiter     = maxit
c
c     == Get input/output vectors selection ==
c
      tag = "wfn1:input vectors"
      if (.not.rtdb_cget(rtdb, tag, 1, wfn1_param%movecs_in)) then
        tag = "scf:output vectors"
        if (.not.rtdb_cget(rtdb, tag, 1, wfn1_param%movecs_in)) then
          tag = "dft:output vectors"
          if (.not.rtdb_cget(rtdb, tag, 1, wfn1_param%movecs_in)) then
            tag = "mcscf:output vectors"
            if (.not.rtdb_cget(rtdb, tag, 1, wfn1_param%movecs_in)) then
              wfn1_param%movecs_in = 'atomic'
            endif
          endif
        endif
      endif
      tag = "wfn1:output vectors"
      if (.not.rtdb_cget(rtdb, tag, 1, wfn1_param%movecs_out)) then
        if (wfn1_param%movecs_in.eq.'atomic') then
          call util_file_name('movecs',.false.,.false.,
     &                        wfn1_param%movecs_out)
        else
          wfn1_param%movecs_out = wfn1_param%movecs_in
        endif
      endif
c
      write(luout,*)
      write(luout,100)nelec,wfn1_param%nea,wfn1_param%neb,int(charge),
     &                nopen+1
  100 format(10x,'No. of electrons :',2x,i4,/,
     &       10x,' Alpha electrons :',2x,i4,/,
     &       10x,'  Beta electrons :',2x,i4,/,
     &       10x,'Charge           :',2x,i4,/,
     &       10x,'Spin multiplicity:',2x,i4)
      write(luout,*)
c
c     == Deal with 1-electron integrals ==
c
      call wfn1_allocate_mat(wfn1_param,wfn1_inst,wfn1_mat)
c
c     == Initialize the integral code ==
c
      call int_init(wfn1_inst%rtdb,1,wfn1_inst%basis)
      call int_acc_high()
      call schwarz_init(wfn1_inst%geom,wfn1_inst%basis)
c     call int_1e_cache_ga(wfn1_inst%basis,oskel)
      call wfn1_gen_mat(wfn1_inst,wfn1_mat)
c
c     == create scratch global arrays that depend on nbf only ==
c
      if (.not.ga_create(MT_DBL,wfn1_param%nbf,wfn1_param%nbf,"g_bfbf1",
     &                   -1,-1,g_bfbf1)) then
        call errquit(pname//" allocate g_bfbf1 failed",
     &               wfn1_ibd_proc(wfn1_param%nbf**2,nproc),GA_ERR)
      endif
      if (.not.ga_create(MT_DBL,wfn1_param%nbf,wfn1_param%nbf,"g_bfbf2",
     &                   -1,-1,g_bfbf2)) then
        call errquit(pname//" allocate g_bfbf2 failed",
     &               wfn1_ibd_proc(wfn1_param%nbf**2,nproc),GA_ERR)
      endif
      if (.not.ga_create(MT_DBL,wfn1_param%nbf,wfn1_param%nbf,"g_bfbf3",
     &                   -1,-1,g_bfbf3)) then
        call errquit(pname//" allocate g_bfbf3 failed",
     &               wfn1_ibd_proc(wfn1_param%nbf**2,nproc),GA_ERR)
      endif
      if (.not.ga_create(MT_DBL,wfn1_param%nbf,wfn1_param%nbf,"g_bfbf4",
     &                   -1,-1,g_bfbf4)) then
        call errquit(pname//" allocate g_bfbf4 failed",
     &               wfn1_ibd_proc(wfn1_param%nbf**2,nproc),GA_ERR)
      endif
c
c     == Find the linearly independent orbitals ==
c     == and orthonormalizing transformation    ==
c
      call ga_sync() ! Overlap matrix has to be complete
      call wfn1_lindep(wfn1_mat%g_s,wfn1_mat%g_s12,g_bfbf1,
     &                 wfn1_param%nmo,ndep_tol)
      if (ndep_nmo.gt.0) wfn1_param%nmo = min(wfn1_param%nbf,ndep_nmo)
      if (wfn1_param%ncorea.gt.wfn1_param%nea) then
        call errquit(pname//" alpha #core orbitals exceeds #electrons",
     &               0,UERR)
      endif
      if (wfn1_param%ncoreb.gt.wfn1_param%neb) then
        call errquit(pname//" beta #core orbitals exceeds #electrons",
     &               0,UERR)
      endif
      if (wfn1_param%nvirta.gt.wfn1_param%nmo-wfn1_param%nea) then
        call errquit(pname//" alpha too many virtual orbitals",
     &               0,UERR)
      endif
      if (wfn1_param%nvirtb.gt.wfn1_param%nmo-wfn1_param%neb) then
        call errquit(pname//" beta too many virtual orbitals",
     &               0,UERR)
      endif
c
c     == create other scratch global arrays that depend on ==
c     == nbf and nmo                                       ==
c
      if (.not.ga_create(MT_DBL,wfn1_param%nbf,wfn1_param%nmo,"g_bfmo",
     &                   -1,-1,g_bfmo)) then
        call errquit(pname//" allocate g_bfmo failed",
     &               wfn1_ibd_proc(wfn1_param%nbf*wfn1_param%nmo,nproc),
     &               GA_ERR)
      endif
      if (.not.ga_create(MT_DBL,wfn1_param%nmo,wfn1_param%nmo,"g_momo1",
     &                   -1,-1,g_momo1)) then
        call errquit(pname//" allocate g_momo1 failed",
     &               wfn1_ibd_proc(wfn1_param%nmo**2,nproc),GA_ERR)
      endif
      if (.not.ga_create(MT_DBL,wfn1_param%nmo,wfn1_param%nmo,"g_momo2",
     &                   -1,-1,g_momo2)) then
        call errquit(pname//" allocate g_momo2 failed",
     &               wfn1_ibd_proc(wfn1_param%nmo**2,nproc),GA_ERR)
      endif
      if (.not.ga_create(MT_DBL,wfn1_param%nmo,wfn1_param%nmo,"g_momo3",
     &                   -1,-1,g_momo3)) then
        call errquit(pname//" allocate g_momo3 failed",
     &               wfn1_ibd_proc(wfn1_param%nmo**2,nproc),GA_ERR)
      endif
      if (.not.ga_create(MT_DBL,wfn1_param%nmo,wfn1_param%nmo,"g_momo4",
     &                   -1,-1,g_momo4)) then
        call errquit(pname//" allocate g_momo4 failed",
     &               wfn1_ibd_proc(wfn1_param%nmo**2,nproc),GA_ERR)
      endif
c
c     == Load any DFT related variables ==
c
      call wfn1_inpdft(wfn1_param,wfn1_inst)
c     call dft_rdinput(rtdb)
c
c     == Allocate memory for the wavefunction ==
c
      call wfn1_allocate_wfn(wfn1_param,wfn1_inst,wfn1_wave)
      call wfn1_allocate_wfn(wfn1_param,wfn1_inst,wfn1_trial_wv)
      allocate(wfn1_t_wv%dcnta(1:wfn1_param%nea))
      allocate(wfn1_t_wv%dcntb(1:wfn1_param%neb))
      allocate(wfn1_t_wv%icnfa(1:wfn1_param%nea))
      allocate(wfn1_t_wv%icnfb(1:wfn1_param%neb))
c
c     == Allocate memory for the natural orbital energies
c
      allocate(v_dna(1:wfn1_param%nmo), v_dnb(1:wfn1_param%nmo))
      allocate(v_dca(1:wfn1_param%nmo), v_dcb(1:wfn1_param%nmo))
      allocate(v_dnna(1:wfn1_param%nmo),v_dnnb(1:wfn1_param%nmo))
      allocate(v_dcna(1:wfn1_param%nmo),v_dcnb(1:wfn1_param%nmo))
      v_dna = 0.0d0
      v_dnb = 0.0d0
      v_dca = 0.0d0
      v_dcb = 0.0d0
      v_dnna = 0.0d0
      v_dnnb = 0.0d0
      v_dcna = 0.0d0
      v_dcnb = 0.0d0
cDEBUG
      if (wfn1_param%corr_expr.eq."wfn1s") then
        read(70,*)v_dna(1:wfn1_param%nmo)
        read(70,*)v_dnb(1:wfn1_param%nmo)
        close(70)
        call wfn1_set_entropy_eorb_hf(v_dna,v_dnb,wfn1_param%nmo)
        call wfn1_set_entropy_eorb(v_dna,v_dnb,wfn1_param%nmo)
        v_dna = 0.0d0
        v_dnb = 0.0d0
      endif
cDEBUG
c
c     == Allocate memory for the derivatives wrt the wavefunctions ==
c
      call wfn1_allocate_drv(wfn1_param,wfn1_deriv)
      call wfn1_allocate_drv(wfn1_param,wfn1_trial_drv)
      call wfn1_allocate_drv(wfn1_param,wfn1_part_drv)
c
c     == Populate the natural orbitals with the eigenfunctions of ==
c     == the S^{-1/2} matrix                                      == 
c
      call ga_copy_patch('N',g_bfbf1,1,wfn1_param%nbf,
     &     wfn1_param%nbf-wfn1_param%nmo+1,wfn1_param%nbf,
     &     wfn1_wave%g_na,1,wfn1_param%nbf,1,wfn1_param%nmo)
      call ga_copy(wfn1_wave%g_na,wfn1_wave%g_nb)
c
c     == Set up electron configuration ==
c
      call wfn1_configuration(wfn1_param,wfn1_wave)
      call wfn1_configuration(wfn1_param,wfn1_trial_wv)
c
c     == Set up initial wave functions ==
c
cDEBUG
      write(*,*)'HVD: wfn1: movecs_in  = ',wfn1_param%movecs_in
      write(*,*)'HVD: wfn1: movecs_out = ',wfn1_param%movecs_out
cDEBUG
      if (wfn1_param%movecs_in == "atomic") then
c
c       == Set up natural orbitals ==
c
        call ga_zero(g_bfbf1)
        call guess_dens(wfn1_inst%rtdb,wfn1_inst%geom,wfn1_inst%basis,
     &                  g_bfbf1)
c
        allocate(eig(1:wfn1_param%nbf))
        call rhf_dens_to_mo(wfn1_inst%rtdb,wfn1_inst%geom,
     &                      wfn1_inst%basis,
     &                      wfn1_param%nea+wfn1_param%neb,
     &                      wfn1_param%nbf,wfn1_param%nmo,1.0d-7,
     &                      g_bfbf1,wfn1_wave%g_na,.false.,eig)
        call ga_copy(wfn1_wave%g_na,wfn1_wave%g_nb)
c
c       == Set up correlation functions == 
c
        call wfn1_guess_corr(wfn1_param,wfn1_wave,eig,wfn1_param%nmo)
        deallocate(eig)
c
c       call wfn1_rot_natorb(wfn1_param,wfn1_inst,wfn1_mat,
c    &                       wfn1_wave,wfn1_deriv,wfn1_energy,
c    &                       wfn1_trial_wv,wfn1_trial_drv,wfn1_part_drv,
c    &                       g_bfbf1,g_bfbf2,g_bfbf3,g_bfbf4,
c    &                       g_bfmo,g_momo1,g_momo2,g_momo3,g_momo4,
c    &                       v_dna,v_dnb,v_dca,v_dcb,
c    &                       v_dnna,v_dnnb,v_dcna,v_dcnb,
c    &                       oprerr,curit)
c
      else
        call wfn1_guess_movecs(wfn1_param,wfn1_wave)
      endif
cDEBUG
c     call wfn1_print_analysis(wfn1_param,wfn1_inst,wfn1_wave,
c    &                         v_dna,v_dnb,v_dca,v_dcb,g_momo1,g_momo2)
cDEBUG
c
c     == Initial wavefunction complete ==
c
      call ga_sync() ! make sure wave function guesses are completed
      if (wfn1_param%solver.eq."linesearch") then
        call errquit("wfn1_task_energy: no linesearch implemented",
     &               0,UERR)
c       oconverged = wfn1_solve_line(wfn1_param,wfn1_inst,wfn1_mat,
c    &                    wfn1_wave,wfn1_deriv,wfn1_energy,
c    &                    wfn1_trial_wv,wfn1_trial_drv,wfn1_part_drv,
c    &                    g_bfbf1,g_bfbf2,g_bfbf3,g_bfbf4,
c    &                    g_bfmo,g_momo1,g_momo2,g_momo3,g_momo4,
c    &                    v_dna,v_dnb,v_dca,v_dcb,
c    &                    v_dnna,v_dnnb,v_dcna,v_dcnb,
c    &                    oprerr,curit)
      else if (wfn1_param%solver.eq."diagonalize") then
        call errquit("wfn1_task_energy: no diagonalize implemented",
     &               0,UERR)
c       oconverged = wfn1_solve_diag(wfn1_param,wfn1_inst,wfn1_mat,
c    &                    wfn1_wave,wfn1_deriv,wfn1_energy,
c    &                    wfn1_trial_wv,wfn1_trial_drv,wfn1_part_drv,
c    &                    g_bfbf1,g_bfbf2,g_bfbf3,g_bfbf4,
c    &                    g_bfmo,g_momo1,g_momo2,g_momo3,g_momo4,
c    &                    v_dna,v_dnb,v_dca,v_dcb,
c    &                    v_dnna,v_dnnb,v_dcna,v_dcnb,
c    &                    oprerr,curit)
      else if (wfn1_param%solver.eq."monte-carlo") then
        oconverged = wfn1_solve_mc(wfn1_param,wfn1_inst,wfn1_mat,
     &                    wfn1_wave,wfn1_deriv,wfn1_energy,
     &                    wfn1_trial_wv,wfn1_trial_drv,wfn1_part_drv,
     &                    g_bfbf1,g_bfbf2,g_bfbf3,g_bfbf4,
     &                    g_bfmo,g_momo1,g_momo2,g_momo3,g_momo4,
     &                    v_dna,v_dnb,v_dca,v_dcb,
     &                    v_dnna,v_dnnb,v_dcna,v_dcnb,
     &                    oprerr,curit)
      else
        if (iproc.eq.0) then
          write(luout,*)pname//" solver = ",wfn1_param%solver
          write(luout,*)pname//" valid options are:"
          write(luout,*)pname//" - linesearch"
          write(luout,*)pname//" - diagonalize"
          write(luout,*)pname//" - monte-carlo"
        endif
        call errquit(pname//" invalid solver",0)
      endif
      wfn1_success = oconverged
      if (iproc.eq.0) then
        if (curit.eq.0) then
c
c         The calculation started on a converged solution,
c         print at least the current energy and error
c
          if (.not.oprerr) then
            write(luout,'(i10,f20.8,f16.6)')curit,wfn1_energy%E_tot,
     &                                      error%total
          else
            write(luout,'(i10,f20.8,5f16.6)')curit,wfn1_energy%E_tot,
     &                                       error%total,
     &                                       error%errna,error%errnb,
     &                                       error%errca,error%errcb
          endif
        endif
        write(LuOut,*)
        if (oconverged) then
          write(LuOut,*)'Reached Convergence'
        else
          write(LuOut,*)'Failed to Converge'
        endif
        write(LuOut,*)
      endif
c
c     == Store results ==
c
      if (.not.rtdb_put(rtdb,'wfn1:energy',MT_DBL,1,wfn1_energy%E_tot))
     &  call errquit(pname//' failed to store energy in rtdb', 
     &               0, RTDB_ERR)
      if (.not.rtdb_put(rtdb,'wfn1:converged',MT_LOG,1,oconverged))
     &  call errquit(pname//' failed to store converged in rtdb', 
     &               0, RTDB_ERR)
      call wfn1_write_movecs(wfn1_param,wfn1_inst,wfn1_wave,
     &                       v_dna,v_dnb,v_dca,v_dcb)
c
c     == Print results ==
c
      call ga_sync() ! for good measure
      if (iproc.eq.0) then
        call wfn1_print_energy(wfn1_energy)
      endif
      call wfn1_print_analysis(wfn1_param,wfn1_inst,wfn1_wave,
     &                         v_dna,v_dnb,v_dca,v_dcb,g_momo1,g_momo2)
cDEBUG
c     call wfn1_print_orb_energies(wfn1_param,
c    &     v_dna,v_dnb,v_dca,v_dcb,
c    &     v_dnna,v_dnnb,v_dcna,v_dcnb)
cDEBUG
c
c     == Clean up ==
c
cDEBUG
c     wfn1_success = .true.
cDEBUG
      call grid_cleanup(.false.)
      call wfn1_deallocate_drv(wfn1_part_drv)
      call wfn1_deallocate_drv(wfn1_trial_drv)
      call wfn1_deallocate_drv(wfn1_deriv)
      deallocate(v_dcna,v_dcnb)
      deallocate(v_dnna,v_dnnb)
      deallocate(v_dca, v_dcb)
      deallocate(v_dna, v_dnb)
      deallocate(wfn1_t_wv%dcntb)
      deallocate(wfn1_t_wv%dcnta)
      deallocate(wfn1_t_wv%icnfb)
      deallocate(wfn1_t_wv%icnfa)
      call wfn1_deallocate_wfn(wfn1_trial_wv)
      call wfn1_deallocate_wfn(wfn1_wave)
      call wfn1_deallocate_mat(wfn1_mat)
c
      if (.not.ga_destroy(g_momo4)) then
        call errquit(pname//" deallocate g_momo4 failed",0,GA_ERR)
      endif
      if (.not.ga_destroy(g_momo3)) then
        call errquit(pname//" deallocate g_momo3 failed",0,GA_ERR)
      endif
      if (.not.ga_destroy(g_momo2)) then
        call errquit(pname//" deallocate g_momo2 failed",0,GA_ERR)
      endif
      if (.not.ga_destroy(g_momo1)) then
        call errquit(pname//" deallocate g_momo1 failed",0,GA_ERR)
      endif
      if (.not.ga_destroy(g_bfmo)) then
        call errquit(pname//" deallocate g_bfmo failed",0,GA_ERR)
      endif
      if (.not.ga_destroy(g_bfbf4)) then
        call errquit(pname//" deallocate g_bfbf4 failed",0,GA_ERR)
      endif
      if (.not.ga_destroy(g_bfbf3)) then
        call errquit(pname//" deallocate g_bfbf3 failed",0,GA_ERR)
      endif
      if (.not.ga_destroy(g_bfbf2)) then
        call errquit(pname//" deallocate g_bfbf2 failed",0,GA_ERR)
      endif
      if (.not.ga_destroy(g_bfbf1)) then
        call errquit(pname//" deallocate g_bfbf1 failed",0,GA_ERR)
      endif
c
c     call int_1e_uncache_ga()
      call schwarz_tidy
      call int_terminate
      call util_print_pop
c
      if (.not.bas_destroy(wfn1_inst%basis))
     &  call errquit(pname//' Failed to destroy basis',0,BASIS_ERR)
      if (.not.geom_destroy(wfn1_inst%geom))
     &  call errquit(pname//' Failed to destroy geometry',0,GEOM_ERR)
c
      end function wfn1_task_energy
!     end submodule wfn1_for_nw
