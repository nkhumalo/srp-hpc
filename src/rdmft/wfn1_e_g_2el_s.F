c-----------------------------------------------------------------------
c
!> \ingroup wfn1
!> @{
!>
!> \file wfn1_e_g_2el_s.F
!> The correlated 2-electron contributions
!>
!> @}
!>
!> \ingroup wfn1_priv
!> @{
!>
!> \brief Evaluate the correlated 2-electron contributions
!>
!> This subroutine evaluates the 2-electron contributions of the
!> N-representable correlation formalism to the energy, the
!> natural orbital Fock matrices, and the correlation function Fock
!> matrices.
!>
!> Note that the Coulomb terms are added on to the v_dc_a and v_dc_b
!> vectors as those contributions are diagonal in the natural orbital
!> basis (just as the 1-electron terms are). The exchange terms do not
!> have this property and are therefore immediately transformed and
!> added onto the correlation function Fock matrices.
!>
!> Within the Hartree-Fock formalism the Coulomb terms can be
!> constructed from regular 1RDMs. The exchange terms require the
!> explicit evaluation of the 2RDM diagonal. The latter form 
!> \f$O(N^2)\f$ data objects and are therefore distributed over 
!> multiple nodes. 
!>
!> Note that naive evaluation of the 2RDM would have a complexity
!> of \f$O(N^6)\f$ as is evident from the expression
!> \f{eqnarray*}{
!>   \Gamma_{abcd} &=& -\sum_{i,j}N_{ai}N^*_{di}N_{bj}N^*_{cj}d_{ij}
!> \f}
!> where all indices \f$a, b, c, d, i\f$ and \f$j\f$ run over all
!> basis functions. This complexity can be reduced to \f$O(N^5)\f$
!> by breaking the equations down into two steps
!> \f{eqnarray*}{
!>   \gamma_{adj} &=& -\sum_{i}N_{ai}N^*_{di}d_{ij} \\\\
!>   \Gamma_{abcd} &=& \sum_{j}\gamma_{adj}N_{bj}N^*_{cj}
!> \f}
!> Note that this approach is still communication intensive as typically
!> the loops over \f$a\f$ and \f$d\f$ would be partitioned into blocks,
!> and for every block of \f$\gamma_{adj}\f$ all \f$O(N^2)\f$ elements
!> of the diagonal \f$d_{ij}\f$ have to be accessed. Hence formally
!> \f$O(N^4)\f$ data needs to be communicated.
!>
!> Ultimately, the effective complexities can be reduced using, for
!> example, the Cauchy-Schwarz inequality [1-4]. For large molecules
!> only \f$O(N^2)\f$ of the \f$a, b, c, d\f$ indeces, only \f$O(N)\f$
!> of the \f$a, d\f$ indeces would survive. This implies that both the
!> computational as well as the communication complexity are reduced to
!> \f$O(N^3)\f$. At present the Cauchy-Schwarz inequality is not used.
!>
!> On the memory management aspects there are costs associated with
!> allocating memory on the heap as in the Fortran ALLOCATE statement.
!> To avoid these overheads, in particular in the inner loops, Fortran
!> POINTERs are used. This means that the actual memory is associated
!> with buffers allocated outside the loops. The pointers are set 
!> inside the loops to provide convenient array views on the memory. 
!> For clarity pairs of variables are used. The ALLOCATABLE arrays
!> having names BUF_<name> are the buffers that hold the actual memory.
!> Associated with each such array is a POINTER <name> that points to
!> BUF_<name> and provides an appropriate number of dimensions and
!> corresponding index limits.
!>
!> In addition this subroutine inherently deals with \f$\alpha\$ and
!> \f$\beta\f$-spin variables. Variables related to \f$\alpha\f$-spin
!> electrons are identified by names with a trailing "_a" whereas 
!> \f$\beta\f$-electron variables are identified by names with a
!> trailing "_b".
!>
!> Another naming convention is that variables include the remaining
!> indeces. For example \f$\gamma_{adj}\f$ given above would have a
!> name that lists the remaining indeces \f$a, d, j\f$ and the order
!> of these indeces corresponds to the order of the dimensions of the
!> 3-dimensional tensor.
!>
!> This subroutine computes the 2-electron energy term in two 
!> different ways. Once in accordance to the N-representable
!> correlation formalism [5], and once, for comparison, according to
!> the conventional Hartree-Fock method [6-8].
!>
!> The electron correlation is modelled as a screened Coulomb
!> interaction where the screening is determined by the entropy.
!> Effectively the 2-electron occupation numbers are scaled by a 
!> factor \f$(1-S(d_i,d_j))\f$ where \f$S\f$ is normalized such that
!> \f$S(1/2,1/2)=1\f$.
!>
!> ### References ###
!>
!> [1] A.-L. Cauchy, "Sur les formules qui résultent de l'emploi du
!>     signe > ou <, et sur les moyennes entre plusieurs quantités",
!>     In Book: Cours d'analyse de l'École Royale Polytechnique, (1821).
!>     Reprinted in Cambridge Library Collection - Mathematics,
!>     Cambridge University Press (2009), pp. 438-459, DOI:
!>     <a href="https://dx.doi.org/10.1017/CBO9780511693328.017>
!>     10.1017/CBO9780511693328.017</a>.
!>
!> [2] V. Bouniakowsky, "Sur quelques inegalités concernant les
!>     intégrales aux différences finies", Memoires de L'Academie
!>     Imperiale des Sciences de St.-Petersbourg, VII, 1 (1859) Nr. 9,
!>     pp. 1-18, URL:
!>     <a href="https://books.google.com/books?id=wZhFAAAAcAAJ>
!>     books.google.com/books?id=wZhFAAAAcAAJ</a>.
!>
!> [3] H. A. Schwarz, "Über ein Flächen kleinsten Flächeninhalts
!>     betreffendes Problem der Variationsrechnung", Acta Societatis
!>     Scientiarum Fennicae, XV: 318 (1885), URL:
!>     <a href="https://books.google.com/books?id=BXtIAQAAMAAJ">
!>     https://books.google.com/books?id=BXtIAQAAMAAJ</a>.
!>
!> [4] J. M. Steele, "The Cauchy-Schwarz Master Class", The Mathematical
!>     Association of America (2004), ISBN: 
!>     <a href="https://isbnsearch.org/isbn/052154677X">
!>     0-521-54677-X</a>.
!>
!> [5] H. J. J. van Dam, "An N-representable Hartree-Fock approach for
!>     non-idempotent density matrixes", in preparation (08/31/2017).
!>
!> [6] C. C. J. Roothaan, "New developments in molecular orbital
!>     theory", Reviews of Modern Physics <b>23</b> (1951) pp. 69-89,
!>     DOI: <a href="https://dx.doi.org/10.1103/RevModPhys.23.69">
!>     10.1103/RevModPhys.23.69</a>.
!>
!> [7] G. G. Hall, "The Molecular Orbital Theory of Chemical Valency.
!>     VIII. A Method of Calculating Ionization Potentials", 
!>     Proceedings of the Royal Society A <b>205</b> (1951) pp. 541-552,
!>     DOI: <a href="https://dx.doi.org/10.1098/rspa.1951.0048">
!>     10.1098/rspa.1951.0048</a>.
!>
!> [8] J. A. Pople, R. K. Nesbet, "Self-consistent orbitals for
!>     radicals", The Journal of Chemical Physics <b>22</b> (1954)
!>     pp. 571-572,  DOI:
!>     <a href="https://dx.doi.org/10.1063/1.1740120">
!>     10.1063/1.1740120</a>.
!>
!> ### TO DO ###
!>
!> 1. Work out whether a given processor needs to execute anything
!>    for a given index of ld, lc, lb, or la and skip empty loops.
!>    At the moment a simple task counter is checked in the inner most
!>    loop. Skipping loop bodies in outer loops can bring further 
!>    efficiency improvements.
!>
!> 2. It seems like in some places consecutive jj loops could be merged
!>    with the potential or reducing the order of intermediate tensors
!>    by 1. If this is actually possible then this could save memory
!>    as well as reduce memory access both of which would be beneficial.
!>
      subroutine wfn1_e_g_2el_s(wfn1_param,wfn1_inst,wfn1_mat,
     &                          wfn1_wave,wfn1_deriv,wfn1_part_drv,
     &                          v_dc_a,v_dc_b,v_dp_a,v_dp_b,
     &                          g_x_aa,g_x_bb,g_dx_aa,g_dx_bb,
     &                          E_2el,E_2el_chf)
      implicit none
c
!include "wfn1_param.fh"
!include "wfn1_matcnst.fh"
!include "wfn1_wfn.fh"
!include "wfn1_drv.fh"
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "bas.fh"
c
      type(wfn1_prm),     intent(in) :: wfn1_param ! the job parameters
      type(wfn1_prminst), intent(in) :: wfn1_inst  ! the job instances
      type(wfn1_matcnst), intent(in) :: wfn1_mat   ! matrix constants
      type(wfn1_wfn),     intent(in) :: wfn1_wave  ! the wavefunction
c
      type(wfn1_drv),  intent(inout) :: wfn1_deriv ! the derivatives
!>    The non-correlation, non-thermal derivatives
      type(wfn1_drv),  intent(inout) :: wfn1_part_drv
c
c     Correlation function Fock matrices in the natural orbital basis
c     alpha-, and beta-spin respectively (excluding exchange)
c
      double precision, intent(inout) :: v_dc_a(1:wfn1_param%nmo)
      double precision, intent(inout) :: v_dc_b(1:wfn1_param%nmo)
c
c     Correlation function Fock matrices in the natural orbital basis
c     alpha-, and beta-spin respectively (excluding exchange) the
c     non-correlation and non-thermal parts
c
      double precision, intent(inout) :: v_dp_a(1:wfn1_param%nmo)
      double precision, intent(inout) :: v_dp_b(1:wfn1_param%nmo)
c
c     This implementation assumes that g_x_aa, g_x_bb, g_dx_aa,
c     and g_dx_bb all have the same data distribution. If the actual
c     data distributions are different then this will not affect the
c     results as these differences will be resolved by the Global
c     Arrays. However, this may generate additional communication and
c     hence reduce the efficiency.
c
c     Nmo x nmo global arrays for the 2RDM exchange diagonals.
c     See exchange terms of Eq.(34) of [5].
c
      integer, intent(in) :: g_x_aa ! exchange alpha-alpha block
      integer, intent(in) :: g_x_bb ! exchange beta-beta block
c
c     Nmo x nmo global arrays for the 2RDM exchange diagonals
c     derivatives. See the diagonal terms of Eq.(119) of [5].
c
      integer, intent(in) :: g_dx_aa ! exchange alpha-alpha block
      integer, intent(in) :: g_dx_bb ! exchange beta-beta block
c
c     The Hartree-Fock 2-electron energies
c
c     The 2-electron energy [5]:
c
      double precision, intent(out) :: E_2el ! The 2-electron energy [5]
c
c     The conventional Hartree-Fock 2-electron energy [6-8]:
c
      double precision, intent(out) :: E_2el_chf 
c
c     Local variables:
c
      integer :: nproc ! the number of processors
      integer :: iproc ! the rank of this processor
      integer :: jproc ! a counter over processors
      integer :: itask ! the current task number
c
      integer :: nbf   ! the number of basis functions
      integer :: nmo   ! the number of MOs
      integer :: ne_a  ! the number of alpha-electrons
      integer :: ne_b  ! the number of beta-electrons
c
      integer :: ii, jj, kk, ll ! natural orbital labels
      integer :: ia, ib, ic, id ! atomic orbital labels
      integer :: la, lb, lc, ld ! atomic orbital shell labels
      integer :: na, nb, nc, nd ! number of shell functions
c
      integer :: ni, nj ! number of natural orbitals
c
      integer :: ia_lo, ia_hi ! a-index atomic orbital limits
      integer :: ib_lo, ib_hi ! b-index atomic orbital limits
      integer :: ic_lo, ic_hi ! c-index atomic orbital limits
      integer :: id_lo, id_hi ! d-index atomic orbital limits
c
      integer :: ii_lo, ii_hi ! i-index natural orbital limits
      integer :: jj_lo, jj_hi ! j-index natural orbital limits
c
      integer :: if_lo, if_hi ! i-index Fock matrix limits
      integer :: jf_lo, jf_hi ! i-index Fock matrix limits
      integer :: ldf          ! leading dimension on Fock matrix
      integer :: idxf         ! index for Fock matrix in local memory
      integer :: ip_lo, ip_hi ! i-index part Fock matrix limits
      integer :: jp_lo, jp_hi ! i-index part Fock matrix limits
      integer :: ldp          ! leading dimension on part Fock matrix
      integer :: idxp         ! index for part Fock matrix in local
                              ! memory
      integer :: ix_lo, ix_hi ! i-index exchange occ. matrix limits
      integer :: jx_lo, jx_hi ! i-index exchange occ. matrix limits
      integer :: ldx          ! leading dimension on exchange occ matrix
      integer :: idxx         ! index for exchange occ. matrix in memory
c
      integer :: max2e  ! quartet maximum no. 2-el integrals
      integer :: mem2   ! amount of scratch space needed
      integer :: maxsbf ! maximum no. of basis functions in a shell
      integer :: nshell ! the number of shells
c
      logical :: ok_dist ! is the data distribution OK?
c
      double precision, allocatable :: occ_a(:) ! a-occupation numbers
      double precision, allocatable :: occ_b(:) ! b-occupation numbers
c
c     Data structure for the electron repulsion integrals, as well as
c     2RDM blocks for the density matrix of [5] and the conventional
c     Hartree-Fock 2RDM.
c
      double precision, allocatable         :: scr(:)
      double precision, allocatable, target :: buf_eri(:)
      double precision, pointer             :: eri(:,:,:,:)
      double precision, allocatable, target :: buf_rdm2(:)
      double precision, pointer             :: rdm2(:,:,:,:)
      double precision, allocatable, target :: buf_rdm2_chf(:)
      double precision, pointer             :: rdm2_chf(:,:,:,:)
c
c     Intermediates in the 2RDM construction
c
      double precision, allocatable, target :: buf_rdm2_acj_aa(:)
      double precision, pointer             :: rdm2_acj_aa(:,:,:)
      double precision, allocatable, target :: buf_rdm2_acj_ab(:)
      double precision, pointer             :: rdm2_acj_ab(:,:,:)
      double precision, allocatable, target :: buf_rdm2_acj_ba(:)
      double precision, pointer             :: rdm2_acj_ba(:,:,:)
      double precision, allocatable, target :: buf_rdm2_acj_bb(:)
      double precision, pointer             :: rdm2_acj_bb(:,:,:)
      double precision, allocatable, target :: buf_rdm2_adj_aa(:)
      double precision, pointer             :: rdm2_adj_aa(:,:,:)
      double precision, allocatable, target :: buf_rdm2_adj_ab(:)
      double precision, pointer             :: rdm2_adj_ab(:,:,:)
      double precision, allocatable, target :: buf_rdm2_adj_ba(:)
      double precision, pointer             :: rdm2_adj_ba(:,:,:)
      double precision, allocatable, target :: buf_rdm2_adj_bb(:)
      double precision, pointer             :: rdm2_adj_bb(:,:,:)
c
c     Same intermediates in the 2RDM construction with entropy
c
      double precision, allocatable, target :: buf_rdm2s_acj_aa(:)
      double precision, pointer             :: rdm2s_acj_aa(:,:,:)
      double precision, allocatable, target :: buf_rdm2s_acj_ab(:)
      double precision, pointer             :: rdm2s_acj_ab(:,:,:)
      double precision, allocatable, target :: buf_rdm2s_acj_ba(:)
      double precision, pointer             :: rdm2s_acj_ba(:,:,:)
      double precision, allocatable, target :: buf_rdm2s_acj_bb(:)
      double precision, pointer             :: rdm2s_acj_bb(:,:,:)
      double precision, allocatable, target :: buf_rdm2s_adj_aa(:)
      double precision, pointer             :: rdm2s_adj_aa(:,:,:)
      double precision, allocatable, target :: buf_rdm2s_adj_ab(:)
      double precision, pointer             :: rdm2s_adj_ab(:,:,:)
      double precision, allocatable, target :: buf_rdm2s_adj_ba(:)
      double precision, pointer             :: rdm2s_adj_ba(:,:,:)
      double precision, allocatable, target :: buf_rdm2s_adj_bb(:)
      double precision, pointer             :: rdm2s_adj_bb(:,:,:)
c
c     Intermediates in the derivative of the 2RDM construction
c
      double precision, allocatable, target :: buf_rdm2_d_acj_a(:)
      double precision, pointer             :: rdm2_d_acj_a(:,:,:)
      double precision, allocatable, target :: buf_rdm2_d_acj_b(:)
      double precision, pointer             :: rdm2_d_acj_b(:,:,:)
      double precision, allocatable, target :: buf_rdm2_d_adj_a(:)
      double precision, pointer             :: rdm2_d_adj_a(:,:,:)
      double precision, allocatable, target :: buf_rdm2_d_adj_b(:)
      double precision, pointer             :: rdm2_d_adj_b(:,:,:)
c
c     Data structures for the 1RDM blocks for the Coulomb terms of 
c     the 2RDM of [5] and the conventional Hartree-Fock 2RDM.
c
c     1RDM factors for Coulomb terms
c
      double precision, allocatable, target :: buf_rdm1_ac_a(:)
      double precision, pointer             :: rdm1_ac_a(:,:)
      double precision, allocatable, target :: buf_rdm1_ac_b(:)
      double precision, pointer             :: rdm1_ac_b(:,:)
      double precision, allocatable, target :: buf_rdm1_bd_a(:)
      double precision, pointer             :: rdm1_bd_a(:,:)
      double precision, allocatable, target :: buf_rdm1_bd_b(:)
      double precision, pointer             :: rdm1_bd_b(:,:)
c
c     1RDM factors for exchange terms of the conventional Hartree-Fock
c
      double precision, allocatable, target :: buf_rdm1_ad_a(:)
      double precision, pointer             :: rdm1_ad_a(:,:)
      double precision, allocatable, target :: buf_rdm1_ad_b(:)
      double precision, pointer             :: rdm1_ad_b(:,:)
      double precision, allocatable, target :: buf_rdm1_bc_a(:)
      double precision, pointer             :: rdm1_bc_a(:,:)
      double precision, allocatable, target :: buf_rdm1_bc_b(:)
      double precision, pointer             :: rdm1_bc_b(:,:)
c
c     Fock matrix blocks. 
c
c     For the Fock matrices for the natural orbitals
c
      double precision, allocatable, target :: buf_fn_ac_a(:)
      double precision, pointer             :: fn_ac_a(:,:)
      double precision, allocatable, target :: buf_fn_ac_b(:)
      double precision, pointer             :: fn_ac_b(:,:)
      double precision, allocatable, target :: buf_fn_ad_a(:)
      double precision, pointer             :: fn_ad_a(:,:)
      double precision, allocatable, target :: buf_fn_ad_b(:)
      double precision, pointer             :: fn_ad_b(:,:)
      double precision, allocatable, target :: buf_fn_bcj_a(:)
      double precision, pointer             :: fn_bcj_a(:,:,:)
      double precision, allocatable, target :: buf_fn_bcj_b(:)
      double precision, pointer             :: fn_bcj_b(:,:,:)
      double precision, allocatable, target :: buf_fn_bdj_a(:)
      double precision, pointer             :: fn_bdj_a(:,:,:)
      double precision, allocatable, target :: buf_fn_bdj_b(:)
      double precision, pointer             :: fn_bdj_b(:,:,:)
c
c     Same as above but including entropy
c
      double precision, allocatable, target :: buf_fns_ac_a(:)
      double precision, pointer             :: fns_ac_a(:,:)
      double precision, allocatable, target :: buf_fns_ac_b(:)
      double precision, pointer             :: fns_ac_b(:,:)
      double precision, allocatable, target :: buf_fns_ad_a(:)
      double precision, pointer             :: fns_ad_a(:,:)
      double precision, allocatable, target :: buf_fns_ad_b(:)
      double precision, pointer             :: fns_ad_b(:,:)
      double precision, allocatable, target :: buf_fns_bcj_a(:)
      double precision, pointer             :: fns_bcj_a(:,:,:)
      double precision, allocatable, target :: buf_fns_bcj_b(:)
      double precision, pointer             :: fns_bcj_b(:,:,:)
      double precision, allocatable, target :: buf_fns_bdj_a(:)
      double precision, pointer             :: fns_bdj_a(:,:,:)
      double precision, allocatable, target :: buf_fns_bdj_b(:)
      double precision, pointer             :: fns_bdj_b(:,:,:)
c
c     - Same as above but partially transformed
c
      double precision, allocatable, target :: buf_fn_ic_a(:)
      double precision, pointer             :: fn_ic_a(:,:)
      double precision, allocatable, target :: buf_fn_ic_b(:)
      double precision, pointer             :: fn_ic_b(:,:)
      double precision, allocatable, target :: buf_fn_id_a(:)
      double precision, pointer             :: fn_id_a(:,:)
      double precision, allocatable, target :: buf_fn_id_b(:)
      double precision, pointer             :: fn_id_b(:,:)
      double precision, allocatable, target :: buf_fn_bj_a(:)
      double precision, pointer             :: fn_bj_a(:,:)
      double precision, allocatable, target :: buf_fn_bj_b(:)
      double precision, pointer             :: fn_bj_b(:,:)
c
c     - With entropy
c
      double precision, allocatable, target :: buf_fns_ic_a(:)
      double precision, pointer             :: fns_ic_a(:,:)
      double precision, allocatable, target :: buf_fns_ic_b(:)
      double precision, pointer             :: fns_ic_b(:,:)
      double precision, allocatable, target :: buf_fns_id_a(:)
      double precision, pointer             :: fns_id_a(:,:)
      double precision, allocatable, target :: buf_fns_id_b(:)
      double precision, pointer             :: fns_id_b(:,:)
      double precision, allocatable, target :: buf_fns_bj_a(:)
      double precision, pointer             :: fns_bj_a(:,:)
      double precision, allocatable, target :: buf_fns_bj_b(:)
      double precision, pointer             :: fns_bj_b(:,:)
c
c     - Same as above but fully transformed
c
      double precision, allocatable, target :: buf_fn_ij_a(:)
      double precision, pointer             :: fn_ij_a(:,:)
      double precision, allocatable, target :: buf_fn_ij_b(:)
      double precision, pointer             :: fn_ij_b(:,:)
c
c     - With entropy
c
      double precision, allocatable, target :: buf_fns_ij_a(:)
      double precision, pointer             :: fns_ij_a(:,:)
      double precision, allocatable, target :: buf_fns_ij_b(:)
      double precision, pointer             :: fns_ij_b(:,:)
c
c     For the Fock matrices for the correlation functions
c
      double precision, allocatable, target :: buf_fc_acj_a(:)
      double precision, pointer             :: fc_acj_a(:,:,:)
      double precision, allocatable, target :: buf_fc_acj_b(:)
      double precision, pointer             :: fc_acj_b(:,:,:)
      double precision, allocatable, target :: buf_fc_adj_a(:)
      double precision, pointer             :: fc_adj_a(:,:,:)
      double precision, allocatable, target :: buf_fc_adj_b(:)
      double precision, pointer             :: fc_adj_b(:,:,:)
c
      complex(kind=8),  allocatable, target :: buf_fc_ij_aa(:)
      complex(kind=8),  pointer             :: fc_ij_aa(:,:)
      complex(kind=8),  allocatable, target :: buf_fc_ij_ab(:)
      complex(kind=8),  pointer             :: fc_ij_ab(:,:)
      complex(kind=8),  allocatable, target :: buf_fc_ij_ba(:)
      complex(kind=8),  pointer             :: fc_ij_ba(:,:)
      complex(kind=8),  allocatable, target :: buf_fc_ij_bb(:)
      complex(kind=8),  pointer             :: fc_ij_bb(:,:)
c     double precision, allocatable, target :: buf_fc_bdj_a(:)
c     double precision, pointer             :: fc_bdj_a(:,:,:)
c     double precision, allocatable, target :: buf_fc_bdj_b(:)
c     double precision, pointer             :: fc_bdj_b(:,:,:)
c
c     Orbital coefficients
c
      double precision, allocatable, target :: buf_N_di_a(:)
      double precision, pointer             :: N_di_a(:,:)
      double precision, allocatable, target :: buf_N_di_b(:)
      double precision, pointer             :: N_di_b(:,:)
      double precision, allocatable, target :: buf_N_ci_a(:)
      double precision, pointer             :: N_ci_a(:,:)
      double precision, allocatable, target :: buf_N_ci_b(:)
      double precision, pointer             :: N_ci_b(:,:)
      double precision, allocatable, target :: buf_N_bi_a(:)
      double precision, pointer             :: N_bi_a(:,:)
      double precision, allocatable, target :: buf_N_bi_b(:)
      double precision, pointer             :: N_bi_b(:,:)
      double precision, allocatable, target :: buf_N_ai_a(:)
      double precision, pointer             :: N_ai_a(:,:)
      double precision, allocatable, target :: buf_N_ai_b(:)
      double precision, pointer             :: N_ai_b(:,:)
c
c     Correlation function coefficients
c
      complex(kind=8),  allocatable, target :: buf_C_is_a(:)
      complex(kind=8),  pointer             :: C_is_a(:,:)
      complex(kind=8),  allocatable, target :: buf_C_is_b(:)
      complex(kind=8),  pointer             :: C_is_b(:,:)
      complex(kind=8),  allocatable, target :: buf_C_jt_a(:)
      complex(kind=8),  pointer             :: C_jt_a(:,:)
      complex(kind=8),  allocatable, target :: buf_C_jt_b(:)
      complex(kind=8),  pointer             :: C_jt_b(:,:)
c
c     Occupation numbers of the 2RDM exchange terms
c
cDEBUG
c     double precision, allocatable, target :: buf_d_aa(:)
c     double precision, pointer             :: d_aa(:,:)
cDEBUG
      double precision, allocatable, target :: buf_x_aa(:)
      double precision, pointer             :: x_aa(:,:)
      double precision, allocatable, target :: buf_x_bb(:)
      double precision, pointer             :: x_bb(:,:)
      double precision, allocatable, target :: buf_dx_aa(:)
      double precision, pointer             :: dx_aa(:,:)
      double precision, allocatable, target :: buf_dx_bb(:)
      double precision, pointer             :: dx_bb(:,:)
#ifdef TEST_RDMS
      double precision, allocatable :: da(:,:)
      double precision, allocatable :: db(:,:)
      double precision, allocatable :: gaa(:,:,:,:)
      double precision, allocatable :: gab(:,:,:,:)
      double precision, allocatable :: gbb(:,:,:,:)
#endif
c
c     Entropy related quantities
c
      double precision :: saaij     ! The entropy for pair i,j
      double precision :: dsaaij_i  ! The entropy derivative wrt d_i
      double precision :: dsaaij_j  ! The entropy derivative wrt d_j
      double precision :: dsaaij_ij ! The entropy derivative wrt d_ij
      double precision :: sabij     ! The entropy for pair i,j
      double precision :: dsabij_i  ! The entropy derivative wrt d_i
      double precision :: dsabij_j  ! The entropy derivative wrt d_j
      double precision :: sbaij     ! The entropy for pair i,j
      double precision :: dsbaij_i  ! The entropy derivative wrt d_i
      double precision :: dsbaij_j  ! The entropy derivative wrt d_j
      double precision :: sbbij     ! The entropy for pair i,j
      double precision :: dsbbij_i  ! The entropy derivative wrt d_i
      double precision :: dsbbij_j  ! The entropy derivative wrt d_j
      double precision :: dsbbij_ij ! The entropy derivative wrt d_ij
c
      double precision, target :: dummy(2)
cDEBUG
      double precision :: ac_aa_mx = 0.0d0
      double precision :: ac_ab_mx = 0.0d0
      double precision :: ac_ba_mx = 0.0d0
      double precision :: ac_bb_mx = 0.0d0
      double precision :: ad_aa_mx = 0.0d0
      double precision :: ad_ab_mx = 0.0d0
      double precision :: ad_ba_mx = 0.0d0
      double precision :: ad_bb_mx = 0.0d0
cDEBUG
c
c     Functions:
c
      double precision, external :: ddot
c
c     Code:
c
      nproc = ga_nnodes()
      iproc = ga_nodeid()
      itask = -1
c
      E_2el     = 0.0d0
      E_2el_chf = 0.0d0
c
      nbf  = wfn1_param%nbf
      nmo  = wfn1_param%nmo
      ne_a = wfn1_param%nea
      ne_b = wfn1_param%neb
c
      if (.not. bas_numcont(wfn1_inst%basis, nshell)) call errquit
     &     ('wfn1_e_g_2el_s: bas_numcont failed for basis',
     &      wfn1_inst%basis, BASIS_ERR)
      if (.not. bas_nbf_cn_max(wfn1_inst%basis, maxsbf)) call errquit
     &     ('wfn1_e_g_2el_s: bas_nbf_cn_max failed for basis',
     &      wfn1_inst%basis, BASIS_ERR)
c
      call ga_distribution(g_x_aa,0,ii_lo,ii_hi,jj_lo,jj_hi)
      ni = ii_hi-ii_lo+1
      nj = jj_hi-jj_lo+1
      call int_mem_2e4c(max2e,mem2)
c
#ifdef TEST_RDMS
      allocate(da(nbf,nbf),db(nbf,nbf))
      allocate(gaa(nbf,nbf,nbf,nbf))
      allocate(gab(nbf,nbf,nbf,nbf))
      allocate(gbb(nbf,nbf,nbf,nbf))
      da  = 0.0d0
      db  = 0.0d0
      gaa = 0.0d0
      gab = 0.0d0
      gbb = 0.0d0
#endif
c
      allocate(occ_a(1:nmo),occ_b(1:nmo))
      allocate(scr(mem2))
      allocate(buf_eri(max2e),buf_rdm2(max2e),buf_rdm2_chf(max2e))
c
      allocate(buf_rdm2_acj_aa(maxsbf*maxsbf*nj))
      allocate(buf_rdm2_acj_ab(maxsbf*maxsbf*nj))
      allocate(buf_rdm2_acj_ba(maxsbf*maxsbf*nj))
      allocate(buf_rdm2_acj_bb(maxsbf*maxsbf*nj))
      allocate(buf_rdm2_adj_aa(maxsbf*maxsbf*nj))
      allocate(buf_rdm2_adj_ab(maxsbf*maxsbf*nj))
      allocate(buf_rdm2_adj_ba(maxsbf*maxsbf*nj))
      allocate(buf_rdm2_adj_bb(maxsbf*maxsbf*nj))
c
      allocate(buf_rdm2s_acj_aa(maxsbf*maxsbf*nj))
      allocate(buf_rdm2s_acj_ab(maxsbf*maxsbf*nj))
      allocate(buf_rdm2s_acj_ba(maxsbf*maxsbf*nj))
      allocate(buf_rdm2s_acj_bb(maxsbf*maxsbf*nj))
      allocate(buf_rdm2s_adj_aa(maxsbf*maxsbf*nj))
      allocate(buf_rdm2s_adj_ab(maxsbf*maxsbf*nj))
      allocate(buf_rdm2s_adj_ba(maxsbf*maxsbf*nj))
      allocate(buf_rdm2s_adj_bb(maxsbf*maxsbf*nj))
c
      allocate(buf_rdm2_d_acj_a(maxsbf*maxsbf*nj))
      allocate(buf_rdm2_d_acj_b(maxsbf*maxsbf*nj))
      allocate(buf_rdm2_d_adj_a(maxsbf*maxsbf*nj))
      allocate(buf_rdm2_d_adj_b(maxsbf*maxsbf*nj))
c
      allocate(buf_N_ai_a(maxsbf*nmo),buf_N_ai_b(maxsbf*nmo))
      allocate(buf_N_bi_a(maxsbf*nmo),buf_N_bi_b(maxsbf*nmo))
      allocate(buf_N_ci_a(maxsbf*nmo),buf_N_ci_b(maxsbf*nmo))
      allocate(buf_N_di_a(maxsbf*nmo),buf_N_di_b(maxsbf*nmo))
c
      allocate(buf_C_is_a(ni*nmo),buf_C_is_b(ni*nmo))
      allocate(buf_C_jt_a(nj*nmo),buf_C_jt_b(nj*nmo))
c
      allocate(buf_rdm1_ac_a(maxsbf*maxsbf))
      allocate(buf_rdm1_ac_b(maxsbf*maxsbf))
      allocate(buf_rdm1_bd_a(maxsbf*maxsbf))
      allocate(buf_rdm1_bd_b(maxsbf*maxsbf))
c
      allocate(buf_rdm1_ad_a(maxsbf*maxsbf))
      allocate(buf_rdm1_ad_b(maxsbf*maxsbf))
      allocate(buf_rdm1_bc_a(maxsbf*maxsbf))
      allocate(buf_rdm1_bc_b(maxsbf*maxsbf))
c
cDEBUG
c     allocate(buf_d_aa(ni*nj))
cDEBUG
      allocate(buf_x_aa(ni*nj), buf_x_bb(ni*nj))
      allocate(buf_dx_aa(ni*nj),buf_dx_bb(ni*nj))
c
      allocate(buf_fn_ac_a(maxsbf*maxsbf),buf_fn_ac_b(maxsbf*maxsbf))
      allocate(buf_fn_ad_a(maxsbf*maxsbf),buf_fn_ad_b(maxsbf*maxsbf))
      allocate(buf_fns_ac_a(maxsbf*maxsbf),buf_fns_ac_b(maxsbf*maxsbf))
      allocate(buf_fns_ad_a(maxsbf*maxsbf),buf_fns_ad_b(maxsbf*maxsbf))
c
      allocate(buf_fn_bcj_a(maxsbf*maxsbf*nj))
      allocate(buf_fn_bcj_b(maxsbf*maxsbf*nj))
      allocate(buf_fn_bdj_a(maxsbf*maxsbf*nj))
      allocate(buf_fn_bdj_b(maxsbf*maxsbf*nj))
      allocate(buf_fns_bcj_a(maxsbf*maxsbf*nj))
      allocate(buf_fns_bcj_b(maxsbf*maxsbf*nj))
      allocate(buf_fns_bdj_a(maxsbf*maxsbf*nj))
      allocate(buf_fns_bdj_b(maxsbf*maxsbf*nj))
c
      allocate(buf_fn_ic_a(nmo*maxsbf),buf_fn_ic_b(nmo*maxsbf))
      allocate(buf_fn_id_a(nmo*maxsbf),buf_fn_id_b(nmo*maxsbf))
      allocate(buf_fn_bj_a(maxsbf*nj),buf_fn_bj_b(maxsbf*nj))
      allocate(buf_fn_ij_a(nmo*nj),buf_fn_ij_b(nmo*nj))
      allocate(buf_fns_ic_a(nmo*maxsbf),buf_fns_ic_b(nmo*maxsbf))
      allocate(buf_fns_id_a(nmo*maxsbf),buf_fns_id_b(nmo*maxsbf))
      allocate(buf_fns_bj_a(maxsbf*nj),buf_fns_bj_b(maxsbf*nj))
      allocate(buf_fns_ij_a(nmo*nj),buf_fns_ij_b(nmo*nj))
c
      allocate(buf_fc_acj_a(maxsbf*maxsbf*nj))
      allocate(buf_fc_acj_b(maxsbf*maxsbf*nj))
      allocate(buf_fc_adj_a(maxsbf*maxsbf*nj))
      allocate(buf_fc_adj_b(maxsbf*maxsbf*nj))
      allocate(buf_fc_ij_aa(ni*nj))
      allocate(buf_fc_ij_ab(ni*nj))
      allocate(buf_fc_ij_ba(ni*nj))
      allocate(buf_fc_ij_bb(ni*nj))
c     allocate(buf_fc_bdj_a(maxsbf*maxsbf*nj))
c     allocate(buf_fc_bdj_b(maxsbf*maxsbf*nj))
c
c     Calculate occupation numbers
c     ----------------------------
c
c     Calculate alpha-occupation numbers
c
      call wfn1_calc_occ(wfn1_wave%g_ca,nmo,ne_a,wfn1_wave%icnfa,
     &                   wfn1_wave%dcnta,occ_a)
c
c     Calculate beta-occupation numbers
c
      call wfn1_calc_occ(wfn1_wave%g_cb,nmo,ne_b,wfn1_wave%icnfb,
     &                   wfn1_wave%dcntb,occ_b)
c
c     Calculate exchange alpha-alpha-occupation numbers
c
      call wfn1_calc_ocx(wfn1_wave%g_ca,nmo,ne_a,wfn1_wave%icnfa,
     &                   wfn1_wave%dcnta,g_x_aa)
c
c     Calculate exchange beta-beta-occupation numbers
c
      call wfn1_calc_ocx(wfn1_wave%g_cb,nmo,ne_b,wfn1_wave%icnfb,
     &                   wfn1_wave%dcntb,g_x_bb)
c
c     Calculate exchange alpha-alpha-occupation number derivative
c
      call wfn1_calc_docx(wfn1_wave%g_ca,nmo,ne_a,wfn1_wave%icnfa,
     &                    wfn1_wave%dcnta,g_dx_aa)
c
c     Calculate exchange beta-beta-occupation number derivative
c
      call wfn1_calc_docx(wfn1_wave%g_cb,nmo,ne_b,wfn1_wave%icnfb,
     &                    wfn1_wave%dcntb,g_dx_bb)
cDEBUG
c     write(*,*)"HVD: occ_a: ",occ_a
c     write(*,*)"HVD: occ_b: ",occ_b
c     write(*,*)"HVD: ocx:"
c     call util_flush(6)
c     call ga_print(g_x_aa)
c     call ga_print(g_x_bb)
c     call util_flush(6)
cDEBUG

c
      call ga_sync()
c
      do ld = 1, nshell
        if (.not. bas_cn2bfr(wfn1_inst%basis, ld, id_lo, id_hi))
     &      call errquit('wfn1_e_g_2el_s: bas_cn2bfr',
     &                   wfn1_inst%basis, BASIS_ERR)
c
        nd = id_hi-id_lo+1
        N_di_a(id_lo:id_hi,1:nmo) => buf_N_di_a
        N_di_b(id_lo:id_hi,1:nmo) => buf_N_di_b
        call ga_get(wfn1_wave%g_na,id_lo,id_hi,1,nmo,N_di_a,nd)
        call ga_get(wfn1_wave%g_nb,id_lo,id_hi,1,nmo,N_di_b,nd)
c
        do lb = 1, nshell
          if (.not. bas_cn2bfr(wfn1_inst%basis, lb, ib_lo, ib_hi))
     &        call errquit('wfn1_e_g_2el_s: bas_cn2bfr',
     &                     wfn1_inst%basis, BASIS_ERR)
c
          nb = ib_hi-ib_lo+1
          N_bi_a(ib_lo:ib_hi,1:nmo) => buf_N_bi_a
          N_bi_b(ib_lo:ib_hi,1:nmo) => buf_N_bi_b
          call ga_get(wfn1_wave%g_na,ib_lo,ib_hi,1,nmo,N_bi_a,nb)
          call ga_get(wfn1_wave%g_nb,ib_lo,ib_hi,1,nmo,N_bi_b,nb)
c
c         Construct the bd 1RDM blocks
c
          rdm1_bd_a(ib_lo:ib_hi,id_lo:id_hi) => buf_rdm1_bd_a
          rdm1_bd_b(ib_lo:ib_hi,id_lo:id_hi) => buf_rdm1_bd_b
          rdm1_bd_a = 0.0d0
          rdm1_bd_b = 0.0d0
          do id = id_lo, id_hi
            do ib = ib_lo, ib_hi
              do ii = 1, nmo
                rdm1_bd_a(ib,id) = rdm1_bd_a(ib,id)
     &                           + N_bi_a(ib,ii)*N_di_a(id,ii)*occ_a(ii)
                rdm1_bd_b(ib,id) = rdm1_bd_b(ib,id)
     &                           + N_bi_b(ib,ii)*N_di_b(id,ii)*occ_b(ii)
              enddo ! ii
            enddo ! ib
          enddo ! id
#ifdef TEST_RDMS
          do id = id_lo, id_hi
            do ib = ib_lo, ib_hi
              da(ib,id) = rdm1_bd_a(ib,id)
              db(ib,id) = rdm1_bd_b(ib,id)
            enddo ! ib
          enddo ! id
#endif
c
          do lc = 1, nshell 
            if (.not. bas_cn2bfr(wfn1_inst%basis, lc, ic_lo, ic_hi))
     &          call errquit('wfn1_e_g_2el_s: bas_cn2bfr',
     &                       wfn1_inst%basis, BASIS_ERR)
c
            nc = ic_hi-ic_lo+1
            N_ci_a(ic_lo:ic_hi,1:nmo) => buf_N_ci_a
            N_ci_b(ic_lo:ic_hi,1:nmo) => buf_N_ci_b
            call ga_get(wfn1_wave%g_na,ic_lo,ic_hi,1,nmo,N_ci_a,nc)
            call ga_get(wfn1_wave%g_nb,ic_lo,ic_hi,1,nmo,N_ci_b,nc)
c
c           Construct the bc 1RDM blocks
c
            rdm1_bc_a(ib_lo:ib_hi,ic_lo:ic_hi) => buf_rdm1_bc_a
            rdm1_bc_b(ib_lo:ib_hi,ic_lo:ic_hi) => buf_rdm1_bc_b
            rdm1_bc_a = 0.0d0
            rdm1_bc_b = 0.0d0
            do ic = ic_lo, ic_hi
              do ib = ib_lo, ib_hi
                do ii = 1, nmo
                  rdm1_bc_a(ib,ic) = rdm1_bc_a(ib,ic)
     &                             + N_bi_a(ib,ii)*N_ci_a(ic,ii)*
     &                               occ_a(ii)
                  rdm1_bc_b(ib,ic) = rdm1_bc_b(ib,ic)
     &                             + N_bi_b(ib,ii)*N_ci_b(ic,ii)*
     &                               occ_b(ii)
                enddo ! ii
              enddo ! ib
            enddo ! ic
c
            do la = 1, nshell
              itask = itask + 1
              if (mod(itask,nproc).ne.iproc) cycle
              if (.not. bas_cn2bfr(wfn1_inst%basis, la, ia_lo, ia_hi))
     &            call errquit('wfn1_e_g_2el_s: bas_cn2bfr',
     &                         wfn1_inst%basis, BASIS_ERR)
c
              na = ia_hi-ia_lo+1
              N_ai_a(ia_lo:ia_hi,1:nmo) => buf_N_ai_a
              N_ai_b(ia_lo:ia_hi,1:nmo) => buf_N_ai_b
              call ga_get(wfn1_wave%g_na,ia_lo,ia_hi,1,nmo,N_ai_a,na)
              call ga_get(wfn1_wave%g_nb,ia_lo,ia_hi,1,nmo,N_ai_b,na)
c
c             Construct the ac 1RDM blocks
c
              rdm1_ac_a(ia_lo:ia_hi,ic_lo:ic_hi) => buf_rdm1_ac_a
              rdm1_ac_b(ia_lo:ia_hi,ic_lo:ic_hi) => buf_rdm1_ac_b
              rdm1_ac_a = 0.0d0
              rdm1_ac_b = 0.0d0
              do ic = ic_lo, ic_hi
                do ia = ia_lo, ia_hi
                  do ii = 1, nmo
                    rdm1_ac_a(ia,ic) = rdm1_ac_a(ia,ic)
     &                               + N_ai_a(ia,ii)*N_ci_a(ic,ii)*
     &                                 occ_a(ii)
                    rdm1_ac_b(ia,ic) = rdm1_ac_b(ia,ic)
     &                               + N_ai_b(ia,ii)*N_ci_b(ic,ii)*
     &                                 occ_b(ii)
                  enddo ! ii
                enddo ! ia
              enddo ! ic
c
c             Construct the ad 1RDM blocks (conventional Hartree-
c             Fock exchange)
c
              rdm1_ad_a(ia_lo:ia_hi,id_lo:id_hi) => buf_rdm1_ad_a
              rdm1_ad_b(ia_lo:ia_hi,id_lo:id_hi) => buf_rdm1_ad_b
              rdm1_ad_a = 0.0d0
              rdm1_ad_b = 0.0d0
              do id = id_lo, id_hi
                do ia = ia_lo, ia_hi
                  do ii = 1, nmo
                    rdm1_ad_a(ia,id) = rdm1_ad_a(ia,id)
     &                               + N_ai_a(ia,ii)*N_di_a(id,ii)*
     &                                 occ_a(ii)
                    rdm1_ad_b(ia,id) = rdm1_ad_b(ia,id)
     &                               + N_ai_b(ia,ii)*N_di_b(id,ii)*
     &                                 occ_b(ii)
                  enddo ! ii
                enddo ! ia
              enddo ! id
c
c             Set up eri, and rdm2 buffers
c
              eri(ia_lo:ia_hi,ic_lo:ic_hi,ib_lo:ib_hi,id_lo:id_hi)
     &        => buf_eri
              rdm2(ia_lo:ia_hi,ic_lo:ic_hi,ib_lo:ib_hi,id_lo:id_hi)
     &        => buf_rdm2
              rdm2_chf(ia_lo:ia_hi,ic_lo:ic_hi,ib_lo:ib_hi,id_lo:id_hi)
     &        => buf_rdm2_chf
c
c             Calculate the 2-electron integrals
c
              call int_2e4c(wfn1_inst%basis, ld, lb,
     &                      wfn1_inst%basis, lc, la,
     &                      mem2, scr, max2e, eri)
c
              rdm2     = 0.0d0
              rdm2_chf = 0.0d0
c
c             Add the 2RDM terms that can be constructed from
c             1RDM factors.
c
              do id = id_lo, id_hi
                do ib = ib_lo, ib_hi
                  do ic = ic_lo, ic_hi
                    do ia = ia_lo, ia_hi
c                     We cannot construct the screened density matrix
c                     this way
c                     rdm2(ia,ic,ib,id) = rdm2(ia,ic,ib,id)
c    &                  + rdm1_ac_a(ia,ic)*rdm1_bd_a(ib,id)
c    &                  - rdm1_ad_a(ia,id)*rdm1_bc_a(ib,ic)
c    &                  + rdm1_ac_b(ia,ic)*rdm1_bd_b(ib,id)
c    &                  - rdm1_ad_b(ia,id)*rdm1_bc_b(ib,ic)
c    &                  + rdm1_ac_a(ia,ic)*rdm1_bd_b(ib,id)
c    &                  + rdm1_ac_b(ia,ic)*rdm1_bd_a(ib,id)
#ifdef TEST_RDMS
c                     This is OK for N-representability testing
                      gaa(ia,ic,ib,id) = gaa(ia,ic,ib,id)
     &                  + rdm1_ac_a(ia,ic)*rdm1_bd_a(ib,id)
     &                  - rdm1_ad_a(ia,id)*rdm1_bc_a(ib,ic)
                      gab(ia,ic,ib,id) = gab(ia,ic,ib,id)
     &                  + rdm1_ac_a(ia,ic)*rdm1_bd_b(ib,id)
                      gbb(ia,ic,ib,id) = gbb(ia,ic,ib,id)
     &                  + rdm1_ac_b(ia,ic)*rdm1_bd_b(ib,id)
     &                  - rdm1_ad_b(ia,id)*rdm1_bc_b(ib,ic)
#endif
c                     This is OK for the conventional Hartree-Fock
c                     method
                      rdm2_chf(ia,ic,ib,id) = rdm2_chf(ia,ic,ib,id)
     &                  + rdm1_ac_a(ia,ic)*rdm1_bd_a(ib,id)
     &                  + rdm1_ac_b(ia,ic)*rdm1_bd_b(ib,id)
     &                  + rdm1_ac_a(ia,ic)*rdm1_bd_b(ib,id)
     &                  + rdm1_ac_b(ia,ic)*rdm1_bd_a(ib,id)
     &                  - rdm1_ad_a(ia,id)*rdm1_bc_a(ib,ic)
     &                  - rdm1_ad_b(ia,id)*rdm1_bc_b(ib,ic)
                    enddo ! ia
                  enddo ! ic
                enddo ! ib
              enddo ! id
c
c             Construct contributions to the natural orbital Fock
c             matrices.
c
              fn_ac_a(ia_lo:ia_hi,ic_lo:ic_hi) => buf_fn_ac_a
              fn_ac_b(ia_lo:ia_hi,ic_lo:ic_hi) => buf_fn_ac_b
              fn_ad_a(ia_lo:ia_hi,id_lo:id_hi) => buf_fn_ad_a
              fn_ad_b(ia_lo:ia_hi,id_lo:id_hi) => buf_fn_ad_b
              fn_ac_a = 0.0d0
              fn_ac_b = 0.0d0
              fn_ad_a = 0.0d0
              fn_ad_b = 0.0d0
c
              fn_ic_a(1:nmo,ic_lo:ic_hi) => buf_fn_ic_a
              fn_ic_b(1:nmo,ic_lo:ic_hi) => buf_fn_ic_b
              fn_id_a(1:nmo,id_lo:id_hi) => buf_fn_id_a
              fn_id_b(1:nmo,id_lo:id_hi) => buf_fn_id_b
              fn_ic_a = 0.0d0
              fn_ic_b = 0.0d0
              fn_id_a = 0.0d0
              fn_id_b = 0.0d0
c
              fns_ac_a(ia_lo:ia_hi,ic_lo:ic_hi) => buf_fns_ac_a
              fns_ac_b(ia_lo:ia_hi,ic_lo:ic_hi) => buf_fns_ac_b
              fns_ad_a(ia_lo:ia_hi,id_lo:id_hi) => buf_fns_ad_a
              fns_ad_b(ia_lo:ia_hi,id_lo:id_hi) => buf_fns_ad_b
              fns_ac_a = 0.0d0
              fns_ac_b = 0.0d0
              fns_ad_a = 0.0d0
              fns_ad_b = 0.0d0
c
              fns_ic_a(1:nmo,ic_lo:ic_hi) => buf_fns_ic_a
              fns_ic_b(1:nmo,ic_lo:ic_hi) => buf_fns_ic_b
              fns_id_a(1:nmo,id_lo:id_hi) => buf_fns_id_a
              fns_id_b(1:nmo,id_lo:id_hi) => buf_fns_id_b
              fns_ic_a = 0.0d0
              fns_ic_b = 0.0d0
              fns_id_a = 0.0d0
              fns_id_b = 0.0d0
c
c             This does not work for the correlated approach
c             do id = id_lo, id_hi
c               do ib = ib_lo, ib_hi
c                 do ic = ic_lo, ic_hi
c                   do ia = ia_lo, ia_hi
c                     fn_ac_a(ia,ic) = fn_ac_a(ia,ic)
c    &                  + eri(ia,ic,ib,id)*rdm1_bd_a(ib,id)
c    &                  + eri(ia,ic,ib,id)*rdm1_bd_b(ib,id)
c                     fn_ac_b(ia,ic) = fn_ac_b(ia,ic)
c    &                  + eri(ia,ic,ib,id)*rdm1_bd_a(ib,id)
c    &                  + eri(ia,ic,ib,id)*rdm1_bd_b(ib,id)
c                     fn_ad_a(ia,id) = fn_ad_a(ia,id)
c    &                  - eri(ia,ic,ib,id)*rdm1_bc_a(ib,ic)
c                     fn_ad_b(ia,id) = fn_ad_b(ia,id)
c    &                  - eri(ia,ic,ib,id)*rdm1_bc_b(ib,ic)
c                   enddo ! ia
c                 enddo ! ic
c               enddo ! ib
c             enddo ! id
c
c             do ic = ic_lo, ic_hi
c               do ia = ia_lo, ia_hi
c                 do ii = 1, nmo
c                   fn_ic_a(ii,ic) = fn_ic_a(ii,ic)
c    &                + fn_ac_a(ia,ic)*occ_a(ii)*N_ai_a(ia,ii)
c                   fn_ic_b(ii,ic) = fn_ic_b(ii,ic)
c    &                + fn_ac_b(ia,ic)*occ_b(ii)*N_ai_b(ia,ii)
c                 enddo ! ii
c               enddo ! ia
c             enddo ! ic
c             do id = id_lo, id_hi
c               do ia = ia_lo, ia_hi
c                 do ii = 1, nmo
c                   fn_id_a(ii,id) = fn_id_a(ii,id)
c    &                + fn_ad_a(ia,id)*occ_a(ii)*N_ai_a(ia,ii)
c                   fn_id_b(ii,id) = fn_id_b(ii,id)
c    &                + fn_ad_b(ia,id)*occ_b(ii)*N_ai_b(ia,ii)
c                 enddo ! ii
c               enddo ! ia
c             enddo ! ic
c
c             do jj = 1, nmo, nj
c               fn_ij_a(1:nmo,jj:min(nmo,jj+nj-1)) => buf_fn_ij_a
c               fn_ij_b(1:nmo,jj:min(nmo,jj+nj-1)) => buf_fn_ij_b
c               fn_ij_a = 0.0d0
c               fn_ij_b = 0.0d0
c               do kk = jj, min(nmo,jj+nj-1)
c                 do ic = ic_lo, ic_hi
c                   do ii = 1, nmo
c                     fn_ij_a(ii,kk) = fn_ij_a(ii,kk)
c    &                  + fn_ic_a(ii,ic)*N_ci_a(ic,kk)
c                     fn_ij_b(ii,kk) = fn_ij_b(ii,kk)
c    &                  + fn_ic_b(ii,ic)*N_ci_b(ic,kk)
c                   enddo ! ii
c                 enddo ! ic
c                 do id = id_lo, id_hi
c                   do ii = 1, nmo
c                     fn_ij_a(ii,kk) = fn_ij_a(ii,kk)
c    &                  + fn_id_a(ii,id)*N_di_a(id,kk)
c                     fn_ij_b(ii,kk) = fn_ij_b(ii,kk)
c    &                  + fn_id_b(ii,id)*N_di_b(id,kk)
c                   enddo ! ii
c                 enddo ! id
c               enddo ! kk
c               call ga_acc(wfn1_deriv%g_dna,1,nmo,jj,min(nmo,jj+nj-1),
c    &                      fn_ij_a,nmo,1.0d0)
c               call ga_acc(wfn1_deriv%g_dnb,1,nmo,jj,min(nmo,jj+nj-1),
c    &                      fn_ij_b,nmo,1.0d0)
c             enddo ! jj
c
c             Reuse the fn_ac_* data to construct the Coulomb terms
c             to the correlation function Fock matrices.
c
c             do ii = 1, nmo
c               do ic = ic_lo, ic_hi
c                 do ia = ia_lo, ia_hi
c                   v_dc_a(ii) = v_dc_a(ii)
c    &                + fn_ac_a(ia,ic)*N_ai_a(ia,ii)*N_ci_a(ic,ii)
c                   v_dc_b(ii) = v_dc_b(ii)
c    &                + fn_ac_b(ia,ic)*N_ai_b(ia,ii)*N_ci_b(ic,ii)
c                 enddo ! ia
c               enddo ! ic
c               do id = id_lo, id_hi
c                 do ia = ia_lo, ia_hi
c                   v_dc_a(ii) = v_dc_a(ii)
c    &                + fn_ad_a(ia,id)*N_ai_a(ia,ii)*N_di_a(id,ii)
c                   v_dc_b(ii) = v_dc_b(ii)
c    &                + fn_ad_b(ia,id)*N_ai_b(ia,ii)*N_di_b(id,ii)
c                 enddo ! ia
c               enddo ! id
c             enddo ! ii
c
c             Add the 2RDM terms that come from the 2-electron
c             occupation number terms of [5].
c
              do jproc = 0, nproc-1
                call ga_distribution(g_x_aa,jproc,ii_lo,ii_hi,
     &                                            jj_lo,jj_hi)
                ni = ii_hi-ii_lo+1
c
c               d_aa(ii_lo:ii_hi,jj_lo:jj_hi)  => buf_d_aa
c               d_bb(ii_lo:ii_hi,jj_lo:jj_hi)  => buf_d_bb
                x_aa(ii_lo:ii_hi,jj_lo:jj_hi)  => buf_x_aa
                x_bb(ii_lo:ii_hi,jj_lo:jj_hi)  => buf_x_bb
                dx_aa(ii_lo:ii_hi,jj_lo:jj_hi) => buf_dx_aa
                dx_bb(ii_lo:ii_hi,jj_lo:jj_hi) => buf_dx_bb
                C_is_a(ii_lo:ii_hi,1:nmo)      => buf_C_is_a
                C_is_b(ii_lo:ii_hi,1:nmo)      => buf_C_is_b
                C_jt_a(jj_lo:jj_hi,1:nmo)      => buf_C_jt_a
                C_jt_b(jj_lo:jj_hi,1:nmo)      => buf_C_jt_b
c
c               call ga_get(g_2aa,ii_lo,ii_hi,jj_lo,jj_hi,d_aa,ni)
c               call ga_get(g_2bb,ii_lo,ii_hi,jj_lo,jj_hi,d_bb,ni)
c
                call ga_get(g_x_aa,ii_lo,ii_hi,jj_lo,jj_hi,x_aa,ni)
                call ga_get(g_x_bb,ii_lo,ii_hi,jj_lo,jj_hi,x_bb,ni)
                call ga_get(g_dx_aa,ii_lo,ii_hi,jj_lo,jj_hi,dx_aa,ni)
                call ga_get(g_dx_bb,ii_lo,ii_hi,jj_lo,jj_hi,dx_bb,ni)
                call ga_get(wfn1_wave%g_ca,ii_lo,ii_hi,1,nmo,C_is_a,ni)
                call ga_get(wfn1_wave%g_cb,ii_lo,ii_hi,1,nmo,C_is_b,ni)
                call ga_get(wfn1_wave%g_ca,jj_lo,jj_hi,1,nmo,C_jt_a,nj)
                call ga_get(wfn1_wave%g_cb,jj_lo,jj_hi,1,nmo,C_jt_b,nj)
c
                rdm2_acj_aa(ia_lo:ia_hi,ic_lo:ic_hi,jj_lo:jj_hi)
     &            => buf_rdm2_acj_aa
                rdm2_acj_ab(ia_lo:ia_hi,ic_lo:ic_hi,jj_lo:jj_hi)
     &            => buf_rdm2_acj_ab
                rdm2_acj_ba(ia_lo:ia_hi,ic_lo:ic_hi,jj_lo:jj_hi)
     &            => buf_rdm2_acj_ba
                rdm2_acj_bb(ia_lo:ia_hi,ic_lo:ic_hi,jj_lo:jj_hi)
     &            => buf_rdm2_acj_bb
                rdm2_adj_aa(ia_lo:ia_hi,id_lo:id_hi,jj_lo:jj_hi)
     &            => buf_rdm2_adj_aa
                rdm2_adj_ab(ia_lo:ia_hi,id_lo:id_hi,jj_lo:jj_hi)
     &            => buf_rdm2_adj_ab
                rdm2_adj_ba(ia_lo:ia_hi,id_lo:id_hi,jj_lo:jj_hi)
     &            => buf_rdm2_adj_ba
                rdm2_adj_bb(ia_lo:ia_hi,id_lo:id_hi,jj_lo:jj_hi)
     &            => buf_rdm2_adj_bb
c
                rdm2_acj_aa = 0.0d0
                rdm2_acj_ab = 0.0d0
                rdm2_acj_ba = 0.0d0
                rdm2_acj_bb = 0.0d0
                rdm2_adj_aa = 0.0d0
                rdm2_adj_ab = 0.0d0
                rdm2_adj_ba = 0.0d0
                rdm2_adj_bb = 0.0d0
c
                rdm2_d_acj_a(ia_lo:ia_hi,ic_lo:ic_hi,jj_lo:jj_hi)
     &            => buf_rdm2_d_acj_a
                rdm2_d_acj_b(ia_lo:ia_hi,ic_lo:ic_hi,jj_lo:jj_hi)
     &            => buf_rdm2_d_acj_b
                rdm2_d_adj_a(ia_lo:ia_hi,id_lo:id_hi,jj_lo:jj_hi)
     &            => buf_rdm2_d_adj_a
                rdm2_d_adj_b(ia_lo:ia_hi,id_lo:id_hi,jj_lo:jj_hi)
     &            => buf_rdm2_d_adj_b
c
                rdm2_d_acj_a = 0.0d0
                rdm2_d_acj_b = 0.0d0
                rdm2_d_adj_a = 0.0d0
                rdm2_d_adj_b = 0.0d0
c
                rdm2s_acj_aa(ia_lo:ia_hi,ic_lo:ic_hi,jj_lo:jj_hi)
     &            => buf_rdm2s_acj_aa
                rdm2s_acj_ab(ia_lo:ia_hi,ic_lo:ic_hi,jj_lo:jj_hi)
     &            => buf_rdm2s_acj_ab
                rdm2s_acj_ba(ia_lo:ia_hi,ic_lo:ic_hi,jj_lo:jj_hi)
     &            => buf_rdm2s_acj_ba
                rdm2s_acj_bb(ia_lo:ia_hi,ic_lo:ic_hi,jj_lo:jj_hi)
     &            => buf_rdm2s_acj_bb
                rdm2s_adj_aa(ia_lo:ia_hi,id_lo:id_hi,jj_lo:jj_hi)
     &            => buf_rdm2s_adj_aa
                rdm2s_adj_ab(ia_lo:ia_hi,id_lo:id_hi,jj_lo:jj_hi)
     &            => buf_rdm2s_adj_ab
                rdm2s_adj_ba(ia_lo:ia_hi,id_lo:id_hi,jj_lo:jj_hi)
     &            => buf_rdm2s_adj_ba
                rdm2s_adj_bb(ia_lo:ia_hi,id_lo:id_hi,jj_lo:jj_hi)
     &            => buf_rdm2s_adj_bb
c
                rdm2s_acj_aa = 0.0d0
                rdm2s_acj_ab = 0.0d0
                rdm2s_acj_ba = 0.0d0
                rdm2s_acj_bb = 0.0d0
                rdm2s_adj_aa = 0.0d0
                rdm2s_adj_ab = 0.0d0
                rdm2s_adj_ba = 0.0d0
                rdm2s_adj_bb = 0.0d0
c
                do jj = jj_lo, jj_hi
                  do ii = ii_lo, ii_hi
                    call wfn1_entropy2_aa(
     &                   occ_a(ii),occ_a(jj),
     &                   x_aa(ii,jj)+occ_a(ii)*occ_a(jj),
     &                   saaij,dsaaij_i,dsaaij_j,dsaaij_ij)
                    sabij    = 0.0d0
                    sbaij    = 0.0d0
                    dsabij_i = 0.0d0
                    dsabij_j = 0.0d0
                    dsbaij_i = 0.0d0
                    dsbaij_j = 0.0d0
c                   if (ii.eq.jj) then
                      call wfn1_entropy2_ab(ii,jj,occ_a(ii),occ_b(jj),
     &                                      sabij,
     &                                      dsabij_i,dsabij_j)
                      call wfn1_entropy2_ab(jj,ii,occ_b(ii),occ_a(jj),
     &                                      sbaij,
     &                                      dsbaij_i,dsbaij_j)
cDEBUG
c      write(*,*)'HVD: ii,jj,sabij = ',ii,jj,occ_a(ii),occ_b(jj),sabij
c      write(*,*)'HVD: ii,jj,sbaij = ',ii,jj,occ_b(ii),occ_a(jj),sbaij
cDEBUG
c                   endif
                    call wfn1_entropy2_aa(
     &                   occ_b(ii),occ_b(jj),
     &                   x_bb(ii,jj)+occ_b(ii)*occ_b(jj),
     &                   sbbij,dsbbij_i,dsbbij_j,dsbbij_ij)
                    do ic = ic_lo, ic_hi
                      do ia = ia_lo, ia_hi
                        rdm2_acj_aa(ia,ic,jj) = rdm2_acj_aa(ia,ic,jj)
     &                    + N_ai_a(ia,ii)*N_ci_a(ic,ii)*
     &                      (x_aa(ii,jj)+occ_a(ii)*occ_a(jj))
                        rdm2_acj_ab(ia,ic,jj) = rdm2_acj_ab(ia,ic,jj)
     &                    + N_ai_a(ia,ii)*N_ci_a(ic,ii)*
     &                      occ_a(ii)*occ_b(jj)
                        rdm2_acj_ba(ia,ic,jj) = rdm2_acj_ba(ia,ic,jj)
     &                    + N_ai_a(ia,ii)*N_ci_a(ic,ii)*
     &                      occ_b(ii)*occ_a(jj)
                        rdm2_acj_bb(ia,ic,jj) = rdm2_acj_bb(ia,ic,jj)
     &                    + N_ai_b(ia,ii)*N_ci_b(ic,ii)*
     &                      (x_bb(ii,jj)+occ_b(ii)*occ_b(jj))
c
                        rdm2s_acj_aa(ia,ic,jj) = rdm2s_acj_aa(ia,ic,jj)
     &                    + N_ai_a(ia,ii)*N_ci_a(ic,ii)*
     &                      (x_aa(ii,jj)+occ_a(ii)*occ_a(jj))*
     &                      (-saaij)
                        rdm2s_acj_ab(ia,ic,jj) = rdm2s_acj_ab(ia,ic,jj)
     &                    + N_ai_a(ia,ii)*N_ci_a(ic,ii)*
     &                      occ_a(ii)*occ_b(jj)*(-sabij)
                        rdm2s_acj_ba(ia,ic,jj) = rdm2s_acj_ba(ia,ic,jj)
     &                    + N_ai_a(ia,ii)*N_ci_a(ic,ii)*
     &                      occ_b(ii)*occ_a(jj)*(-sbaij)
                        rdm2s_acj_bb(ia,ic,jj) = rdm2s_acj_bb(ia,ic,jj)
     &                    + N_ai_b(ia,ii)*N_ci_b(ic,ii)*
     &                      (x_bb(ii,jj)+occ_b(ii)*occ_b(jj))*
     &                      (-sbbij)
                      enddo ! ia
                    enddo ! ic
                    do id = id_lo, id_hi
                      do ia = ia_lo, ia_hi
                        rdm2_adj_aa(ia,id,jj) = rdm2_adj_aa(ia,id,jj)
     &                    + N_ai_a(ia,ii)*N_di_a(id,ii)*
     &                      (x_aa(ii,jj)+occ_a(ii)*occ_a(jj))
c                       rdm2_adj_ab(ia,id,jj) = rdm2_adj_ab(ia,id,jj)
c    &                    + N_ai_a(ia,ii)*N_di_a(id,ii)*
c    &                      occ_a(ii)*occ_b(jj)*(1.0d0-sabij)
c                       rdm2_adj_ba(ia,id,jj) = rdm2_adj_ba(ia,id,jj)
c    &                    + N_ai_a(ia,ii)*N_di_a(id,ii)*
c    &                      occ_b(ii)*occ_a(jj)*(1.0d0-sbaij)
                        rdm2_adj_bb(ia,id,jj) = rdm2_adj_bb(ia,id,jj)
     &                    + N_ai_b(ia,ii)*N_di_b(id,ii)*
     &                      (x_bb(ii,jj)+occ_b(ii)*occ_b(jj))
c
                        rdm2s_adj_aa(ia,id,jj) = rdm2s_adj_aa(ia,id,jj)
     &                    + N_ai_a(ia,ii)*N_di_a(id,ii)*
     &                      (x_aa(ii,jj)+occ_a(ii)*occ_a(jj))*
     &                      (-saaij)
                        rdm2s_adj_bb(ia,id,jj) = rdm2s_adj_bb(ia,id,jj)
     &                    + N_ai_b(ia,ii)*N_di_b(id,ii)*
     &                      (x_bb(ii,jj)+occ_b(ii)*occ_b(jj))*
     &                      (-sbbij)
                      enddo ! ia
                    enddo ! id
                  enddo ! ii
                enddo ! jj
cDEBUG
                ac_aa_mx = max(ac_aa_mx,maxval(rdm2_acj_aa))
                ac_ab_mx = max(ac_ab_mx,maxval(rdm2_acj_ab))
                ac_ba_mx = max(ac_ba_mx,maxval(rdm2_acj_ba))
                ac_bb_mx = max(ac_bb_mx,maxval(rdm2_acj_bb))
                ad_aa_mx = max(ad_aa_mx,maxval(rdm2_adj_aa))
                ad_ab_mx = max(ad_ab_mx,maxval(rdm2_adj_ab))
                ad_ba_mx = max(ad_ba_mx,maxval(rdm2_adj_ba))
                ad_bb_mx = max(ad_bb_mx,maxval(rdm2_adj_bb))
cDEBUG
c
cDEBUG
c               do jj = jj_lo, jj_hi
c                 do ii = ii_lo, ii_hi
c                   do id = id_lo, id_hi
c                     do ib = ib_lo, ib_hi
c                       do ic = ic_lo, ic_hi
c                         do ia = ia_lo, ia_hi
c                           gaa(ia,ic,ib,id) = gaa(ia,ic,ib,id)
c    &                        + d_aa(ii,jj)*N_ai_a(ia,ii)*N_ci_a(ic,ii)*
c    &                          N_bi_a(ib,jj)*N_di_a(id,jj)
c    &                        + d_aa(ii,jj)*N_ai_a(ia,jj)*N_ci_a(ic,jj)*
c    &                          N_bi_a(ib,ii)*N_di_a(id,ii)
c    &                        - d_aa(ii,jj)*N_ai_a(ia,ii)*N_di_a(id,ii)*
c    &                          N_bi_a(ib,jj)*N_ci_a(ic,jj)
c    &                        - d_aa(ii,jj)*N_bi_a(ib,ii)*N_ci_a(ic,ii)*
c    &                          N_ai_a(ia,jj)*N_di_a(id,jj)
c                         enddo ! ia
c                       enddo ! ic
c                     enddo ! ib
c                   enddo ! id
c                 enddo ! ii
c               enddo ! jj
cDEBUG
c
                do jj = jj_lo, jj_hi
                  do id = id_lo, id_hi
                    do ib = ib_lo, ib_hi
                      do ic = ic_lo, ic_hi
                        do ia = ia_lo, ia_hi
c
                         rdm2(ia,ic,ib,id) = rdm2(ia,ic,ib,id)
     &                   -(rdm2_adj_aa(ia,id,jj)+rdm2s_adj_aa(ia,id,jj))
     &                    *N_bi_a(ib,jj)*N_ci_a(ic,jj)
     &                   +(rdm2_acj_aa(ia,ic,jj)+rdm2s_acj_aa(ia,ic,jj))
     &                    *N_bi_a(ib,jj)*N_di_a(id,jj)
     &                   +(rdm2_acj_ab(ia,ic,jj)+rdm2s_acj_ab(ia,ic,jj))
     &                    *N_bi_b(ib,jj)*N_di_b(id,jj)
     &                   +(rdm2_acj_ba(ia,ic,jj)+rdm2s_acj_ba(ia,ic,jj))
     &                    *N_bi_a(ib,jj)*N_di_a(id,jj)
     &                   -(rdm2_adj_bb(ia,id,jj)+rdm2s_adj_bb(ia,id,jj))
     &                    *N_bi_b(ib,jj)*N_ci_b(ic,jj)
     &                   +(rdm2_acj_bb(ia,ic,jj)+rdm2s_acj_bb(ia,ic,jj))
     &                    *N_bi_b(ib,jj)*N_di_b(id,jj)
#ifdef TEST_RDMS
c This is non-sense now
c                         gaa(ia,ic,ib,id) = gaa(ia,ic,ib,id)
c    &                      - rdm2_adj_aa(ia,id,jj)*N_bi_a(ib,jj)*
c    &                        N_ci_a(ic,jj)
c    &                      + rdm2_acj_aa(ia,ic,jj)*N_bi_a(ib,jj)*
c    &                        N_di_a(id,jj)
c                         gbb(ia,ic,ib,id) = gbb(ia,ic,ib,id)
c    &                      - rdm2_adj_bb(ia,id,jj)*N_bi_b(ib,jj)*
c    &                        N_ci_b(ic,jj)
c    &                      + rdm2_acj_bb(ia,ic,jj)*N_bi_b(ib,jj)*
c    &                        N_di_b(id,jj)
#endif
                        enddo ! ia
                      enddo ! ic
                    enddo ! ib
                  enddo ! id
                enddo ! jj
c
c               Calculate the natural orbital Fock matrix contributions
c               (Exchange term from the 2RDM occupation numbers)
c
                fn_bcj_a(ib_lo:ib_hi,ic_lo:ic_hi,jj_lo:jj_hi)
     &            => buf_fn_bcj_a
                fn_bcj_b(ib_lo:ib_hi,ic_lo:ic_hi,jj_lo:jj_hi)
     &            => buf_fn_bcj_b
                fn_bdj_a(ib_lo:ib_hi,id_lo:id_hi,jj_lo:jj_hi)
     &            => buf_fn_bdj_a
                fn_bdj_b(ib_lo:ib_hi,id_lo:id_hi,jj_lo:jj_hi)
     &            => buf_fn_bdj_b
                fn_bcj_a = 0.0d0
                fn_bcj_b = 0.0d0
                fn_bdj_a = 0.0d0
                fn_bdj_b = 0.0d0
c
                fns_bcj_a(ib_lo:ib_hi,ic_lo:ic_hi,jj_lo:jj_hi)
     &            => buf_fns_bcj_a
                fns_bcj_b(ib_lo:ib_hi,ic_lo:ic_hi,jj_lo:jj_hi)
     &            => buf_fns_bcj_b
                fns_bdj_a(ib_lo:ib_hi,id_lo:id_hi,jj_lo:jj_hi)
     &            => buf_fns_bdj_a
                fns_bdj_b(ib_lo:ib_hi,id_lo:id_hi,jj_lo:jj_hi)
     &            => buf_fns_bdj_b
                fns_bcj_a = 0.0d0
                fns_bcj_b = 0.0d0
                fns_bdj_a = 0.0d0
                fns_bdj_b = 0.0d0
c
                do jj = jj_lo, jj_hi
                  do id = id_lo, id_hi
                    do ib = ib_lo, ib_hi
                      do ic = ic_lo, ic_hi
                        do ia = ia_lo, ia_hi
                          fn_bcj_a(ib,ic,jj) = fn_bcj_a(ib,ic,jj)
     &                      + eri(ia,ic,ib,id)*rdm2_adj_aa(ia,id,jj)
!    &                      + eri(ia,ic,ib,id)*rdm2_adj_ba(ia,id,jj) !hvd
                          fn_bcj_b(ib,ic,jj) = fn_bcj_b(ib,ic,jj)
     &                      + eri(ia,ic,ib,id)*rdm2_adj_bb(ia,id,jj)
!    &                      + eri(ia,ic,ib,id)*rdm2_adj_ab(ia,id,jj)
                          fn_bdj_a(ib,id,jj) = fn_bdj_a(ib,id,jj)
     &                      + eri(ia,ic,ib,id)*rdm2_acj_aa(ia,ic,jj)
     &                      + eri(ia,ic,ib,id)*rdm2_acj_ba(ia,ic,jj)
                          fn_bdj_b(ib,id,jj) = fn_bdj_b(ib,id,jj)
     &                      + eri(ia,ic,ib,id)*rdm2_acj_bb(ia,ic,jj)
     &                      + eri(ia,ic,ib,id)*rdm2_acj_ab(ia,ic,jj)
c
                          fns_bcj_a(ib,ic,jj) = fns_bcj_a(ib,ic,jj)
     &                      + eri(ia,ic,ib,id)*rdm2s_adj_aa(ia,id,jj)
!    &                      + eri(ia,ic,ib,id)*rdm2s_adj_ba(ia,id,jj) !hvd
                          fns_bcj_b(ib,ic,jj) = fns_bcj_b(ib,ic,jj)
     &                      + eri(ia,ic,ib,id)*rdm2s_adj_bb(ia,id,jj)
!    &                      + eri(ia,ic,ib,id)*rdm2s_adj_ab(ia,id,jj)
                          fns_bdj_a(ib,id,jj) = fns_bdj_a(ib,id,jj)
     &                      + eri(ia,ic,ib,id)*rdm2s_acj_aa(ia,ic,jj)
     &                      + eri(ia,ic,ib,id)*rdm2s_acj_ba(ia,ic,jj)
                          fns_bdj_b(ib,id,jj) = fns_bdj_b(ib,id,jj)
     &                      + eri(ia,ic,ib,id)*rdm2s_acj_bb(ia,ic,jj)
     &                      + eri(ia,ic,ib,id)*rdm2s_acj_ab(ia,ic,jj)
                        enddo ! ia
                      enddo ! ic
                    enddo ! ib
                  enddo ! id
                enddo ! jj
c
                fn_bj_a(ib_lo:ib_hi,jj_lo:jj_hi) => buf_fn_bj_a
                fn_bj_b(ib_lo:ib_hi,jj_lo:jj_hi) => buf_fn_bj_b
                fn_bj_a = 0.0d0
                fn_bj_b = 0.0d0
c
                fns_bj_a(ib_lo:ib_hi,jj_lo:jj_hi) => buf_fns_bj_a
                fns_bj_b(ib_lo:ib_hi,jj_lo:jj_hi) => buf_fns_bj_b
                fns_bj_a = 0.0d0
                fns_bj_b = 0.0d0
c
                do jj = jj_lo, jj_hi
                  do ic = ic_lo, ic_hi
                    do ib = ib_lo, ib_hi
                      fn_bj_a(ib,jj) = fn_bj_a(ib,jj)
     &                  - fn_bcj_a(ib,ic,jj)*N_ci_a(ic,jj)
                      fn_bj_b(ib,jj) = fn_bj_b(ib,jj)
     &                  - fn_bcj_b(ib,ic,jj)*N_ci_b(ic,jj)
c
                      fns_bj_a(ib,jj) = fns_bj_a(ib,jj)
     &                  - fns_bcj_a(ib,ic,jj)*N_ci_a(ic,jj)
                      fns_bj_b(ib,jj) = fns_bj_b(ib,jj)
     &                  - fns_bcj_b(ib,ic,jj)*N_ci_b(ic,jj)
                    enddo ! ib
                  enddo ! ic
                  do id = id_lo, id_hi
                    do ib = ib_lo, ib_hi
                      fn_bj_a(ib,jj) = fn_bj_a(ib,jj)
     &                  + fn_bdj_a(ib,id,jj)*N_di_a(id,jj)
                      fn_bj_b(ib,jj) = fn_bj_b(ib,jj)
     &                  + fn_bdj_b(ib,id,jj)*N_di_b(id,jj)
c
                      fns_bj_a(ib,jj) = fns_bj_a(ib,jj)
     &                  + fns_bdj_a(ib,id,jj)*N_di_a(id,jj)
                      fns_bj_b(ib,jj) = fns_bj_b(ib,jj)
     &                  + fns_bdj_b(ib,id,jj)*N_di_b(id,jj)
                    enddo ! ib
                  enddo ! ic
                enddo ! jj
c
                fn_ij_a(1:nmo,jj_lo:jj_hi) => buf_fn_ij_a
                fn_ij_b(1:nmo,jj_lo:jj_hi) => buf_fn_ij_b
                fn_ij_a = 0.0d0
                fn_ij_b = 0.0d0
c
                fns_ij_a(1:nmo,jj_lo:jj_hi) => buf_fns_ij_a
                fns_ij_b(1:nmo,jj_lo:jj_hi) => buf_fns_ij_b
                fns_ij_a = 0.0d0
                fns_ij_b = 0.0d0
c
                do jj = jj_lo, jj_hi
                  do ii = 1, nmo
                    do ib = ib_lo, ib_hi
                      fn_ij_a(ii,jj) = fn_ij_a(ii,jj)
     &                  + fn_bj_a(ib,jj)*N_bi_a(ib,ii)
                      fn_ij_b(ii,jj) = fn_ij_b(ii,jj)
     &                  + fn_bj_b(ib,jj)*N_bi_b(ib,ii)
c
                      fns_ij_a(ii,jj) = fns_ij_a(ii,jj)
     &                  + fns_bj_a(ib,jj)*N_bi_a(ib,ii)
                      fns_ij_b(ii,jj) = fns_ij_b(ii,jj)
     &                  + fns_bj_b(ib,jj)*N_bi_b(ib,ii)
                    enddo ! ib
                  enddo ! ii
                enddo ! jj
cDEBUG
c               fn_ij_a = 0.5d0*fn_ij_a
c               fn_ij_b = 0.5d0*fn_ij_b
cDEBUG
c
                call ga_acc(wfn1_deriv%g_dna,1,nmo,jj_lo,jj_hi,
     &                      fn_ij_a,nmo,1.0d0)
                call ga_acc(wfn1_deriv%g_dnb,1,nmo,jj_lo,jj_hi,
     &                      fn_ij_b,nmo,1.0d0)
c
                call ga_acc(wfn1_part_drv%g_dna,1,nmo,jj_lo,jj_hi,
     &                      fns_ij_a,nmo,1.0d0)
                call ga_acc(wfn1_part_drv%g_dnb,1,nmo,jj_lo,jj_hi,
     &                      fns_ij_b,nmo,1.0d0)
c
c               Calculate the correlation function Fock matrix
c               contributions (Exchange term from the 2RDM occupation
c               numbers)
c
c We need a loop over all jj here
c Taken care of by the loop over jproc
                fc_acj_a(ia_lo:ia_hi,ic_lo:ic_hi,jj_lo:jj_hi)
     &            => buf_fc_acj_a
                fc_acj_b(ia_lo:ia_hi,ic_lo:ic_hi,jj_lo:jj_hi)
     &            => buf_fc_acj_b
                fc_adj_a(ia_lo:ia_hi,id_lo:id_hi,jj_lo:jj_hi)
     &            => buf_fc_adj_a
                fc_adj_b(ia_lo:ia_hi,id_lo:id_hi,jj_lo:jj_hi)
     &            => buf_fc_adj_b
                fc_acj_a = 0.0d0
                fc_acj_b = 0.0d0
                fc_adj_a = 0.0d0
                fc_adj_b = 0.0d0
c
c We need a loop over all ii here
c Taken care of by the loop over jproc
                fc_ij_aa(ii_lo:ii_hi,jj_lo:jj_hi)
     &            => buf_fc_ij_aa
                fc_ij_ab(ii_lo:ii_hi,jj_lo:jj_hi)
     &            => buf_fc_ij_ab
                fc_ij_ba(ii_lo:ii_hi,jj_lo:jj_hi)
     &            => buf_fc_ij_ba
                fc_ij_bb(ii_lo:ii_hi,jj_lo:jj_hi)
     &            => buf_fc_ij_bb
                fc_ij_aa = 0.0d0
                fc_ij_ab = 0.0d0
                fc_ij_ba = 0.0d0
                fc_ij_bb = 0.0d0
c               fc_bdj_a = 0.0d0
c               fc_bdj_b = 0.0d0
c
                do jj = jj_lo, jj_hi
                  do id = id_lo, id_hi
                    do ib = ib_lo, ib_hi
                      do ic = ic_lo, ic_hi
                        do ia = ia_lo, ia_hi
                          fc_acj_a(ia,ic,jj) = fc_acj_a(ia,ic,jj)
     &                    + eri(ia,ic,ib,id)*N_bi_a(ib,jj)*N_di_a(id,jj)
                          fc_acj_b(ia,ic,jj) = fc_acj_b(ia,ic,jj)
     &                    + eri(ia,ic,ib,id)*N_bi_b(ib,jj)*N_di_b(id,jj)
                          fc_adj_a(ia,id,jj) = fc_adj_a(ia,id,jj)
     &                    + eri(ia,ic,ib,id)*N_bi_a(ib,jj)*N_ci_a(ic,jj)
                          fc_adj_b(ia,id,jj) = fc_adj_b(ia,id,jj)
     &                    + eri(ia,ic,ib,id)*N_bi_b(ib,jj)*N_ci_b(ic,jj)
                        enddo ! ia
                      enddo ! ic
                    enddo ! ib
                  enddo ! id
                enddo ! jj
c
                do jj = jj_lo, jj_hi
                  do ii = ii_lo, ii_hi
                    do ic = ic_lo, ic_hi
                      do ia = ia_lo, ia_hi
                        fc_ij_aa(ii,jj) = fc_ij_aa(ii,jj)
     &                  + fc_acj_a(ia,ic,jj)*N_ai_a(ia,ii)*N_ci_a(ic,ii)
                        fc_ij_ab(ii,jj) = fc_ij_ab(ii,jj)
     &                  + fc_acj_b(ia,ic,jj)*N_ai_a(ia,ii)*N_ci_a(ic,ii)
                        fc_ij_ba(ii,jj) = fc_ij_ba(ii,jj)
     &                  + fc_acj_a(ia,ic,jj)*N_ai_b(ia,ii)*N_ci_b(ic,ii)
                        fc_ij_bb(ii,jj) = fc_ij_bb(ii,jj)
     &                  + fc_acj_b(ia,ic,jj)*N_ai_b(ia,ii)*N_ci_b(ic,ii)
                      enddo ! ia
                    enddo ! ic
c We need to compute all ii and jj
c Taken care of by the loop over jproc
                    do id = id_lo, id_hi
                      do ia = ia_lo, ia_hi
                        fc_ij_aa(ii,jj) = fc_ij_aa(ii,jj)
     &                  - fc_adj_a(ia,id,jj)*N_ai_a(ia,ii)*N_di_a(id,ii)
!                       fc_ij_ab(ii,jj) = fc_ij_ab(ii,jj) !hvd
!    &                  - fc_adj_b(ia,id,jj)*N_ai_a(ia,ii)*N_di_a(id,ii)
!                       fc_ij_ba(ii,jj) = fc_ij_ba(ii,jj)
!    &                  - fc_adj_a(ia,id,jj)*N_ai_b(ia,ii)*N_di_b(id,ii)
                        fc_ij_bb(ii,jj) = fc_ij_bb(ii,jj)
     &                  - fc_adj_b(ia,id,jj)*N_ai_b(ia,ii)*N_di_b(id,ii)
                      enddo ! ia
                    enddo ! id
                  enddo ! ii
                enddo ! jj
c
c               All non-exchange derivative terms (the derivatives
c               of the exchange occupation numbers are dealt with
c               towards the end of the subroutine
c
c               1. Terms differentiated occupation numbers in
c                  position ii
c                  Alpha result:
c                  1.1 occ_a(ii)*occ_a(jj)
c                  1.2 occ_a(ii)*occ_b(jj)
c                  Beta result:
c                  1.3 occ_b(ii)*occ_a(jj)
c                  1.4 occ_b(ii)*occ_b(jj)
c               2. Term differentiated entropy in position ii
c                  Alpha result:
c                  2.1 x_aa(ii,jj)+occ_a(ii)*occ_a(jj)
c                  2.2 occ_a(ii)*occ_b(jj)
c                  Beta result:
c                  2.3 occ_b(ii)*occ_a(jj)
c                  2.4 x_bb(ii,jj)+occ_b(ii)*occ_b(jj)
c
                do jj = jj_lo, jj_hi
                  do ii = ii_lo, ii_hi
                    call wfn1_entropy2_aa(
     &                   occ_a(ii),occ_a(jj),
     &                   x_aa(ii,jj)+occ_a(ii)*occ_a(jj),
     &                   saaij,dsaaij_i,dsaaij_j,dsaaij_ij)
                    sabij    = 0.0d0
                    sbaij    = 0.0d0
                    dsabij_i = 0.0d0
                    dsabij_j = 0.0d0
                    dsbaij_i = 0.0d0
                    dsbaij_j = 0.0d0
c                   if (ii.eq.jj) then
                      call wfn1_entropy2_ab(ii,jj,occ_a(ii),occ_b(jj),
     &                                      sabij,
     &                                      dsabij_i,dsabij_j)
                      call wfn1_entropy2_ab(jj,ii,occ_b(ii),occ_a(jj),
     &                                      sbaij,
     &                                      dsbaij_i,dsbaij_j)
c                   endif
                    call wfn1_entropy2_aa(
     &                   occ_b(ii),occ_b(jj),
     &                   x_bb(ii,jj)+occ_b(ii)*occ_b(jj),
     &                   sbbij,dsbbij_i,dsbbij_j,dsbbij_ij)
                    v_dc_a(ii) = v_dc_a(ii)
!                   1.1 occ_a(ii)*occ_a(jj)
c    &              + 0.25d0*(1.0d0-saaij)*occ_a(jj)*fc_ij_aa(ii,jj)
     &              + 1.00d0*(-saaij)*occ_a(jj)*fc_ij_aa(ii,jj)
!                   1.2 occ_a(ii)*occ_b(jj)
c    &              + 0.25d0*(1.0d0-sabij)*occ_b(jj)*fc_ij_ab(ii,jj)
     &              + 1.00d0*(-sabij)*occ_b(jj)*fc_ij_ab(ii,jj)
                    v_dc_b(ii) = v_dc_b(ii)
!                   1.3 occ_b(ii)*occ_a(jj)
c    &              + 0.25d0*(1.0d0-sbaij)*occ_a(jj)*fc_ij_ba(ii,jj)
     &              + 1.00d0*(-sbaij)*occ_a(jj)*fc_ij_ba(ii,jj)
!                   1.4 occ_b(ii)*occ_b(jj)
c    &              + 0.25d0*(1.0d0-sbbij)*occ_b(jj)*fc_ij_bb(ii,jj)
     &              + 1.00d0*(-sbbij)*occ_b(jj)*fc_ij_bb(ii,jj)
                    v_dc_a(ii) = v_dc_a(ii)
!                   2.1 x_aa(ii,jj)+occ_a(ii)*occ_a(jj)
     &              - 1.0d0*(dsaaij_i+dsaaij_ij*occ_a(jj))*
     &                (x_aa(ii,jj)+occ_a(ii)*occ_a(jj))*
     &                fc_ij_aa(ii,jj)
!                   2.2 occ_a(ii)*occ_b(jj)
     &              - 1.0d0*dsabij_i*occ_a(ii)*occ_b(jj)*fc_ij_ab(ii,jj)
                    v_dc_b(ii) = v_dc_b(ii)
!                   2.3 occ_b(ii)*occ_a(jj)
     &              - 1.0d0*dsbaij_i*occ_b(ii)*occ_a(jj)*fc_ij_ba(ii,jj)
!                   2.4 x_bb(ii,jj)+occ_b(ii)*occ_b(jj)
     &              - 1.0d0*(dsbbij_i+dsbbij_ij*occ_b(jj))*
     &                (x_bb(ii,jj)+occ_b(ii)*occ_b(jj))*
     &                fc_ij_bb(ii,jj)
c
                    v_dp_a(ii) = v_dp_a(ii)
!                   1.1 occ_a(ii)*occ_a(jj)
c    &              + 0.25d0*(1.0d0-saaij)*occ_a(jj)*fc_ij_aa(ii,jj)
     &              + 1.00d0*occ_a(jj)*fc_ij_aa(ii,jj)
!                   1.2 occ_a(ii)*occ_b(jj)
c    &              + 0.25d0*(1.0d0-sabij)*occ_b(jj)*fc_ij_ab(ii,jj)
     &              + 1.00d0*occ_b(jj)*fc_ij_ab(ii,jj)
                    v_dp_b(ii) = v_dp_b(ii)
!                   1.3 occ_b(ii)*occ_a(jj)
c    &              + 0.25d0*(1.0d0-sbaij)*occ_a(jj)*fc_ij_ba(ii,jj)
     &              + 1.00d0*occ_a(jj)*fc_ij_ba(ii,jj)
!                   1.4 occ_b(ii)*occ_b(jj)
c    &              + 0.25d0*(1.0d0-sbbij)*occ_b(jj)*fc_ij_bb(ii,jj)
     &              + 1.00d0*occ_b(jj)*fc_ij_bb(ii,jj)
                  enddo ! ii
                enddo ! jj
cDEBUG
c               fc_ij_a = 0.5d0*fc_ij_a
c               fc_ij_b = 0.5d0*fc_ij_b
cDEBUG
c
                call ga_acc(wfn1_deriv%g_dca,ii_lo,ii_hi,jj_lo,jj_hi,
     &                      fc_ij_aa,ni,1.0d0)
                call ga_acc(wfn1_deriv%g_dcb,ii_lo,ii_hi,jj_lo,jj_hi,
     &                      fc_ij_bb,ni,1.0d0)
c
                call ga_acc(wfn1_part_drv%g_dca,ii_lo,ii_hi,jj_lo,jj_hi,
     &                      fc_ij_aa,ni,1.0d0)
                call ga_acc(wfn1_part_drv%g_dcb,ii_lo,ii_hi,jj_lo,jj_hi,
     &                      fc_ij_bb,ni,1.0d0)
c
              enddo ! jproc
c
c             Add the 2-electron energy terms
c
              E_2el     = E_2el     + ddot(na*nb*nc*nd,eri,1,rdm2,1)
              E_2el_chf = E_2el_chf + ddot(na*nb*nc*nd,eri,1,rdm2_chf,1)
c
            enddo ! la
          enddo ! lc
        enddo ! lb
      enddo ! ld
c
c     At this point g_dca and g_dcb contain exchange part of the 
c     2-electron occupation numbers. We have not taken into account
c     the derivatives of the 2-electron exchange occupation numbers.
c     Those factor are introduced here.
c
c     It is assumed that the data distribution of the correlation 
c     function Fock matrices and the exchange occupation number 
c     derivatives are the same.
c
      call ga_sync()
c
c     do alpha-alpha-spin
c
      call ga_distribution(wfn1_deriv%g_dca,iproc,if_lo,if_hi,
     &                     jf_lo,jf_hi)
      call ga_distribution(wfn1_part_drv%g_dca,iproc,ip_lo,ip_hi,
     &                     jp_lo,jp_hi)
      call ga_distribution(g_dx_aa,iproc,ix_lo,ix_hi,jx_lo,jx_hi)
      ok_dist =             if_lo.eq.ip_lo.and.if_hi.eq.ip_hi
      ok_dist = ok_dist.and.jf_lo.eq.jp_lo.and.jf_hi.eq.jp_hi
      ok_dist = ok_dist.and.if_lo.eq.ix_lo.and.if_hi.eq.ix_hi
      ok_dist = ok_dist.and.jf_lo.eq.jx_lo.and.jf_hi.eq.jx_hi
      if (.not.ok_dist) then
        call errquit("Wfn1_e_g_2el: Data distribution mismatch",0,UERR)
      endif
      call ga_access(wfn1_deriv%g_dca,if_lo,if_hi,jf_lo,jf_hi,idxf,ldf)
      call ga_access(wfn1_part_drv%g_dca,ip_lo,ip_hi,jp_lo,jp_hi,
     &               idxp,ldp)
      call ga_access(g_dx_aa,ix_lo,ix_hi,jx_lo,jx_hi,idxx,ldx)
      do jj=0,jf_hi-jf_lo
        do ii=0,if_hi-if_lo
          call wfn1_entropy2_aa(occ_a(ii),occ_a(jj),
     &                          x_aa(ii,jj)+occ_a(ii)*occ_a(jj),saaij,
     &                          dsaaij_i,dsaaij_j,dsaaij_ij)
          dbl_mb(idxf+ii+jj*ldf) = dbl_mb(idxf+ii+jj*ldf)
     &                           * ((     -saaij)-dsaaij_ij)
     &                           * dbl_mb(idxx+ii+jj*ldx)
          dbl_mb(idxp+ii+jj*ldp) = dbl_mb(idxp+ii+jj*ldp)
     &                           * ((1.0d0      )          )
     &                           * dbl_mb(idxx+ii+jj*ldx)
        enddo
      enddo
c
c     do beta-beta-spin
c
      call ga_distribution(wfn1_deriv%g_dcb,iproc,if_lo,if_hi,
     &                     jf_lo,jf_hi)
      call ga_distribution(wfn1_part_drv%g_dcb,iproc,ip_lo,ip_hi,
     &                     jp_lo,jp_hi)
      call ga_distribution(g_dx_bb,iproc,ix_lo,ix_hi,jx_lo,jx_hi)
      ok_dist =             if_lo.eq.ip_lo.and.if_hi.eq.ip_hi
      ok_dist = ok_dist.and.jf_lo.eq.jp_lo.and.jf_hi.eq.jp_hi
      ok_dist = ok_dist.and.if_lo.eq.ix_lo.and.if_hi.eq.ix_hi
      ok_dist = ok_dist.and.jf_lo.eq.jx_lo.and.jf_hi.eq.jx_hi
      if (.not.ok_dist) then
        call errquit("Wfn1_e_g_2el: Data distribution mismatch",0,UERR)
      endif
      call ga_access(wfn1_deriv%g_dcb,if_lo,if_hi,jf_lo,jf_hi,idxf,ldf)
      call ga_access(wfn1_part_drv%g_dcb,ip_lo,ip_hi,jp_lo,jp_hi,
     &               idxp,ldp)
      call ga_access(g_dx_bb,ix_lo,ix_hi,jx_lo,jx_hi,idxx,ldx)
      do jj=0,jf_hi-jf_lo
        do ii=0,if_hi-if_lo
          call wfn1_entropy2_aa(occ_b(ii),occ_b(jj),
     &                          x_bb(ii,jj)+occ_b(ii)*occ_b(jj),sbbij,
     &                          dsbbij_i,dsbbij_j,dsbbij_ij)
          dbl_mb(idxf+ii+jj*ldf) = dbl_mb(idxf+ii+jj*ldf)
     &                           * ((     -sbbij)-dsbbij_ij)
     &                           * dbl_mb(idxx+ii+jj*ldx)
          dbl_mb(idxp+ii+jj*ldp) = dbl_mb(idxp+ii+jj*ldp)
     &                           * ((1.0d0      )          )
     &                           * dbl_mb(idxx+ii+jj*ldx)
        enddo
      enddo
c
      E_2el     = 0.5d0*E_2el
      E_2el_chf = 0.5d0*E_2el_chf
c
      call ga_dgop(1300,E_2el,1,"+")
      call ga_dgop(1301,E_2el_chf,1,"+")
cDEBUG
c     write(*,*)'HVD: ac_aa_mx = ',ac_aa_mx
c     write(*,*)'HVD: ac_ab_mx = ',ac_ab_mx
c     write(*,*)'HVD: ac_ba_mx = ',ac_ba_mx
c     write(*,*)'HVD: ac_bb_mx = ',ac_bb_mx
c     write(*,*)'HVD: ad_aa_mx = ',ad_aa_mx
c     write(*,*)'HVD: ad_ab_mx = ',ad_ab_mx
c     write(*,*)'HVD: ad_ba_mx = ',ad_ba_mx
c     write(*,*)'HVD: ad_bb_mx = ',ad_bb_mx
cDEBUG
c
#ifdef TEST_RDMS
      flush(6)
c     gaa = 0.25d0*gaa
      call wfn1_test_rdms(wfn1_param,wfn1_mat,da,db,gaa,gab,gbb,
     &                    1.0d-7)
#endif
c
c     fc_bcj_a(0:0,0:0,0:0) => dummy
c     fc_bcj_b(0:0,0:0,0:0) => dummy
c     fc_bdj_a(0:0,0:0,0:0) => dummy
c     fc_bdj_b(0:0,0:0,0:0) => dummy
c     fn_ij_a(0:0,0:0) => dummy
c     fn_ij_b(0:0,0:0) => dummy
c     fn_bj_a(0:0,0:0) => dummy
c     fn_bj_b(0:0,0:0) => dummy
c     fn_ic_a(0:0,0:0) => dummy
c     fn_ic_b(0:0,0:0) => dummy
c     fn_id_a(0:0,0:0) => dummy
c     fn_id_b(0:0,0:0) => dummy
c
      deallocate(buf_fc_acj_a,buf_fc_acj_b)
      deallocate(buf_fc_adj_a,buf_fc_adj_b)
      deallocate(buf_fc_ij_aa,buf_fc_ij_bb)
      deallocate(buf_fc_ij_ab,buf_fc_ij_ba)
c
      deallocate(buf_fns_ij_a,buf_fns_ij_b)
      deallocate(buf_fns_bj_a,buf_fns_bj_b)
      deallocate(buf_fns_ic_a,buf_fns_ic_b)
      deallocate(buf_fns_id_a,buf_fns_id_b)
      deallocate(buf_fns_bcj_a,buf_fns_bcj_b)
      deallocate(buf_fns_bdj_a,buf_fns_bdj_b)
      deallocate(buf_fns_ac_a,buf_fns_ac_b)
      deallocate(buf_fns_ad_a,buf_fns_ad_b)
c
      deallocate(buf_fn_ij_a,buf_fn_ij_b)
      deallocate(buf_fn_bj_a,buf_fn_bj_b)
      deallocate(buf_fn_ic_a,buf_fn_ic_b)
      deallocate(buf_fn_id_a,buf_fn_id_b)
      deallocate(buf_fn_bcj_a,buf_fn_bcj_b)
      deallocate(buf_fn_bdj_a,buf_fn_bdj_b)
      deallocate(buf_fn_ac_a,buf_fn_ac_b)
      deallocate(buf_fn_ad_a,buf_fn_ad_b)
c
      deallocate(buf_dx_aa,buf_dx_bb)
      deallocate(buf_x_aa, buf_x_bb)
cDEBUG
c     deallocate(buf_d_aa)
cDEBUG
c
      deallocate(buf_rdm1_bc_a,buf_rdm1_bc_b)
      deallocate(buf_rdm1_ad_a,buf_rdm1_ad_b)
c
      deallocate(buf_rdm1_bd_a,buf_rdm1_bd_b)
      deallocate(buf_rdm1_ac_a,buf_rdm1_ac_b)
c
      deallocate(buf_N_di_a,buf_N_di_b)
      deallocate(buf_N_ci_a,buf_N_ci_b)
      deallocate(buf_N_bi_a,buf_N_bi_b)
      deallocate(buf_N_ai_a,buf_N_ai_b)
c
      deallocate(buf_rdm2_d_acj_a,buf_rdm2_d_acj_b)
      deallocate(buf_rdm2_d_adj_a,buf_rdm2_d_adj_b)
c
      deallocate(buf_rdm2s_acj_aa,buf_rdm2s_acj_bb)
      deallocate(buf_rdm2s_acj_ab,buf_rdm2s_acj_ba)
      deallocate(buf_rdm2s_adj_aa,buf_rdm2s_adj_bb)
      deallocate(buf_rdm2s_adj_ab,buf_rdm2s_adj_ba)
c
      deallocate(buf_rdm2_acj_aa,buf_rdm2_acj_bb)
      deallocate(buf_rdm2_acj_ab,buf_rdm2_acj_ba)
      deallocate(buf_rdm2_adj_aa,buf_rdm2_adj_bb)
      deallocate(buf_rdm2_adj_ab,buf_rdm2_adj_ba)
c
      deallocate(buf_eri,buf_rdm2,buf_rdm2_chf)
      deallocate(scr)
      deallocate(occ_a,occ_b)
#ifdef TEST_RDMS
      deallocate(da,db)
      deallocate(gaa,gab,gbb)
#endif
c
      end subroutine wfn1_e_g_2el_s
!>
!> @}
c
c-----------------------------------------------------------------------
