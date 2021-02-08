#$Id: GNUmakefile 25209 2014-02-02 18:41:05Z niri $
# makefile,v 1.8 1994/12/05 20:37:08 og845 Exp

   OBJ = wfn1_module.o \
         wfn1_2el_hf.o \
         wfn1_2el_wfn1.o \
         wfn1_closest_three.o \
         wfn1_df1df2df3.o \
         wfn1_dlasyf.o \
         wfn1_dnrm2.o \
         wfn1_dsytf2.o \
         wfn1_dsytrf.o \
         wfn1_dsytri.o \
         wfn1_element.o \
         wfn1_entropy.o \
         wfn1_exp.o \
         wfn1_f0df0f1.o \
         wfn1_f0f1df1.o \
         wfn1_f1f2f3.o \
         wfn1_fix_phase.o \
         wfn1_fix_xm.o \
         wfn1_gen_trial_wfn.o \
         wfn1_gsl_wrap.o \
         wfn1_idmin.o \
         wfn1_inpdft.o \
         wfn1_insert.o \
         wfn1_lindep.o \
         wfn1_linesearch.o \
         wfn1_print.o \
         wfn1_print_occ.o \
         wfn1_rot_natorb.o \
         wfn1_size.o \
         wfn1_skew.o \
         wfn1_symm.o \
         wfn1_test_rdms.o \
         wfn1_thimble.o \
         wfn1_transform.o \
         wfn1_write_movecs.o \
         wfn1_zero_vv.o 

#        wfn1.o \
#        wfn1_4index.o \
#        wfn1_allocate_drv.o \
#        wfn1_allocate_mat.o \
#        wfn1_allocate_wfn.o \
#        wfn1_backorder.o \
#        wfn1_calc_occ.o \
#        wfn1_column_max.o \
#        wfn1_configuration.o \
#        wfn1_diag.o \
#        wfn1_e_g.o \
#        wfn1_e_g_2el_corr.o \
#        wfn1_e_g_2el_hf.o \
#        wfn1_e_g_2el_mx.o \
#        wfn1_e_g_2el_vx.o \
#        wfn1_e_g_2el_s.o \
#        wfn1_e_g_dft.o \
#        wfn1_e_g_redist.o \
#        wfn1_gemm.o \
#        wfn1_gen_mat.o \
#        wfn1_guess_corr.o \
#        wfn1_guess_movecs.o \
#        wfn1_mcsearch.o \
#        wfn1_ortho.o \
#        wfn1_print_analysis.o \
#        wfn1_scale_offdiag.o \
#        wfn1_solve_diag.o \
#        wfn1_solve_line.o \
#        wfn1_solve_mc.o \
#        wfn1_trace.o \
#        wfn1_vectors.o \

   LIBRARY = librdmft.a

 LIB_DEFINES =
ifdef TEST_RDMS
 LIB_DEFINES += -DTEST_RDMS
endif
ifdef USE_GSL
 LIB_DEFINES += -DUSE_GSL
endif

   SUBDIRS = 
LIB_INCLUDES = 

HEADERS = 

include ../config/makefile.h
include ../config/makelib.h

