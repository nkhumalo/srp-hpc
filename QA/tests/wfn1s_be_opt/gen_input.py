#!/bin/env python
laa_fac_hh=[0.00,0.25,0.50,0.75,1.00]
laa_fac_hl=[0.00,0.25,0.50,0.75,1.00]
laa_fac_ll=[-0.25,0.00,0.25,0.50,0.75,1.00]
lab_fac_hh=[0.00]
lab_fac_hl=[0.00]
lab_fac_ll=[0.00]

for aa_fac_hh in laa_fac_hh:
  for aa_fac_hl in laa_fac_hl:
    for aa_fac_ll in laa_fac_ll:
      for ab_fac_hh in lab_fac_hh:
        for ab_fac_hl in lab_fac_hl:
          for ab_fac_ll in lab_fac_ll:
            test_name="wfn1s_be_aa_{0:.2f}_{1:.2f}_{2:.2f}_ab_{3:.2f}_{4:.2f}_{5:.2f}".format(
                      aa_fac_hh,aa_fac_hl,aa_fac_ll,ab_fac_hh,ab_fac_hl,ab_fac_ll)
            input_name=test_name+".nw"
            data_name=test_name+"_dat"
            vecs_name=data_name+".movecs"
            f = open(input_name,"w")
            f.write("echo\n")
            f.write("start {0}\n".format(data_name))
            f.write("geometry\n")
            f.write("  be 0 0 0\n")
            f.write("end\n")
            f.write("basis\n")
            f.write("  * library 6-31g*\n")
            f.write("end\n\n")
            f.write("set wfn1:corr_expr \"wfn1_mx\"\n")
            f.write("set wfn1:t_bath 0.05\n")
            f.write("set wfn1:maxit 3000\n")
            f.write("set wfn1:print_error T\n")
            f.write("task wfn1 energy\n\n")
            f.write("set \"wfn1:input vectors\" \"{0}\"\n".format(vecs_name))
            f.write("set wfn1:corr_expr \"wfn1s\"\n")
            f.write("set wfn1:t_bath 0.0\n")
            f.write("set wfn1:maxit 500\n")
            f.write("set wfn1:print_error T\n")
            f.write("set wfn1:aa_hh {0:.2f}\n".format(aa_fac_hh))
            f.write("set wfn1:aa_hl {0:.2f}\n".format(aa_fac_hl))
            f.write("set wfn1:aa_ll {0:.2f}\n".format(aa_fac_ll))
            f.write("set wfn1:ab_hh {0:.2f}\n".format(ab_fac_hh))
            f.write("set wfn1:ab_hl {0:.2f}\n".format(ab_fac_hl))
            f.write("set wfn1:ab_ll {0:.2f}\n".format(ab_fac_ll))
            f.write("task wfn1 energy\n")
            f.close()
