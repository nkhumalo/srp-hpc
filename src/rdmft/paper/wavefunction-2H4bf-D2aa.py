'''
Generate a wave function for 2 electrons in 4 orbitals

The chemical system is 2 H atoms in a row with equal bond distances.
The H2 molecules is aligned with Z-axis.
'''

def symmetry(i):
     '''
     Return the factor that the orbital 'i' gets under the mirror
     operation in the xy-plane.
     '''
     if i == 1:
         return 1
     elif i == 2:
         return -1
     elif i == 3:
         return 1
     elif i == 4:
         return -1
     else:
         print("symmetry invalid input %d" % i)


def multiply(i1,i2):
     '''
     Multiply the symmetry operations to get the overall symmetry of the determinant.
     '''
     return symmetry(i1)*symmetry(i2)

def iterate():
     '''
     Loop over all possible pairs of determinants and print out
     the contributions to the 1-electron density matrix.
     '''
     files={}
     for ia1 in range(1,5):
         for ib1 in range(1,5):
             for ja1 in range(1,5):
                 for jb1 in range(1,5):
                     t4 = (ia1,ja1,ib1,jb1)
                     files[t4]=""
     for ia1 in range(1,5):
         for ib1 in range(1,5):
             for ja1 in range (1,5):
                 for jb1 in range (1,5):
                     faci=multiply(ia1,ib1)
                     facj=multiply(ja1,jb1)
                     if (faci==1 and facj==1):
                         if (ia1==ja1 and ib1==jb1):
                             t4 = (ia1,ja1,ib1,jb1)
                             files[t4] = files[t4] + format("\Ctwo{%d}{%d}^2" % (ia1,ib1))
                         else:
                             t4 = (ia1,ja1,ib1,jb1)
                             files[t4] = files[t4] + format("\Ctwo{%d}{%d}\Ctwo{%d}{%d}^{*}" % (ia1,ib1,ja1,jb1))
     print("\\begin{eqnarray}")
     print("   D &=&") 
     print("   \\begin{pmatrix}") 
     for ia1 in range(1,5):
         for ib1 in range(1,5):
             for ja1 in range(1,5):
                 for jb1 in range(1,5):
                     t4 = (ia1,ja1,ib1,jb1)
                     if (len(files[t4]) == 0):
                       files[t4] = "0"
                     if (ja1 == 4 and jb1 == 4):
                       print("    %s \\\\" % files[t4])
                     else:
                       print("    %s &" % files[t4])
     print("   \\end{pmatrix}") 
     print("\\end{eqnarray}")

if __name__ == "__main__":
    iterate()
