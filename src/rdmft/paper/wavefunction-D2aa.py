'''
Generate a wave function for 4 electrons in 4 orbitals

The chemical system is 4 H atoms in a row with equal bond distances.
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


def multiply(i1,i2,i3,i4):
     '''
     Multiply the symmetry operations to get the overall symmetry.
     '''
     return symmetry(i1)*symmetry(i2)*symmetry(i3)*symmetry(i4)

def iterate():
     '''
     Loop over all possible pairs of determinants and print out
     the contributions to the 1-electron density matrix.
     '''
     files={}
     for ia1 in range(1,5):
         for ia2 in range(ia1+1,5):
             for ja1 in range(1,5):
                 for ja2 in range(ja1+1,5):
                     t4 = (ia1,ia2,ja1,ja2)
                     files[t4]=""
     for ia1 in range(1,5):
         for ja1 in range(1,5):
             for ia2 in range (ia1+1,5):
                 for ja2 in range (ja1+1,5):
                     for ib1 in range(1,5):
                         for ib2 in range (ib1+1,5):
                             faci=multiply(ia1,ia2,ib1,ib2)
                             facj=multiply(ja1,ja2,ib1,ib2)
                             if (faci==1 and facj==1):
                               if (ia1==ja1 and ia2==ja2):
                                 t4 = (ia1,ia2,ja1,ja2)
                                 files[t4] = files[t4] + format(" + \C{%d%d}{%d%d}^2" % (ia1,ia2,ib1,ib2))
                               else:
                                 t4 = (ia1,ia2,ja1,ja2)
                                 files[t4] = files[t4] + format(" + \C{%d%d}{%d%d}\C{%d%d}{%d%d}" % (ia1,ia2,ib1,ib2,ja1,ja2,ib1,ib2))
     print("\\begin{eqnarray}")
     print("   \\left(\\begin{array}{cccccc}") 
     for ia1 in range(1,5):
         for ia2 in range(ia1+1,5):
             for ja1 in range(1,5):
                 for ja2 in range(ja1+1,5):
                     t4 = (ia1,ia2,ja1,ja2)
                     if (len(files[t4]) == 0):
                       files[t4] = "0"
                     if (ja1 == 3 and ja2 == 4):
                       print("    %s \\\\" % files[t4])
                     else:
                       print("    %s &" % files[t4])
     print("   \\end{array}\\right)") 
     print("\\end{eqnarray}")

if __name__ == "__main__":
    iterate()
