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
       for ib1 in range(1,5):
         for ja1 in range(1,5):
           for jb1 in range(1,5):
             t4 = (ia1,ib1,ja1,jb1)
             files[t4]=""
     for ia1 in range(1,5):
       for ja1 in range(1,5):
         for ib1 in range (1,5):
           for jb1 in range (1,5):
             for ka2 in range(1,5):
               if (ka2 != ia1 and ka2 != ja1):
                 for kb2 in range (1,5):
                   if (kb2 != ib1 and kb2 != jb1):
                     faci=multiply(ia1,ka2,ib1,kb2)
                     facj=multiply(ja1,ka2,jb1,kb2)
                     if (faci==1 and facj==1):
                       if (ia1==ja1 and ib1==jb1):
                         t4 = (ia1,ib1,ja1,jb1)
                         files[t4] = files[t4] + format(" + \C{%d%d}{%d%d}^2" % (min(ia1,ka2),max(ia1,ka2),min(ib1,kb2),max(ib1,kb2)))
                       else:
                         t4 = (ia1,ib1,ja1,jb1)
                         fc = 1
                         if (ia1 > ka2):
                           fc = -fc
                           i1 = ka2
                           i2 = ia1
                         else:
                           i1 = ia1
                           i2 = ka2
                         if (ib1 > kb2):
                           fc = -fc
                           i3 = kb2
                           i4 = ib1
                         else:
                           i3 = ib1
                           i4 = kb2
                         if (ja1 > ka2):
                           fc = -fc
                           i5 = ka2
                           i6 = ja1
                         else:
                           i5 = ja1
                           i6 = ka2
                         if (jb1 > kb2):
                           fc = -fc
                           i7 = kb2
                           i8 = jb1
                         else:
                           i7 = jb1
                           i8 = kb2
                         if (fc > 0):
                           files[t4] = files[t4] + format(" + \C{%d%d}{%d%d}\C{%d%d}{%d%d}" % (i1,i2,i3,i4,i5,i6,i7,i8))
                         else:
                           files[t4] = files[t4] + format(" - \C{%d%d}{%d%d}\C{%d%d}{%d%d}" % (i1,i2,i3,i4,i5,i6,i7,i8))
     print("\\begin{eqnarray}")
     print("   \\left(\\begin{array}{cccccccccccccccc}") 
     for ia1 in range(1,5):
       for ib1 in range(1,5):
         for ja1 in range(1,5):
           for jb1 in range(1,5):
             t4 = (ia1,ib1,ja1,jb1)
             if (len(files[t4]) == 0):
               files[t4] = "0"
             if (ja1 == 4 and jb1 == 4):
               print("    %s \\\\" % files[t4])
             else:
               print("    %s &" % files[t4])
     print("   \\end{array}\\right)") 
     print("\\end{eqnarray}")

if __name__ == "__main__":
    iterate()
