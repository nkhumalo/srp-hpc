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
     Loop over all possible determinants and print out any that
     match the symmetry criteria.
     '''
     for ia1 in range(1,5):
         for ia2 in range (ia1+1,5):
             for ib1 in range(1,5):
                 for ib2 in range (ib1+1,5):
                     fac=multiply(ia1,ia2,ib1,ib2)
                     if (fac==1):
                       print("+ \C{%d%d}{%d%d}\XXa{%d}{%d}\XXb{%d}{%d}" % (ia1,ia2,ib1,ib2,ia1,ia2,ib1,ib2))

if __name__ == "__main__":
    iterate()
