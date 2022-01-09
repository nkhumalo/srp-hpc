#!/bin/env python
#
# Translate a subset of all atoms in a PDB file.
#
def read_pdb_file(filename):
    '''
    Read the contents of the PDB file and return a list of lines.
    '''
    fp = open(filename,'r')
    pdb = fp.readlines()
    fp.close()
    return pdb

def get_cell_dimensions(pdb):
    '''
    Take the PDB file as a list of lines. 
    Read the cell dimensions from the first line (starting with "CRYST1")
    and return this as a tuple.
    '''
    line = pdb[0]
    elements = line.split()
    x = float(elements[1])
    y = float(elements[2])
    z = float(elements[3])
    return (x,y,z)

def gen_translation(direction,number,pdb):
    '''
    Generate the translation vector based on the choice of direction (given as a string),
    and the cell dimension of the PDB file.
    Return this a tuple.
    '''
    (x,y,z) = get_cell_dimensions(pdb)
    trans = []
    if direction == "x":
        for n in range(1,number+1):
            trans.append((n*x,0.0,0.0))
    elif direction == "y":
        for n in range(1,number+1):
            trans.append((0.0,n*y,0.0))
    elif direction == "z":
        for n in range(1,number+1):
            trans.append((0.0,0.0,n*z))
    else:
        print("Interesting coordinate system you are using. Your choice of direction is %s\n"%direction)
    return trans

def new_cell_dimensions(direction,number,pdb):
    '''
    Given the PDB data and the direction compute the new cell dimensions
    after extending the system in the given dimension. 
    Return the new dimensions as a tuple.
    '''
    (x,y,z) = get_cell_dimensions(pdb)
    if direction == "x":
        cell = (number*x,y,z)
    elif direction == "y":
        cell = (x,number*y,z)
    elif direction == "z":
        cell = (x,y,number*z)
    else:
        print("Interesting coordinate system you are using. Your choice of direction is %s\n"%direction)
    return cell

def split_pdb(pdb):
    '''
    Given the PDB file, separate the contents into three parts:
    - the cell definition
    - the solute atoms
    - the solvent atoms
    Return these results as a tuple of lists of lines.
    '''
    cryst   = []
    solute  = []
    solvent = []
    for line in pdb:
        elements = line.split()
        if elements[0] == "CRYST1":
            cryst.append(line)
        elif elements[0] == "ATOM":
            if elements[3] == "HOH":
                solvent.append(line)
            else:
                solute.append(line)
        elif elements[0] == "TER":
            solute.append(line)
    return (cryst,solute,solvent)

def translate_atoms(atomsin,target_atoms,translation):
    '''
    Take a given list of atoms and create a new list of atoms
    where all atoms in the range of target_atos are translated by
    the given translations.
    Return the new list of atoms. The input list remains unchanged.
    '''
    atomsout = []
    #DEBUG
    print(target_atoms)
    #DEBUG
    (dx,dy,dz) = translation
    for atom in atomsin:
        elements = atom.split()
        if elements[0] == "ATOM" or elements[0] == "HETATM":
            if int(atom[6:12]) in target_atoms:
                sx = atom[30:38]
                sy = atom[38:46]
                sz = atom[46:54]
                x  = float(sx)+float(dx)
                y  = float(sy)+float(dy)
                z  = float(sz)+float(dz)
                new_atom = "%s%8.3f%8.3f%8.3f%s"%(atom[:30],x,y,z,atom[54:])
            else:
                new_atom = atom
        else:
            new_atom = atom
        atomsout.append(new_atom)
    return atomsout

def new_pdb(cryst,oldsolute,oldsolvent,newsolute,newsolvent):
    '''
    Given the different sections of the old PDB and the translated parts of the 
    solute and solvent lists construct a new PDB file. 
    Return the list of lines of the new PDB file.
    '''
    end = ["END"]
    newpdb = cryst+oldsolute+newsolute+oldsolvent+newsolvent+end
    return newpdb

def new_cell(oldcryst,newcell):
    '''
    Given the old crystal line and the new cell dimensions construct 
    a new crystal line. 
    Return the list of new crystal lines.
    '''
    inline = oldcryst[0]
    (x,y,z) = newcell
    outline = "%s%9.3f%9.3f%9.3f%s"%(inline[:6],x,y,z,inline[33:])
    outlist = []
    outlist.append(outline)
    return outlist

def write_pdb(filename,pdb):
    '''
    Given the name of the output file and the list of PDB lines 
    write the data to the PDB file.
    '''
    fp = open(filename,'w')
    for line in pdb:
        fp.write(line)
    fp.close()

def parse_arguments():
    '''
    Parse command line arguments.
    '''
    from argparse import ArgumentParser
    prs = ArgumentParser(description='''
    A script to translate some atoms in a PDB file.

    You can provide a range of atom numbers following the usual
    Python convention (i.e. "10,16" produces a range including
    atoms 10, 11, 12, 13, 14, and 15). The script will iterate
    over all atoms in the PDB file and translate the atoms with
    serial numbers within the specified range.
    ''')
    prs.add_argument("infile",help="the input PDB file")
    prs.add_argument("outfile",help="the output PDB file")
    prs.add_argument("vector",help="the translation vector, e.g. \"1.0,0.5,0.3\"")
    prs.add_argument("atoms",help="the range of atom numbers to translate, e.g. \"5118,5175\"")
    args = prs.parse_args()
    return args

def execute_with_arguments(args):
    inputfile  = args.infile
    outputfile = args.outfile
    vec1 = args.vector.split(",")
    vector  = (float(vec1[0]),float(vec1[1]),float(vec1[2]))
    atm1 = args.atoms.split(",")
    atoms  = range(int(atm1[0]),int(atm1[1]))
    inputdata = read_pdb_file(inputfile)
    #trans = gen_translation(direction,number,inputdata)
    #newcell = new_cell_dimensions(direction,number,inputdata)
    #(cryst,solute,solvent) = split_pdb(inputdata)
    #newcryst = new_cell(cryst,newcell)
    #newsolute = translate_atoms(solute,trans)
    #newsolvent = translate_atoms(solvent,trans)
    outputdata = translate_atoms(inputdata,atoms,vector)
    write_pdb(outputfile,outputdata)

def main():
    execute_with_arguments(parse_arguments())

if __name__ == "__main__":
    main()
