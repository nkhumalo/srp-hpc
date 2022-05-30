#!/bin/env python
#
# Take a PDB file and a trajectory file in Amber's CRD format.
# Then write every frame from the trajectory file as a PDB file.
# The PDB files are numbered to distinguish them from each other.
#
def read_pdb_file(filename):
    '''
    Read the contents of the PDB file and return a list of lines.
    '''
    fp = open(filename,'r')
    pdb = fp.readlines()
    fp.close()
    return pdb

def move_atoms(atomsin,newpos):
    '''
    Take the positions of atoms specified in a PDB and a list of 
    new atom positions. For each atom replace the old position with
    the new position and add the line to the output structure.
    Return the output structure, the input structure remains unchanged.
    '''
    atomsout = []
    ii = -1
    for atom in atomsin:
        elements = atom.split()
        if elements[0] == "ATOM" or elements[0] == "HETATM":
            ii += 1
            newcrd = newpos[ii]
            (x,y,z) = newcrd
            new_atom = "%s%8.3f%8.3f%8.3f%s"%(atom[:30],x,y,z,atom[54:])
        else:
            new_atom = atom
        atomsout.append(new_atom)
    return atomsout

def num_atoms(pdbin):
    '''
    Take the input PDB file and count all atoms.
    '''
    ii = 0
    for atom in pdbin:
        elements = atom.split()
        if elements[0] == "ATOM" or elements[0] == "HETATM":
            ii += 1
    return ii

def write_pdb(filename,pdb):
    '''
    Given the name of the output file and the list of PDB lines 
    write the data to the PDB file.
    '''
    fp = open(filename,'w')
    for line in pdb:
        fp.write(line)
    fp.close()

def read_crdframe(fp,numatm):
    '''
    Given the current position in the CRD file and the number of 
    atoms in the molecule, read the next frame from the file.
    '''
    numcrd = 3*numatm
    # in a CRD file every line holds 10 coordinates
    numline = numcrd//10
    if numcrd%10 != 0:
        numline += 1
    crdstring = ""
    ii = 0
    while ii < numline:
        ii += 1
        line = fp.readline()
        if not line:
            # reached end-of-file
            raise StopIteration("Reached end of CRD file")
        crdstring = crdstring + " " + line
    # after a frame is read there is one more line specifying the box
    line = fp.readline()
    # now repackage the coordinates as a list of floating point triples
    crdlist = crdstring.split()
    ii = 0
    outcrd = []
    while ii < numatm:
        iix = 3*ii
        iiy = 3*ii+1
        iiz = 3*ii+2
        x   = float(crdlist[iix])
        y   = float(crdlist[iiy])
        z   = float(crdlist[iiz])
        outcrd.append((x,y,z))
        ii += 1
    return outcrd

def open_crdfile(filename):
    '''
    Given the filename open the CRD file and return the filepointer

    The CRD file has a comment on the first line. So we need to open
    the file, skip the first line and then return the filepointer.
    '''
    fp = open(filename)
    line = fp.readline()
    return fp

def close_crdfile(fp):
    '''
    Close the CRD file
    '''
    fp.close()

def parse_arguments():
    '''
    Parse command line arguments.
    '''
    from argparse import ArgumentParser
    prs = ArgumentParser(description='''
    A script to write a CRD file as a sequence of PDB files.

    The output PDB files will be numbered. In each one the coordinates
    are the ones from the corresponding frame in the CRD file.
    ''')
    prs.add_argument("infile",help="the input PDB file")
    prs.add_argument("crdfile",help="the input CRD file")
    prs.add_argument("outfile",help="the prefix for the output PDB files")
    args = prs.parse_args()
    return args

def execute_with_arguments(args):
    inputfile  = args.infile
    crdfile    = args.crdfile
    outputfile = args.outfile
    inputpdb   = read_pdb_file(inputfile)
    numatm     = num_atoms(inputpdb)
    crdfp      = open_crdfile(crdfile)
    ii         = 1
    try:
        while True:
            crdframe  = read_crdframe(crdfp,numatm)
            outputpdb = move_atoms(inputpdb,crdframe)
            outfile   = outputfile + "_" + f'{ii:05d}' + ".pdb"
            write_pdb(outfile,outputpdb)
            ii += 1
    except StopIteration:
        pass
    close_crdfile(crdfp)

def main():
    execute_with_arguments(parse_arguments())

if __name__ == "__main__":
    main()
