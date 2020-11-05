#!/bin/env python
"""
This script takes a PDB file and generates the corresponding
AIMD QM/MM input file. This requires some special steps:

1. The PDB atom names need to be converted into atom names 
   that make sense to an QM code. I.e. the element needs to
   be clear.
2. The coordinates as given in the PDB file and translated
   atom names need to be output in an XYZ format.
3. A topology file needs to be generated that specifies all
   the force field parameters.
4. Topology file needs to be output in an input block that
   the PSPW code can understand. 

The overall execution is:

1. Take PDB file and write an input file for the NWChem
   prepare module.
2. Run the NWChem prepare module to generate a topology
   file.
3. Use the temporary PDB file to produce a geometry block.
4. Use the topology file to define:
   - masses
   - charges
   - bond springs for each bond
   - bond angle springs for all bond angles
   - torsion angle springs for all torsions

Note that we need to know which atoms are in the QM region
to avoid defining force field parameters between them.
"""

t_description="""
Generate an input file for the AIMD QM/MM capability for a biological
system specified by a PDB file.
"""
t_epilog="""
The list of QM atoms is specified as a comma separated list of
ranges. Each range may either be single atom number or the start
and end number separated by a colon. Examples:
  "1"
  "1,4,6,8"
  "1:10"
  "5:10,15:20"
"""

import os

def parse_arguments():
    """
    Parse the command line arguments and return them as
    fields in the namespace.
    """
    from argparse import ArgumentParser
    prs = ArgumentParser(prog="pdb_gen_aimd.py",description=t_description,
                         epilog=t_epilog)
    prs.add_argument("-i","--input",dest="input",
                     help="The PDB file of chemical system")
    prs.add_argument("-o","--output",dest="output",
                     help="The NW file for the AIMD run")
    prs.add_argument("-q","--qm",dest="qm_atoms",help="The QM atoms")
    prs.add_argument("-n","--nwchem",dest="nwchem_exe",
                     help="The NWChem executable")
    args = prs.parse_args()
    return args

def complement_arguments(args):
    """
    The command line arguments may not provide all necessary 
    information. Additional information might be obtained from
    the environment variables. Here we complement the command
    line arguments with this additional information.
    The complemented arguments namespace is returned.
    """
    if not args.nwchem_exe:
        args.nwchem_exe = os.environ.get("NWCHEM_EXECUTABLE")
    if not args.nwchem_exe:
        nwchem_top = os.environ.get("NWCHEM_TOP")
        nwchem_target = os.environ.get("NWCHEM_TARGET")
        if nwchem_top and nwchem_target:
            args.nwchem_exe = f"{nwchem_top}/bin/{nwchem_target}/nwchem"
    return args

def get_element(name):
    """
    Given an atom label determine the chemical element. The atom label
    is expected to correspond to an atom name that occurs in a PDB file.
    An atom label consists of 4 characters. The element information
    is contained in the first 2 characters. If the first character is a
    number this can be replace with a space. After that the first 2
    characters correspond to the chemical element.
    
    For example:
    
    Atom name     Chemical element
    " CA "        "C"
    "3HB "        "H"
    "Ca  "        "Ca"
    """
    el1=name[0:2]
    char0=el1[0:1]
    try:
        int(char0)
        el2=" "+el1[1:2]
    except ValueError:
        el2=el1
    if el2[0:1]==" ":
        el3=el2[1:2]
    else:
        el3=el2
    return el3
    

class base_atom:
    """
    A class for general properties of atoms.
    """
    def __init__(self,name,number,resname,resnumber):
        self.name=None
        self.number=None
        self.element=None
        self.residue_name=None
        self.residue_number=None
        if name:
            self.name=name
            self.element=get_element(name)
        if number:
            self.number=int(number)
        if resname:
            self.residue_name=resname
        if resnumber:
            self.residue_number=int(resnumber)

    def set_name(self,name):
        """
        Sets the name and chemical element of the atom.
        """
        self.name = name
        self.element = get_element(name)

    def set_residue_name(self,name):
        """
        Sets the residue name of the atom.
        """
        self.residue_name=name

    def set_number(self,number):
        """
        Sets the atom number.
        """
        self.number=int(number)

    def set_residue_number(self,number):
        """
        Sets the residue number.
        """
        self.residue_number=int(number)

class pdb_atom(base_atom):
    """
    A class for atoms from a PDB file.
    In addition to a base atom this atom also has spatial
    coordinates. Also it might an MM or QM atom.
    """
    def __init__(self,line,flag_large):
        """
        Instantiate a pdb_atom instance from a line of a
        PDB file. NWChem supports regular PDB files and
        large PDB files. The flag_large argument is True
        when we have a large PDB file format and False
        otherwise.
        The PDB file format is detailed at:
        https://www.wwpdb.org/documentation/file-format-content/format33/v3.3.html [accessed 11/04/2020]
        """
        record=line[0:6]
        if record != "ATOM  " and record != "HETATM":
            raise ValueError("Not an atom record")
        atmnum=line[6:11]
        atmnam=line[12:16]
        resnam=line[17:20]
        if flag_large:
            resnum=line[20:26]
        else:
            resnum=line[22:26]
        coords=line[30:54]
        base_atom.__init__(atmnam,atmnum,resnam,resnum)
        self.coordinates=coords
        self.qm=False

    def set_qm(self):
        """
        Make this atom a QM atom.
        """
        self.qm=True

    def set_mm(self):
        """
        Make this atom an MM atom.
        """
        self.qm=False
        

def write_prepare_input(pdb_filename,prepare_input_filename):
    """
    Write the input for the NWChem prepare module to generate the
    topology file. 
    """
    fp = open(prepare_input_filename,"w")
    fp.write("echo\n")
    fp.write("start pdb-gen-aimd-prepare-dat\n")
    fp.write("prepare\n")
    fp.write("  system pdb-gen-aimd-prepare-sys\n")
    fp.write(f"  source {pdb_filename}\n")
    fp.write("  new_top new_seq\n")
    fp.write("  new_rst\n")
    fp.write("end\n")
    fp.write("task prepare")
    fp.close()

def delete_prepare_input(prepare_input_filename):
    """
    Delete the input file for the NWChem prepare module.
    """
    if os.path.exists(prepare_input_filename):
        os.remove(prepare_input_filename)

def run_nwchem(nwchem_exe,input_filename):
    """
    Run NWChem in serial on the given input file.
    """
    os.system(f"{nwchem_exe} {input_filename} > {input_filename}.out")

def execute_with_arguments(args):
    """
    Execute the whole pipeline with the command line arguments given.
    """
    print(args)
    prepare_input = "pdb-gen-aimd-prepare.nw"
    write_prepare_input(args.input,prepare_input)
    run_nwchem(args.nwchem_exe,prepare_input)


if __name__ == "__main__":
    execute_with_arguments(complement_arguments(parse_arguments()))
