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
            args.nwchem_exe = str(nwchem_top)+"/bin/"+str(nwchem_target)+"/nwchem"
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


def get_lattice(line):
    """
    Extract the lattice size a, b, and c. We expect all lattice
    angle to be 90.0 degrees. Return the lattice as a tuple.
    """
    a = line[6:15]
    b = line[15:24]
    c = line[24:33]
    alpha = line[33:40]
    beta  = line[40:47]
    gamma = line[47:54]
    a = float(a)
    b = float(b)
    c = float(c)
    if (abs(float(alpha)-90.0)>1.0e-2):
        raise ValueError("Alpha not equal 90")
    if (abs(float(beta)-90.0)>1.0e-2):
        raise ValueError("Beta not equal 90")
    if (abs(float(gamma)-90.0)>1.0e-2):
        raise ValueError("Gamma not equal 90")
    return (a,b,c)


class base_atom:
    """
    A class for general properties of atoms.
    """
    def __init__(self,name,number,resname,resnumber):
        self.name=None
        self.number=None
        self.element=None
        self.atom_label = None # Name to be used in QM code
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

    def set_atom_label(self,atom_label):
        """
        The atom label is a unique identifier corresponding to
        an atom_type that can also be used as an atom name in
        the input of a QM code. Typically such a name would be
        <element><number>.
        """
        self.atom_label = atom_label

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
        base_atom.__init__(self,atmnam,atmnum,resnam,resnum)
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
        
class top_atom(base_atom):
    """
    A class for atoms from a topology file.
    In addition to a base atom this atom also has an atom type, 
    bonds, Lennard-Jones potentials, mass, and partial charges.
    """
    def __init__(self,atmnum):
        """
        Initialize a topology atom. Because of the structure
        of the topology file we cannot set all properties at
        once. In fact we initially just have the number of the
        atom in the structure.
        """
        base_atom.__init__(self,None,atmnum,None,None)
        self.atom_type = None  # MD atom type
        self.atom_type_num = None  # MD atom type number
        self.bound_atoms = []
        self.lj_c6 = None
        self.lj_c12 = None
        self.mass = None
        self.charge = None

    def set_atom_type(self,atom_type):
        """
        Set the atom type (string).
        """
        self.atom_type = atom_type
        #self.set_name(atom_type)

    def set_atom_type_num(self,atom_type):
        """
        Set the atom type number (non-negative integer).
        """
        num = int(atom_type)
        if num < 0:
            raise ValueError(f"invalid atom_type_number {num}")
        self.atom_type_num = num

    def add_bound_atom(self,bound_atom):
        """
        Add the number of a bound atom (integer).
        """
        ii = int(bound_atom)
        if self.number:
            if self.number == ii:
                raise ValueError("Atom cannot be bound to itself")
        if ii in self.bound_atom:
            raise ValueError("Atom cannot bind to the same atom twice")
        self.bound_atoms.append(ii)

    def set_lennard_jones(self,c6,c12):
        """
        Set the C_6 and C_12 coefficients of the Lennard Jones
        potential.
        """
        self.lj_c6 = float(c6)
        self.lj_c12 = float(c12)

    def set_mass(self,mass):
        """
        Set the atom mass.
        """
        self.mass = float(mass)

    def set_charge(self,charge):
        """
        Set the atom charge.
        """
        self.charge = float(charge)


def write_prepare_input(pdb_filename,prepare_input_filename):
    """
    Write the input for the NWChem prepare module to generate the
    topology file. 
    """
    newline = "\n"
    fp = open(prepare_input_filename,"w")
    fp.write("echo\n")
    fp.write("start pdb-gen-aimd-prepare-dat\n")
    fp.write("prepare\n")
    fp.write("  system pdb-gen-aimd-prepare-sys\n")
    fp.write(f"  source {pdb_filename}{newline}")
    fp.write("  new_top new_seq\n")
    fp.write("  new_rst\n")
    fp.write("end\n")
    fp.write("task prepare")
    fp.close()

def read_pdb_file(pdb_filename):
    """
    Read the PDB file and return a list of pdb_atoms,
    and the lattice dimensions.
    Note we need to check whether this PDB file is
    written in NWChem's large PDB format.
    """
    atom_list = []
    lattice   = (0.0,0.0,0.0)
    large_pdb = False
    fp = open(pdb_filename,"r")
    line = fp.readline()
    while line:
        if   line[0:4] == "ATOM":
            atom_list.append(pdb_atom(line,large_pdb))
        elif line[0:6] == "HETATM":
            atom_list.append(pdb_atom(line,large_pdb))
        elif line[0:6] == "LRGPDB":
            large_pdb = True
        elif line[0:6] == "CRYST1":
            lattice = get_lattice(line)
        elif line[0:3] == "END":
            break
        line = fp.readline()
    fp.close()
    return (atom_list,lattice)

def read_top_file(top_filename):
    """
    Read the topology file and return:
    - a list of top_atoms
    - a list of bond distance parameters
    - a list of bond angle parameters
    - a list of torsion angle parameters
    All this data will be returned as a tuple.
    """
    nwl = "\n"
    charge_list = []
    atom_types = []
    final_atom_types = []
    element_count = {}
    atom_labels = {}
    fp = open(top_filename,"r")
    for ii in range(4):
        next(fp)
    ipnum = int(fp.readline())
    nats  = int(fp.readline()) # Number of atom types
    nqu   = int(fp.readline()) # Number of solute atoms
    nseq  = int(fp.readline()) # Number of residues
    next(fp)
    #
    # Read atom types
    #
    for kk in range(nats):
        line = fp.readline()
        atmnum = int(line[5:10])
        atmtyp = line[10:14]
        atmmss = float(line[17:30])
        atom = top_atom(atmnum)
        atom.set_atom_type(atmtyp)
        atom.set_mass(atmmss)
        atom_types.append(atom)
    #
    # Read Lennard-Jones coefficients
    #
    numlines = int(nats*(nats+1)/2)
    for kk in range(numlines):
        line = fp.readline()
        ii = int(line[0:5])
        jj = int(line[5:10])
        c6 = float(line[10:22])
        c6_2 = float(line[22:34])
        c12 = float(line[34:46])
        c12_2 = float(line[46:58])
        if ii == jj:
            atom_types[ii-1].set_lennard_jones(c6,c12)
    #
    # Read partial charges
    # - These are partial charges for a kind of atom
    #
    for kk in range(nqu):
        line = fp.readline()
        qq = float(line[5:17])
        charge_list.append(qq)
    #
    # Skip residue info
    #
    for kk in range(nseq):
        line = fp.readline()
    #
    # Solvent number of parameters?
    # - Solvent is water and the structure seems to be modeled
    #   with three bond lengths
    #
    line = fp.readline()
    solvent_num_atom = int(line[0:7])
    solvent_num_bond = int(line[7:14])
    solvent_num_angl = int(line[14:21])
    if solvent_num_angl > 0:
        raise ValueError(f"Unexpected number of solvent angles: {solvent_num_angl}")
    #
    # Solute number of parameters?
    #
    line = fp.readline()
    solute_num_atom = int(line[0:7])
    solute_num_bond = int(line[7:14])
    solute_num_angl = int(line[14:21])
    solute_num_tors = int(line[21:28])
    #
    # Read solvent data
    #
    # - Read solvent atom names
    #
    solvent_atom_list = []
    for kk in range(solvent_num_atom):
        line = fp.readline()
        atom = top_atom(kk+1)
        atom.set_name(line[10:14])
        charge_num = int(line[46:51])-1
        type_num = int(line[41:46])-1
        atom.set_charge(charge_list[charge_num])
        atom.set_atom_type(atom_types[type_num].atom_type)
        atom.set_atom_type_num(type_num)
        solvent_atom_list.append(atom)
    #
    # - Read bond parameters (bond length and force constant)
    #
    solvent_bond_parameters = []
    for kk in range(solvent_num_bond):
        line = fp.readline()
        iatm = int(line[0:7])
        jatm = int(line[7:14])
        line = fp.readline()
        r_eq = float(line[0:12])
        fc   = float(line[12:24])
        solvent_bond_parameters.append((iatm,jatm,r_eq,fc))
    #
    # Skip next 2 line
    #
    next(fp)
    next(fp)
    #
    # Read solute data
    #
    # - Read solute atom names
    #
    solute_atom_list = []
    for kk in range(solute_num_atom):
        line = fp.readline()
        atom = top_atom(kk+1)
        atom.set_name(line[10:14])
        charge_num = int(line[53:58])-1
        type_num = int(line[47:52])-1
        atom.set_charge(charge_list[charge_num])
        atom.set_atom_type(atom_types[type_num].atom_type)
        atom.set_atom_type_num(type_num)
        solute_atom_list.append(atom)
    #
    # - Read bond parameters (bond length and force constant)
    #
    solute_bond_parameters = []
    for kk in range(solute_num_bond):
        line = fp.readline()
        iatm = int(line[0:7])
        jatm = int(line[7:14])
        line = fp.readline()
        r_eq = float(line[0:12])
        fc   = float(line[12:24])
        solute_bond_parameters.append((iatm,jatm,r_eq,fc))
    #
    # - Read bond angle parameters (bond angle and force constant)
    #
    solute_angle_parameters = []
    for kk in range(solute_num_angl):
        line = fp.readline()
        iatm = int(line[0:7])
        jatm = int(line[7:14])
        katm = int(line[14:21])
        line = fp.readline()
        a_eq = float(line[0:10])
        fc   = float(line[10:22])
        solute_angle_parameters.append((iatm,jatm,katm,a_eq,fc))
    #
    # - Read dihedral angle parameters (angle and force constant)
    #
    solute_torsion_parameters = []
    for kk in range(solute_num_tors):
        line = fp.readline()
        iatm = int(line[0:7])
        jatm = int(line[7:14])
        katm = int(line[14:21])
        latm = int(line[21:28])
        line = fp.readline()
        kdih = int(line[0:3])
        d_eq = float(line[3:13])
        fc   = float(line[13:25])
        solute_torsion_parameters.append((iatm,jatm,katm,latm,kdih,d_eq,fc))
    #
    fp.close()
    #
    # Patchwork: We have a list of atom names as they appear in the PDB
    # file. We have a list of charges. We have a list of atom types.
    # In a QM code we need a unique name for every center that has a
    # unique combination of atom type and charge. (Associating Lennard-
    # Jones parameters and charges goes by atom label). In the MD code
    # atom types and charges are separately managed properties.
    # Hence at this point we need to create a unique atom label for
    # every unique atom type - charge combination. Down stream these
    # labels need to replace the atom names in the data from the PDB
    # file.
    #
    for atom in solute_atom_list:
        element = atom.element
        charge = atom.charge
        type = atom.atom_type
        type_num = atom.atom_type_num
        number = atom.number
        name = atom.name
        c6 = atom_types[type_num].lj_c6
        c12 = atom_types[type_num].lj_c12
        atom_type = (type,charge)
        if atom_type in atom_labels:
            atom_label = atom_labels[atom_type]
        else:
            if element in element_count:
                count = element_count[element]
                count += 1
            else: 
                count = 1
            element_count[element] = count
            atom_label = f"{element}{count}"
            atom_labels[atom_type] = atom_label
            new_atom = top_atom(number)
            new_atom.set_atom_type(type)
            new_atom.set_charge(charge)
            new_atom.set_lennard_jones(c6,c12)
            new_atom.set_atom_label(atom_label)
            new_atom.set_name(name)
            final_atom_types.append(new_atom)
        atom.set_atom_label(atom_label)
    for atom in solvent_atom_list:
        element = atom.element
        charge = atom.charge
        type = atom.atom_type
        type_num = atom.atom_type_num
        number = atom.number
        name = atom.name
        c6 = atom_types[type_num].lj_c6
        c12 = atom_types[type_num].lj_c12
        atom_type = (type,charge)
        if atom_type in atom_labels:
            atom_label = atom_labels[atom_type]
        else:
            if element in element_count:
                count = element_count[element]
                count += 1
            else: 
                count = 1
            element_count[element] = count
            atom_label = f"{element}{count}"
            atom_labels[atom_type] = atom_label
            new_atom = top_atom(number)
            new_atom.set_atom_type(type)
            new_atom.set_charge(charge)
            new_atom.set_lennard_jones(c6,c12)
            new_atom.set_atom_label(atom_label)
            new_atom.set_name(name)
            final_atom_types.append(new_atom)
        atom.set_atom_label(atom_label)
    #
    return (final_atom_types,
            solvent_atom_list,
            solvent_bond_parameters,
            solute_atom_list,solute_bond_parameters,
            solute_angle_parameters,solute_torsion_parameters)
    
def label_pdb_atom(solvent_atom_list,solute_atom_list,atoms_of_pdb):
    """
    For the input file we need atom labels that can be used to
    associate the right charges and Lennard-Jones parameters with
    the right atoms. The names in the PDB file do not allow for that
    mapping. Hence we use the labels that were added to the
    the solvent_atom_list and solute_atom_list as obtained from the
    topology file.
    Return the update atoms_of_pdb list.
    """
    num_solute = len(solute_atom_list)
    num_pdb = len(atoms_of_pdb)
    num_solvent = num_pdb - num_solute
    num_water = int(num_solvent/3)
    if num_water*3 != num_solvent:
        raise ValueError(f"Number of atoms mismatch: {num_pdb} != {num_solute}+{num_water}*3")
    for ii in range(num_solute):
        solute_name = solute_atom_list[ii].name
        pdb_name = atoms_of_pdb[ii].name
        if pdb_name != solute_name:
            raise ValueError(f"PDB - solute atoms: name mismatch: *{pdb_name}-{solute_name}*")
        atoms_of_pdb[ii].set_atom_label(solute_atom_list[ii].atom_label)
    for ii in range(num_solute,num_pdb):
        solvent_name = solvent_atom_list[ii%3].name
        pdb_name = atoms_of_pdb[ii].name
        if pdb_name != solvent_name:
            raise ValueError(f"PDB - solvent atoms: name mismatch: *{pdb_name}-{solvent_name}*")
        atoms_of_pdb[ii].set_label(solvent_atom_list[ii%3].atom_label)
    return atoms_of_pdb

def write_structure(fileptr,pdb_atoms):
    """
    Write the PDB atoms to the file given by the file pointer
    (which is the NWChem PSPW QMMM input file).
    """
    newline = "\n"
    for atom in pdb_atoms:
        qm_flag = " "
        if not atom.qm:
            qm_flag = "^"
        label = atom.atom_label+qm_flag
        coords = atom.coordinates
        line = f"  {label:<10s} {coords}{newline}"
        fileptr.write(line)

def write_qmmm_input(qmmm_input_filename,atoms_pdb,lattice,
                     atom_types,
                     solvent_atom_list,
                     solvent_bond_parameters,
                     solute_atom_list,
                     solute_bond_parameters,solute_angle_parameters,
                     solute_torsion_parameters):
    """
    Write the QMMM input file. This input file can be thought
    of as consisting of multiple sections:
    - Preamble
    - Geometry
    - PSPW input
    - PSPW QMMM block
    - Task
    """
    (aa,bb,cc) = lattice
    nwl = "\n"
    fp = open(qmmm_input_filename,"w")
    fp.write("echo\n")
    fp.write(f"title {qmmm_input_filename}\n")
    fp.write(f"start {qmmm_input_filename}-dat\n") 
    fp.write("geometry units angstroms center autosym autoz print xyz\n")
    write_structure(fp,atoms_pdb)
    fp.write("end\n")
    fp.write("pswp\n")
    fp.write("  qmmm\n")
    for atom in atom_types:
        label = atom.atom_label
        c6 = atom.lj_c6
        c12 = atom.lj_c12
        fp.write(f"    lj_mm_parameters {label:<10s} {c6:12.6e} {c12:12.6e}{nwl}")
    for atom in atom_types:
        label = atom.atom_label
        charge = atom.charge
        fp.write(f"    mm_psp {label:<10s} {charge:10.5f} 4 0.55{nwl}")
    fp.write("  end\n")
    fp.write("  simulation_cell units angstrom\n")
    fp.write("    boundry_conditions: aperiodic\n")
    fp.write("    cell_name: cell\n")
    fp.write("    lattice\n")
    fp.write(f"      lat_a {aa}{nwl}")
    fp.write(f"      lat_b {bb}{nwl}")
    fp.write(f"      lat_c {cc}{nwl}")
    fp.write("    end\n")
    fp.write("    ngrid 16 16 16\n")
    fp.write("  end\n")
    fp.write("  wavefunction_initializer\n")
    fp.write("    restricted\n")
    fp.write("    restricted_electrons: nn\n")
    fp.write("    cell_name: cell\n")
    fp.write("  end\n")
    fp.write("  steepest_descent\n")
    fp.write("    cell_name: cell\n")
    fp.write("    time_step: 5\n")
    fp.write("    loop: 10 10\n")
    fp.write("    tolerances: 1.0d-9 1.0d-9 1.0d-4\n")
    fp.write("    energy_cutoff:       15.0d0\n")
    fp.write("    wavefunction_cutoff: 15.0d0\n")
    fp.write("    exchange_correlation: pbe96\n")
    fp.write("  end\n")
    fp.write("  conjugate_gradient\n")
    fp.write("    cell_name: cell\n")
    fp.write("    loop: 1 10\n")
    fp.write("    tolerances: 1.0d-9 1.0d-9 1.0d-4\n")
    fp.write("    energy_cutoff:       15.0d0\n")
    fp.write("    wavefunction_cutoff: 15.0d0\n")
    fp.write("    exchange_correlation: pbe96\n")
    fp.write("  end\n")
    fp.write("end\n")
    fp.write("task pspw wavefunction_initializer\n")
    fp.write("task pspw steepest_descent\n")
    fp.write("task cg_pspw energy\n")
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
    topology_file = "pdb-gen-aimd-prepare-sys.top"
    qmmm_input    = args.output
    write_prepare_input(args.input,prepare_input)
    run_nwchem(args.nwchem_exe,prepare_input)
    (atoms_of_pdb,lattice) = read_pdb_file(args.input)
    (atom_types,
     solvent_atom_list,
     solvent_bond_parameters,
     solute_atom_list,solute_bond_parameters,
     solute_angle_parameters,solute_torsion_parameters) = read_top_file(
        topology_file)
    atoms_of_pdb = label_pdb_atom(solvent_atom_list,solute_atom_list,
        atoms_of_pdb)
    write_qmmm_input(qmmm_input,atoms_of_pdb,lattice,atom_types,
        solvent_atom_list,
        solvent_bond_parameters,
        solute_atom_list,solute_bond_parameters,
        solute_angle_parameters,solute_torsion_parameters)


if __name__ == "__main__":
    execute_with_arguments(complement_arguments(parse_arguments()))
