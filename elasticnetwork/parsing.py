import networkx as nx
import os
from os import path
import urllib
import subprocess
from collections import defaultdict
import gzip
import string
import numpy as np
import csv
import elasticnetwork.settings
import elasticnetwork.molecules
import warnings
import itertools

#####################################
#                                   #
# The following functions deal with #
#      parsing of the PDB file      #
#                                   #
#####################################

class PDBParser(object):
    """Class for loading information on atoms and residues 
    from specified PDB file. 

    Parameters
    -----------
    pdb_filename : str 
      PDB ID or file name
    """

    def __init__(self, pdb_filename):
        self.pdb = pdb_filename
        # filename after adding hydrogens/stripping unwanted atoms etc.
        self.pdb_final = self.pdb         
        
    def parse(self, protein, bio=False, model=None, chain='all', 
                strip='default', strip_ANISOU=True, remove_LINK=False,
              add_H=None):
        """ Takes a proteingraph3 Protein object and loads it with data 
        from the given PDB file.

        Parameters
        ----------
        bio : bool
          Indicates whether the file is a biological pdb file.
          If True, multiple models in the same file will be combined.
        model : str
          Model to be loaded.  Should be specified if there
          are multiple models in the same file (e.g. NMR structures)
        chain : tuple of str
          Tuple of chains to load.  
        strip : tuple of str 
          name of objects to be removed from PDB file
        strip_ANISOU : bool
          specifies whether remove ANISOU entries should be
          removed from PDB file
        remove_LINK : bool 
          specifies whether LINK entries should be removed 
          from the PDB file
        """
    
        if strip_ANISOU:
            print("Removing ANISOU entries")
            self.strip_ANISOU()

        # Symmetric subunits are often stored as separate models in bio files         
        if self.check_models():
            print("Combining models")
            self.combine_models()

        if remove_LINK:
            print("Removing LINK entries")
            self.remove_LINK_entries()

        print("Stripping unwanted atom types from the PDB file")
        if not (strip=='default' or isinstance(strip, dict)):
            raise TypeError("'strip' should be a dict")
        self.strip_atoms(strip)

        self.renumber_atoms()

        if add_H or not self.has_hydrogens():
            print("Adding hydrogens to PDB file")
            self.add_hydrogens()
        
        # Parse individual lines of the pdb file
        print("Loading atoms into protein object")
        (protein.atoms,
         protein.residues,
         protein.chains,
         protein.pdbNum_to_atomID) = self.parse_pdb_lines(model, chain)
        
        protein.pdb_id = self.pdb_final


    def strip_ANISOU(self):
        """Strips ANISOU records from the PDB file by calling 
        grep and cat bash commands

        """
        subprocess.call("grep -v ANISOU " + self.pdb_final + " > " + 
                        self.pdb_final[0:-4] + "_temp.pdb",
                        shell = True)
        subprocess.call("cat " + self.pdb_final[0:-4] + "_temp.pdb > " + 
                        self.pdb_final, 
                        shell = True) 
        os.remove(self.pdb_final[0:-4] + "_temp.pdb")

    def remove_LINK_entries(self):
        subprocess.call("grep -v LINK " + self.pdb_final + " > " + 
                        self.pdb_final[0:-4] + "_temp.pdb", shell = True)
        subprocess.call("cat " + self.pdb[0:-4] + "_temp.pdb > " + 
                        self.pdb_final, shell = True)
        os.remove(self.pdb_final[0:-4] + "_temp.pdb")

    def strip_atoms(self, strip):
        """Creates a new PDB file with atoms specified in the dictionary
        'strip' removed.
        
        Parameters
        ----------
        strip : dict
          Dictionary 
        """
        aa_to_eliminate = elasticnetwork.settings.aa_to_eliminate

        if strip == 'default':
            strip = {}
            strip['res_name'] = aa_to_eliminate
        else:
            strip['res_name'] = strip.get('res_name', []) + aa_to_eliminate
            
        old_f = open(self.pdb_final, 'r')
        new_f = open(self.pdb_final[0:-4] + '_stripped.pdb', 'w')
        print (strip['res_name'])
        for line in old_f:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                atom = parse_atom_line(line)
                if (
                        atom.PDBnum not in strip.get('PDBnum', []) and
                        atom.name not in strip.get('name', []) and
                        atom.element not in strip.get('element', []) and
                        atom.chain not in strip.get('chain', []) and
                        atom.res_num not in strip.get('res_num', []) and
                        atom.res_name not in strip.get('res_name', []) and
                        [atom.res_num, atom.chain] not in strip.get('residues', [])
                ):
                    new_f.write(line)
            elif line.startswith('CONECT'):
                continue
            else:
                new_f.write(line)

        self.pdb_final = self.pdb_final[0:-4] + '_stripped.pdb'

    
    def renumber_atoms(self):
        subprocess.call('python3 ' + elasticnetwork.settings.DEPENDENCIES_ROOT +
                        '/atom_renumbering3.py ' + self.pdb_final,
                        shell=True)



    def has_hydrogens(self):
        """Search PDB file for hydrogen atoms. If one is found the 
        file is assumed to have had hydrogens previously added.
        """

        try:
            pdbf = open(self.pdb_final,'r')
        except:
            raise IOError("Couldn't open PDB file " + self.pdb_final + 
                            "to check for hydrogens")

        for line in pdbf:
            if line[0:4] == 'ATOM' and line[13] == 'H':
                return True
        return False


    def add_hydrogens(self):
        """ Runs the command-line program Reduce to add hydrogens
        to the specified pdb file
        """

        subprocess.call(elasticnetwork.settings.REDUCE + ' -Quiet -BUILD -DB ' +
                        elasticnetwork.settings.FIRST_ROOT + '/lib/het_dict.txt ' 
                        + self.pdb_final + ' > ' + self.pdb_final[0:-4] 
                        + '_H_temp.pdb', shell=True)
        
        subprocess.call('python3 ' + elasticnetwork.settings.DEPENDENCIES_ROOT +
                        '/atom_renumbering3.py ' + self.pdb_final[0:-4] + 
                        '_H_temp.pdb',
                        shell=True)
        
        subprocess.call('mv ' + self.pdb_final[0:-4] + '_H_temp.pdb ' + self.pdb_final[0:-4] + '_H.pdb', 
                        shell=True)

        
        self.pdb_final = self.pdb_final[0:-4] + '_H.pdb'
        

    def parse_pdb_lines(self, model, chain):
        """Parses the details of the atoms from a pdb file.

        Returns
        -------
        pdb_data : tuple
          tuple of data for the atoms in the protein containing
            1. AtomList of the atoms in the PDB file 
            2. Dict mapping residues to lists of atoms
            3. Dict mapping chains to atom IDs
            4. Dict mapping PDB numbers to unique atom ids
        """

        # Containers for data from PDB files
        residues = defaultdict(list)
        chains = defaultdict(list)
        pdbNum_to_atomID = defaultdict(list)
                        
        #TO DO : use a generator here
        atoms = load_atoms(self.pdb_final, model=model, chain=chain)
        print (atoms)
        for atom in atoms:
            residues[(atom.res_num, atom.chain)].append(atom.id)
            # Add atom id to the chain -> atom map
            chains[atom.chain].append(atom.id)
            # Map between atom number in PDB file and its atom ID in protein
            pdbNum_to_atomID[atom.PDBnum] = atom.id

        # Turn the defaultdict automatic adding of keys off, so
        # strange things do not happen when the user accesses the dicts
        residues.default_factory = None
        chains.default_factory = None
        pdbNum_to_atomID.default_factory = None
        pdb_data = (atoms, residues, chains, pdbNum_to_atomID)
        
        return pdb_data

    def check_models(self):
        """ Checks if there are multiple models in the PDB file
        """
    
        pdbf = open(self.pdb_final, 'r')
        models = False
        for line in pdbf:
            if line.startswith('MODEL'):
                models = True
                break
        return models


    def combine_models(self):
        """ Combine multiple models into a single model. Chains and atoms
        are renamed so that they do not clash.
        """

        pdbf = open(self.pdb_final, 'r')
        lines = pdbf.readlines()
        lines_new = []
        for i, line in enumerate(lines):
            lines_new.append(list(line.rstrip('\r\n')))

        # Find MODEL and ENDMDL lines in file
        MODEL_lines = [(i, line) for i, line
                    in enumerate(lines)
                    if line.startswith('MODEL')]
        ENDMDL_lines = [(i, line) for i, line
                        in enumerate(lines)
                        if line.startswith('ENDMDL')]

        # Check there are the same number of MODEL and ENDMDL lines
        if not (len(MODEL_lines) == len(ENDMDL_lines)):
            raise Exception('The number of MODEL and ENDMDL lines do not match')

        # Get models
        models = [int(x[1][13]) for i, x in enumerate(MODEL_lines)]
        if not (models == range(1, len((models))+1)):
            raise Exception('The model numbers are not consecutive')
        nModels = len(models)

        # Get start and end lines for each model
        lineNumbers = dict((i+1, (MODEL_lines[i][0], ENDMDL_lines[i][0]))
                        for i in range(0, nModels))

        # Check that all entries between the start and end
        # of the MODEL are ATOM/HETATM/TER
        for model in models:
            for line in lines[lineNumbers[model][0]+1:lineNumbers[model][1]-1]:
                if not (line.startswith('ATOM')
                        or line.startswith('HETATM')
                        or line.startswith('TER')):
                    raise Exception('There are non-atom lines inside MODEL ' +
                                    ' environment at line: ' + str(i))

        # Get the chain IDs for each model
        # (these will probably be duplicate for each model)
        chainIDs_old = dict((model, []) for model in models)
        for model in models:
            for line in lines[lineNumbers[model][0]:lineNumbers[model][1]]:
                if line.startswith('ATOM') or line.startswith('HETATM'):
                    chainIDs_old[model].append(line[21])
            chainIDs_old[model] = sorted(list(set(chainIDs_old[model])))

        # Create new unique chain IDs
        chainIDs_new = dict((model, []) for model in models)
        chainCount = 0
        for model in models:
            for chain in chainIDs_old[model]:
                chainIDs_new[model].append(string.ascii_uppercase[chainCount])
                chainCount = chainCount+1

        # Map old chain id to new chain id
        chainIDmapping = dict((el, {}) for el in models)
        for model in models:
            for i, chain in enumerate(chainIDs_old[model]):
                chainIDmapping[model][chain] = chainIDs_new[model][i]

        for model in models:
            for i, line in enumerate(lines_new[lineNumbers[model][0]:lineNumbers[model][1]]):
                if ''.join(line[0:4]) == 'ATOM' or "".join(line[0:6]) == 'HETATM':
                    line[21] = chainIDmapping[model][line[21]]
                    lines_new[lineNumbers[model][0]+i] = line
        for i, line in enumerate(lines_new):
            lines_new[i] = "".join(line)

        i = 0
        for model in models:
            lines_new.pop(lineNumbers[model][0]-i)
            i += 1
            lines_new.pop(lineNumbers[model][1]-i)
            i += 1

        f = open(self.pdb_final, 'w')
        for line in lines_new:
            f.write(line)
        f.close()

        return True


def load_atoms(pdb, PDBnum='all', name='all', res_name='all', chain = 'all',
               res_num='all', model=None):
    """ Load a list of atoms from a pdb file
    
    Returns
    --------
    atoms : :proteingraph3:AtomList
      AtomList of atoms with given PDBnums/name/residue/model
    """
    
    atoms = elasticnetwork.molecules.AtomList()
    currentModel = False
    id_counter = 0
    with open(pdb, 'r') as f:
        for line in f:
            if not currentModel:
                if ((line.startswith('MODEL') and int(line[10:14]) == model)
                    or model is None):
                    currentModel = True
            if line.startswith('ATOM') or line.startswith('HETATM'):
                atom = parse_atom_line(line)
                if (
                    currentModel and 
                    (PDBnum == 'all' or atom.PDBnum in PDBnum) and
                    (name == 'all' or atom.name in name) and
                    (res_name == 'all' or atom.res_name in res_name) and
                    (chain == 'all' or atom.chain in chain) and
                    (res_num == 'all' or atom.res_num in res_num)
                ):
                    atom.id = id_counter
                    atoms.append(atom)
                    id_counter += 1
                    if line.startswith('ENDMDL'):
                        currentModel = False
    return atoms
        
def parse_atom_line(line):
    """ Extract data from a single ATOM line in a pdb file

    Returns
    -------
    atom : :proteingraph3:Atom
      atom object representing the atom in the input line
    """

    element = line[76:78].strip()
    name = line[12:16].strip()
    PDBnum = int(line[6:11])
    chainID = line[21]
    res_num = line[22:27].strip()
    res_name = line[17:20].strip() + line[26].strip()
    try:
        bfactor = float(line[60:66])
    except ValueError:
        bfactor = None

    coordinates = np.array([float(line[30:38]),
                            float(line[38:46]),
                            float(line[46:54])])
                
    # Create Atom object with this information
    atom = elasticnetwork.molecules.Atom(0, PDBnum, element, name,
                                        chainID, res_num, res_name,
                                        bfactor, coordinates)
    return atom


def convert_atom_name_to_standard(atom_name):
        """ PDB files have alternative ways of naming 
        atoms (Int-ELEMENT-Int or ELEMENT-Int-Int).
        For consistency with the way we store atoms in 
        the energy dictionaries, this function converts
        atom names into Int-ELEMENT-Int format.
        """

        if atom_name[-2:].isdigit():
            atom_name = atom_name[-1:] + atom_name[:-1]
        return atom_name




"""
Can keep all of the code above as it is.  Loads all the atom info into 
Protein() object
"""


#####################################
#                                   #
# The following functions deal with #
# generation of the protein network #
# including angles and dihedrals if #
# specified.                        #
#                                   #
#####################################

class FIRST_Network_Generator(object):

    def __init__(self):
        self.bonds = elasticnetwork.molecules.BondList()
        #added here
        self.angles = elasticnetwork.molecules.AngleList()
        self.dihedrals = elasticnetwork.molecules.DihedralList()
        #
        self.atom_list = elasticnetwork.molecules.AtomList()
        self.pdbNum_to_atomID = {}
        self.pdb = ''

    def generate_elastic_network(self, protein, hbond_threshold=-0.01, bond_files_path=None, remove_data_files=True, hphobes=True,  
                              angles=True, dihedrals=True):
        """ Adds list of bonds (central forces) to construct the given protein network
            and angle + dihedral constraints if specified

            Maybe add van der Waals too?

        Parameters
        ----------
        pdb : str
          4 digit PDB ID code
        hbond_threshold : float
          Energy cut-off for inclusion of hydrogen bonds in graph. 
          Default is -0.01 kcal/mol.
        atom_list : :proteingraph3:AtomList
          List of atoms in protein
        pdbNum_to_atomID : dict
          map between atom PDB number and Atom object id

        """

        self.pdb = protein.pdb_id
        self.atom_list = protein.atoms
        self.pdbNum_to_atomID = protein.pdbNum_to_atomID
        if bond_files_path is None:
            self.bond_files_path = os.path.dirname(os.path.abspath(protein.pdb_id))
        else:
            self.bond_files_path = bond_files_path

    
        # Call FIRST to generate info about covalent/hydrogen/hydrophobic bonds
        subprocess.call(elasticnetwork.settings.FIRST_ROOT + '/FIRST ' + self.pdb +
                        ' -c 2.0 -covout -E ' + str(hbond_threshold) +
                        ' -H 3 -hbout -o 3 -phout -v 3 -non -L ' +
                        elasticnetwork.settings.FIRST_ROOT, shell=True)
        subprocess.call('mv ' + self.bond_files_path + '/hphobes.out ' + 
                        self.bond_files_path + '/hphobes5A.out', shell=True)
        subprocess.call(elasticnetwork.settings.FIRST_ROOT + '/FIRST ' + self.pdb +
                        ' -c 4.5 -covout -E ' + str(hbond_threshold) +
                        ' -H 3 -hbout -o 3 -phout -v 3 -non -L '
                        + elasticnetwork.settings.FIRST_ROOT, shell=True)
        subprocess.call('mv ' + self.bond_files_path + '/hphobes.out ' + 
                        self.bond_files_path + '/hphobes8A.out', shell=True)
        
        

        """
        2-centre interactions (bonds or central forces)
        """

        self.parse_covalent(protein, self.bond_files_path + '/cov.out') 
        #maybe add double/peptide info here?
        self.parse_hbonds(protein, self.bond_files_path + '/hbonds.out', -0.01)
        if hphobes:
            self.parse_hphobes(protein, self.bond_files_path + '/hphobes5A.out', 1)
            self.parse_hphobes(protein, self.bond_files_path + '/hphobes8A.out', 1) #will have to just play around with this
        # self.parse_electrostatics(protein, self.bond_files_path + '/cov.out')
        # self.parse_stacked(protein)

        # Make list of bonds unique in the sense that there is only
        # one bond per pair of atoms (that can include multiple
        # contributions from hydrophobic, hydrogen, electrostatics etc.)
        self.bonds = uniquify(self.bonds)
        protein.bonds = self.bonds

        # Give each bond an ID
        id_counter = 0
        for bond in self.bonds:
            bond.id = id_counter
            id_counter += 1
        
        """
        3-centre interactions (angles)
        """
        if angles:
            self.parse_angles(protein)
            protein.angles = self.angles
            
            # Give each angle an ID
            id_counter = 0
            for angle in self.angles:
                angle.id = id_counter
                id_counter += 1 

            self.add_angle_connected_atoms(protein)

        """
        4-centre interactions (dihedrals)
        """
        if dihedrals:
            self.parse_dihedrals(protein)
            protein.dihedrals = self.dihedrals

            # Give each dihedral an ID
            id_counter = 0
            for dihedral in self.dihedrals:
                dihedral.id = id_counter
                id_counter += 1 
            
            self.add_dihedral_connected_atoms(protein)
        
        
        # Remove files generated by FIRST
        if remove_data_files:
            os.remove(self.pdb[:-4] + '_results.txt')
            os.remove(self.pdb[:-4] + '_bond.txt')
            os.remove(self.pdb[:-4] + '_data.txt')
            os.remove(self.pdb[:-4] + '_RCD.pdb')
            os.remove(self.pdb[:-4] + '_RCD.pml')
            os.remove(self.bond_files_path + '/cov.out')
            os.remove(self.bond_files_path + '/hbonds.out')
            os.remove(self.bond_files_path + '/hphobes5A.out')
            os.remove(self.bond_files_path + '/hphobes8A.out')
     


    def parse_hbonds(self, protein, hbonds_file, threshold):
        """ Load hydrogen bonds from FIRST output file into a BondList.

        Parameters
        ----------
        hbonds_file : str
          FIRST output file with hydrogen bond info
        threshold : float 
          Energy cut-off for inclusion of bond in the graph
        atom_list : AtomList 
          List of atoms in the protein
        pdbNum_to_atomID : dict 
          map from atom PDB number to unique atom ID

        """

        bonds = elasticnetwork.molecules.BondList()
        hbonds_file = open(hbonds_file, 'r')
        for line in csv.reader(hbonds_file, delimiter=' ', skipinitialspace=True):
            if float(line[2]) < threshold:
                atom1_id = self.pdbNum_to_atomID[int(line[0])]
                atom2_id = self.pdbNum_to_atomID[int(line[1])]
                atom1 = self.atom_list[atom1_id]
                atom2 = self.atom_list[atom2_id]
                force_constant = 1 #this may change
                self.bonds.append(
                    elasticnetwork.molecules.Bond([], atom1, atom2, 
                                                  force_constant, 'HYDROGEN', 0)  #added 0 bars here
                    )
                if atom2_id not in protein.atoms[atom1_id].connected_atoms:
                    protein.atoms[atom1_id].connected_atoms[atom2_id] = 0
                else:
                    if protein.atoms[atom1_id].connected_atoms[atom2_id] < 5:
                        protein.atoms[atom1_id].connected_atoms[atom2_id] = 0
            
                if atom1_id not in protein.atoms[atom2_id].connected_atoms:
                    protein.atoms[atom2_id].connected_atoms[atom1_id] = 0
                else:
                    if protein.atoms[atom2_id].connected_atoms[atom1_id] < 5:
                        protein.atoms[atom2_id].connected_atoms[atom1_id] = 0



    def parse_covalent(self, protein, cov_file):
        """ Load covalent bonds from FIRST output file into a BondList.

        Parameters
        ----------
        cov_file : str
          FIRST output file with covalent bond info
        """

        bonds = elasticnetwork.molecules.BondList()
        cov_file = open(cov_file, 'r')
        for line in csv.reader(cov_file, delimiter=' ', skipinitialspace=True):
            atom1_id = self.pdbNum_to_atomID[int(line[0])]
            atom2_id = self.pdbNum_to_atomID[int(line[1])]
            #added here
            bars = int(line[2])
            atom1 = self.atom_list[atom1_id]
            atom2 = self.atom_list[atom2_id]
            # If in different amino acids
            force_constant = 1
            self.bonds.append(elasticnetwork.molecules.Bond([], atom1, atom2, 
                                                            force_constant,
                                                            'COVALENT', bars))
            if atom2_id not in protein.atoms[atom1_id].connected_atoms:
                protein.atoms[atom1_id].connected_atoms[atom2_id] = bars
            else:
                if protein.atoms[atom1_id].connected_atoms[atom2_id] < 5:
                    protein.atoms[atom1_id].connected_atoms[atom2_id] = bars
            
            if atom1_id not in protein.atoms[atom2_id].connected_atoms:
                protein.atoms[atom2_id].connected_atoms[atom1_id] = bars
            else:
                if protein.atoms[atom2_id].connected_atoms[atom1_id] < 5:
                    protein.atoms[atom2_id].connected_atoms[atom1_id] = bars
                
        

    def parse_hphobes(self, protein, hphobes_file, force_constant):
        """ Load hydrophobic bonds from FIRST output file into a BondList.
    
        Parameters
        ----------
        cov_file : str
          FIRST output file with hydrophobic bond info
        bond_strength : float 
          strength of hydrophobic interactions

        """

        hphobes_file = open(hphobes_file, 'r')
        bonds = elasticnetwork.molecules.BondList()
        for line in csv.reader(hphobes_file, delimiter=' ', skipinitialspace=True):
            atom1_id = self.pdbNum_to_atomID[int(line[0])]
            atom2_id = self.pdbNum_to_atomID[int(line[1])]
            atom1 = self.atom_list[atom1_id]
            atom2 = self.atom_list[atom2_id]
            self.bonds.append(elasticnetwork.molecules.Bond([], atom1, atom2, 
                                                            force_constant,
                                                            'HYDROPHOBIC', 0))  #added 0 bars here
            if atom2_id not in protein.atoms[atom1_id].connected_atoms:
                protein.atoms[atom1_id].connected_atoms[atom2_id] = 0
            else:
                if protein.atoms[atom1_id].connected_atoms[atom2_id] < 5:
                    protein.atoms[atom1_id].connected_atoms[atom2_id] = 0
            
            if atom1_id not in protein.atoms[atom2_id].connected_atoms:
                protein.atoms[atom2_id].connected_atoms[atom1_id] = 0
            else:
                if protein.atoms[atom2_id].connected_atoms[atom1_id] < 5:
                    protein.atoms[atom2_id].connected_atoms[atom1_id] = 0     
             
    # #might just get rid of this...
    # def parse_stacked(self, protein):
    #     # Default pi-pi energy is -10 kcal/mol
    #     electrostatic_energy = float(-10)
    #     # Conversion to J/A^2
    #     electrostatic_energy = -120*(electrostatic_energy)*4.184/6.022

    #     stacked = []
    #     FIRST_results_file = self.pdb[0:-4] + '_results.txt'
    #     try:
    #         with open(FIRST_results_file, 'r') as f:
    #             w = False
    #             for line in f:
    #                 if w and line.isspace():
    #                     w = False
    #                     break
    #                 elif w:
    #                     stacked.append(line.strip().split())
    #                 if not w and line.startswith(' Stacked'):
    #                     w = True
    #     except IOError:
    #         print ('No file named ' + FIRST_results_file + ' in current' +
    #                 'directory. Cannot load stacked interactions')

    #     for bond in stacked:
    #         atom1_id = self.pdbNum_to_atomID[int(bond[0])]
    #         atom2_id = self.pdbNum_to_atomID[int(bond[1])]
    #         atom1 = self.atom_list[atom1_id]
    #         atom2 = self.atom_list[atom2_id]
    #         self.bonds.append(proteingraph3.molecules3D.Bond([], atom1, atom2, 
    #                                                         electrostatic_energy, 
    #                                                         'STACKED', 0))  #added 0 bars here
    #         if atom2_id not in protein.atoms[atom1_id].connected_atoms:
    #             protein.atoms[atom1_id].connected_atoms[atom2_id] = 0
    #         else:
    #             if protein.atoms[atom1_id].connected_atoms[atom2_id] < 5:
    #                 protein.atoms[atom1_id].connected_atoms[atom2_id] = 0
            
    #         if atom1_id not in protein.atoms[atom2_id].connected_atoms:
    #             protein.atoms[atom2_id].connected_atoms[atom1_id] = 0
    #         else:
    #             if protein.atoms[atom2_id].connected_atoms[atom1_id] < 5:
    #                 protein.atoms[atom2_id].connected_atoms[atom1_id] = 0 

    #and this...
    # def parse_electrostatics(self, protein, cov_file):
    #     """ Load electrostatic bonds using LINK entries in the PDB file.
    
    #     cov_file : str
    #       FIRST output file with covalent bond info

    #     """
    
    #     # Get atom numbers and bond length of LINK entries in PDB file
    #     LINK_entries = self.parse_LINK_entries()
    #     LINK_atoms_PDBnum = LINK_entries[0]
    #     LINK_bond_distance = LINK_entries[1]
    
    #     f = open(cov_file, 'r')
    #     cov_bonds = []
    
    #     # Load covalent bonds from cov.out file to check whether LINK entries
    #     # have already been included
    #     for line in csv.reader(f, delimiter=' ', skipinitialspace=True):   
    #         atom1_id = int(line[0])
    #         atom2_id = int(line[1])
    #         cov_bonds.append((atom1_id, atom2_id))

    #     electro_bonds = proteingraph3.molecules3D.BondList()
    #     # epsilon is the dielectric constant
    #     epsilon = 4
    #     for i in range(0, len(LINK_atoms_PDBnum), 2):
    #         atom1 = self.atom_list[self.pdbNum_to_atomID[LINK_atoms_PDBnum[i]]]
    #         atom2 = self.atom_list[self.pdbNum_to_atomID[LINK_atoms_PDBnum[i+1]]]        
    #         bond_length = LINK_bond_distance[i//2]
    #         # Treat as covalent bond if length is less than 1.8A
    #         if (bond_length < 1.8 and 
    #             (LINK_atoms_PDBnum[i], LINK_atoms_PDBnum[i+1]) not in cov_bonds):
    #             bond_strength = self._covalent_bond_strength(atom1, atom2)
    #             bond = proteingraph3.molecules3D.Bond([], atom1, atom2, 
    #                                                 bond_strength, 'COVALENT')
    #             self.bonds.append(bond)

    #             if atom2_id not in protein.atoms[atom1_id].connected_atoms:
    #                 protein.atoms[atom1_id].connected_atoms[atom2_id] = 0
    #             else:
    #                 if protein.atoms[atom1_id].connected_atoms[atom2_id] < 5:
    #                     protein.atoms[atom1_id].connected_atoms[atom2_id] = 0
                
    #             if atom1_id not in protein.atoms[atom2_id].connected_atoms:
    #                 protein.atoms[atom2_id].connected_atoms[atom1_id] = 0
    #             else:
    #                 if protein.atoms[atom2_id].connected_atoms[atom1_id] < 5:
    #                     protein.atoms[atom2_id].connected_atoms[atom1_id] = 0 

    #         else:
    #             try:
    #                 # Find charges of the two atoms using charge database 
    #                 # in proteingraph3.energies file
    #                 atom1_name_standard = convert_atom_name_to_standard(atom1.name)
    #                 atom2_name_standard = convert_atom_name_to_standard(atom2.name)
    #                 charge1 = proteingraph3.energies.charges[atom1.res_name][atom1_name_standard]
    #                 charge2 = proteingraph3.energies.charges[atom2.res_name][atom2_name_standard]
    #                 # Calculate bond strength using Coulomb's law
    #                 bond_strength = (1389*charge1*charge2)/(epsilon*bond_length)
    #                 bond_strength *= (-120/6.022)
    #                 if bond_strength > 0:
    #                     bond = proteingraph3.molecules3D.Bond([], atom1, atom2, 
    #                                                         bond_strength, 
    #                                                         'ELECTROSTATIC', 0)  #added 0 bars here
    #                     self.bonds.append(bond)

    #                     if atom2_id not in protein.atoms[atom1_id].connected_atoms:
    #                         protein.atoms[atom1_id].connected_atoms[atom2_id] = 0
    #                     else:
    #                         if protein.atoms[atom1_id].connected_atoms[atom2_id] < 5:
    #                             protein.atoms[atom1_id].connected_atoms[atom2_id] = 0
                        
    #                     if atom1_id not in protein.atoms[atom2_id].connected_atoms:
    #                         protein.atoms[atom2_id].connected_atoms[atom1_id] = 0
    #                     else:
    #                         if protein.atoms[atom2_id].connected_atoms[atom1_id] < 5:
    #                             protein.atoms[atom2_id].connected_atoms[atom1_id] = 0 

    #             except(KeyError):
    #                 KeyError("Charge of atom {0} from residue {1} ".format(atom1.name, 
    #                                                                     atom1.res_name) + 
    #                         "or atom {0} from residue {1} ".format(atom2.name, 
    #                                                                 atom2.res_name) + 
    #                         "could not be found in the charge database, please update")

    # def parse_LINK_entries(self):
    #     """Find the PDB numbers of the atoms in the LINK entries of
    #     a pdb file and the bond length of the specified interactions
    #     """

    #     pdb = open(self.pdb, 'r')
    #     lines = pdb.readlines()
    #     LINK_bonds = []
    #     LINK_bond_distance = []
    #     num_link_entries = 0
    #     for line in lines:
    #         if line.startswith('LINK'):
    #             atom1 = {'name':line[12:16].strip(), 
    #                      'res_name': line[17:20].strip(), 
    #                      'res_num': line[22:26].strip(),
    #                      'chain': line[21]}
    #             atom2 = {'name': line[42:46].strip(), 
    #                      'res_name': line[47:50].strip(), 
    #                      'res_num': line[52:56].strip(),
    #                      'chain': line[51]}
    #             LINK_bonds.append((atom1, atom2))
    #             LINK_bond_distance.append(float(line[74:78]))
        
    #     LINK_atoms_PDBnum = []
    #     for LINK_bond in LINK_bonds:
    #         found_atoms = [None, None]
    #         for i in range(2):
    #             found = False
    #             LINK_atom = LINK_bond[i]
    #             for line in lines:
    #                 if (line.startswith('ATOM') or line.startswith('HETATM')):
    #                     atom = parse_atom_line(line)
    #                     if ( LINK_atom['name'] == atom.name and
    #                          LINK_atom['res_name'] == atom.res_name and
    #                          LINK_atom['res_num'] == atom.res_num and
    #                          LINK_atom['chain'] == atom.chain ):
    #                         found = True
    #                         found_atoms[i] = atom
    #                         break
    #             if not found:
    #                 warnings.warn(
    #                     ('LINK atom {0} {1} {2} {3} not found'.format(LINK_atom['name'], 
    #                                                                   LINK_atom['res_name'],
    #                                                                   LINK_atom['res_num'],
    #                                                                   LINK_atom['chain'])
    #                      + ' ignoring this LINK entry.')
    #                 )
    #                 found = False
    #         if found_atoms[0] and found_atoms[1]:
    #             LINK_atoms_PDBnum.extend([found_atoms[0].PDBnum, found_atoms[1].PDBnum])

    #     assert len(LINK_atoms_PDBnum) % 2 == 0, "Odd number of LINK atoms"

    #     return (LINK_atoms_PDBnum, LINK_bond_distance)

    
    
    """Angles"""

    def parse_angles(self, protein):
        """ Load angle interactions from the Protein object into an AngleList.
        Loop over each Atom object and identify the pairs of connected atoms that 
        have either 5 or 6 bars (single or double covalent bond)"""


        for atom in protein.atoms:
            covalently_connected_atoms = []
            for connected_atom, bars in atom.connected_atoms.items():
                if bars >= 5:
                    covalently_connected_atoms.append(connected_atom)
            
            connected_atom_pairs = [comb for comb in itertools.combinations(covalently_connected_atoms,2)]

            angle_interactions = []
            for pair in connected_atom_pairs:
                angle = elasticnetwork.molecules.Angle([], protein.atoms[pair[0]], protein.atoms[atom.id], protein.atoms[pair[1]], 1)  #value subject to change...
                self.angles.append(angle)

            
    def add_angle_connected_atoms(self, protein):
        """This has to be run after parse_angles, maybe put a check in here"""

        for angle in protein.angles:
            atom1, atom3 = angle.atom1, angle.atom3
            
            if atom3.id not in protein.atoms[atom1.id].connected_atoms:
                protein.atoms[atom1.id].connected_atoms[atom3.id] = 0
            else:
                if protein.atoms[atom1.id].connected_atoms[atom3.id] < 5:
                    protein.atoms[atom1.id].connected_atoms[atom3.id] = 0
            
            if atom1.id not in protein.atoms[atom3.id].connected_atoms:
                protein.atoms[atom3.id].connected_atoms[atom1.id] = 0
            else:
                if protein.atoms[atom3.id].connected_atoms[atom1.id] < 5:
                    protein.atoms[atom3.id].connected_atoms[atom1.id] = 0 


    
    """Dihedrals"""

    def parse_dihedrals(self, protein):
        """ Load diehdral interactions from the Protein object into a DihedralList.
        Loop over each Bond object and identify the double and peptide bonds i.e. 
        those that have 6 bars in the FIRST formalism """

        #May need to think about how this currently handles aromatic rings
        #each bond in the ring is given 6 bars

        double_and_peptide = []        
        for bond in protein.bonds:
            if bond.bars == 6:
                double_and_peptide.append(bond)
        
        for bond in double_and_peptide:
            atom1 = bond.atom1
            atom2 = bond.atom2

            atom1_connected = list(atom1.connected_atoms.keys())
            atom1_connected.remove(atom2.id)

            atom2_connected = list(atom2.connected_atoms.keys())
            atom2_connected.remove(atom1.id)
            

            joined_pairs = list(itertools.product(atom1_connected, atom2_connected))
            for pair in joined_pairs:
                dihedral = elasticnetwork.molecules.Dihedral([], protein.atoms[pair[0]], atom1, protein.atoms[pair[1]], atom2, 1) #value subject to change...
                self.dihedrals.append(dihedral)



    def add_dihedral_connected_atoms(self, protein):
        """This has to be run after parse_dihedrals, maybe put a check in here"""

        for dihedral in protein.dihedrals:
            atom1, atom4 = dihedral.atom1, dihedral.atom4
            
            if atom4.id not in protein.atoms[atom1.id].connected_atoms:
                protein.atoms[atom1.id].connected_atoms[atom4.id] = 0
            else:
                if protein.atoms[atom1.id].connected_atoms[atom4.id] < 5:
                    protein.atoms[atom1.id].connected_atoms[atom4.id] = 0
            
            if atom1.id not in protein.atoms[atom4.id].connected_atoms:
                protein.atoms[atom4.id].connected_atoms[atom1.id] = 0
            else:
                if protein.atoms[atom4.id].connected_atoms[atom1.id] < 5:
                    protein.atoms[atom4.id].connected_atoms[atom1.id] = 0 



"""Helper functions """

def uniquify(bond_list):
    """Takes a list of bonds which may contain multiple bonds between
    the same two atoms and returns a list where these have been combined 
    into a single bond with the interaction energies summed.

    Parameters
    ----------
    bond_list : BondList 
      List of non-unique bonds

    Returns
    -------
    bonds_unique : BondList 
      List of bonds where separate Bonds between the same two atoms are
      combined into a single Bond with the interaction energies summed
    """

    bonds_unique = elasticnetwork.molecules.BondList()
    count = 0
    seen = {}
    for bond in bond_list:
        try: 
            bonds_unique[seen[bond]].add_interaction(bond.weight, bond.bond_type)
        except:
            seen[bond] = count*1
            bonds_unique.append(bond)
            count += 1
    return bonds_unique
