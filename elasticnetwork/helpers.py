import scipy.spatial
import numpy as np
import csv
import datetime
import urllib
import gzip
import os
import sys
from elasticnetwork.molecules import Atom, AtomList, BondList

def distance_between(list1, list2):
    """ Returns an M * N matrix of pairwise distances between
        two sets of atoms or bonds.  The coordinates of a bond
        is defined as half-way along the line connecting its two atoms.
        
        Parameters
        ----------
        list1 : :proteingraph:AtomList/:proteingraph:BondList/ 
                :proteingraph:Atom
          List of M atoms or bonds
        list2 : :proteingraph:AtomList or :proteingraph:BondList
                :proteingraph:Atom
          List of N atoms or bonds

        Returns
        -------
        pairwise_distances : :scipy:
          M * N scipy matrix of pairwise distances
    """
    def get_coordinates(list_):
        if isinstance(list_, BondList):
            coordinates = np.zeros((len(list_),3))
            for i, bond in enumerate(list_):
                midpoint = (bond.atom1.coordinates() + bond.atom2.coordinates())/2
                coordinates[i,:] = midpoint
        elif isinstance(list_, AtomList):
            coordinates = list_.coordinates()
        elif isinstance(list_, Atom):
            coordinates = np.ndarray((1, 3))
            coordinates[0] = list_.xyz
        else:
            raise TypeError("Input should be AtomLists or BondList")
        return coordinates

    coordinates1 = get_coordinates(list1)
    coordinates2 = get_coordinates(list2)
    pairwise_distances = scipy.spatial.distance.cdist(coordinates1, coordinates2)
    return pairwise_distances

def dist(coords1, coords2):
    """ Finds distance between two sets of coords, where the residue
    is a coarse-grained GNM and so the coordinates are
    the c-alpha atoms

    Parameters
    ----------
    resi1 : :proteingraph:Residue
     First residue
    resi2 : :proteingraph:Residue
     Second residue

    Returns
    -------
    dist : float
      distance between the two residues
    """

    coords1 = np.array(coords1)
    coords2 = np.array(coords2)
    coord_diff = coords1 - coords2
    dist = np.sqrt(sum(coord_diff**2))
    return dist

def find_atoms(atom_list, PDBnum='any', element='any', name='any', res_name='any',
                  chain = 'any', res_num='any', model=None):
    """Find atoms in an AtomList which meet certain requirements:
        TODO: add params list
    """

    atoms_out = AtomList()
    for atom in atom_list:
        if ((PDBnum == 'any' or PDBnum == atom.PDBnum) and
            (name == 'any' or name == atom.name) and
            (element == 'any' or element == atom.element) and
            (res_name == 'any' or atom.res_name) and
            (chain == 'any' or atom.chain) and
            (res_num == 'any' or atom.res_num)
        ):
            atoms_out.append(atom)
    return atoms_out

def atoms_within(cut_off, atom_site_1, atom_site_2):
    """ Find atoms in site 2 which are within a given cut-off distance
    of site 1.
    
    :param cut_off: radius around site 1 to search (in angstroms)
    :type cut_off: float
    :param atom_site_1: site to search around
    :type atom_site_1: AtomList
    :param atom_site_2: atoms to search for
    :type atom_site_2: AtomList

    :returns: sublist of atom_site_2 of all atoms less than the cutoff
    distance away from atom_site_1
    """

    distance_matrix = distance_between(atom_site_1, atom_site_2)
    within_cut_off = np.where(distance_matrix < cut_off)[1]
    within_cut_off = list(set(within_cut_off))
    return atom_site_2[within_cut_off]

def weak_bonds(atoms, bond_list):
    """For a given set of bonds and a given set of atoms, identify
    the sublist of bonds containing (at least one of) those atoms. 

    Parameters
    -----------
    atoms : :proteingraph:AtomList
      List of atoms whose bonds are to be identified
    bond_list: :proteingraph:BondList 
      List of bonds to search in

    Returns
    -------
    bond : :proteingraph:BondList 
      sublist of bonds
    """

    atom_ids = atoms.id()
    # The criterion used here for deciding if the bond is covalent
    # or not is an arbritrary energy cut-off (all covalent bonds
    # so far encountered are well above the threshold of 3000)
    bond_sublist = BondList([bond for bond in bond_list
                            if (bond.atom1.id in atom_ids or bond.atom2.id in atom_ids) 
                            and bond.force_constant < 10e5]) #this can be adjusted
    
    return bond_sublist

def write_atom_id(atom_list, filename):
    """ Write the PDB numbers of a set of atoms to a text file.
    
    Parameters
    ----------
    atom_list : :proteingraph:AtomList 
       AtomList of atoms to write
    filename : str 
      name of file to be written
    """

    with open(filename, 'w') as f:
        for atom in atom_list:
            f.write(str(atom.PDBnum) + '\n')

def write_atom_data(atom_list, atom_data, filename):
    """ Write the PDB numbers of a set of atoms to a text file
    which can be read by the pymol function datab2_atom.

    Parameters
    -----------
    atom_list : :proteingraph:AtomList
      list of atoms to write
    filename : str 
      name of file to be written
    """

    with open(filename, 'w') as f:
        for atom in atom_list:
            f.write(atom.chain + ' ' +
                    atom.res_num + ' ' +
                    atom.res_name + ' ' + 
                    atom.name + ' ' + 
                    str(atom_data[atom.id]) + '\n')


def write_bonds(bond_list, filename, data_vector={}):
    """ Write the bonds in a given list to a text file in the format:
        bond_id, atom1_PDBnumber, atom2_PDBnumber

    Parameters
    ----------
    bond_list : :ProteinGraph:
      List of bonds to be written
    filename : str
      Name of file to be written
    data_vector : dict
      dictionary of (bond id, values)
    """

    if (data_vector and
        not len(bond_list) == len(data_vector)):
        raise ValueError("Data vector and list of bonds must have same length")

    with open(filename, 'w') as f:
        wtr = csv.writer(f)
        for i, bond in enumerate(bond_list):
            if data_vector:
                wtr.writerow([bond.id, bond.atom1.PDBnum, 
                                bond.atom2.PDBnum, data_vector[bond.id]])
            else:
                wtr.writerow([bond.id, bond.atom1.PDBnum, 
                                bond.atom2.PDBnum])


def write_bonds_gnm(gnm_protein, file_name):
    
    with open(file_name, 'w') as f:
        for i, bond in enumerate(gnm_protein.graph.edges()):
            resi1 = gnm_protein.residues[bond[0]]
            resi2 = gnm_protein.residues[bond[1]]
            out = "dist bond{0}, resi {1} and chain {2} and name ca, resi {3} and chain {4} and name ca".format(i, resi1.res_num, resi1.chain, resi2.res_num, resi2.chain)
            f.write(out + '\n')

def write_results(output_name, pdb_id,
                  source_residues, target_residues, 
                  nTarget_atoms, nTarget_bonds, nSamples,
                  target_atom_hits, target_bond_hits,
                  atom_p_values, bond_p_values):
    """Write the results from a combined Markov transient and edge-edge analysis
    to a results file.
    """
    
    now = datetime.datetime.now()
    datestring = str(now.day) + '/' + str(now.month) + '/' + str(now.year)
    header_length = len(pdb_id) + 22
    with open(output_name + '.results', 'w') as f:
        # Write header        
        f.write('#' * header_length + '\n')
        f.write('#' + ' ' * (header_length-2) + '#\n')
        f.write('#' + ' ' * 5 + pdb_id + ' ' + datestring + ' ' * 5 + '#\n')
        f.write('#' + ' ' * (header_length-2) + '#\n')
        f.write('#' * header_length + '\n\n')

        # Write source residues
        f.write('Source residues\n')
        f.write('#' * 15 + '\n')
        f.write(str(source_residues) + '\n\n')

        # Write target residues
        f.write('Target residues\n')
        f.write('#' * 15 + '\n')
        f.write(str(target_residues) + '\n\n')

        # Write transient results
        f.write('Transient results\n')
        f.write('#' * 17 + '\n')
        for cut_off, hits in sorted(target_atom_hits.iteritems()):
            f.write('Target atom hits @ {0}: {1}/{2}.  P = {3}\n'
                    .format(cut_off, hits, nTarget_atoms, 
                            atom_p_values[cut_off]))
        f.write('\n')

        # Write perturbation results
        f.write('Perturbation results\n')
        f.write('#' * 20 + '\n')
        for cut_off, hits in sorted(target_bond_hits.iteritems()):
            f.write('Target bond hits @ {0}: {1}/{2}.  P = {3}\n'
                    .format(cut_off, hits, nTarget_bonds, 
                            bond_p_values[cut_off]))
        f.write('\n')

        f.write('P-values estimated from {0} samples'.format(str(nSamples)))

    return None
 
def retrieve_pdb(pdb_id, target_folder=None, bio=True):
    """ Pulls a pdb file with the specified PDB ID from
    the RCSB protein data bank.

    Parameters
    -----------
    pdb_id : str 
      4 letter PDB identifier
    bio : bool
      If True, biological assembly file will be retrieved
    """
    if target_folder:
        target_folder = target_folder + '/'
    else:
        target_folder = ''

    assert len(pdb_id) == 4, 'PDB ID should be 4 characters'
    if bio:
        URL = "http://www.rcsb.org/pdb/files/" + pdb_id + ".pdb1.gz"
        urllib.urlretrieve(URL, target_folder + pdb_id + ".pdb1.gz")
        in_f = gzip.open(target_folder + pdb_id + ".pdb1.gz")
        out_f = open(target_folder + pdb_id + "_bio.pdb",'wb')
        out_f.write(in_f.read())
        in_f.close()
        out_f.close()
        os.remove(target_folder + pdb_id + ".pdb1.gz")
    else:
        URL = "http://www.rcsb.org/pdb/files/" + pdb_id + ".pdb"
        urllib.urlretrieve(URL, target_folder + pdb_id + ".pdb")
    return None

