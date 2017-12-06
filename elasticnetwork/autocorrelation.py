import numpy as np
import scipy
import scipy.sparse
import scipy.linalg
import pandas as pd

class Autocorrelation_run(object):
    """ Class which runs and holds the results of a 3D autocorrelation
    analysis of a protein.
    
    Parameters
    ----------
    protein : Protein
      Protein object to be analysed
    angles:
      Specifies whether angles are included in the stiffness matrix
    dihedrals:
      Specifies whether dihedrals are included in the stiffness matrix 

    """

    def __init__(self, protein, angles=True, dihedrals=True):
        """
       
        """

        self.protein = protein
        self.angles = angles
        self.dihedrals = dihedrals

        self.bond_results = []
        self.angle_results = []

    
    def calculate_autocorrelation(self, bonds_out=True, angles_out=False):
        """
        bonds_out: calculate autocorrelation for bonds
        angles_out: calculate autocorrelation for angles
        """

        natoms = len(self.protein.atoms) 
        nbonds = len(self.protein.bonds)
        nangles = len(self.protein.angles)
        ndihedrals = len(self.protein.dihedrals)

        #two-centre (bond) adjacency matrix
        A_bonds, G_bonds = generate_bond_matrices(self.protein)
        K_bonds = (A_bonds.T).dot(G_bonds).dot(A_bonds)

        #three-centre (angle) adjacency matrix
        if self.angles:
            A_angles, G_angles = generate_angle_matrices(self.protein)
            K_angles = (A_angles.T).dot(G_angles).dot(A_angles)
        else:
            #just set angle stiffness matrix to zero
            K_angles = scipy.sparse.csr_matrix((3*natoms, 3*natoms))
        
        #four-centre (dihedral) adjacency matrix
        if self.dihedrals:
            A_dihedrals, G_dihedrals = generate_dihedral_matrices(self.protein)
            K_dihedrals = (A_dihedrals.T).dot(G_dihedrals).dot(A_dihedrals)
        else:
            #just set dihedral stiffness matrix to zero
            K_dihedrals = scipy.sparse.csr_matrix((3*natoms, 3*natoms))

        K_total = K_bonds + K_angles + K_dihedrals
    
        print('Calculating pseudoinverse of the stiffness matrix...')
        K_pinv = np.linalg.pinv(K_total.toarray())
        print('..done.')


        if bonds_out:

            #'trick' for getting the diagonal elements only of the autocorrelation matrix
            autocorrelation = np.multiply(((A_bonds.todense()).dot(K_pinv)), A_bonds.todense()).sum(1)
            autocorrelation = [element[0] for element in autocorrelation.tolist()]

            bond_ids = [bond.id for bond in self.protein.bonds]

            bond_names = [make_bond_name(bond) for bond in self.protein.bonds]
            atom1_pdbnum = [bond.atom1.PDBnum for bond in self.protein.bonds]
            atom2_pdbnum = [bond.atom2.PDBnum for bond in self.protein.bonds]

            bond_results_frame = pd.DataFrame(
                {
                'bond_name' : bond_names,
                'atom1_name': atom1_pdbnum,
                'atom2_name': atom2_pdbnum, 
                'autocorrelation' : autocorrelation,
                }, 
                index=bond_ids, 
            )
            self.bond_results = bond_results_frame

        if angles_out:
            
            #'trick' for getting the diagonal elements only of the autocorrelation matrix
            autocorrelation = result = np.multiply(((A_angles.todense()).dot(K_pinv)), A_angles.todense()).sum(1)
            autocorrelation = [element[0] for element in autocorrelation.tolist()]            
            
            angle_ids = [angle.id for angle in self.protein.angles]

            angle_names = [make_angle_name(angle) for angle in self.protein.angles]
            atom1_pdbnum = [angle.atom1.PDBnum for angle in self.protein.angles]
            atom2_pdbnum = [angle.atom2.PDBnum for angle in self.protein.angles]
            atom3_pdbnum = [angle.atom3.PDBnum for angle in self.protein.angles]
            

            angle_results_frame = pd.DataFrame(
                {
                'angle_name' : angle_names,
                'atom1_name': atom1_pdbnum,
                'atom2_name': atom2_pdbnum, 
                'atom3_name': atom3_pdbnum,
                'autocorrelation' : autocorrelation,
                }, 
                index=angle_ids, 
            )
            self.angle_results = angle_results_frame
    
    


def generate_bond_matrices(protein):
    """
    Returns the (sparse) m by n incidence matrix A for the two-centre
    interaction (i.e. the bonds) and the m by m diagonal matrix of
    force constants.
    """
    natoms = len(protein.atoms)
    nbonds = len(protein.bonds)

    A = np.zeros([nbonds, 3*natoms])
    force_constants = np.zeros(nbonds)
    for bond in protein.bonds:
        
        atom1_id = bond.atom1.id
        atom2_id = bond.atom2.id

        atom1_xyz = bond.atom1.xyz
        atom2_xyz = bond.atom2.xyz

        bond_length = np.linalg.norm(atom1_xyz - atom2_xyz)

        row = A[bond.id]
        row[[3*atom1_id, (3*atom1_id)+1, (3*atom1_id)+2]] = (atom1_xyz - atom2_xyz)/bond_length
        row[[3*atom2_id, (3*atom2_id)+1, (3*atom2_id)+2]] = (atom2_xyz - atom1_xyz)/bond_length

        force_constant = bond.force_constant # this needs to change when force constants figured out
        force_constants[bond.id] = force_constant

    A = scipy.sparse.csr_matrix(A)
    G = scipy.sparse.diags(force_constants) 

    return (A, G)


def generate_angle_matrices(protein):
    """
    Returns the (sparse) m by n incidence matrix for the three-centre
    interaction (i.e. the angles) and the m by m diagonal matrix of
    force constants.
    """

    #double check maths for this to be safe (particularly signs)

    natoms = len(protein.atoms)
    nangles = len(protein.angles)

    A = np.zeros([nangles, 3*natoms])
    force_constants = np.zeros(nangles)
    for angle in protein.angles:

        atom1_id = angle.atom1.id
        atom2_id = angle.atom2.id
        atom3_id = angle.atom3.id

        atom1_xyz = angle.atom1.xyz
        atom2_xyz = angle.atom2.xyz
        atom3_xyz = angle.atom3.xyz

        three_centre_length = np.linalg.norm(atom1_xyz - atom3_xyz)

        row = A[angle.id]
        row[[3*atom1_id, (3*atom1_id)+1, (3*atom1_id)+2]] = (atom2_xyz - atom3_xyz)/three_centre_length
        row[[3*atom2_id, (3*atom2_id)+1, (3*atom2_id)+2]] = -((atom2_xyz - atom1_xyz) + (atom2_xyz - atom3_xyz))/three_centre_length
        row[[3*atom3_id, (3*atom3_id)+1, (3*atom3_id)+2]] = (atom2_xyz - atom1_xyz)/three_centre_length

        force_constant = angle.force_constant
        force_constants[angle.id] = force_constant
    
    A = scipy.sparse.csr_matrix(A)
    G = scipy.sparse.diags(force_constants)

    return (A,G)

def generate_dihedral_matrices(protein):
    """
    Returns the (sparse) m by n incidence matrix for the four-centre
    interaction (i.e. the dihedrals) and the m by m diagonal matrix of
    force constants.
    """

    #double check maths for this to be safe (particularly signs)

    natoms = len(protein.atoms)
    ndihedrals = len(protein.dihedrals)

    A = np.zeros([ndihedrals, 3*natoms])
    force_constants = np.zeros(ndihedrals)
    for dihedral in protein.dihedrals:
        
        atom1_id = dihedral.atom1.id
        atom2_id = dihedral.atom2.id
        atom3_id = dihedral.atom3.id
        atom4_id = dihedral.atom4.id

        atom1_xyz = dihedral.atom1.xyz
        atom2_xyz = dihedral.atom2.xyz
        atom3_xyz = dihedral.atom3.xyz
        atom4_xyz = dihedral.atom4.xyz

        four_centre_length = np.linalg.norm(atom1_xyz - atom4_xyz)

        row = A[dihedral.id]
        row[[3*atom1_id, (3*atom1_id)+1, (3*atom1_id)+2]] = -((atom1_xyz - atom3_xyz) + (atom4_xyz - atom2_xyz))/four_centre_length 
        row[[3*atom2_id, (3*atom2_id)+1, (3*atom2_id)+2]] = -((atom2_xyz - atom1_xyz) + (atom2_xyz - atom3_xyz) + (atom2_xyz - atom4_xyz))/four_centre_length
        row[[3*atom3_id, (3*atom3_id)+1, (3*atom3_id)+2]] = -((atom3_xyz - atom4_xyz) + (atom3_xyz - atom1_xyz) + (atom3_xyz - atom2_xyz))/four_centre_length
        row[[3*atom4_id, (3*atom4_id)+1, (3*atom4_id)+2]] = -((atom4_xyz - atom2_xyz) + (atom1_xyz - atom3_xyz))/four_centre_length

        force_constant = dihedral.force_constant
        force_constants[dihedral.id] = force_constant

    A = scipy.sparse.csr_matrix(A)
    G = scipy.sparse.diags(force_constants)

    return (A, G)


def make_bond_name(bond):
    """ Returns a string giving the name of the two atoms in 
    a given bond.
    
    Parameters
    ----------
    bond : :class:'Bond'
      bond object for which string is to be returned
    
    Returns
    -------
    result : string
      Name of the bond in the form 'Atom1 Name : Atom 2 Name'
    
    """
    return (bond.atom1.res_name +
            bond.atom1.res_num + ' ' + 
            bond.atom1.chain + ' ' + 
            bond.atom1.name + ' : ' + 
            bond.atom2.res_name +
            bond.atom2.res_num + ' ' + 
            bond.atom2.chain + ' ' + 
            bond.atom2.name)       


def make_angle_name(angle):
    """ Returns a string giving the name of the two atoms in 
    a given bond.
    
    Parameters
    ----------
    bond : :class:'Bond'
        bond object for which string is to be returned
    
    Returns
    -------
    result : string
        Name of the bond in the form 'Atom1 Name : Atom 2 Name'
    
    """

    return (angle.atom1.res_name +
            angle.atom1.res_num + ' ' + 
            angle.atom1.chain + ' ' + 
            angle.atom1.name + ' : ' + 
            angle.atom2.res_name +
            angle.atom2.res_num + ' ' + 
            angle.atom2.chain + ' ' + 
            angle.atom2.name + ' : ' + 
            angle.atom3.res_name +
            angle.atom3.res_num + ' ' + 
            angle.atom3.chain + ' ' + 
            angle.atom3.name) 
