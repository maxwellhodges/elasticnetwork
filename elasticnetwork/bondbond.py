import numpy as np
import networkx as nx
import scipy.sparse
import pandas as pd
import elasticnetwork.helpers
from elasticnetwork.molecules import BondList
import elasticnetwork.statistics
#import proteingraph.null_models
#import matplotlib
#matplotlib.use('Agg')
#from pylab import plt
import csv
from collections import defaultdict

class BondBond_run(object):
    """ Class which runs and holds the results of a 3D bond-bond propensity
    analysis of a protein.
    
    Parameters
    ----------
    protein : Protein
      Protein object to be analysed
    source_residues : list of tuples
      Residues to be used as the source site in the analysis

    """

    def __init__(self, protein, source_residues):
        """


        """

        self.protein = protein
        self.source_residues = source_residues
        self.source_atoms = protein.get_atoms(source_residues)
        self.source_bonds = elasticnetwork.helpers.weak_bonds(
            self.source_atoms, self.protein.bonds
        )

        self.nbonds = None
        self.bond_results = None
        self.angle_results = None
        self.dihedral_results = None


                


    def calculate_propensities(self, bonds=True, angles=False):
       
 
        natoms = len(self.protein.atoms)
        nbonds = len(self.protein.bonds)

        #two-centre (bond) adjacency matrix
        A_bonds, G_bonds = generate_bond_matrices(self.protein)
        K_bonds = (A_bonds.T).dot(G_bonds).dot(A_bonds)

        #check for whether K is sparse
        if self.protein.angles:
            A_angles, G_angles = generate_angle_matrices(self.protein)
            K_angles = (A_angles.T).dot(G_angles).dot(A_angles)
        else:
            #just set angle stiffness matrix to zero
            K_angles = scipy.sparse.csr_matrix((3*natoms, 3*natoms))
        
        if self.protein.dihedrals:
            A_dihedrals, G_dihedrals = generate_dihedral_matrices(self.protein)
            K_dihedrals = (A_dihedrals.T).dot(G_dihedrals).dot(A_dihedrals)
        else:
            #just set dihedral stiffness matrix to zero
            K_dihedrals = scipy.sparse.csr_matrix((3*natoms, 3*natoms))


        K_total = K_bonds + K_angles + K_dihedrals

        source_bonds_ids = np.array(self.source_bonds.id())

        
        f_0 = np.zeros(nbonds)
        f_0[source_bonds_ids] = -1

        #currently, input is limited to two-centre stretches/compression
        f_ext = (A_bonds.T).dot(f_0.T)

        #sparse least squares method
        print('Solving least squares equation ...')
        #atom displacements
        u = scipy.sparse.linalg.lsqr(K_total, f_ext)
        print('..done.')


        if bonds:
            #bond stretches/compressions    
            bond_stretches = A_bonds.dot(u[0])

            #set source bonds to zero
            bond_stretches[source_bonds_ids] = 0
                   
            bond_ids = [bond.id for bond in self.protein.bonds]
            bond_names = [make_bond_name(bond) for bond in self.protein.bonds]

            atom1_pdbnum = [bond.atom1.PDBnum for bond in self.protein.bonds]
            atom2_pdbnum = [bond.atom2.PDBnum for bond in self.protein.bonds]

            atom1_chain = [bond.atom1.chain for bond in self.protein.bonds]
            atom2_chain = [bond.atom2.chain for bond in self.protein.bonds]

            bond_force_constants = [bond.force_constant for bond in self.protein.bonds]

            abs_change = np.abs(bond_stretches)

            abs_change_squared = np.multiply(abs_change, abs_change)
            strain_energy = np.multiply(bond_force_constants, abs_change_squared)

            bond_distance = elasticnetwork.helpers.distance_between(
                self.source_bonds, self.protein.bonds
            ).min(0)

            bond_results_frame = pd.DataFrame(
                {
                'bond_name' : bond_names, 
                'atom1' : atom1_pdbnum,
                'atom2' : atom2_pdbnum, 
                'distance': bond_distance,
                'stretch' : bond_stretches.tolist(),
                'abs_change' : abs_change.tolist(),
                'force_constants' : bond_force_constants,
                'strain_energy' : strain_energy,
                }, 
                index=bond_ids, 
            )
            self.bond_results = bond_results_frame

        
        """Need to write a distance_between function for angles - 'centre' of angle 
           is just the centre of mass of the three atoms """
        if angles:
            #angle stretches/compressions
            angle_stretches = A_angles.dot(u[0])
       

            angle_ids = [angle.id for angle in self.protein.angles]

            angle_names = [make_angle_name(angle) for angle in self.protein.angles]

            atom1_pdbnum = [angle.atom1.PDBnum for angle in self.protein.angles]
            atom2_pdbnum = [angle.atom2.PDBnum for angle in self.protein.angles]
            atom3_pdbnum = [angle.atom3.PDBnum for angle in self.protein.angles]

            atom1_chain = [angle.atom1.chain for angle in self.protein.angles]
            atom2_chain = [angle.atom2.chain for angle in self.protein.angles]
            atom3_chain = [angle.atom3.chain for angle in self.protein.angles]

            angle_force_constants = [angle.force_constant for angle in self.protein.angles]

            abs_change = np.abs(angle_stretches)
            abs_change_squared = np.multiply(abs_change, abs_change)

            strain_energy = np.multiply(angle_force_constants, abs_change_squared)
            angle_results_frame = pd.DataFrame(
                {
                'angle_name' : angle_names, 
                'atom1' : atom1_pdbnum,
                'atom2' : atom2_pdbnum,
                'atom3' : atom3_pdbnum, 
                'stretch' : angle_stretches.tolist(),
                'abs_change' : abs_change.tolist(),
                'force_constants' : angle_force_constants,
                'strain_energy' : strain_energy,
                }, 
                index=angle_ids, 
            )
            self.angle_results = angle_results_frame



    def calculate_bond_quantile_scores(self, meshsize=10):
        """ Calculates the quantile scores for each bond from 
        the perturbation propensity and distance from the source
        by using quantile regression to fit quantiles to the 
        results.  Stores the results in the bond_results dataframe.

        """
        if self.bond_results is None or self.bond_results.empty:
            raise ValueError('No bond results to calculate quantile scores for.' + 
                            'Run calculate_bond_propensities() first.')
        else:
            bond_results_nonzero = self.bond_results[self.bond_results['abs_change']!=0]
            nbonds_nonzero = len(bond_results_nonzero)
            n_quantiles = int(nbonds_nonzero/meshsize)
            quantile_fit_matrix = np.ndarray((nbonds_nonzero, n_quantiles))
            for i in range(n_quantiles):
                p = float(i)/n_quantiles
                log_fit = elasticnetwork.statistics.quant_reg(
                    bond_results_nonzero['distance'], 
                    np.log(bond_results_nonzero['abs_change']),
                    p
                )
                quantile_fit_matrix[:,i] = np.exp(log_fit)
        
            #: Use searchsorted to find the first fit line which is above
            #: each data point - this gives the quantile score for this bond
            qs = np.zeros(nbonds_nonzero)
            for i, bond_idx in enumerate(bond_results_nonzero.index):
                qs[i] = np.searchsorted(
                    quantile_fit_matrix[i,:], bond_results_nonzero['abs_change'][bond_idx]
                )
            
            bond_results_nonzero['qs'] = qs/n_quantiles
            self.bond_results = self.bond_results.join(bond_results_nonzero['qs'])
            self.bond_results.fillna(0, inplace=True)

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


def make_dihedral_name(dihedral):
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

    return (dihedral.atom1.res_name +
            dihedral.atom1.res_num + ' ' + 
            dihedral.atom1.chain + ' ' + 
            dihedral.atom1.name + ' : ' + 
            dihedral.atom2.res_name +
            dihedral.atom2.res_num + ' ' + 
            dihedral.atom2.chain + ' ' + 
            dihedral.atom2.name + ' : ' +
            dihedral.atom3.res_name +
            dihedral.atom3.res_num + ' ' + 
            dihedral.atom3.chain + ' ' + 
            dihedral.atom3.name + ' : ' +  
            dihedral.atom4.res_name +
            dihedral.atom4.res_num + ' ' + 
            dihedral.atom4.chain + ' ' + 
            dihedral.atom4.name)  

    
    

def generate_bond_matrices(protein):
    """
    Returns the (sparse) m by n incidence matrix A for the two-centre
    interaction (i.e. the bonds) and the m by m diagonal matrix of
    force constants.
    """
    natoms = len(protein.atoms)
    nbonds = len(protein.bonds)

    #A = np.zeros([nbonds, 3*natoms])
    A = scipy.sparse.lil_matrix((nbonds, 3*natoms))        

    force_constants = np.zeros(nbonds)
    for bond in protein.bonds:
        
        atom1_id = bond.atom1.id
        atom2_id = bond.atom2.id

        atom1_xyz = bond.atom1.xyz
        atom2_xyz = bond.atom2.xyz

        bond_length = np.linalg.norm(atom1_xyz - atom2_xyz)

        #row = A[bond.id]
        A[bond.id, [3*atom1_id, (3*atom1_id)+1, (3*atom1_id)+2]] = (atom1_xyz - atom2_xyz)/bond_length
        A[bond.id, [3*atom2_id, (3*atom2_id)+1, (3*atom2_id)+2]] = (atom2_xyz - atom1_xyz)/bond_length

        force_constant = bond.force_constant
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

    #A = np.zeros([nangles, 3*natoms])
    A = scipy.sparse.lil_matrix((nangles, 3*natoms))    

    force_constants = np.zeros(nangles)
    for angle in protein.angles:

        atom1_id = angle.atom1.id
        atom2_id = angle.atom2.id
        atom3_id = angle.atom3.id

        atom1_xyz = angle.atom1.xyz
        atom2_xyz = angle.atom2.xyz
        atom3_xyz = angle.atom3.xyz

        three_centre_length = np.linalg.norm(atom1_xyz - atom3_xyz)

        #row = A[angle.id]
        A[angle.id ,[3*atom1_id, (3*atom1_id)+1, (3*atom1_id)+2]] = (atom2_xyz - atom3_xyz)/three_centre_length
        A[angle.id ,[3*atom2_id, (3*atom2_id)+1, (3*atom2_id)+2]] = -((atom2_xyz - atom1_xyz) + (atom2_xyz - atom3_xyz))/three_centre_length
        A[angle.id ,[3*atom3_id, (3*atom3_id)+1, (3*atom3_id)+2]] = (atom2_xyz - atom1_xyz)/three_centre_length

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

    #A = np.zeros([ndihedrals, 3*natoms])
    A = scipy.sparse.lil_matrix((ndihedrals, 3*natoms))
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

        #row = A[dihedral.id]
        A[dihedral.id, [3*atom1_id, (3*atom1_id)+1, (3*atom1_id)+2]] = -((atom1_xyz - atom3_xyz) + (atom4_xyz - atom2_xyz))/four_centre_length 
        A[dihedral.id, [3*atom2_id, (3*atom2_id)+1, (3*atom2_id)+2]] = -((atom2_xyz - atom1_xyz) + (atom2_xyz - atom3_xyz) + (atom2_xyz - atom4_xyz))/four_centre_length
        A[dihedral.id, [3*atom3_id, (3*atom3_id)+1, (3*atom3_id)+2]] = -((atom3_xyz - atom4_xyz) + (atom3_xyz - atom1_xyz) + (atom3_xyz - atom2_xyz))/four_centre_length
        A[dihedral.id, [3*atom4_id, (3*atom4_id)+1, (3*atom4_id)+2]] = -((atom4_xyz - atom2_xyz) + (atom1_xyz - atom3_xyz))/four_centre_length

        force_constant = dihedral.force_constant
        force_constants[dihedral.id] = force_constant

    A = scipy.sparse.csr_matrix(A)
    G = scipy.sparse.diags(force_constants)

    return (A, G)

    
        

        










        # if dihedrals:
        #     #dihedral stretches/compressions
        #     dihedral_stretches = A_dihedrals.dot(u[0])
       

        #     dihedral_ids = [dihedral.id for dihedral in self.protein.dihedral.dihedral]

        #     dihedral_names = [make_dihedral_name(dihedral) for dihedral in self.protein.dihedral.dihedral]

        #     atom1_pdbnum = [dihedral.atom1.PDBnum for dihedral in self.protein.dihedral.dihedral]
        #     atom2_pdbnum = [dihedral.atom2.PDBnum for dihedral in self.protein.dihedral.dihedral]
        #     atom3_pdbnum = [dihedral.atom3.PDBnum for dihedral in self.protein.dihedral.dihedral]
        #     atom4_pdbnum = [dihedral.atom4.PDBnum for dihedral in self.protein.dihedral.dihedral]

        #     atom1_chain = [dihedral.atom1.chain for angle in self.protein.dihedral.dihedral]
        #     atom2_chain = [dihedral.atom2.chain for angle in self.protein.dihedral.dihedral]
        #     atom3_chain = [dihedral.atom3.chain for angle in self.protein.dihedral.dihedral]
        #     atom4_chain = [dihedral.atom4.chain for angle in self.protein.dihedral.dihedral]

        #     dihedral_force_constants = [dihedral.force_constant for dihedral in self.protein.dihedral.dihedral]

        #     abs_change = np.abs(dihedral_stretches)
        #     abs_change_squared = np.multiply(abs_change, abs_change)

        #     strain_energy = np.multiply(dihedral_force_constants, abs_change_squared)
        #     angle_results_frame = pd.DataFrame(
        #         {
        #         'dihedral_name' : angle_names, 
        #         'atom1' : atom1_pdbnum,
        #         'atom2' : atom2_pdbnum,
        #         'atom3' : atom3_pdbnum, 
        #         'atom4' : atom4_pdbnum,
        #         'stretch' : dihedral_stretches.tolist(),
        #         'abs_change' : abs_change.tolist(),
        #         'force_constants' : dihedral_force_constants,
        #         'strain_energy' : strain_energy,
        #         }, 
        #         index=angle_ids, 
        #     )
        #     self.dihedral_results = dihedral_results_frame



    



  


