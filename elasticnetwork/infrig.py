import numpy as np
import scipy 
import scipy.sparse
import elasticnetwork.molecules
from itertools import combinations
import pandas as pd


class Infrig_run(object):

    def __init__(self, protein):

        self.protein = protein
        self.atom_list = []
        self.flexible_atoms = []
        #Rigidity matrix
        self.clusters = dict()
        self.R = np.array([])
        self.results = None


    
    
    def calculate_infrig(self, angles=True, dihedrals = True):
        """ This is the main function that calculates rigidity and returns a "tbc" 
        containing the infomation about the rigid clusters.
        """
        
        # """Keep track of atoms that still need to be idenitified as rigid or flexible"""
        # self.atom_list = []
        # self.flexible_atoms = []

        """Generate the rigidity matrix """
        R_bonds = generate_bond_rigidity(self.protein)
        self.R = R_bonds

        if angles and dihedrals:
            R_angles = generate_angle_rigidity(self.protein)
            R_dihedrals = generate_dihedral_rigidity(self.protein)
            self.R = np.vstack((R_bonds, R_angles, R_dihedrals))

        elif angles:
            R_angles = generate_angle_rigidity(self.protein)
            self.R = np.vstack((R_bonds, R_angles))    

        """ 
        Find the infinitesimal motions of the protein from the overall Rigidity matrix.
        Use SVD - only want the vectors corresponding to 0 singular values i.e. bottom n - m rows of V^T and any extra rows
        that go with 0 singular values in S.
        """        
        print("Performing SVD of the Rigidity Matrix...")
        u,s,vt = np.linalg.svd(self.R)
        print("SVD calculation completed.")
        singular_vals, = np.where(s < 10e-15)
        #this code only currently works for when m > 3n
        #so will need to add more logic for case where only bonds are included
        #keep it simpler for now whilst testing...
        if self.R.shape[0] > self.R.shape[1]:
            inf_motions = vt[singular_vals,:]
            #scale the inf motions by the size of the protein for numerical stability
            inf_motions = (len(self.protein.atoms))*inf_motions

        else:
            inf_motions = np.vstack((vt[singular_vals,:], vt[self.R.shape[0]:, :]))

        
        """only want loop over atoms that form a tetrahedron"""
        for atom in self.protein.atoms:
            if len(atom.connected_atoms) > 2:
                self.atom_list.append(atom)
            else:
                self.flexible_atoms.append(atom.id)
        

        #think I want the while loop here - while the atom_list still has members, continue looping over
        #the list, which should get smaller every loop as atoms are added either to the flexible or to 
        #a rigid cluster.  i.e. while len(atom_list) != 0.  Need to make sure there is no way to get an
        #infinite loop here...
        cluster_index = 0
        while len(self.atom_list) != 0:
            
            
            """code for finding rigid tetrahedron - adds atoms to flexible_atoms and removes from
            atom_list if they are found to not be part of a rigid tetrahedron"""
            tet = self.find_rigid_tetrahedron()
            if len(tet) == 0:
                for atom in self.atom_list:
                    self.flexible_atoms.append(atom.id)
                    self.atom_list.remove(atom)
                break 
            #create dictionary to allow easy assignment of R_test below
            tet_dict = dict(zip(tet, [0,1,2,3]))

            """After previous iteration of while loop and after searching for tetrahedrons - create a dictionary that maps atom id to index within the atom_list """
            atomid_to_index = dict()
            for i, atom in enumerate(self.atom_list):
                atomid_to_index[atom.id] = i

            
            """Now tet is our tetrahedron that we use as the basis of the infinitesimal motions. """
            tet_edges = [edge for edge in combinations(tet,2)]
            R_test = np.zeros([6,12])
            for i, edge in enumerate(tet_edges):

                atom1_id = edge[0]
                atom2_id = edge[1]

                atom1_xyz = self.protein.atoms[edge[0]].xyz
                atom2_xyz = self.protein.atoms[edge[1]].xyz

                row = R_test[i]
                row[[3*tet_dict[atom1_id], (3*tet_dict[atom1_id])+1, (3*tet_dict[atom1_id])+2]] = (atom1_xyz - atom2_xyz)
                row[[3*tet_dict[atom2_id], (3*tet_dict[atom2_id])+1, (3*tet_dict[atom2_id])+2]] = (atom2_xyz - atom1_xyz)

            
          
            rigid_atoms = []
            for i, inf_motion in enumerate(inf_motions):
                #this then needs to be part of a for loop over all of the infinitesimal motions
                #should stick in some progress info to keep track of where the calculation is
                if i%10 == 0:
                    print("Loop " + str(i) + " of " + str(len(inf_motions)) )

                """Need to change original coordinate system to the principal axes of the tetrahedron so the 6 rigid motions in the 12-space are orthogonal """
                p0, p1, p2, p3 = (self.protein.atoms[tet[0]]).xyz, (self.protein.atoms[tet[1]]).xyz, (self.protein.atoms[tet[2]]).xyz, (self.protein.atoms[tet[3]]).xyz
                

                """Get the infinitesimal motions for all atoms left in atom_list in the reference frame of the tetrahedron."""
                coords = []
                for atom in self.atom_list:
                    coord = atom.xyz
                    coords.append(coord)

                new_coords, K = self.new_coord_frame(coords, p0,p1,p2,p3)
                
                tet_indices = [atomid_to_index[tet[0]], atomid_to_index[tet[1]], atomid_to_index[tet[2]], atomid_to_index[tet[3]]]
                    
                
                """Generate the 6 orthogonal 12-vectors describing the rigid motion of the tetrahedron """
                inf_tx, inf_ty, inf_tz, t_x_t, t_y_t, t_z_t = generate_translation_vectors(new_coords, tet_indices , K)
                inf_rx, inf_ry, inf_rz, r_x_t, r_y_t, r_z_t = generate_rotation_vectors(new_coords, tet_indices , K)

                """ Project the infinitesimal motions of the tetrahedron onto each of the six orthogonal infinitesimal motions.
                From each infinitesimal motion of all the atoms in the protein, subtract the component in each of the 6 directions.
                Any atom that is also part of the rigid cluster with the tetrahedron should go back to its original position in all infinitesimal motions """

                #only want to calculate those atoms that are still in self.atom_list
                #split into arrays of 3 (xyz of each atom) and convert to np array then slice
                inf_motion = np.array(np.split(inf_motion, len(inf_motion)/3))[np.array([atom.id for atom in self.atom_list])]
                #need get the tet inf motions from here, leave in list as will be checked by abs < 10e-4 anyway
                tet_motions = (inf_motion[np.array(tet_indices)]).flatten()

                #then flatten again
                inf_motion = inf_motion.flatten()

                final_motions = inf_motion - ( (np.dot(tet_motions, t_x_t)*inf_tx +
                                                np.dot(tet_motions, t_y_t)*inf_ty +
                                                np.dot(tet_motions, t_z_t)*inf_tz + 
                                                np.dot(tet_motions, r_x_t)*inf_rx +
                                                np.dot(tet_motions, r_y_t)*inf_ry +
                                                np.dot(tet_motions, r_z_t)*inf_rz ) )
                
                #finally, need some logic here to identify those atoms that move back to their original positions
                #with the tetrahedron.  Have to go back over all infinitesimal motions so want intersection of results
                #from each inf motion.  Will need a for loop over each of the atom xyz coords from the atom_list.  If the 
                #the euclidean norm of the atom's xyz coords is < 10e-4 then add it to a list within this loop.

                #turn 3N array into N array - makes it easier to index
                final_motions_split = np.split(final_motions, len(inf_motion)/3)
                rigid_atoms_iteration = []
                for i, atom in enumerate(self.atom_list):
                    if np.linalg.norm(final_motions_split[i]) < 10e-4:
                        rigid_atoms_iteration.append(atom.id)

                rigid_atoms.append(rigid_atoms_iteration)

            #now need to only keep those atoms that appear in every iteration
            rigid_atoms_final = set(rigid_atoms[0]).intersection(*rigid_atoms)

            self.clusters[cluster_index] = rigid_atoms_final

            #need to remove rigid atoms from atom_list
            for i in rigid_atoms_final:
                (self.atom_list).remove(self.protein.atoms[i])
            
            #also need to remove atoms that now have connected < 3 once rigid atoms removed
            atom_list_ids = [atom.id for atom in self.atom_list]
            for atom in self.atom_list:
                if len(set(atom.connected_atoms).intersection(atom_list_ids)) < 3:
                    self.flexible_atoms.append(atom.id)
                    self.atom_list.remove(atom)


            print ('Cluster ' + str(cluster_index) + ' completed...')
            print (str(len(self.atom_list)) + ' atoms left in atom_list...')
            cluster_index += 1

        print ('0 atoms left in atom_list.')

        atom_names = [make_atom_name(atom) for atom in self.protein.atoms]
        cluster_ids = np.zeros(len(self.protein.atoms))
        for i in self.clusters:                              
            for j in self.clusters[i]:         
                cluster_ids[j] = i+1 #so that 'cluster 0' are the flexible atoms

        cluster_ids = [int(i) for i in cluster_ids.tolist()]

        atom_results_frame = pd.DataFrame(
            {
            'atom_name' : atom_names, 
            'cluster_id' : cluster_ids,
            }, 
            index=self.protein.atoms.id(), 
        )
        self.results = atom_results_frame






    def find_rigid_tetrahedron(self):

        #need to add in line here to add atom to 'flexible' list if rigid tet
        #cannot be found for it as it must by definition be flexible at that point
        atom_list_ids = [atom.id for atom in self.atom_list]
        #on last loop, need to return something or tet is not assigned
        tet = []
        for atom in self.atom_list:
            connected = set(atom.connected_atoms).intersection(atom_list_ids)
            combs = [i for i in combinations(connected,3)]
            for comb in combs:
                test_tet_atoms = np.hstack([atom.id, comb])
                indices = np.hstack((3*test_tet_atoms, (3*test_tet_atoms)+1, (3*test_tet_atoms)+2)) 
                #need to sort indices or mix up columns of R
                indices = np.sort(indices)
                #filter out only those columns corresponding to chosen atoms
                R_tet = self.R[:, indices]
                #filter out rows that contain fewer than 6 entries - other rows are 'dangling edges' 
                remove_dangling = np.sum(np.abs(R_tet) > 10e-15, axis=1) == 6
                R_tet_stripped = R_tet[remove_dangling]
                #now only those edges between the set of 4 atoms are present.
                #Do SVD to check for rigidity (can't just count edges as may not be independent).
                u,s,vh = np.linalg.svd(R_tet_stripped)
                if (s.shape[0] - vh.shape[0]) + np.sum(s < 10e-15) == 6:
                    tet = [atom.id, comb[0], comb[1], comb[2]]
                    break
            else:
                self.flexible_atoms.append(atom.id)
                self.atom_list.remove(atom)
                continue
            
            break

        return tet


    def new_coord_frame(self, coords, p0, p1, p2, p3):

        centre_of_mass = (p0 + p1 + p2 + p3)/4
        (p0, p1, p2, p3) = (p0, p1, p2, p3) - centre_of_mass

        I = np.zeros([3,3])
        for i in [0,1,2]:
            for j in [0,1,2]:
                I[i,j] = ((np.linalg.norm(p0)*(i == j) - p0[i]*p0[j]) +
                        (np.linalg.norm(p1)*(i == j) - p1[i]*p1[j]) +
                        (np.linalg.norm(p2)*(i == j) - p2[i]*p2[j]) +
                        (np.linalg.norm(p3)*(i == j) - p3[i]*p3[j]) )

        #K has rows that are the eigenvectors of I
        w, v = np.linalg.eig(I)

        eig0 = np.array([v[0][0], v[1][0], v[2][0]])
        eig1 = np.array([v[0][1], v[1][1], v[2][1]])
        eig2 = np.array([v[0][2], v[1][2], v[2][2]])

        K = np.array([eig0,eig1,eig2])

        new_coord_list = []
        for i, coord in enumerate(coords):
            new_coord = coord - centre_of_mass
            new_coord = np.dot(K,new_coord)
            new_coord_list.append(new_coord)


        return (new_coord_list, K)



def generate_rotation_vectors(coords, tet_indices, K):
    """ 
    coords: 
    tet_indices:
    K: the rotation matrix that converts the coordinates to the principal axes frame
    """
    e_x, e_y, e_z = np.array([1,0,0]), np.array([0,1,0]), np.array([0,0,1])
    """Generate infinitesimal rigid motions in the new coordinate system """
    inf_rig_rxs_new, inf_rig_rys_new, inf_rig_rzs_new = [],[],[]
    for coord in coords:
        inf_rig_rx = np.cross(e_x, coord)
        inf_rig_rxs_new.append(inf_rig_rx)
        inf_rig_ry = np.cross(e_y, coord)
        inf_rig_rys_new.append(inf_rig_ry)
        inf_rig_rz = np.cross(e_z, coord)
        inf_rig_rzs_new.append(inf_rig_rz)
    
    """Convert back to the original coordinate frame """
    inf_rig_rxs_orig, inf_rig_rys_orig, inf_rig_rzs_orig = [],[],[]
    for xrot in inf_rig_rxs_new:
        xrot_orig = np.dot(K.T, xrot)
        inf_rig_rxs_orig.append(xrot_orig)
    for yrot in inf_rig_rys_new:
        yrot_orig = np.dot(K.T, yrot)
        inf_rig_rys_orig.append(yrot_orig)
    for zrot in inf_rig_rzs_new:
        zrot_orig = np.dot(K.T, zrot)
        inf_rig_rzs_orig.append(zrot_orig)

    """Extract the tet inf rigid motions """
    r_x_t = (np.array(inf_rig_rxs_orig)[np.array(tet_indices)]).flatten()
    r_y_t = (np.array(inf_rig_rys_orig)[np.array(tet_indices)]).flatten()
    r_z_t = (np.array(inf_rig_rzs_orig)[np.array(tet_indices)]).flatten()

    """Normalise over the tetrahedron motions """
    inf_rig_rxs_orig = (np.array(inf_rig_rxs_orig).flatten())/np.linalg.norm(r_x_t)
    inf_rig_rys_orig = (np.array(inf_rig_rys_orig).flatten())/np.linalg.norm(r_y_t)
    inf_rig_rzs_orig = (np.array(inf_rig_rzs_orig).flatten())/np.linalg.norm(r_z_t)

    r_x_t = r_x_t/np.linalg.norm(r_x_t)
    r_y_t = r_y_t/np.linalg.norm(r_y_t)
    r_z_t = r_z_t/np.linalg.norm(r_z_t)

    return (inf_rig_rxs_orig, inf_rig_rys_orig, inf_rig_rzs_orig,
            r_x_t, r_y_t, r_z_t)



def generate_translation_vectors(coords, tet_indices, K):
    """
    node_coords: a 1 by 12 numpy array of the xyz coords of the tetrahedron
    K: the rotation matrix that converts the coordinates to the principal axes frame
    """
    e_x, e_y, e_z = np.array([1,0,0]), np.array([0,1,0]), np.array([0,0,1])

    """Generate infinitesimal rigid motions in the new coordinate system """
    t_x, t_y, t_z = e_x, e_y, e_z

    """Convert back to the original coordinate frame """
    t_x, t_y, t_z = np.dot(K.T, t_x), np.dot(K.T, t_y), np.dot(K.T, t_z)
    
    inf_rig_txs_orig = np.tile(t_x, len(coords))
    inf_rig_tys_orig = np.tile(t_y, len(coords))
    inf_rig_tzs_orig = np.tile(t_z, len(coords))
    
    """Get the tetrahedron motions"""
    t_x_t = np.hstack((t_x, t_x, t_x, t_x))
    t_y_t = np.hstack((t_y, t_y, t_y, t_y))
    t_z_t = np.hstack((t_z, t_z, t_z, t_z))

    """Normalise over the tetrahedron motions """
    inf_rig_txs_orig = inf_rig_txs_orig/np.linalg.norm(t_x_t)
    inf_rig_tys_orig = inf_rig_tys_orig/np.linalg.norm(t_y_t)
    inf_rig_tzs_orig = inf_rig_tzs_orig/np.linalg.norm(t_z_t) 

    t_x_t = t_x_t/np.linalg.norm(t_x_t)
    t_y_t = t_y_t/np.linalg.norm(t_y_t)
    t_z_t = t_z_t/np.linalg.norm(t_z_t)    


    return (inf_rig_txs_orig, inf_rig_tys_orig, inf_rig_tzs_orig,
            t_x_t, t_y_t, t_z_t)

    
    


def generate_bond_rigidity(protein):
    """
    Returns the (sparse) m by n rigidity matrix A for the two-centre
    interaction (i.e. the bonds) and the m by m diagonal matrix of
    force constants.
    """
    natoms = len(protein.atoms)
    nbonds = len(protein.bonds)

    R = np.zeros([nbonds, 3*natoms])
    for bond in protein.bonds:
        
        atom1_id = bond.atom1.id
        atom2_id = bond.atom2.id

        atom1_xyz = bond.atom1.xyz
        atom2_xyz = bond.atom2.xyz


        row = R[bond.id]
        row[[3*atom1_id, (3*atom1_id)+1, (3*atom1_id)+2]] = (atom1_xyz - atom2_xyz)
        row[[3*atom2_id, (3*atom2_id)+1, (3*atom2_id)+2]] = (atom2_xyz - atom1_xyz)


    return R


def generate_angle_rigidity(protein):
    """
    Returns the (sparse) m by n incidence matrix for the three-centre
    interaction (i.e. the angles) and the m by m diagonal matrix of
    force constants.
    """

    #double check maths for this to be safe (particularly signs)

    natoms = len(protein.atoms)
    nangles = len(protein.angles)

    R = np.zeros([nangles, 3*natoms])
    for angle in protein.angles:

        atom1_id = angle.atom1.id
        atom3_id = angle.atom3.id

        atom1_xyz = angle.atom1.xyz
        atom3_xyz = angle.atom3.xyz

        row = R[angle.id]
        row[[3*atom1_id, (3*atom1_id)+1, (3*atom1_id)+2]] = (atom1_xyz - atom3_xyz)
        row[[3*atom3_id, (3*atom3_id)+1, (3*atom3_id)+2]] = (atom3_xyz - atom1_xyz)
        

    return R


def generate_dihedral_rigidity(protein):
    """
    Returns the (sparse) m by n incidence matrix for the four-centre
    interaction (i.e. the dihedrals) and the m by m diagonal matrix of
    force constants.
    """

    #double check maths for this to be safe (particularly signs)

    natoms = len(protein.atoms)
    ndihedrals = len(protein.dihedrals)

    R = np.zeros([ndihedrals, 3*natoms])
    for dihedral in protein.dihedrals:
        
        atom1_id = dihedral.atom1.id
        atom4_id = dihedral.atom4.id

        atom1_xyz = dihedral.atom1.xyz
        atom4_xyz = dihedral.atom4.xyz


        row = R[dihedral.id]
        row[[3*atom1_id, (3*atom1_id)+1, (3*atom1_id)+2]] = (atom1_xyz - atom4_xyz)
        row[[3*atom4_id, (3*atom4_id)+1, (3*atom4_id)+2]] = (atom4_xyz - atom1_xyz)


    return R



def make_atom_name(atom):
        """ Returns a string giving the name of the two atoms in 
        a given atom.
        
        Parameters
        ----------
        atom : :class:'Atom'
          atom object for which string is to be returned
        
        Returns
        -------
        result : string
          Name of the atom in the form 'ResNameResNum Chain AtomName'
        
        """

        return (atom.res_name +
                atom.res_num + ' ' +
                atom.chain + ' ' +
                atom.name)



