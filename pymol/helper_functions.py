import colorsys,sys,re
from pymol import cmd
import numpy as np
import csv

def color_atoms(filename, color_by='thalf', power=1, gradient="bwr", lower=None, upper=None):

    power = float(power)

    atoms = []
    f = open(filename, "r")
    reader = csv.reader(f)
    header = reader.next()

    column = [i for i, col in enumerate(header)
              if col == color_by]

    if len(column) == 0:
        raise ValueError('No %s column found.  Check your input file' % (color_by))
    elif len(column) > 1:
        raise ValueError('More than one %s column' % (color_by))
    else:
        column = column[0]

    for row in reader:
        resi, chain, atom = row[1].split(" ")
        res_name = resi[0:3]
        res_num = int(resi[3:])
        chain = chain.strip()
        atom = atom.strip()
        result = row[column]
        atoms.append([chain, res_num, res_name, atom, float(result)])
    
    f.close()
    num_atoms = len(atoms)
        
    scores = [atom[4] for atom in atoms]
    if lower:
        score_min = float(lower)
    else:
        score_min = min(scores)    
    if upper:
        score_max = float(upper)
    else:
        score_max = max(scores)
    
    for atom in atoms:
        chain = atom[0]
        res_num = atom[1]
        res_name = atom[2]
        name = atom[3]
        score = atom[4]

        scaled_score = (score - score_min)/(score_max - score_min)
        
        curve1 = round((2**power)*(scaled_score**power), 5)
        curve2 = round((2**power)*((1-scaled_score)**power), 5)
        
        rgb = [min(1, curve1), 
               min(curve1, curve2),
               min(curve2, 1)]
        
        cmd.set_color('r_color' + chain + str(res_num) + res_name + name,rgb)
        cmd.color('r_color' + chain + str(res_num) + res_name + name, 
                    'chain ' + chain + ' AND resi ' + str(res_num) + ' AND resn ' + res_name + ' AND name ' + name)

cmd.extend('color_atoms',color_atoms)

def show_atoms(filename, threshold, show_by='thalf', below=False, type_ = 'spheres', hide = 'no'):

    threshold = float(threshold)

    atoms = []
    f = open(filename, "r")
    reader = csv.reader(f)
    header = reader.next()

    column = [i for i, col in enumerate(header)
              if col == show_by]

    if len(column) == 0:
        raise ValueError('No %s column found.  Check your input file' % (color_by))
    elif len(column) > 1:
        raise ValueError('More than one %s column' % (color_by))
    else:
        column = column[0]

    for row in reader:
        resi, chain, atom = row[1].split(" ")
        res_name = resi[0:3]
        res_num = int(resi[3:])
        chain = chain.strip()
        atom = atom.strip()
        result = row[column]
        atoms.append([chain, res_num, res_name, atom, float(result)])

    if hide == 'yes':
        cmd.hide(type_)
    
    for atom in atoms:
        score = float(atom[4])
        res_num = atom[1]
        atom_name = atom[3]
        chain = atom[0]
        if not below:
            if score > threshold:
                cmd.show(type_,'resi ' + str(res_num) + ' and name ' + atom_name + ' and chain ' + chain)
        else:
            if score < threshold:
                cmd.show(type_,'resi ' + str(res_num) + ' and name ' + atom_name + ' and chain ' + chain)
    f.close()


def color_residues(filename, color_by='propensity', power=1, gradient="bwr", scale=None):

    power = float(power)

    f = open(filename, "r")
    reader = csv.reader(f)

    residues = []
    
    header = reader.next()
    column = [i for i, col in enumerate(header)
              if col == color_by]

    if len(column) == 0:
        raise ValueError('No %s column found.  Check your input file' % (color_by))
    elif len(column) > 1:
        raise ValueError('More than one %s column' % (color_by))
    else:
        column = column[0]
    
    for line in reader:
        residues.append((line[0], line[1], float(line[column])))   #changed from line[1] and line[2]
    f.close()
    num_residues = len(residues)
    
    scores = [float(residue[2]) for residue in residues] 
    if scale:
        min_value = scale[0]
        max_value = scale[1]
    else:
        max_value = max(scores)
        min_value = min(scores)
    
    if gradient=='bwr':
        for residue in residues:
            res_num, chain, score = residue
            
            scaled_score =  (score - min_value)/(max_value - min_value)
            
            curve1 = round((2**power)*(scaled_score**power), 5)
            curve2 = round((2**power)*((1-scaled_score)**power), 5)
        
            rgb = [min(1, curve1), 
                    min(curve1, curve2),
                    min(curve2, 1)]
        
            cmd.set_color('r_color' + res_num + chain, rgb)
            cmd.color('r_color' + res_num + chain, 'resi '+ res_num +' AND chain '+ chain)
            print 'colored residue' + res_num + ' chain '+ chain

def show_residues(filename, threshold, show_by='pp', rep='sticks', mol=None):

    threshold = float(threshold)

    f = open(filename, "r")
    reader = csv.reader(f)

    residues = []
    
    header = reader.next()
    column = [i for i, col in enumerate(header)
              if col == show_by]

    if len(column) == 0:
        raise ValueError('No %s column found.  Check your input file' % (color_by))  #should say show_by ?
    elif len(column) > 1:
        raise ValueError('More than one %s column' % (color_by))  #should say show_by ?
    else:
        column = column[0]
    
    for line in reader:
        residues.append((line[0], line[1], float(line[column]))) #changed from 1,2 to 0,1
    f.close()
    num_residues = len(residues)
    
    for residue in residues:
        res_num, chain, score = residue
        if score > threshold:
            if mol:
                cmd.show(rep, mol + ' and resi ' + str(res_num) + ' and chain ' + str(chain))
            else:
                cmd.show(rep, 'resi ' + str(res_num) + ' and chain ' + str(chain))

    return True

def load_bonds(filename, prefix='bond', mol=None):

    f = open(filename, 'r')
    reader = csv.reader(f)
    header = reader.next()
    
    # if header != ['bond_id', 'atom1_id', 'atom2_id',
    #               'bond_name', 'weight', 'distance',
    #               'propensity', 'qs', 'qs_test_set']:
    #     raise ValueError('header is not correct format')
    
    for row in reader:
        bond_id, atom1_id, atom2_id = row[0], row[1], row[2]
        if mol:
            cmd.distance(prefix + bond_id,
                        'id ' + atom1_id + ' and ' + mol,   #changed AND to and
                        'id ' + atom2_id + ' and ' + mol)   #changed AND to and
        else:
            cmd.distance(prefix + bond_id,
                        'id ' + atom1_id,
                        'id ' + atom2_id)
    cmd.hide('labels')
    cmd.set('dash_gap',0)

cmd.extend('load_bonds', load_bonds)

def color_bonds(filename, color_by='pp', name='bond', 
                power=1, gradient="bwr", scale=None):

    power = float(power)

    bonds = []

    f = open(filename, "r")
    reader = csv.reader(f)

    header = reader.next()
    
    column = [i for i, col in enumerate(header)
              if col == color_by]
    if len(column) > 1:
        raise ValueError('More than one {0} column'.format(color_by))
    else:
        column = column[0]

    for line in reader:
        bonds.append([line[0], float(line[column])])
    f.close()
    
    num_bonds = len(bonds)
        
    scores = [bond[1] for bond in bonds]
    if scale:
        score_min = float(scale[0])
        score_max = float(scale[1])
    else:
        score_max = max(scores)
        score_min = min(scores)

    for bond in bonds:
        id_ = bond[0]
        score = bond[1]

        scaled_score = (score - score_min)/(score_max - score_min)
    
        curve1 = round((2**power)*(scaled_score**power), 5)
        curve2 = round((2**power)*((1-scaled_score)**power), 5)
        
        rgb = [min(1, curve1), 
               min(curve1, curve2),
               min(curve2, 1)]
    
        cmd.set_color('b_color' + name + id_, rgb)
        cmd.color('b_color' + name + id_, name + id_)


def show_bonds(filename, threshold, show_by='pp', name='bond', 
               mol=None, hide = 'no', residues='no'):

    threshold = float(threshold)

    if hide == 'yes':
        cmd.hide('dashes')
        cmd.hide('sticks')

    f =  open(filename, "r")
    reader = csv.reader(f)
    
    header = reader.next()
    
    score_column = [i for i, col in enumerate(header)
              if col == show_by]
    if len(score_column) > 1:
        raise ValueError('More than one {0} column'.format(show_by))
    else:
        score_column = score_column[0]
    
    if residues == "yes": #need to add scenario where residues != "yes" or atom1_column is not assigned
        atom1_column = [i for i, col in enumerate(header)
                        if col == 'atom1_id']
        if len(atom1_column) > 1:
            raise ValueError('More than one {0} column'.format('atom1_id'))
        else:
            atom1_column = atom1_column[0]
                
        atom2_column = [i for i, col in enumerate(header)
                        if col == 'atom2_id']
        if len(atom2_column) > 1:
            raise ValueError('More than one {0} column'.format('atom2_id'))
        else:
            atom2_column = atom2_column[0]
    
    bond_list = []
    for line in reader:
        bond_list.append(
            [line[0], line[atom1_column], line[atom2_column],
             float(line[score_column])]
        )
    f.close()

    for bond in bond_list:
        if bond[3] >= threshold:
            cmd.show('dashes', name + str(bond[0]))
            if residues == 'yes' and mol:
                cmd.show('sticks','br. id '+ str(bond[1]) + ' AND ' + mol)
                cmd.show('sticks','br. id '+ str(bond[2]) + ' AND ' + mol)
            elif residues == 'yes':
                cmd.show('sticks','br. id '+ str(bond[1]))
                cmd.show('sticks','br. id '+ str(bond[2]))

                                        
cmd.extend('color_residues',color_residues)
cmd.extend('show_residues', show_residues)
cmd.extend('show_atoms',show_atoms)
cmd.extend('color_atoms',color_atoms) 
cmd.extend('load_bonds', load_bonds)
cmd.extend('show_bonds',show_bonds)
cmd.extend('color_bonds',color_bonds)


