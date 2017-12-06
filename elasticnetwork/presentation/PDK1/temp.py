# """PDK1 3D """

# import proteingraph3 as pg
# import pandas as pd

# pdk1 = pg.molecules3D.Protein()

# pdk1_parser = pg.parsing3D.PDBParser('3orz_mm1.pdb')
# pdk1_parser.parse(pdk1, strip={'res_name': ['HOH']})

# """Bonds"""
# ggenerator = pg.parsing3D.FIRST_Network_Generator()
# ggenerator.generate_elastic_network(pdk1, angles=False, dihedrals=False)


# source_residues = [('1','A')] 

# propensity_calcs = pg.bondbond3D.BondBond3D_run(pdk1, source_residues)
# propensity_calcs.calculate_3D_propensities()
# (propensity_calcs.bond_results).to_csv('pdk1_bondbond.csv')



# infrig_run = pg.infrig.Infrig_run(pdk1)
# infrig_run.calculate_infrig(angles=False, dihedrals=False)
# (infrig_run.results).to_csv('pdk1_infrig.csv')


# """Angles"""
# ggenerator_angles = pg.parsing3D.FIRST_Network_Generator()
# ggenerator_angles.generate_elastic_network(pdk1, angles=True, dihedrals=False)

# propensity_calcs_angles = pg.bondbond3D.BondBond3D_run(pdk1, source_residues)
# propensity_calcs_angles.calculate_3D_propensities()
# (propensity_calcs_angles.bond_results).to_csv('pdk1_bondbond_angles.csv')

# infrig_run_angles = pg.infrig.Infrig_run(pdk1)
# infrig_run_angles.calculate_infrig(angles=True, dihedrals=False)
# (infrig_run_angles.results).to_csv('pdk1_infrig_angles.csv')


import proteingraph as pg
pdk1 = pg.molecules.Protein()
pdk1_parser = pg.parsing.PDBParser('3orz_mm1.pdb')
pdk1_parser.parse(pdk1, strip={'res_name': ['HOH']})
ggenerator = pg.parsing.FIRST_Graph_Generator()
ggenerator.generate_graph(pdk1)
source_residues = [('1','A')] 
propensity_calcs = pg.bondbond.BondBond_run(pdk1, source_residues, bond_type='all')
propensity_calcs.calculate_bond_propensities()
propensity_calcs.bond_results_to_csv('/home/maxhodges/proteinnetwork/proteingraph3/presentation/PDK1/pdk1_proteingraph.csv')
