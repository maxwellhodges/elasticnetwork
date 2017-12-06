import proteingraph3 as pg
import numpy as np
import pickle
pdk1 = pg.molecules3D.Protein()

pdk1_parser = pg.parsing3D.PDBParser('3orz_mm1.pdb')
pdk1_parser.parse(pdk1, strip={'res_name': ['HOH']})

ggenerator = pg.parsing3D.FIRST_Network_Generator()
ggenerator.generate_elastic_network(pdk1)

autocorrelation_run = pg.autocorrelation3D.Autocorrelation3D_run(pdk1)
autocorrelation_run.calculate_interaction_autocorrelation()

pickle.dump(pdk1, open('pdk1.p', 'wb'))
pickle.dump(autocorrelation_run.results, open('pdk1_autocorr_results.p', 'wb'))
