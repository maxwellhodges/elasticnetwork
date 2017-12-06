import proteingraph3 as pg
import numpy as np
import pickle
adk = pg.molecules3D.Protein()

adk_parser = pg.parsing3D.PDBParser('2rgx.pdb')
adk_parser.parse(adk, strip={'res_name': ['HOH']})

ggenerator = pg.parsing3D.FIRST_Network_Generator()
ggenerator.generate_elastic_network(adk)

autocorrelation_run = pg.autocorrelation3D.Autocorrelation3D_run(adk)
autocorrelation_run.calculate_interaction_autocorrelation()

pickle.dump(adk, open('adk.p', 'wb'))
pickle.dump(autocorrelation_run.results, open('adk_autocorr_results.p', 'wb'))
