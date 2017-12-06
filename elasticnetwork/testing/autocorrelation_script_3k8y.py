import proteingraph3 as pg
import numpy as np
import pickle

hras = pg.molecules3D.Protein()

hras_parser = pg.parsing3D.PDBParser('3k8y.pdb')
hras_parser.parse(hras)

ggenerator = pg.parsing3D.FIRST_Network_Generator()
ggenerator.generate_elastic_network(hras)

autocorrelation_run = pg.autocorrelation3D.Autocorrelation3D_run(hras)
autocorrelation_run.calculate_interaction_autocorrelation()


pickle.dump(autocorrelation_run.results, open('hras_water_autocorr_results.p', 'wb'))
