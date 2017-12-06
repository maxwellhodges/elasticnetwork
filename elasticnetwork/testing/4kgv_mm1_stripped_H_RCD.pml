# Rigid cluster decomposition coloring script for PyMol
#
# Created by Dan Farrell, Brandon Hespenheide.
# Department of Physics and Astronomy
# Biophysics Theory Group
# Arizona State University
############################################################

from pymol import cmd
from pymol.cgo import *

bg_color white
load 4kgv_mm1_stripped_H_RCD.pdb

set line_width = 1
color black

color 0x0000b2, ( b > 0.99 and b < 1.01)
# Rigid Cluster 1 has 41446 atoms.
create RC1, ( b > 0.99 and b < 1.01)
show lines, RC1
set line_width = 3, RC1
color 0x0000b2, RC1

# Rigid Cluster BIN2
create BIN2, ( b > 1.99 and b < 2.01)
show lines, BIN2
set line_width = 3, BIN2
color gray, BIN2
disable BIN2

# Rigid Cluster BIN1
create BIN1, ( b > 2.99 and b < 490.01)
show lines, BIN1
set line_width = 3, BIN1
color gray, BIN1
disable BIN1

