import os
import sys

import manipulate_xyz as m

save='shifted.xyz'
read='tetra.xyz'
shiftx=1.0
shifty=1.0
shiftz=1.0

m.shift_xyz(save,read,shiftx,shifty,shiftz)

