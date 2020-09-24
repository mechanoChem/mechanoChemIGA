#!/usr/bin/env python

import numpy as np
from time import time
import ctypes
ctypes.CDLL("libmpi.so", mode=ctypes.RTLD_GLOBAL)


from PYmechanoChem import PYmechanoChem
problem = PYmechanoChem()
problem.setup_mechanoChemIGA()

begin = time()
for i in range(10):
    problem.simulate()
print(time() - begin)
