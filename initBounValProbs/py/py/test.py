#!/usr/bin/env python


# https://stackoverflow.com/questions/3761391/boostpython-python-list-to-stdvector
# https://wiki.python.org/moin/boost.python/extract
# https://www.boost.org/doc/libs/1_39_0/libs/python/test/vector_indexing_suite.cpp

import numpy as np
from time import time

from PYmechanoChem import PYmechanoChem
problem = PYmechanoChem()
problem.setup_mechanoChemIGA()

begin = time()
for i in range(10):
    problem.simulate()
print(time() - begin)
