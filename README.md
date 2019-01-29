# ase_dcdftbmd
ASE Calculator for DC-DFTB-MD

Example:
```python
from dcdftbmd import DCDFTBMD

from ase import Atoms
from ase.optimize import BFGS
from ase.io.xyz import write_xyz

import numpy as np

dcdftb_keywords = """
SCC=(DAMPXH=TRUE DAMPXHZETA=4.00 THIRDFULL=TRUE)
DC=FALSE
MISC=(FORCE=TRUE)
"""


dcdftb_skinfo = """
2
O 2 -0.1575
     O-O.skf O-H.skf
H 1 -0.1857
     H-O.skf H-H.skf
"""


d = 0.9575
t = np.pi / 180 * 104.51
water = Atoms('H2O',
              positions=[(d, 0, 0),
                         (d * np.cos(t), d * np.sin(t), 0),
                         (0, 0, 0)],
              calculator=DCDFTBMD(dcdftb_keywords, dcdftb_skinfo))


dyn = BFGS(water)
dyn.run(fmax=0.0001)

with open('final.xyz', 'w') as f:
    write_xyz(f, water)

```
