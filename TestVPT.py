#!/usr/bin/env python3

from pyIMSRG import *
from math import pi

emax=3
hw=20
ref='O16'
val='sd-shell'

ms=ModelSpace(emax,ref,val)
ms.SetHbarOmega(hw)


#### Syntax for this operator is
#### VPT_alpha0_alpha1_alpha2
mpi = (2*M_PION_CHARGED + M_PION_NEUTRAL)/3
a0 = mpi/(8*pi*M_NUCLEON)  #g0=1
a1 = -a0/2
a2 = a0


VPT = OperatorFromString( ms, "VPT_{}_0_0".format(a0))
rw = ReadWrite()

VPT.PrintTwoBody()

# Output format for two body matrix elements is
### J_bra parity_bra Tz_bra   J_ket parity_ket Tz_ket   a  b  c  d  <ab|VPT|cd>
rw.WriteOperatorHuman(VPT,'VMTg0_1.op')

