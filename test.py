
from pyIMSRG import *
from math import pi

emax=2
hw=16
ref='O16'
val='sd-shell'

ms=ModelSpace(emax,ref,val)
ms.SetHbarOmega(hw)


#### Syntax for this operator is
#### VPT_alpha0_alpha1_alpha2
mpi = (2*M_PION_CHARGED + M_PION_NEUTRAL)/3
a0 = mpi/(8*pi*M_NUCLEON)  #g0=1
a1 = -a0/2
a2 = a0a1 = -a0/2

print('ind\tn\tl\tj\ttz')
for i in ms.all_orbits:
   oi = ms.GetOrbit(i)
   print('{}\t{}\t{}\t{}\t{}'.format(i,oi.n,oi.l,oi.j2,oi.tz2))
ch=6
tbc = ms.GetTwoBodyChannel(ch)
nkets = tbc.GetNumberKets()
for iket in range(nkets):
   ket = tbc.GetKet(iket)
   print(iket,' a b = ',ket.p, ket.q)

VPT = OperatorFromString( ms, "VPT_{}_0_0".format(1))  
VPT.PrintTwoBody()
#VPT.PrintOneBody()
print(VPT.TwoBody.GetTBME_norm(6,7, 0,8,2,6))

print(HBARC)
