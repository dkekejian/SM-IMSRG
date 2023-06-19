#!/usr/bin/env python3


from pyIMSRG import *



emax = 3       # maximum number of oscillator quanta in the model space
ref = 'O17'     # reference used for normal ordering
val = 'sd-shell' # valence space
A = 17

hw = 20    # harmonic oscillator basis frequency
core_generator = 'white'   # definition of generator eta for decoupling the core (could also use 'white')
valence_generator = 'shell-model-white'  # definition of generator for decoupling the valence space (could also use 'shell-model-white'

ms = ModelSpace(emax,ref,val)
ms.SetHbarOmega(hw)

### Otherwise, we use the Minnesota NN potential
LECs = 'Minnesota'
f3b = 'none'

H0 = Operator(ms, 0, 0, 0, 2)

if LECs == 'Minnesota':
    H0 += OperatorFromString(ms,'VMinnesota')


H0 += OperatorFromString(ms,'Trel')

hf = HartreeFock(H0)
hf.Solve()
hf.PrintSPEandWF()



H0NO = hf.GetNormalOrderedH(2)


RMS = 1.2*pow(A,1/3)
Schiff = ISDipoleOp(ms,3,1,RMS)/10

Schiffpp = UnitTest(ms).RandomOp(ms, 1, 0, 0, 2, 1 )


delta_s = 1
H0 = H0NO
VPTs = hf.TransformToHFBasis(UnitTest(ms).RandomOp(ms, 0,0,1,2,-1))


#eta = Operator(ms, 0, 0, 0, 2)
#etapv = Operator(ms, 0, 0, 1, 2)
eta = UnitTest(ms).RandomOp(ms, 0, 0, 0, 2, -1 )
etapv = UnitTest(ms).RandomOp(ms, 0, 0, 1, 2, -1 )

#eta.SetAntiHermitian()
#etapv.SetAntiHermitian()
#generator =  GeneratorPV()


for s in range(0,10,delta_s):
 # generator.Update(H0,VPTs,eta,etapv)
 # H0 +=  delta_s*Commutator.Commutator(eta,H0)
 # VPTs +=  delta_s*Commutator.Commutator(etapv,H0) + delta_s*Commutator.Commutator(eta,VPTs)
  Schiff +=  delta_s*Commutator.Commutator(eta,Schiff)
 # Schiffpp +=  delta_s*Commutator.Commutator(etapv,Schiff) + delta_s*Commutator.Commutator(eta,Schiffpp)
  Schiffpp += delta_s*Commutator.Commutator(eta,Schiffpp)
print('OK')

Schiff.PrintTwoBody()
Schiffpp.PrintTwoBody()
