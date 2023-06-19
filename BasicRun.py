#!/usr/bin/env python3

from pyIMSRG import *
from math import pi


emax = 3       # maximum number of oscillator quanta in the model space
ref = 'O17'     # reference used for normal ordering
val = 'psd-shell' # valence space
A = 17

hw = 20    # harmonic oscillator basis frequency
core_generator = 'white'   # definition of generator eta for decoupling the core (could also use 'white')
valence_generator = 'shell-model-white'  # definition of generator for decoupling the valence space (could also use 'shell-model-white'
smax_core = 5       # limit of integration in flow parameter s for first stage of decoupling
smax_valence = 10   # limit of s for second stage of decoupling


### Otherwise, we use the Minnesota NN potential
LECs = 'Minnesota'
f3b = 'none'


##########################################################################
###  END PARAMETER SETTING. BEGIN ACTUALLY DOING STUFF ##################
##########################################################################


### Create an instance of the ModelSpace class
ms = ModelSpace(emax,ref,val)
ms.SetHbarOmega(hw)

### the ReadWrite object handles reading and writing of files
rw = ReadWrite()

rank_j, parity, rank_Tz, particle_rank = 0,0,0,2
if f3b != 'none':
   particle_rank = 3

### Create an instance of the Operator class, representing the Hamiltonian
H0 = Operator(ms,rank_j, rank_Tz, parity, particle_rank)


VPT = OperatorFromString( ms, "VPT_a0_0_0")

### Either generate the matrix elements of the Minnesota potential, or read in matrix elements from file
if LECs == 'Minnesota':
    H0 += OperatorFromString(ms,'VMinnesota')

else:
 rw.ReadBareTBME_Darmstadt(f2b,H,f2e1,f2e2,f2l)
 if f3b != 'none':
     rw.Read_Darmstadt_3body(f3b,H,f3e1,f3e2,f3e3)


### Add the relative kinetic energy, so H = Trel + V
H0 += OperatorFromString(ms,'Trel')


### Create an instance of the HartreeFock class, used for solving the Hartree-Fock equations
hf = HartreeFock(H0)
hf.Solve()
hf.PrintSPEandWF()
VPT = hf.TransformToHFBasis(VPT)
 
### Do normal ordering with respect to the HF basis
H0NO = hf.GetNormalOrderedH(2)
ueta = UnitTest(ms).RandomOp(ms, 0, 0, 0, 1, -1 )
uetapt=hf.TransformToHFBasis(UnitTest(ms).RandomOp(ms, 0,0,1,2,1))
print('rank of ueta=',uetapt.GetParticleRank())

RMS = 1.2*pow(A,1/3)
Schiff = ISDipoleOp(ms,3,1,RMS)/10
#print('parity=',Schiff.GetParity())
#print('particlerank=',Schiff.GetParticleRank())
#print('Jrank=',Schiff.GetJRank())

Schiffpp = Operator(ms,1,0,0,2)

delta_s = 1
H0 = H0NO 
VPTs = uetapt #VPT 

eta = Operator(ms, 0, 0, 0, 2)
etapv = Operator(ms, 0, 0, 1, 2)
eta.SetAntiHermitian()
etapv.SetAntiHermitian()
generator =  GeneratorPV()

for s in range(0,10,delta_s):
  generator.Update(H0,VPTs,eta,etapv) 
  print('OK')
  print('norm of etapt1body=',etapv.OneBodyNorm())
  print('norm of etapt2body=',etapv.TwoBodyNorm())
  H0 +=  delta_s*Commutator.Commutator(eta,H0)
  VPTs +=  delta_s*Commutator.Commutator(etapv,H0) + delta_s*Commutator.Commutator(eta,VPTs)
  Schiff +=  delta_s*Commutator.Commutator(eta,Schiff)
#  Schiffpp +=  delta_s*Commutator.Commutator(etapv,Schiff) + delta_s*Commutator.Commutator(eta,Schiffpp)
#  Schiffpp += delta_s*Commutator.Commutator(eta,Schiffpp)
  Schiffpp += delta_s*Commutator.Commutator(etapv,Schiff)
#  Schiffpp += delta_s*Commutator.Commutator(etapv,Schiff)
  Schiffpp += delta_s*Commutator.Commutator(eta,Schiffpp)
  print('parity=',Schiffpp.GetParity())
  print('particlerank=',Schiffpp.GetParticleRank())
  print('Jrank=',Schiffpp.GetJRank())
  Schiffpp.PrintOneBody()
print('OK')

#valence_fnameH0 = 'output/{}_{}_{}_e{}_hw{}_H0.snt'.format(val,ref,LECs,emax,hw)
#rw.WriteTokyo( H0, valence_fnameH0 ,'H0')
#valence_fnameVPT = 'output/{}_{}_{}_e{}_hw{}_VPT.snt'.format(val,ref,LECs,emax,hw)
#rw.WriteTensorTokyo(valence_fnameVPT ,Schiffpp)
#rw.WriteOperatorHuman(VPTs,'VPT.op')
#rw.WriteTokyo(H0,'H0.snt','')
#rw.WriteTokyo(Schiffpp,'Schiffpp.snt','')
### Finally, print out profiling information so we know why this took so dang long to run...
#prof = IMSRGProfiler()
#prof.PrintAll()




