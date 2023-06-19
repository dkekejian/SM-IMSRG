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

#### Example format of how to read input interaction matrix elements from file (these are not included with the code)
#f2b = 'input/chi2b_srg0800_eMax16_EMax16_hwHO020.me2j.gz'
#f2e1,f2e2,f2l = 16,16,16
#f3b = 'input/chi2b3b400cD-02cE0098_srg0800ho40C_eMax12_EMax12_hwHO020.me3j.gz'
#f3e1,f3e2,f3e3 = 12,24,12
#LECs = 'srg0800'

### Otherwise, we use the Minnesota NN potential
LECs = 'Minnesota'
f3b = 'none'

### name of file to write resulting shell model effective interaction.
### *.snt is the exension used with KSHELL

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

mpi = (2*M_PION_CHARGED + M_PION_NEUTRAL)/3
a0 = mpi/(8*pi*M_NUCLEON)  #g0=1
a1 = -a0/2
a2 = a0

VPT = OperatorFromString( ms, "VPT_a0_0_0")
#VPT = UnitTest(ms).RandomOp(ms, 0, 0, 1, 2, 1 )
#print(VPT.Norm())
#VPT.PrintOneBody()
#VPT.PrintTwoBody()
#print(VPT.TwoBodyNorm())

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
#print(VPT.Norm())
#print(VPT.OneBodyNorm())
#print(VPT.TwoBodyNorm())
#VPT.PrintTwoBody()

### Do normal ordering with respect to the HF basis
H0NO = hf.GetNormalOrderedH(2)
#VPTNO = hfvpt.GetNormalOrderedH(2)
ueta = UnitTest(ms).RandomOp(ms, 0, 0, 0, 1, -1 )
uetapt=hf.TransformToHFBasis(UnitTest(ms).RandomOp(ms, 0,0,1,2,1))
print('rank of ueta=',uetapt.GetParticleRank())

RMS = 1.2*pow(A,1/3)
Schiff = ISDipoleOp(ms,3,1,RMS)/10
#Schiff = hf.TransformToHFBasis(Schiff)
#Schiff = UnitTest(ms).RandomOp(ms, 1, 0, 1, 2, 1 )
#Schiff.PrintOneBody()
#Schiff.PrintTwoBody()
print('parity=',Schiff.GetParity())
print('particlerank=',Schiff.GetParticleRank())
print('Jrank=',Schiff.GetJRank())

Schiffpp = Operator(ms,1,0,0,2)
#Schiffpp = UnitTest(ms).RandomOp(ms, 1, 0, 0, 2, 1 ) #Operator(ms,1,0,0,1)
#Schiffpp.SetHermitian()


#print(Schiffpp.TwoBody.GetTBME_norm(6,6, 0,8,0,8))
#Schiffpp.TwoBody.SetTBME(6,6, 0,8,0,8,5)
#print(Schiffpp.TwoBody.GetTBME_norm(6,6, 0,8,0,8))

#print(Schiff.TwoBody.GetTBME_norm(40,41, 13,13,7,13))
#Schiff.TwoBody.SetTBME(40,41, 13,13,7,13,5)
#print(Schiff.TwoBody.GetTBME_norm(40,41, 13,13,7,13))

#Schiffpp.SetTwoBody(1,0,0,0,0,0,0,1,0,0,5)

#Schiffpp.PrintTwoBody()

#Schiff.PrintTwoBody()

#print('after')

#GeneratorPV().Update(H0,VPT,ueta,uetapt)
#VPT.PrintOneBody()
#VPT.PrintTwoBody()
#print('after')



delta_s = 1
H0 = H0NO 
VPTs = uetapt #VPT 
#imsrgsolverH0 = IMSRGSolver(H0NO)
#imsrgsolverH0.SetGenerator('white')
#imsrgsolverVPT = IMSRGSolver(VPT)
#imsrgsolverVPT.SetGenerator('white')

#print(VPTs.Norm())
#print(VPTs.OneBodyNorm())
#print(VPTs.TwoBodyNorm())

eta = Operator(ms, 0, 0, 0, 2)
etapv = Operator(ms, 0, 0, 1, 2)
eta.SetAntiHermitian()
etapv.SetAntiHermitian()
generator =  GeneratorPV()

#H0.PrintOneBody()

#prina( VPTs.Norm() )
for s in range(0,10,delta_s):
#  imsrgsolverH0.SetH_s(H0)
#  imsrgsolverH0.UpdateEta()
#  eta = imsrgsolverH0.Eta
#  imsrgsolverVPT.SetH_s(VPTs)
#  imsrgsolverVPT.UpdateEta()
#  etapv = imsrgsolverVPT.Eta
  print('OK')
  generator.Update(H0,VPTs,eta,etapv) 
  print('OK')
#  Generator().Update(H0,eta)
 # print('norm of etapt=',etapt.Norm())
  print('norm of etapt1body=',etapv.OneBodyNorm())
  print('norm of etapt2body=',etapv.TwoBodyNorm())
  H0 +=  delta_s*Commutator.Commutator(eta,H0)
 # VPTs +=   delta_s*Commutator.Commutator(eta,VPTs)
  VPTs +=  delta_s*Commutator.Commutator(etapv,H0) + delta_s*Commutator.Commutator(eta,VPTs)
  Schiff +=  delta_s*Commutator.Commutator(eta,Schiff)
#  Schiffpp +=  delta_s*Commutator.Commutator(etapv,Schiff) + delta_s*Commutator.Commutator(eta,Schiffpp)
  Schiffpp += delta_s*Commutator.Commutator(eta,Schiffpp)
#  Schiffpp += delta_s*Commutator.Commutator(etapv,Schiff)
  print('parity=',Schiffpp.GetParity())
  print('particlerank=',Schiffpp.GetParticleRank())
  print('Jrank=',Schiffpp.GetJRank())
 # Schiff.PrintOneBody()
  Schiffpp += delta_s*Commutator.Commutator(etapv,Schiff)
#  Schiffpp += delta_s*Commutator.Commutator(eta,Schiffpp)
  print('parity=',Schiffpp.GetParity())
  print('particlerank=',Schiffpp.GetParticleRank())
  print('Jrank=',Schiffpp.GetJRank())
  Schiffpp.PrintOneBody()
print('OK')

#Schiff.PrintOneBody()
#Schiffpp.PrintOneBody()
print('parity=',Schiffpp.GetParity())
print('particlerank=',Schiffpp.GetParticleRank())
print('Jrank=',Schiffpp.GetJRank())

#eta.PrintTwoBody()
#print(VPTs.TwoBodyNorm())
#print("VPT")
#VPTs.PrintTwoBody()


### Create an instance of the IMSRGSolver class, used for solving the IMSRG flow equations
#imsrgsolver = IMSRGSolver(H0NO)
#imsrgsolver.SetMethod('magnus')  # Solve using the Magnus formulation. Could also be 'flow_RK4'

#imsrgsolver.SetGenerator(core_generator)
#imsrgsolver.SetSmax(smax_core)

### Do the first stage of integration to decouple the core
#imsrgsolver.Solve()

### Now set the generator for the second stage to decouple the valence space
#imsrgsolver.SetGenerator(valence_generator)
#imsrgsolver.SetSmax(smax_valence)
#imsrgsolver.Solve()


### Hs is the IMSRG-evolved Hamiltonian
#Hs = imsrgsolver.GetH_s()

### Shell model codes assume the interaction is normal ordered with respect to the core
### and typically we choose a reference different from the core, so we need to re-normal-order
### with respect to the core
#Hs = Hs.UndoNormalOrdering()

#Hs = Hs.DoNormalOrderingCore()

### Write out the effective valence space interaction

#valence_fnameH0 = 'output/{}_{}_{}_e{}_hw{}_H0.snt'.format(val,ref,LECs,emax,hw)
#rw.WriteTokyo( H0, valence_fnameH0 ,'H0')
valence_fnameVPT = 'output/{}_{}_{}_e{}_hw{}_VPT.snt'.format(val,ref,LECs,emax,hw)
#rw.WriteTensorTokyo(valence_fnameVPT ,Schiffpp)
#rw.WriteOperatorHuman(VPTs,'VPT.op')
#rw.WriteTokyo(H0,'H0.snt','')
#rw.WriteTokyo(Schiffpp,'Schiffpp.snt','')
### Finally, print out profiling information so we know why this took so dang long to run...
#prof = IMSRGProfiler()
#prof.PrintAll()




