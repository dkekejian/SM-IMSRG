#!/usr/bin/env python3

from pyIMSRG import *
emax = 2         # maximum number of oscillator quanta in the model space
ref = 'O16'     # reference used for normal ordering
val = 'psd-shell' # valence space
#ms, jrank, tz, parity, particle_rank, hermitian
ms = ModelSpace(emax,ref,val)
Xs = UnitTest(ms).RandomOp(ms, 0, 0, 1, 2, 1 )
Ys = UnitTest(ms).RandomOp(ms, 0, 0, 1, 2, 1 )
Zs1 = Operator(ms, 0, 0, 0, 2)
Zs2 = Operator(ms, 0, 0, 0, 2)
Commutator.comm122ss(Xs,Ys,Zs1)
Commutator.comm122st(Xs,Ys,Zs2)
#Zs = Commutator.Commutator(Xs,Ys)
#print('parity=',Zs.GetParity())
#print('particlerank=',Zs.GetParticleRank())
#print('Jrank=',Zs.GetJRank())

#Ys.PrintOneBody()

#hf = HartreeFock(Xs)
#hf.Solve()
#Xs = hf.TransformToHFBasis(Xs)
Zs1.PrintTwoBody()
Zs2.PrintTwoBody()
#if Zs.PrintTwoBody() == Zt.PrintTwoBody():
# print(True)


#rw = ReadWrite()
#rw.WriteTokyo(Xs,'Xs.snt','')
#print(Xs.GetParity())


