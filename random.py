from pyIMSRG import *
emax = 2    # maximum number of oscillator quanta in the model space
ref = 'O17'     # reference used for normal ordering
val = 'psd-shell' # valence space
#ms, jrank, tz, parity, particle_rank, hermitian
ms = ModelSpace(emax,ref,val)
Xs = UnitTest(ms).RandomOp(ms, 0, 0, 0, 2, -1 ) #etapv
Yt = UnitTest(ms).RandomOp(ms, 1, 0, 0, 2, 1 )  #shiff
#Yt = UnitTest(ms).RandomOp(ms, 0, 0, 0, 2, 1)
#Yt = 1*Ys
#Yt.MakeReduced()
Zt = Operator(ms, 1, 0, 0, 2)
#Zt = Operator(ms, 1, 0, 1, 2)
#Commutator.comm111st(Xs,Yt,Zt)
#Zt=Commutator.Commutator(Xs,Yt)
#Commutator.comm121st(Xs,Yt,Zt)
#Zt.MakeNotReduced()
Commutator.comm122st(Xs,Yt,Zt)
#Commutator.comm222_pp_hh_221st(Xs,Yt,Zt)
#Commutator.comm222_phst(Xs,Yt,Zt)

#Zss.PrintOneBody()
#Xs.PrintOneBody()
#Ys.PrintOneBody()
print('space')
Zt.PrintOneBody()
Zt.PrintTwoBody()
print('space')
#Zs.PrintOneBody()
#Zt.PrintOneBody()
#print(Zs.GetParity())
#print(Zt.GetParity())
#print(Zss.GetParity())

#if (ZsOneBody() == Zt.OneBody()):
# print(True)




#UnitTest(ms).TestCommutators()
