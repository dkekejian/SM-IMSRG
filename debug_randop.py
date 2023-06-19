#!/usr/bin/env python3


from pyIMSRG import *



emax=3
hw=16
ref = 'C12'
val = 'p-shell'

ms=  ModelSpace(emax,ref,val)
ms.SetHbarOmega(hw)


ut = UnitTest(ms)
ut.RandomOp(ms, 0, 0, 1, 2, -1)

UnitTest(ms).RandomOp(ms, 0,0,1,2,-1)

