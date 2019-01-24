import sys,os,math
from util import *

# assumes 3 identical objs called r1,r2,r3
# places r2 an r3 based on r1

axis = Vec(1,1,1).normalized()
# cen  = Vec(22.487567, 67.710289, -0.116232 )
cen  = Vec(19.855878, 64.881564, -2.624470 )
print axis,cen
rots = {"obj02":120.0,"obj03":240.0}

for n in ("obj02","obj03"):
   R = rotation_matrix(axis,rots[n])
   m = cmd.get_model("obj01")
   print m.atom[0].coord
   for i in range(len(m.atom)):
       tmp = ( R * (Vec(m.atom[i].coord)-cen) ) + cen
       m.atom[i].coord = [tmp.x,tmp.y,tmp.z]
   print m.atom[0].coord
   cmd.load_model(m,n,1,discrete=1)

