from pymol.cgo import *
from pymol import cmd
from pymol.vfont import plain

# create the axes object, draw axes with cylinders coloured red, green,
#blue for X, Y and Z

obj = [
   BEGIN, LINES,
   COLOR, 1.0, 0.0, 0.0, 
   VERTEX, 0.0, 0.0, 0.0,
   VERTEX, 20.0, 0.0, 0.0,
   COLOR, 0.0, 1.0, 0.0, 
   VERTEX, 0.0, 0.0, 0.0,
   VERTEX, 0.0, 20.0, 0.0, 
   COLOR, 0.0, 0.0, 1.0, 
   VERTEX, 0.0, 0.0, 0.0,
   VERTEX, 00, 0.0, 20.0,
 
   END
]                                                                                            


# obj = [
#    CYLINDER, 0., 0., 0., 10., 0., 0., 0.2, 1.0, 1.0, 1.0, 1.0, 0.0, 0.,
#    CYLINDER, 0., 0., 0., 0., 10., 0., 0.2, 1.0, 1.0, 1.0, 0., 1.0, 0.,
#    CYLINDER, 0., 0., 0., 0., 0., 10., 0.2, 1.0, 1.0, 1.0, 0., 0.0, 1.0,
# 
#    ]
# 
# # add labels to axes object
# 
# cyl_text(obj,plain,[-5.,-5.,-1],'Origin',0.20,axes=[[3.0,0.0,0.0],[0.0,3.0,0.0],[0.0,0.0,3.0]])
# cyl_text(obj,plain,[10.,0.,0.],'X',0.20,axes=[[3.0,0.0,0.0],[0.0,3.0,0.0],[0.0,0.0,3.0]])
# cyl_text(obj,plain,[0.,10.,0.],'Y',0.20,axes=[[3.0,0.0,0.0],[0.0,3.0,0.0],[0.0,0.0,3.0]])
# cyl_text(obj,plain,[0.,0.,10.],'Z',0.20,axes=[[3.0,0.0,0.0],[0.0,3.0,0.0],[0.0,0.0,3.0]])

# then we load it into PyMOL
cmd.load_cgo(obj,'axes')


