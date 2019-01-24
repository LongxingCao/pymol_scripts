from pymol.cgo import *
from pymol import cmd
from pymol.vfont import plain


#this is a plugin version of the axes_cyl scripts by Dr. Robert L. Campbell (in fact I only added
#the __init__ function and englobed the resto of the code in a main function)
#As a very minor change I replaced the "ORIGIN" label for the origin with an "O"

def main():
	# create the axes object, draw axes with cylinders coloured red, green,
	#blue for X, Y and Z

	obj = [
	   CYLINDER, 0., 0., 0., 10., 0., 0., 0.2, 1.0, 1.0, 1.0, 1.0, 0.0, 0.,
	   CYLINDER, 0., 0., 0., 0., 10., 0., 0.2, 1.0, 1.0, 1.0, 0., 1.0, 0.,
	   CYLINDER, 0., 0., 0., 0., 0., 10., 0.2, 1.0, 1.0, 1.0, 0., 0.0, 1.0,

	   ]

	# add labels to axes object

	cyl_text(obj,plain,[-5.,-5.,-1],'O',0.20,axes=[[3.0,0.0,0.0],[0.0,3.0,0.0],[0.0,0.0,3.0]])
	cyl_text(obj,plain,[10.,0.,0.],'X',0.20,axes=[[3.0,0.0,0.0],[0.0,3.0,0.0],[0.0,0.0,3.0]])
	cyl_text(obj,plain,[0.,10.,0.],'Y',0.20,axes=[[3.0,0.0,0.0],[0.0,3.0,0.0],[0.0,0.0,3.0]])
	cyl_text(obj,plain,[0.,0.,10.],'Z',0.20,axes=[[3.0,0.0,0.0],[0.0,3.0,0.0],[0.0,0.0,3.0]])

	# then we load it into PyMOL
	cmd.load_cgo(obj,'axes')

def __init__(self):
	self.menuBar.addmenuitem('Plugin', 'command',
                             'Showaxes',
                             label = 'Showaxes',
	command = lambda s=self : main())


