
import pymol
from pymol import *

from rosetta_vdw import useRosettaRadii
from PackingMeasureUtils import useTempRadii
from GenUtils import zload

def viewCav(prot = 'all' , lig = 'hetatm', niram = 'cbas'):
	cmd.hide('everything')
	cmd.show('stick', prot)
	cmd.hide('everything', 'elem h')
	cmd.show('spheres', 'resn cav')
	cmd.hide('sticks', 'resn cav')
	useTempRadii()
	cmd.center('resn cav')

def loadCavPDB(fname,name=None):

	if name is None:
		name = name = os.path.basename(fname)
	if name.endswith('.gz'):
		name = name[:-3] 
	if name.endswith('.pdb'):
		name = name[:-4] 
	if name.endswith('.'):
		name = name[:-1] 

	zload(fname,name)
	viewCav()

cmd.extend('viewCav',viewCav)
cmd.extend('loadCavPDB',loadCavPDB)
