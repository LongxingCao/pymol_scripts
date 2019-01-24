
import pymol
from pymol import *

from rosetta_vdw import useRosettaRadii
from PackingMeasureUtils import useTempRadii
from GenUtils import zload

def viewLig(prot = 'all' , lig = 'hetatm', niram = 'cbas'):

	cmd.show('spheres', prot)
	useRosettaRadii()
	cmd.hide('everything', lig)
	cmd.show('sticks', lig)
	cmd.hide('sticks', 'resn cav')
	cmd.show('spheres', 'resn cav')
	cmd.hide('everything', 'elem h and (neighbor elem c)')
	useTempRadii('resn cav')
	util.cbas()
	cmd.color("green", "elem c and hetatm and not resn cav")
	cmd.center('hetatm and not resn cav')

def loadLigPDB(fname,name=None):

	if name is None:
		name = name = os.path.basename(fname)
	if name.endswith('.gz'):
		name = name[:-3] 
	if name.endswith('.pdb'):
		name = name[:-4] 
	if name.endswith('.'):
		name = name[:-1] 

	zload(fname,name)
	viewLig()

cmd.extend('viewLig',viewLig)
cmd.extend('loadLigPDB',loadLigPDB)
