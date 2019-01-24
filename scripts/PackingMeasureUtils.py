import pymol
from pymol import cmd
import sys,os,random
from math import floor

sys.path.append("/Users/sheffler/svn/trunk/python_scripts/pymol_scripts")
from GenUtils import zload

rainbow = ['0xCCFF00', '0xFF0000', '0xCC00FF', '0x00FF66', '0x80FF00', '0x0019FF', '0x00FF19', '0x33FF00', '0xFF0099', '0xFF004D', '0xFF00E6', '0x00FFFF', '0x0066FF', '0x8000FF', '0x00B3FF', '0xFFE500', '0x00FFB2', '0xFF4C00', '0x3300FF', '0xFF9900']


def useRosettaRadii():
	cmd.alter("element C", "vdw=2.00")
	cmd.alter("element N", "vdw=1.75")
	cmd.alter("element O", "vdw=1.55")
	cmd.alter("element H", "vdw=1.00")
	cmd.alter("element P", "vdw=1.90")
	cmd.set("sphere_scale", 1.0)
	
def expandRadii(delta=1.0, sel='all'):
	for a in cmd.get_model(sel).atom:	
		r = float(a.vdw) + float(delta)
		cmd.alter("index "+`a.index`,'vdw='+`r`)
	cmd.rebuild(sel,"spheres")

def contractRadii(delta=1.0, sel='all'):
	for a in cmd.get_model(sel).atom:	
		r = float(a.vdw) - float(delta)
		cmd.alter("index "+`a.index`,'vdw='+`r`)
	cmd.rebuild(sel,"spheres")

def useOccColors(sel="all"):	
	d = {}
	for a in cmd.get_model().atom:
		d[a.q] = True
	colors = rainbow
   # random.shuffle(colors)
   # colors *= len(d)/len(colors)+1
	for ii in range(len(d.keys())):
		cmd.color( colors[ii] ,"%s and q=%i"%(sel,d.keys()[ii]))

def useTempColors(sel="all"):
	for a in cmd.get_model(sel).atom:
		q = a.b
		c = intcolors[ int(floor(q))%len(intcolors) ]
		cmd.color( c ,"%s and resi %s and name %s"%(sel,a.resi,a.name))


def useOccRadii(sel="all"):
	for a in cmd.get_model(sel).atom:
		q = a.q
		if q >= 3:
			print "shrik radius"
			q <- 0.1
		cmd.alter("%s and resi %s and name %s"%(sel,a.resi,a.name),"vdw=%f"%(q))
	cmd.rebuild()

def useTempRadii(sel="all"):
	for ii in range(30):
		radius = "%0.1f"%(float(ii+1)/10)
		vdw = "%0.1f"%(float(ii+1)/5)
		cmd.alter(sel+" and b="+radius,"vdw="+vdw)
	cmd.rebuild()

def loadPackingPDB(file,name=None,native=None):
	"""
	usage: loadPackingPDB <file> , [<name for object>]
	loads a foo_packing.pdb file and colors it all pretty-like
	creates two selections along with the loaded object called
	NAMEcavities and NAMEprotein which are the heteratoms representing
	holes and everything else, respectively. Names can get pretty long,
	by pymol lets you do good stuff like "select NA*cav*", which will
	match a selection called NAMEISREALLYLONGcavities.
	"""

	if name is None:
		name = name = os.path.basename(file)
	if name.endswith('.gz'):
		name = name[:-3] 
	if name.endswith('.pdb'):
		name = name[:-4] 
	if name.endswith('.'):
		name = name[:-1] 
	if name.endswith("_packing"):
		name = name[:-8]

	zload(file,name)
	cmd.hide('everything',name)
	
	if native is not None:
		cmd.align(name,native)	
		cmd.zoom(native)

	useRosettaRadii()

	cavselname  = name+"cavities"
	protselname = name+"protein"

	cmd.select(cavselname, "resn CAV and b > 0.1 and %s"%(name) )
	cmd.select(protselname,"(not resn CAV) and %s"%(name) )

	useTempRadii(cavselname)
	useOccColors(cavselname)
	cmd.color("white",protselname)

	cmd.show('spheres', cavselname )
	cmd.show("cartoon",protselname)
	cmd.show("lines",protselname)

	cmd.select("none")
	cmd.delete("sele*")
	cmd.move('z',-50)

	return name

cmd.extend("loadPackingPDB",loadPackingPDB)		   
cmd.extend("useOccRadii",useOccRadii)
cmd.extend("useOccColors",useOccColors)
cmd.extend("useTempRadii",useTempRadii)
cmd.extend("useTempColors",useTempColors)
cmd.extend('useRosettaRadii', useRosettaRadii)
cmd.extend('expandRadii',expandRadii)
cmd.extend('contractRadii',contractRadii)


