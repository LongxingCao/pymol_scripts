
import pymol
from pymol import cmd

# Select protons which can participate in Hbonds (note: Rosetta names expected)
def selectPolarProtons( name = None, hideother = False ):
	if name is None: name = "polar_protons"
	# protein
	cmd.select( name, 'not hetatm and name h') # backbone hn
	# iterative reselection to avoid large single selection (PyMOL can choke)
	cmd.select( name, '%s or (resn arg and name 1hh1+2hh1+1hh2+2hh2+he)' % name )
	cmd.select( name, '%s or (resn asn and name 1hd2+2hd2)' % name )
	cmd.select( name, '%s or (resn gln and name 1he2+2he2)' % name )
	cmd.select( name, '%s or (resn his and name he2+hd1)' % name )
	cmd.select( name, '%s or (resn lys and name 1hz+2hz+3hz)' % name )
	cmd.select( name, '%s or (resn ser and name hg)' % name )
	cmd.select( name, '%s or (resn thr and name hg1)' % name )
	cmd.select( name, '%s or (resn trp and name he1)' % name )
	cmd.select( name, '%s or (resn tyr and name hh)' % name )
	cmd.select( name, '%s or (resn a and name 1h6+2h6)' % name ) # DNA
	cmd.select( name, '%s or (resn c and name 1h4+2h4)' % name ) # DNA
	cmd.select( name, '%s or (resn g and name h1+1h2+2h2)' % name ) # DNA
	cmd.select( name, '%s or (resn t and name h3)' % name ) # DNA
	cmd.disable( name ) # hide selection visuals
	if hideother: cmd.hide( 'everything', 'e. h and not %s' % name )
	return name

# Select protons which cannot participate in Hbonds
def selectApolarProtons( name = None, hideother = False ):
	if name is None: name = "apolar_protons"
	polh = selectPolarProtons()
	cmd.select( name, 'symbol h and not %s' % polh )
#	cmd.delete( polh ) # but why not keep it?
	if hideother: cmd.hide( 'everything', 'e. h and not %s' % name )
	return name

# Color by atom, with support for alternate C/H colors
def colorCPK(selection="all",carbon='grey',hydrogen='white'):
	cmd.color(carbon,"(%s) and symbol C"%(selection))
	cmd.color(hydrogen,"(%s) and symbol H"%(selection))
	cmd.color('blue',"(%s) and symbol N"%(selection))
	cmd.color('red',"(%s) and symbol O"%(selection))
	cmd.color('yellow',"(%s) and symbol S"%(selection))

cmd.extend("selectPolarProtons",selectPolarProtons)
cmd.extend("selectApolarProtons",selectApolarProtons)
cmd.extend("colorCPK",colorCPK)

