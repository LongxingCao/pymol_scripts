from pymol import cmd
import sys,os,glob,random


POINTER  = -1
LOADLIST = []



def loadTag( tag , directory , name=None, native=None, increment=True, extra="" ):
	global POINTER,LOADLIST
	if name is None:
		name = tag
	if tag.endswith("_packing.pdb"):
		if increment:
			LOADLIST.append( (tag,directory,name,native) )
			POINTER = len(LOADLIST)-1
			print POINTER
		if native is None:
			cmd.do('loadPackingPDB '+directory+'/'+tag+","+name )
		else:
			cmd.do('loadPackingPDB '+directory+'/'+tag+","+name+','+native)
	elif tag.endswith('.pdb'):
		cmd.load(directory+'/'+tag)
	else:
		pattern = directory+'/'+tag+'*.pdb'
		print pattern
		g = glob.glob( pattern )
		#print g
		cmd.load( g[0], name )
		cmd.show('cartoon')
	print "DO:",extra #TODO FIXME
	cmd.do(extra)



def loadFromGlob(pattern, name=None, native=None, delete=True,extra="", pickrandom=False):
	if delete:
		cmd.delete('all')
	if not ( pattern.endswith(".pdb") or pattern.endswith(".pdb.gz") ):
		pattern += "*.pdb"
		pattern += "*.pdb.gz"
	print pattern
	g = glob.glob( pattern )
	if len(g) < 1:
		print "CAN'T FIND ANY FILES MATCHING:",pattern
		return
	if pickrandom:
		random.shuffle(g)
	#print g
	directory = os.path.dirname(g[0])
	tag = os.path.basename(g[0])
	print directory,tag
	loadTag( tag , directory, name, native, extra=extra )

def loadprev():
	global POINTER, LOADLIST
	if POINTER < 1:
		print "loadprev CAN'T GO BACK ANY FURTHER"
		return
	cmd.delete( LOADLIST[POINTER][2] ) # delete name
	POINTER -= 1
	tag,directory,name,native = LOADLIST[POINTER]
	loadTag(tag,directory,name,native,increment=False)
	print POINTER

def browseReset():
	POINTER = -1
	LOADLIST = []

def loadnext():
	global POINTER, LOADLIST
	if POINTER >= len(LOADLIST)-1:
		print "loadnext CAN'T GO FORWARD ANY FURTHER"
		return
	print "delete", LOADLIST[POINTER][2]
	cmd.delete( LOADLIST[POINTER][2] ) # delete name
	POINTER += 1
	tag,directory,name,native = LOADLIST[POINTER]
	loadTag(tag,directory,name,native,increment=False)
	print POINTER

cmd.extend("browseReset",browseReset)
cmd.extend('loadTag',loadTag)
cmd.extend('loadFromGlob',loadFromGlob)
cmd.extend('loadprev',loadprev)
cmd.extend('loadnext',loadnext)

