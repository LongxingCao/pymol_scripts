from pymol import cmd,viewing
import string,os,re

# ashworth
# minimal general classes to support convenient "list mode" behavior in pymol (left/right arrows cycle through list)
# the 'Lite' classes use the LoadDeleteCycler instead of the EnableCycler, in order that only a single pdb from the list is loaded into memory at any given time. These 'Lite' versions are preferable for large numbers of pdbs that would exceed system memory if loaded all at once.

####################################################################################################
# relates paths to object names in a way that matches the result of cmd.load
def objname(path):
	return re.sub( '\.pdb$', '', path.split('/')[-1] )

# base class cycler for enable/disable behavior, with all objects (pdbs) preloaded
class EnableCycler:
	def iter(self,by=1):
		enabled = cmd.get_names('objects',enabled_only=1)[0]
		choices = self.choices()
		l = len(choices)
		next = 0
		for i in range(l):
			if objname(choices[i]) == enabled:
				next = objname(choices[ (i+by) % l ])
				break
		cmd.disable(enabled)
		cmd.enable(next)
		cmd.zoom(next)
		cmd.replace_wizard('message',next)
	def choices(self):
		print 'ERROR unimplemented base class method called'

# base class cycler for load/delete behavior (to be employed when there are too many pdbs to hold in memory all at once)
class LoadDeleteCycler:
	def iter(self,by=1):
		loaded = cmd.get_names('objects')[0]
		choices = self.choices()
		l = len(choices)
		next = 0
		for i in range(l):
			if objname(choices[i]) == loaded:
				next = choices[ (i+by) % l ]
				break
		cmd.delete('all')
		if not os.path.exists(next):
			print 'ERROR: can\'t locate file %s!' %next
			return
		cmd.load(next)
		cmd.zoom()
		cmd.replace_wizard('message',next)

####################################################################################################
# cycler over all pdbs in directory
class PDBDirCycler(EnableCycler):
	def __init__(self,dir='.'):
		self.pdbs = [ f for f in os.listdir(dir) if re.search('.pdb$',f) ]
		self.loadpdbs()
	def loadpdbs(self):
		for pdb in self.pdbs:
			cmd.load(pdb)
			cmd.disable(objname(pdb))
		cmd.enable(objname(self.pdbs[0]))
	def choices(self):
		return self.pdbs

class PDBDirCyclerLite(LoadDeleteCycler):
	def __init__(self,dir='.'):
		self.pdbs = [ f for f in os.listdir(dir) if re.search('.pdb$',f) ]
		pdb = self.pdbs[0]
		cmd.load(pdb)
	def choices(self):
		return self.pdbs

####################################################################################################
# cycler over pdbs in list file
class PDBListFileCycler(EnableCycler):
	def __init__(self,file):
		self.pdbs = [ l.strip() for l in open(file) ]
		self.loadpdbs()
	def loadpdbs(self):
		for pdb in self.pdbs:
			cmd.load(pdb)
			cmd.disable(objname(pdb))
		cmd.enable(objname(self.pdbs[0]))
	def choices(self):
		return self.pdbs

class PDBListFileCyclerLite(LoadDeleteCycler):
	def __init__(self,file):
		self.pdbs = [ l.strip() for l in open(file) ]
		pdb = self.pdbs[0]
		cmd.load(pdb)
	def choices(self):
		return self.pdbs

####################################################################################################
# cycler over all PyMOL objects (including new ones)
class ObjectCycler(EnableCycler):
	def __init__(self):
		cmd.disable('all')
		cmd.enable( cmd.get_names('objects')[0] )
	def choices(self):
		return cmd.get_names('objects')

####################################################################################################
def prev_pdb():
	viewing.cycler.iter(-1)
def next_pdb():
	viewing.cycler.iter(1)

def spawnPDBDirCycler(dir='.',lite=False):
	if lite: viewing.cycler = PDBDirCyclerLite(dir)
	else: viewing.cycler = PDBDirCycler(dir)
	cmd.set_key('left',prev_pdb)
	cmd.set_key('right',next_pdb)

def spawnPDBListFileCycler(file,lite=False):
	if lite: viewing.cycler = PDBListFileCyclerLite(file)
	else: viewing.cycler = PDBListFileCycler(file)
	cmd.set_key('left',prev_pdb)
	cmd.set_key('right',next_pdb)

def spawnObjectCycler():
	viewing.cycler = ObjectCycler()
	cmd.set_key('left',prev_pdb)
	cmd.set_key('right',next_pdb)

cmd.extend( 'pdbdircycler', lambda dir='.': spawnPDBDirCycler('.',False) )
cmd.extend( 'pdbdircyclerlite', lambda dir='.': spawnPDBDirCycler('.',True) )
cmd.extend( 'pdblistfilecycler', lambda file: spawnPDBListFileCycler(file,False) )
cmd.extend( 'pdblistfilecyclerlite', lambda file: spawnPDBListFileCycler(file,True) )
cmd.extend( 'objcycler', spawnObjectCycler )
