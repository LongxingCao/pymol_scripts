"""Utility commands to cycle active objects in pymol.

Adds cycler types to pymol:
    pdbdircycler          - Cycles through all pdbs in a given directory.
    pdbdircyclerlite      - Cycles through all pdbs in a given directory, loading one object at a time.
    pdblistfilecycler     - Cycles through all pdbs listed in a file.
    pdblistfilecyclerlite - Cycles through all pdbs listed in a file, loading one object at a time.
    objcycler             - Cycles though all objects.

Adds cycler commands:
    set_cycler_command    - Sets command run on each cycler iteration. Use to init object representation in lite cyclers.
"""

import logging
logger = logging.getLogger("Cycler")

from pymol import cmd,viewing
import os,re

from glob import glob
from os import path

# ashworth
# minimal general classes to support convenient "list mode" behavior in pymol (left/right arrows cycle through list)
# the 'Lite' classes use the LoadDeleteCycler instead of the EnableCycler, in order that only a single pdb from the list is loaded into memory at any given time. These 'Lite' versions are preferable for large numbers of pdbs that would exceed system memory if loaded all at once.

####################################################################################################
# relates paths to object names in a way that matches the result of cmd.load
def objname(objpath):
    return re.sub( r'(\.pdb|\.pdb.gz)$', '', path.basename(objpath))

# base class cycler for enable/disable behavior, with all objects (pdbs) preloaded
class EnableCycler(object):
    def __init__(self):
        self.current_index = 0
        self.auto_zoom = False
        self.onload_command = None

    def iter(self,by=1):
        #enabled = cmd.get_names('objects',enabled_only=1)[0]

        choices = self.choices()
        l = len(choices)

        assert self.current_index < l
        next_object = (self.current_index + by) % l

        cmd.disable(objname(choices[self.current_index]))

        self.current_index = next_object
        cmd.enable(objname(choices[self.current_index]))
        if self.auto_zoom:
            cmd.zoom(objname(choices[self.current_index]))

        if self.onload_command:
            logging.debug("onload_command: %s", self.onload_command)
            cmd.do(self.onload_command)

        cmd.replace_wizard('message',choices[self.current_index])

    def choices(self):
        raise NotImplementedError("EnableCycler.choices")

# base class cycler for load/delete behavior (to be employed when there are too many pdbs to hold in memory all at once)
class LoadDeleteCycler(object):
    def __init__(self):
        self.auto_zoom = False
        self.onload_command = None

    def iter(self,by=1):
        loaded = cmd.get_names('objects')[0]
        choices = self.choices()
        l = len(choices)
        next_file = 0
        for i in range(l):
            if objname(choices[i]) == loaded:
                next_file = choices[ (i+by) % l ]
                break
        cmd.delete('all')
        if not os.path.exists(next_file):
            raise ValueError("Can not locate file: %s" % next)
        cmd.load(next_file)
        if self.auto_zoom:
            cmd.zoom()

        if self.onload_command:
            logging.debug("onload_command: %s", self.onload_command)
            cmd.do(self.onload_command)

        cmd.replace_wizard('message',next_file)

    def choices(self):
        raise NotImplementedError("EnableCycler.choices")

####################################################################################################
# cycler over all pdbs in directory
class PDBDirCycler(EnableCycler):
    def __init__(self,target_dir='.'):
        super(PDBDirCycler, self).__init__()
        self.pdbs = glob(path.join(target_dir, "*.pdb*"))
        self.loadpdbs()
    def loadpdbs(self):
        for pdb in self.pdbs:
            cmd.load(pdb)
            cmd.disable(objname(pdb))
        cmd.enable(objname(self.pdbs[0]))
    def choices(self):
        return self.pdbs

class PDBDirCyclerLite(LoadDeleteCycler):
    def __init__(self,target_dir='.'):
        super(PDBDirCyclerLite, self).__init__()
        self.pdbs = [ f for f in os.listdir(target_dir) if re.search('.pdb.*$',f) ]
        pdb = self.pdbs[0]
        cmd.load(pdb)
    def choices(self):
        return self.pdbs

####################################################################################################
# cycler over pdbs in list file
class PDBListFileCycler(EnableCycler):
    def __init__(self,list_file):
        super(PDBListFileCycler, self).__init__()
        self.pdbs = [ l.strip() for l in open(list_file) ]
        self.loadpdbs()
    def loadpdbs(self):
        for pdb in self.pdbs:
            cmd.load(pdb)
            cmd.disable(objname(pdb))
        cmd.enable(objname(self.pdbs[0]))
    def choices(self):
        return self.pdbs

class PDBListFileCyclerLite(LoadDeleteCycler):
    def __init__(self,list_file):
        super(PDBListFileCyclerLite, self).__init__()
        self.pdbs = [ l.strip() for l in open(list_file) ]
        pdb = self.pdbs[0]
        cmd.load(pdb)
    def choices(self):
        return self.pdbs

####################################################################################################
# cycler over all PyMOL objects (including new ones)
class ObjectCycler(EnableCycler):
    def __init__(self):
        super(ObjectCycler, self).__init__()
        cmd.disable('all')
        cmd.enable( cmd.get_names('objects')[0] )
    def choices(self):
        return cmd.get_names('objects')

####################################################################################################
def prev_pdb():
    viewing.cycler.iter(-1)
def next_pdb():
    viewing.cycler.iter(1)

def spawnPDBDirCycler(target_dir='.',lite=False):
    if lite: viewing.cycler = PDBDirCyclerLite(target_dir)
    else: viewing.cycler = PDBDirCycler(target_dir)
    cmd.set_key('left',prev_pdb)
    cmd.set_key('right',next_pdb)

def spawnPDBListFileCycler(list_file,lite=False):
    if lite: viewing.cycler = PDBListFileCyclerLite(list_file)
    else: viewing.cycler = PDBListFileCycler(list_file)
    cmd.set_key('left',prev_pdb)
    cmd.set_key('right',next_pdb)

def spawnObjectCycler():
    viewing.cycler = ObjectCycler()
    cmd.set_key('left',prev_pdb)
    cmd.set_key('right',next_pdb)

def setCyclerOnloadCommand(command_string):
    if command_string[0] == '"' or command_string[0] == "'":
        command_string = command_string[1:-1]
    if viewing.cycler:
        logging.debug("Setting cycler onload_command: %s", command_string)
        viewing.cycler.onload_command = command_string

cmd.extend( 'pdbdircycler', lambda dir='.': spawnPDBDirCycler('.',False) )
cmd.extend( 'pdbdircyclerlite', lambda dir='.': spawnPDBDirCycler('.',True) )
cmd.extend( 'pdblistfilecycler', lambda file: spawnPDBListFileCycler(file,False) )
cmd.extend( 'pdblistfilecyclerlite', lambda file: spawnPDBListFileCycler(file,True) )

cmd.extend( 'objcycler', spawnObjectCycler )
cmd.extend( 'set_cycler_command', setCyclerOnloadCommand )

# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
