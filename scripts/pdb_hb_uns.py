# ashworth

# this is a script for reading the source file (in pwd) for a pdb that is already loaded, and displaying the hbonds and buried unsatisfieds contained in their respective tables (if found).  Similar functionality was recently introduced into the Rosetta Pymol Plugin by Ron Jacak, and this script borrows a couple of improvements from it.  I hope to keep this script as a simple standalone (it's not meant to compete with the Rosetta Pymol Plugin).

import re, string, gzip
from pymol import cmd, cgo

class UnsAtom:
	def __init__( self, uns_info ):
		info = uns_info.split()
		self.resi = info[6]
		self.chain = info[7]
		self.atom = info[10]

def hbond(don,acc,red,green,blue):

	obj = []
	rad = 0.08
	num_dash = 4
	steps = num_dash*3+1

	xyz = []
	step = []
	for dim in range(3):
		xyz.append( don.coord[dim] )
		dis = acc.coord[dim]-don.coord[dim]
		step.append( dis/steps )

	obj = [ cgo.LINEWIDTH, 3.0, cgo.BEGIN, cgo.LINES, cgo.COLOR, red, green, blue ]

	for dash in range(num_dash): # number of dashes
		for dim in range(3): xyz[dim] += 2*step[dim]
		obj.append( cgo.VERTEX )
		obj.extend( xyz )
		for dim in range(3): xyz[dim] += step[dim]
		obj.append( cgo.VERTEX )
		obj.extend( xyz )
	obj.append( cgo.END )
	return obj

def create_hbonds( lines, name ):

	model = cmd.get_model(name)

	# Ron Jacak's cool monster regex, borrowed from the Rosetta Pymol Plugin
	hb_re = re.compile("(?:PROT|BASE) \s*[A-Z]+ \s*\d+ \s*(\d+) ([ A-Z]) \s*([A-Z0-9]+) \s*[A-Z]+ \s*\d+ \s*(\d+) ([ A-Z]) \s*([A-Z0-9]+)\s* (-?[\d\.]*)")
	hbonds = []
	for line in lines:
		match = hb_re.search(line)
		if match == None: continue
		(d_resi, d_chain, d_atom, a_resi, a_chain, a_atom, energy) = match.groups()
		energy = float(energy)

		if energy < -0.05: # ingores very weak "hydrogen bonds"

			d_addr = '/%s//%s/%s/%s' % ( name, d_chain, d_resi, d_atom )
			a_addr = '/%s//%s/%s/%s' % ( name, a_chain, a_resi, a_atom )

			d_atm = model.atom[ cmd.index(d_addr)[0][1] - 1 ]
			a_atm = model.atom[ cmd.index(a_addr)[0][1] - 1 ]

			if energy <= -0.9: colorscale = 1.0
			else: colorscale = -1 * energy + 0.1 # ratio to strong, with offset
			hbonds.extend( hbond( d_atm, a_atm, colorscale, colorscale, 0.0 ) )
	cmd.load_cgo( hbonds, 'hb_%s' % name )

def create_ds_uns( lines, name ):
	uns = { 'uns_sc':[], 'uns_bb':[] }
	for line in lines:
		type = line[3:8]
		if type == 'SCACC' or type == 'SCDON': uns['uns_sc'].append( UnsAtom(line) )
		elif type == 'BBACC' or type == 'BBDON': uns['uns_bb'].append( UnsAtom(line) )

	for type,list in uns.items():
		if list == []: continue
		selstr = string.join( [ '/%s//%s/%s/%s' % ( name, atom.chain, atom.resi, atom.atom ) for atom in list ], ' or ' )
		typename = '%s_%s' % ( type, name )
		cmd.select( typename, selstr )

	# copy unsatisfieds into a separate object (allows unique transparency, among other things)
	for type in uns:
		typename = '%s_%s' % ( type, name )
		cmd.disable( typename )
		obj = '%s_obj' % typename
		cmd.create( obj, typename )
		cmd.show( 'spheres', obj )
		cmd.set( 'sphere_scale' , '0.75', obj )
		cmd.set( 'sphere_transparency', '0.5', obj )

###########################
### begin main function ###
###########################

def pdb_hbonds_and_uns():

	# look in object list for rosetta pdb's with source files in pwd
	for name in cmd.get_names():

		source = None
		if os.path.exists( '%s.pdb' % name ): source = file( name + '.pdb', 'r' )
		elif os.path.exists( '%s.pdb.gz' % name ): source = gzip.open( '%s.pdb.gz' % name, 'r' )
		else: print 'cannot find source pdb file for', name; continue

		hb_lines = []; ds_lines = []

		hb_start = False
		for line in source:
			# hbond line collection uses Ron Jacak's method, for improved generality
			if line == '\n': hb_start = False
			if hb_start == True: hb_lines.append(line)
			if line.startswith("Loc, res, pos, pdb"): hb_start = True
			if line.startswith("DS "): ds_lines.append( line )

		if hb_lines == []: print 'no hbond lines found for %s' % name
		else:
			print 'Showing Rosetta hbonds for %s...' % name
			create_hbonds( hb_lines, name )
		if ds_lines == []: print 'no decoystats lines found for %s' % name
		else:
			print 'Showing Rosetta buried unsatisfieds for %s...' % name
			create_ds_uns( ds_lines, name )

cmd.extend('pdb_hb_uns',pdb_hbonds_and_uns)
