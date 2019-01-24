#!/usr/bin/env python
#-----------------------------------------------------------
# -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:
#
# rosetta_tools.py
#
# This PyMOL plugin performs a number of tasks using the modeling program Rosetta.
# As of Jan 2010, the plugin can be used to visualize hydrogen bonds according to their Rosetta energy
# function energy. It colors the hbonds by their score. The plugin can also be used to display the
# largest hydrophobic patches on the surface of a protein. It uses the program QUILT to find these 
# patches.
#
# Authors:  Ron Jacak <ronj@email.unc.edu>
#
# Creation Date:  27 July 2006
# Revisions: 1.0   RJ Stripped all of the design functionality from the RosettaDesign.py plugin
#                     and made this file work.
#            1.1   RJ Added the ability to display hydrophobic patches and renamed the file to 
#                     rosetta_tools.py
#
#------------------------------------------------------------------------------------------
# This file is made available under the Rosetta Commons license.
# See http:# www.rosettacommons.org/license
# (C) 199x-2006 University of Washington
# (C) 199x-2006 University of California Santa Cruz
# (C) 199x-2006 University of California San Francisco
# (C) 199x-2006 Johns Hopkins University
# (C) 199x-2006 University of North Carolina, Chapel Hill
# (C) 199x-2006 Vanderbilt University
#
#
import re   # for breaking up the hbond data
import os   # for os.path.join
from subprocess import Popen, PIPE
from time import sleep
import math       # floor, range, among others
#import operator  # used once to sort a list of tuples on a particular element
import colorsys  # to convert hsv values to rgb
from tkMessageBox import showerror  # to get showerror method

# test import
#from pmg_tk.skins import PMGSkin
#from pmg_tk.skins.normal import Normal

try:
	import pymol
	from pymol.cgo import *
	from pymol import cmd
	from pymol import _cmd  # used by the hbond function
	import chempy
except ImportError:
	print "Redefining pymol functions for testing..."
	class pymol:
		class cmd:
			pass

#
# GLOBAL variables
#

#rj Rosetta unsatisfied-hb stats list the group that's unsatisfied and the res-type.
#rj This dictionary will convert Rosetta groups into pymol atom selection macros
#group2atoms = {}
#group2atoms = { 
#		("COO","ASP"): ["OD1", "OD2"],
#		("COO-","ASP"): ["OD1", "OD2"],
#		("COO","GLU"): ["OE1", "OE2"],
#		("COO-","GLU"): ["OE1", "OE2"],
#		("CO","ASN"): ["OD1"],
#		("NH2","ASN"): ["ND2"],
#		("CO","GLN"): ["OE1"],
#		("NH2","GLN"): ["NE2"],
#		("N","HIS"): ["ND1"],
#		("NH","HIS"): ["NE2"],
#		("OH","SER"): ["OG"],
#		("OH","THR"): ["OG1"],
#		("OH","TYR"): ["OH"],
#		("NH","TRP"): ["NE1"],
#		("NH","ARG"): ["NE"],
#		("NH2","ARG"): ["NH1", "NH2"],
#		("NH2+","ARG"): ["NH1"],
#		("NHarg","ARG"): ["NE"],
#		("NH3","LYS"): ["NZ"],
#		("NH3+","LYS"): ["NZ"],
#		("NHaro", "TRP"): ["NE1"],
#		("CO", "A"): ["O1P", "O2P"],
#		("CO", "C"): ["O1P", "O2P"],
#		("CO", "G"): ["O1P", "O2P"],
#		("CO", "T"): ["O1P", "O2P"],
#		("NH2", "C"): ["N4"],
#		("NH2", "G"): ["N4"],
#		("N", "A"): ["N6"],
#		("N", "T"): ["N6"]
#		}

def __init__( self ):
	""" init method for the entire plugin """

	#rj addcascademenu() can be used to create a cascading menu to which items can then
	#rj be added using addmenuitem().  For instance, we may want to make the bonds visualized
	#rj by color rather than thickness, which could be toggled using options in a cascading menu.
	#rj for now, make the plugin a single item.
	#rj
	self.menuBar.addcascademenu('Plugin', 'RosettaToolsCascade',
				label='Rosetta Tools',
				statusHelp='A suite of tools using Rosetta' )

	#rj addmenuitem() type can be one of command, separator, checkbutton, radiobutton, or cascade
	self.menuBar.addmenuitem('RosettaToolsCascade', 'command',
				label='display hydrogen bonds',
				statusHelp = 'Displays hydrogen bonds colored by energy',
				command = lambda s=self : RDGroup(s).createHBonds() )

	self.menuBar.addmenuitem('RosettaToolsCascade', 'command',
				label='display hydrogen bonds, interface only',
				statusHelp = 'Displays interface hydrogen bonds colored by energy',
				command = lambda s=self : RDGroup(s).createHBonds(True) )

	self.menuBar.addmenuitem('RosettaToolsCascade', 'command',
				label='display largest hydrophobic patch',
				statusHelp = 'Displays largest hydrophobic patch by color',
				command = lambda s=self : RDGroup(s).displayLargestHydrophobicPatch() )

	#self.menuBar.addmenuitem('RosettaToolsCascade', 'command',
	#			label = 'color residues by sasapack score',
	#			command = lambda s=self: RDGroup(s).ColorSasaMenu() )


#
# The following function was written by Gareth Stockwell and downloaded from
# his EBI homepage: http://www.ebi.ac.uk/~gareth/pymol/downloads/scripts/hbond.py.
#
# Parts of the function have been modified to fit with the needs of this plugin, but
# the core functionality is the same. Rather than calling this function in a loop
# making new selection each time, the code has been altered to accept atoms between
# which to draw hydrogen bonds. By not having to cmd.select new selection each iteration
# the CGO objects can be drawn much faster.
#

#from cmd import lock,unlock,_cmd # for dist counter, as needed by the hbond function
#
#-------------------------------------------------
# Draw a hydrogen bond
#---------------------------------

def hbond(name=None,a1=None,a2=None,r=1.0,g=1.0,b=0.2,weight=0.03,dash_gap=0.20,dash_length=0.10,transparency=0.4):

	"""
	DESCRIPTION

	"hbond" creates a dashed line between two selections, marked with
	an arrow.

	USAGE

	hbond
	hbond (selection1), (selection2)
	hbond name
	hbond name, (selection1), (selection2) , [, r, [, g, [, b] ] ]
	hbond r, g, b
	hbond name, r, g, b
	hbond (selection1), (selection2), r, g, b

	name = name of hbond object
	selection1, selection2 = atom selections
	r, g, b = colour
	"""

	"""
	# Deal with 'optional' arguments
	if name is not None:
		try:
			# if name is a float, it's really r.
			# this isn't *quite* true, because someone could name something '1'
			# but we'll assume that they don't do that.
			float(name)
			name,r,g,b = None,name,s1,s2
			s1,s2 = '(lb)','(rb)'

		except ValueError:
			pass
		try:
			# if s1,s2,r are floats, s1,s2,r meant to be r,g,b
			float(s1),float(s2),float(r)
			r,g,b = s1,s2,r
			s1,s2 = '(lb)','(rb)'
		except ValueError:
				pass

		if name is not None and len(name):
			if name[0]=='(' or ' ' in name or '/' in name: # 'name' is really s1.
				try:
					# if s2 is a float, it meant to be r.
					float(s2)
					r,g,b = s2,r,g
				except ValueError:
					pass
				name,s1,s2 = None,name,s1

		if name is None:
			name = 'hbond'
			try:
				cmd.lock()
				# I'm re-using dist_counter here because
				# 1) hbond() is very similar to distance()
				# 2) making a new thing that you can set is a pain and
				#    involves changes to modules/pymol/settings.py
				#    (I think .. I don't totally grok settings.py)
				cnt = _cmd.get('dist_counter') + 1.0
				_cmd.legacy_set('dist_counter','%10f' % cnt)
				name = 'hbond%02.0f' % cnt
			finally:
				cmd.unlock()
	"""

	if (a1 == None or a2 == None):
		print "a1 and a2 cannot be None. Please specify values for a1 and a2."
		return

	# Convert arguments into floating point values
	rr = float(r)
	gg = float(g)
	bb = float(b)

	# Get dash length, gap length and dash radius from PyMOL
	# settings
	#dl = float(cmd.get_setting_tuple("dash_length")[1][0])
	#gl = float(cmd.get_setting_tuple("dash_gap")[1][0])
	#dr = float(cmd.get_setting_tuple("dash_radius")[1][0])

	# added by rj - weight is a passed parameter
	dl = float(dash_length)
	gl = float(dash_gap)
	dr = float(weight)

	"""
	# Get tuple containing object and index of atoms in these
	# selections
	x1 = cmd.index(s1,1)
	x2 = cmd.index(s2,1)

	# Get number of atoms in each selection
	n1 = len(x1)
	n2 = len(x2)
	if(n1 < 1):
		print "Error: selection " + s1 + " has no atoms"
		return
	if(n2 < 1):
		print "Error: selection " + s2 + " has no atoms"
		return

	# Get objects and atom indices
	o1 = x1[0][0]
	i1 = x1[0][1]
	o2 = x2[0][0]
	i2 = x2[0][1]

	# Get ChemPy models
	m1 = cmd.get_model(o1)
	m2 = cmd.get_model(o2)

	# Get atoms
	a1 = m1.atom[i1-1]
	a2 = m2.atom[i2-1]
	"""

	# Use the atoms that were passed in
	# Get coords
	x1 = a1.coord[0]
	y1 = a1.coord[1]
	z1 = a1.coord[2]
	x2 = a2.coord[0]
	y2 = a2.coord[1]
	z2 = a2.coord[2]

	# Make some nice strings for user feedback
	#s1 = o1 + "/" + a1.chain + "/" + a1.resn + "." + a1.resi + "/" + a1.name
	#print s1 + "(" + str(x1) + "," + str(y1) + "," + str(z1) + ")"
	#s2 = o2 + "/" + a2.chain + "/" + a2.resn + "." + a2.resi + "/" + a2.name
	#print s2 + "(" + str(x2) + "," + str(y2) + "," + str(z2) + ")"

	# Calculate distances
	dx = x2 - x1
	dy = y2 - y1
	dz = z2 - z1
	d  = math.sqrt((dx*dx) + (dy*dy) + (dz*dz))
	#print "distance = " + str(d) + "A"

	# Work out how many times (dash_len + gap_len) fits into d
	dash_tot = dl + gl
	n_dash = int(math.floor(d / dash_tot))

	# Work out step lengths
	dx1 = (dl / dash_tot) * (dx / n_dash)
	dy1 = (dl / dash_tot) * (dy / n_dash)
	dz1 = (dl / dash_tot) * (dz / n_dash)
	dx2 = (dx / n_dash)
	dy2 = (dy / n_dash)
	dz2 = (dz / n_dash)


	# Empty CGO object
	obj = []

	# Generate dashes
	x = x1
	y = y1
	z = z1
	for i in range(n_dash):
		# Generate a dash cylinder
		alpha = float(transparency)
		obj.extend( [ ALPHA, alpha, CYLINDER, x, y, z, x+dx1, y+dy1, z+dz1, dr, rr, gg, bb, rr, gg, bb ] )

		# Move to start of next dash
		x = x + dx2
		y = y + dy2
		z = z + dz2

	"""
	# Add an arrow half way along
	# Calculate midpoint
	xm = (x1 + x2) / 2
	ym = (y1 + y2) / 2
	zm = (z1 + z2) / 2

	# Vector pointing along the bond
	xv = (3.5 * dr / d) * dx
	yv = (3.5 * dr / d) * dy
	zv = (3.5 * dr / d) * dz

	# Rotate step vector 90deg around X axis
	xxv = xv
	yyv = -1 * zv
	zzv = yv

	# Add lines
	obj.extend( [ CYLINDER, xm-xv-xxv, ym-yv-yyv, zm-zv-zzv, xm+xv, ym+yv, zm+zv, dr, rr, gg, bb, rr, gg, bb ] )
	obj.extend( [ CYLINDER, xm-xv+xxv, ym-yv+yyv, zm-zv+zzv, xm+xv, ym+yv, zm+zv, dr, rr, gg, bb, rr, gg, bb ] )
	"""

	# Load the object into PyMOL
	#cmd.load_cgo(obj, name)
	return obj

# Add to PyMOL API
#cmd.extend("hbond",hbond)


class RDGroup:
	""" The RDGroup class implements all of the functionality for the plugin."""

	def __init__(self, tkApp):
		""" initialization stuff for the RDGroup class """

		#rj define parent Tk application, used to diplay dialog box error messages
		self.parent = tkApp.root

		#rj flag used for determining when the parent thread can go on; used by the createHBonds method
		#self.childFinished = 0
		

	"""
	def colorSasa(self,outputpdb):
		#yl run the perl script to get sasa engergy for each residue, then set col
		sasa_command = "/usr/local/bin/sasapick.pl " + outputpdb
		commands.getoutput(sasa_command)
		color_object_name = outputpdb.split('.')[0]
		color_object_name = color_object_name.split('/')[-1]
		if(os.path.exists("/tmp/sasapick.txt") and (os.path.getsize("/tmp/sasapick.txt")!=0)):
			find_sasa = open("/tmp/sasapick.txt","r")
			sasa_lines = find_sasa.readlines()
			color_step = 1.0/20				# get the number for RGB color increase/decrease. there are total 20 selections
			R = 0.0
			G = 0.0
			B = 1.0
			for single_line in sasa_lines:
				print "single line: ", single_line
				sasa_group = single_line.split(',')
				sasa_color_group=''
				for aa_info in sasa_group:
					if(aa_info.count(' ')==1):
						if(aa_info.split(' ')[1]):
							sasa_resi = aa_info.split(' ')[0]
							sasa_chain = aa_info.split(' ')[1]
							sasa_color_group = sasa_color_group + "(resi "+ str(sasa_resi) + ")and(chain " + sasa_chain + ")and(byobj " + color_object_name + ")or"
						else:
							sasa_resi = aa_info.split(' ')[0]
							sasa_color_group = sasa_color_group + "(resi "+ str(sasa_resi) + ")and(byobj " + color_object_name + ")or"

				sasa_color_group = sasa_color_group[0:-2]
				this_color = str(R) + "," + str(G)+","+ str(B)
				color_name = str(R)
				color_select_name = "color_select"
				cmd.select(color_select_name,sasa_color_group)
				cmd.set_color(color_name,this_color)
				cmd.color(color_name,color_select_name)
				R += color_step
				B -= color_step
			#disable the selections made for color the residues
			cmd.delete(color_select_name)
		else:
			showinfo('Information','There is no sasapack score found in the pdb file.')
			return

	def makeSphereObject(self, sele_exp):
		# Return a string containing a CGO sphere object(s) on the atom of the given selection.
		# get_model returns ChemPy atom objects corresponding to the passed in selection
		list = cmd.get_model(sele_exp)
		atoms = list.atom

		r = 0.25
		try:
			x = atoms[0].coord[0]
			y = atoms[0].coord[1]
			z = atoms[0].coord[2]
		except:
			print "No coordinates available for atoms in " + sele_exp
			return None


		return [ SPHERE, float(x), float(y), float(z), float(r) ]


	def getUnsatisfiedHydrogenBondData(self, unsHbDataLines, object=None, model=None):
		# Given a list of strings with the unsatisfied hydrogen bond data, return a list of
		# tuples have the relevant information, exactly like getHydrogenBondData().

		# the line consists of (GU)(space)(3-letter-AA-code)(space)(1-5letter-space-padded explanation of
		# what group is unsatisfied)(space)(4-char-padded-int-res-num)(space)(4-char-padded-int-res-num-pdb)
		# GU ALA   BBH   16   32
		# GU ASP   BBH   17   33
		# GU TYR    OH   22   38
		# GU ARG   BBH   24   40
		#
		# dna-protein lines look as follows
		# GU   C    CO  163  502
		# GU   C    CO  163  502
		# GU   C   NH2  163  502
		# GU   A    CO  164  503
		# GU   A    CO  164  503
		# GU   A     N  164  503
		# GU   G    CO  165  504
		# GU   G IGNOR  165  504

		#
		# (This comment made more sense with the old way of printing unsatisfied hydrogen bonds. It's left
		# in for explanatory purposes only.)
		# the DS lines we can ignore. the GU are the ones we really care about. the difference is
		# that the DS lines print info on all hydrogens not in a hydrogen bond. the GU lines print
		# only the groups that have not a single hydrogen bond. since a group can contain more than
		# one hydrogen, we're interested in groups that have 0 hbonds.

		# old style of printing unsatisfied
		# GU T BBH 2
		# GU A BBH 16

		unsHbDataTuples = []
		pattObject = re.compile("GU \s*([A-Z]+) \s*(BBO|COO|CO|N|OH|BBH|NH3|NH2|NH|NHarg|NHaro|COO-|NH2\+|NH3\+|IGNOR) \s*\d+ \s*(\d+)")
		pattObject2 = re.compile("GU ([A-Z]) (BBO|COO|COO-|CO|N|OH|BBH|NH3|NH2|NH) (\d+)")
		for line in unsHbDataLines:
			if ( line.startswith("DS ") or (line == "\n") ):
				continue
			matchObject = pattObject.search(line)
			if ( matchObject == None ):
				matchObject = pattObject2.search(line)
				if ( matchObject == None ):
					print "Rosetta Tools Plugin: Unable to parse unsatisfied hydrogen bond info. Format not understood: " + line,
					continue
				(res_type, group_type, rosetta_res_num) = matchObject.groups()
				# allow for compatibility with older versions of hbond info printing
				if ( group_type == "COO-" ):
					group_type = "COO"
				res_num = self.rosettaNumToPDBNum(rosetta_res_num, object, model)
				if ( res_num == None ):
					break
			else:
				(res_type, group_type, res_num) = matchObject.groups()
				group_type = group_type.lstrip()

			unsHbDataTuples.append( (group_type, res_type, res_num) )

		return unsHbDataTuples

	def rosettaNumToPDBNum( self, rosetta_num, object, modelList=None ):
		# Translate rosetta nums to PDB nums.

		try:
			return self.unique_residues[int(rosetta_num)-1]
		except:
			self.unique_residues = []
			residues = []
			#list = cmd.get_model(object)
			for atom in modelList.atom:
				residues.append( atom.resi )
			unique_residues = []
			for resi in residues:
				if not resi in unique_residues:
					unique_residues.append(resi)
			self.unique_residues = unique_residues

		try:
			return self.unique_residues[int(rosetta_num)-1]
		except:
			print self.unique_residues
			print sys.exc_info()
			print "Rosetta Tools Plugin: Failed to translate rosetta num " + rosetta_num + " to a PDB number."
			return None
	"""

	def getHydrogenBondData(self, raw_hb_data):
		"""Given a list of strings containing the hydrogen bond information as given by Rosetta, returns a list of tuples containing the relevant info."""

		# the lines of the PDB file we're interested in look like below.
		# this is the output when -output_hbond_info is passed on the rosetta command line.
		# note: nucleic acids use 1-letter residue codes
		#
		#(1to3-letter-don-aa-code)(space)(4character-padded-int-don-res-num)(space)(4char-leftpadded-int-don-res-num-pdb)(space)(1char-don-chain)(space)(4char-atom-name)(space)(1to3-letter-acc-aa-code)(space)
		#(4character-leftpadded-int-acc-res-num)(space)(4char-leftpadded-int-acc-res-num-pdb)(space)(1char-acc-chain)(space)(4char-atom-name)(space)(6char-padded-float-hbEnergy)(space)(6char-padded-float-distance)(space)
		#(specificity:none|ALL|something else)
		# don't worry about matching the distance or specificity, though

		# don_chain don_resname don_resnum don_pdbresnum don_atomname - acc_chain acc_resname acc_resnum acc_pdbresnum acc_atomname energy weight don_nb don_sasa don_bfactor donSS - acc_nb acc_sasa acc_bfactor accSS 
		# A_LYS_   6_   6A_ NZ  - A GLU    3    3  OE2; hbE: -0.0926, weight: 0.624; 16 23.77  2.5 L - 16 21.90 15.6 L 
		# A ASN  272  272_ ND2_- A GLU    3    3  O  ; hbE: -0.1073, weight: 0.765; 22  0.00  0.0 L - 16  0.00  0.0 L 
		# A GLY    5    5A_ N   - A ASN  272  272  OD1; hbE: -0.4400, weight: 0.718; 14  0.00  0.0 L - 22  1.88  1.9 L 

		hb_data_tuples = []
		pattern_object = re.compile("([ A-Z]) ([A-Z]+) \s*\d+ \s*(\d+)([ A-Z]) \s*([A-Z0-9]+)\s*\- ([ A-Z]) ([A-Z]+) \s*\d+ \s*(\d+)([ A-Z]) \s*([A-Z0-9]+)")

		for line in raw_hb_data:
			# first split off the containing the atom/chain/resnum info from the rest of the string
			(bond_info, energy_info, extra_info) = line.split(';')
			
			match_object = pattern_object.search( bond_info )
			if ( match_object == None ):
				print line
				print "No match found"
				continue
			(don_chain, don_resn, don_resi, don_icode, don_res_atom, acc_chain, acc_resn, acc_resi, acc_icode, acc_res_atom) = match_object.groups()

			(energy_field, weight_field) = energy_info.split(", ")
			str_energy = energy_field.split(": ")[1]
			str_weight = weight_field.split(": ")[1]
			energy = float(str_energy)  # groups in match objects are 1-based
			weight = float(str_weight)

			# make default values for the chain if it wasn't in the rosetta output
			if ( don_chain == " " ):
				don_chain = "*"
			if ( acc_chain == " " ):
				acc_chain = "*"

			if ( don_icode != " " ):
				don_resi = don_resi + don_icode
			if ( acc_icode != " " ):
				acc_resi = acc_resi + acc_icode
			
			hb_data_tuples.append( (don_chain, don_resi, don_resn, don_res_atom, acc_chain, acc_resi, acc_resn, acc_res_atom, energy, weight) )

		return hb_data_tuples

	def getPatchData(self, raw_patch_data):
		"""Given a list of strings containing the nonpolar patch information as given by quilt, returns a list of tuples containing the relevant info."""

		# the lines of the output we're interested in look like below.

		# C 1
		# A 17 @ CB =1855-;P 19 @ CB =220-;P 19 @ CG =3023-;T 21 @ CG2=2495-;I 73 @ CG1=1514-;
		# I 73 @ CG2=1402-;I 73 @ CD1=773-;K 137 @ CG =56-;K 137 @ CD =2175-;L 140 @ CD1=334-;

		patch_data_tuples = []
		
		for line in raw_patch_data:
			#print line
			atom_data = line.split(";")
			for atom_datum in atom_data:
				if ( atom_datum == "" ): # skip the field after the line-terminating ';'
					continue
				#print atom_datum
				# A 17 @ CB =1855-
				# T A21 @ CG2=2495-
				pattern_object = re.compile("\s*[A-Z] ([A-Z]?)(\d*) @ ([A-Z0-9]*)\s?=.*")
				
				match_object = pattern_object.search( atom_datum )
				if ( match_object == None ):
					print line
					print "No match found"
					continue
				
				(chain, res_num, atom_name) = match_object.groups()

				patch_data_tuples.append( (res_num, chain, atom_name) )

		return patch_data_tuples


	def runRosettaForHbonds( self, pathToDatabaseFiles, pathToPDBFile ):
		""" Run a Rosetta executable which will output hbond information. The intention behind this plugin is that users 
			can load any structure and display it's Hbonds."""
		command = None
		try:
			command = "%s -database %s -s %s -ignore_unrecognized_res -no_output" % ( self.execpath, pathToDatabaseFiles, pathToPDBFile )
		except:
			self.childFinished = -1
			print "Rosetta Tools Plugin: Path to executable not set. Please specify the path to a Rosetta executable in the .rosettaplugin file."
			return

		self.commandRun = command  # save this so we can output the command to the PyMOL window
		output = Popen([command], stdout=PIPE).communicate()[0]

		#rj not sure how to check the return status and capture the output with this new Popen/subprocess python module
		#rj so I'm just going to assume that if there's output the command finished successfully.
		if ( output == None ):
			self.childFinished = -1

		self.childFinished = 1
		
		return output
		
			
	def createHBonds(self, interface_only=False):
		""" Method which gets executed when the user selected "Quick Visualize Hydrogen Bonds" or "Apply"
		on the hydrogen bonds tab of the notebook in the main Rosetta Tools GUI."""

		self.interface_hbonds_only = interface_only

		# save the view so we can restore it later
		my_view = cmd.get_view()

		# we have to pass the path to the rosetta executable and the database files to this method somehow
		# the easiest way is to have a config file in the users home directory
		init_file_name = os.environ['HOME'] + '/.rosettatoolsplugin'

		#rj Done with paths file. Now read the att file. All the att file contains is the path to the executable, right Yi?
		if ( os.path.exists(init_file_name) ):
			init_file = open( init_file_name, "r" )
			for line in init_file.xreadlines():
				if ( line.find("rosetta_executable") != -1 ):
					self.path_to_executable = line.split()[1]
				if ( line.find("rosetta_database") != -1 ):
					self.path_to_database_files = line.split()[1]
			init_file.close()

		else:
			showerror('Rosetta Tools Plugin: ERROR', "You need to specify a path to the Rosetta executable " \
						"'report_hbonds_for_plugin', and to the Rosetta database in a file named '.rosettatoolsplugin' " \
						"in your home directory.\nrosetta_executable path/to/exec\nrosetta_database path/to/database", parent=self.parent )
			return

		print "Rosetta Tools Plugin: Scoring structure to find hydrogen bonds..."

		# set values for some of the cgo object parameters
		bond_width = 0.04
		gap_length = 0.20
		dash_length = 0.10

		# if no objects are loaded, do nothing except print to the log
		if ( len(cmd.get_object_list()) == 0 ):   # cmd.get_names could also be used
			print "Rosetta Tools Plugin: ERROR 001: No structure found. Please load a structure before using this plugin."
			return

		# first check for files in the current working directory
		# if not there, skip that structure
		# at some point in the future, we should pop up a small window that asks the user to input the path to the file
		path_to_pdb_files = os.getcwd()

		# if multiple objects have been loaded, draw the hydrogen bonds for all of them
		# for each structure....
		object_list = cmd.get_names("objects", enabled_only=1)
		selection_list = cmd.get_names("all")
		for object in object_list:

			# don't try to draw Hbonds for CGO objects that are hbond - this is a consequence of calling get_names() for
			# all enabled objects
			if ( object.endswith("-hb") ):
				continue
			
			# don't redraw H-bonds if they're already there. this allows one to load a structure,
			# visualize the H-bonds in that structures, load another structure and visualize again
			# but bonds aren't redrawn for the first object loaded. 
			# only checks for the -hb. if that's been deleted, then redraw for the object
			# I could check to see if both the hb's and unsats have been deleted, but there's no good reason to.
			
			if ( self.interface_hbonds_only ):
				if ( ( object + "-iface-hb" ) in selection_list ):
					print "Rosetta Tools Plugin: Interface H-bonds for object '" + object + "' already exist. Delete CGO object to recreate bonds."
					continue
			else:
				if ( ( object + "-hb" ) in selection_list ):
					print "Rosetta Tools Plugin: Bonds for object '" + object + "' already exist. Delete CGO object to recreate bonds."
					continue

			# concatenate the path to the filename and check if the file exists
			filename = object + ".pdb"
			absolute_file_path = os.path.join( path_to_pdb_files, filename )

			if ( os.path.exists(path_to_pdb_files) != True ):
				print "Rosetta Tools Plugin: ERROR 002: Path doesn't exist. A valid path to the PDB file must be specified. Skipping this structure."
				continue
			if ( os.path.isfile(absolute_file_path) != True ):
				print "Rosetta Tools Plugin: ERROR 003: Structure file " + filename + " not found. Issue fetch PDBID and retry."
				continue

			# Need to start a new thread/process to run rosetta with the right flags and create the hb output.
			# There are two ways to parallelize in python: forks and threads.  Threads are more lightweight
			# but created threads can't run on their own - the main thread must persist or the child threads
			# are killed. Also, PyMOL doesn't appear to run plugin code in a separate thread/process. When
			# a plugins code is executed, the main PyMOL GUI window stop responding until the plugin code 
			# finishes.  What this means is that it would do no good to make the Rosetta call a thread because 
			# the main plugin code would still have to return before control is returned to the main GUI window.
			# So, I believe a fork is the way to go. Can a fork exist on its own?  Presumably, the parent could
			# exit and the child process will be orphaned. Since orphaned processes are reaped somehow by the 
			# kernel, this should be ok. The child process is what will go on to display the CGO objects since
			# the hbond information will be created by it.  The parent process will not be able to do anything
			# really.

			command = "%s -database \"%s\" -s \"%s\" -ignore_unrecognized_res -no_output" % (self.path_to_executable, self.path_to_database_files, absolute_file_path)
			print "Rosetta Tools Plugin: Running command: '" + command + "'"
			try:
				p = Popen(command, stdout=PIPE, stderr=PIPE, shell=True)
				output,error = p.communicate()
			except:
				showerror('Rosetta Tools Plugin: ERROR', 'An error occurred while trying to run Rosetta to generate hydrogen bond information.', parent=self.parent)
				return

			#rj not sure how to check the return status and capture the output with this new Popen/subprocess python module
			#rj so I'm just going to assume that if there's output the command finished successfully.

			# old way of running Rosetta
			#thread.start_new_thread( self.runRosettaForHbonds, (object, pathToDatabaseFiles, path_to_pdb_file))
			#while (self.childFinished == 0):
			#	time.sleep(.1)

			# an error occurred (-1 return value) while trying to run Rosetta to generate the hbond information
			#if (self.childFinished == -1):
			#	showerror('Rosetta Tools Plugin: ERROR', 'An error occurred while trying to run Rosetta to generate hydrogen bond information.' \
			#		'The command run was as follows:\n\n' + self.commandRun + "\n\n" + 
			#		'Please try this command manually to determine why the command failed.', parent=self.parent )
			#	return
				

			# parse the output of the command
			hbDataStart = 0
			hbData = []
			#unsHbDataStart = 0
			#unsHbData = []
			for line in output.split("\n"):
				if ( string.find( line, "DONE" ) != -1 ):
					break
				if ( hbDataStart ):
					hbData.append(line)
				#elif ( unsHbDataStart ):
				#	unsHbData.append(line)

				if ( string.find( line,"don_resname") != -1 ):
					hbDataStart = 1
				#elif (line.startswith("GU ")):
				#	unsHbDataStart = 1

			if ( len(hbData) == 0 ):
				print "Rosetta Tools Plugin: Hydrogen Bond info not found in output from Rosetta command."
			else:
				print "Rosetta Tools Plugin: The following hydrogen bonds were found:"
				print "\n".join( hbData )

				# parse up the lines into 9-tuples containing the residue numbers in the hydrogen bond, the energy and the distance
				# don_chain, don_resi, don_resn, don_res_atom, acc_chain, acc_resi, acc_resn, acc_atom, energy, distance) )
				hbDataTuples = self.getHydrogenBondData( hbData )

				hbDataTuples.sort( lambda x, y: cmp(float(x[8]),float(y[8])) )
				self.sortedHbDataTuples = hbDataTuples

				try:
					best_energy = float(self.sortedHbDataTuples[0][8])
					worst_energy = float(self.sortedHbDataTuples[-1][8]) # will be zero, since Ehb shouldn't give positive values
				except:
					print "Error in hydrogen bond data."
					print self.sortedHbDataTuples
					return


				# let's start with 50 bins for the energy values, giving us 50 different shades of color possible
				# in the future, this can be expanded to be continuous, instead of discretized.
				nbins = 50
				bin_width = (worst_energy - best_energy) / nbins
				
				# the lowest value for transparency that can be seen in the viewer is 0.4
				# the highest value possible is 1.0
				# we want the low energy bonds to have a transparency of 1.0, and higher energies to be more see thru
				trans_bin_width = (1.0 - 0.3) / nbins

				# use default saturation and value(brightness) values
				hue = 0.16667
				sat = 1.0
				value = 1.0
				
				# create colors in pymol
				# use the colorsys module to make a gradient of colors in hsv color space and then convert
				# it into rgb values.  then save them in a dictionary for lookup later.
				self.colorDict = {}
				self.transparencyDict = {}
				for binNum in range(nbins):

					# create colors using hsv scale (fractional)
					# use a gradient of red to white. going through the color spectrum just causes confusion when
					# looking at structures. keep the lower limit a light red - sat of 0.2 - not gray.
					hsv = (hue, sat * float(nbins-binNum)/nbins, value)
					#convert to rgb and append to color list
					rgb = colorsys.hsv_to_rgb(hsv[0],hsv[1],hsv[2])
					color_name = str(best_energy + binNum*bin_width)
					self.colorDict[ color_name ] = rgb

					# the dictionary looks something like below
					#{ '-1.4829': (1.0, 1.0, 0.0), 
					#  '-1.305216': (1.0, 1.0, 0.12),
					#  '-1.157146': (1.0, 1.0, 0.22),
					# note, keys are not sorted in a dictionary
					
					# we want to do something similar with the transparency parameter
					transparency = 1.0 - trans_bin_width*binNum
					self.transparencyDict[ color_name ] = transparency
				
				# need the color keys to be in sorted order, as floats. unfortunately, no way to convert a list
				# of string to a list of floats easily. needed because this list is what we'll go through to figure
				# out what to color each hbond.
				colorKeys = self.colorDict.keys()
				fColorKeys = []
				for each in colorKeys:
					fColorKeys.append(float(each))
				fColorKeys.sort()
				colorKeys = fColorKeys
				# the list looks as [-1.4829, -1.453286, -1.423672, ... ]

				totalCount = len(self.sortedHbDataTuples)

				# partition the bonds into strong, mid, and weak
				bondsObjects = []


				# create a dictionary for holding labels. Since labels have to go on the donor atom and since a 
				# donor atom can participate in more than one h-bond, we need to save the labels in a dictionary
				# keyed on the donor atom and then go through this dict after all h-bonds have been drawn.
				labelList = {}

				# before we start creating labels, tell PyMOL to color them yellow and size them nicely
				cmd.set("label_color", "yellow")
				try:
					cmd.set("label_size", "12")
				except:
					print "Rosetta Tools Plugin: Old version of PyMOL. Default sized labels being used."

				# last thing before we start going through the tuples is to print out a line to the background
				# to explain what the lines being printed after this represent (interface bonds)
				print "Rosetta Tools Plugin: The following interface spanning hydrogen bonds were found:"
				
				# now go through the data tuples and call the hbond function
				model = cmd.get_model(object)
				for tuple in self.sortedHbDataTuples:

					# don_chain, don_resi, don_resn, don_res_atom, acc_chain, acc_resi, acc_resn, acc_atom, energy, distance
					don_chain = tuple[0]
					don_resi = tuple[1]
					don_resn = tuple[2]
					don_res_atom = tuple[3]
					acc_chain = tuple[4]
					acc_resi = tuple[5]
					acc_resn = tuple[6]
					acc_res_atom = tuple[7]
					energy = tuple[8]
					distance = tuple[9]
					
					name = "hb." + don_resi + '-' + don_res_atom + '.' + acc_resi + '-' + acc_res_atom

					# figure out which atoms to draw lines between
					# tuple: (donor_chain, donor res, donor atom, acceptor_chain, acceptor res, acceptor atom, energy)

					sele1_exp = '/' + '/'.join([ object, '', don_chain, don_resi, don_res_atom ])
					sele2_exp = '/' + '/'.join([ object, '', acc_chain, acc_resi, acc_res_atom ])

					# index returns tuples containing object name and atom object index
					x1 = cmd.index( sele1_exp, 1 )
					x2 = cmd.index( sele2_exp, 1 )

					# Check to make sure we got something out of index
					if( len(x1) < 1):
						print "Rosetta Tools Plugin: Selection " + sele1_exp + " has no atoms."
						continue
					if( len(x2) < 1):
						print "Rosetta Tools Plugin: Selection " + sele2_exp + " has no atoms."
						continue

					a1 = model.atom[ x1[0][1] - 1 ]
					a2 = model.atom[ x2[0][1] - 1 ]


					# figure out the rgb code to use
					# traverse the list of keys of the color dictionary and check if the energy is greater(worse)
					# than the current key and less(better) than the next key. If it is, the value in the color
					# dictionary for that key is what we want to use
					
					if ( energy <= float(colorKeys[0]) ):
						rgbTuple = self.colorDict[ str(colorKeys[0]) ]
						transparency = self.transparencyDict[ str(colorKeys[0]) ]
					else:
						for index in range( len(colorKeys) -1):
							if ( (energy > float(colorKeys[index])) and (energy <= float(colorKeys[index+1])) ):
								rgbTuple = self.colorDict[ str(colorKeys[index]) ]
								transparency = self.transparencyDict[ str(colorKeys[index]) ]
								break

					if ( len(rgbTuple) == 0 ):
						rgbTuple = self.colorDict[ str(colorKeys[-1]) ]
						transparency = self.transparencyDict[ str(colorKeys[-1]) ]

					if ( self.interface_hbonds_only ):
						if ( don_chain == acc_chain ):
							continue
						
					# make the call to Gareth's function to construct the CGO
					try:
						cgoObject = hbond(name=name, a1=a1, a2=a2, r=rgbTuple[0], g=rgbTuple[1], b=rgbTuple[2],
											weight=bond_width, dash_gap=gap_length, dash_length=dash_length, transparency=transparency )
					except:
						print sele1_exp
						print sele2_exp
						print "Rosetta Tools Plugin: ERROR 007: Error in creating CGO object for bond."
						# to next bond
						continue
					
					# now create a label next to this bond containing the distance and Rosetta energy for this H-bond
					#distanceAndEnergyText = '"%sA E:%s"' % (str(round(distance, 2)), str(round(energy, 2)))
					energyText = '"E: %s"' % (str(round(energy, 2)))
					if (don_chain != acc_chain):
						if (labelList.has_key(sele1_exp)):
							labelList[sele1_exp] = "%s, %s" % (labelList[sele1_exp], energyText)
						else:
							labelList[sele1_exp] = energyText
					
						# print out a line containing all the info in the PDB file so people can check the energy easily
						print "Rosetta Tools Plugin: %s %s-%s %s - %s-%s %s: distance: %2.2f, hbE: %2.2f" % \
							(object, don_resn, don_resi, don_res_atom, acc_resn, acc_resi, acc_res_atom, distance, energy)

					# save into our arrays for "load"ing later
					bondsObjects.extend( cgoObject )


				# now that the CGOs are all ready, make them visible in pymol
				if ( self.interface_hbonds_only ):
					cgo_objects_display_name = object + "-iface-hb"
				else:
					cgo_objects_display_name = object + "-hb"
					
				cmd.load_cgo( bondsObjects, cgo_objects_display_name )

				# go through the list of labels and display them
				for sele_exp in labelList.keys():
					try:
						cmd.label(sele_exp, labelList[sele_exp])
					except:
						print sele_exp
						print "Rosetta Tools Plugin: ERROR 008: Error in creating bond label."
						# to next bond
					
			
			# Now do the UNSATISFIED hydrogen bonds
			# make another list to hold the unsatisfied selections
			#unsatisfied = []
			#sele_exp = None

			#if ( len(unsHbData) == 0 ):
			#	print "Rosetta Tools Plugin: Unsatisfied hydrogen Bond info not found in file " + filename + "."
			#else:
			#	self.unsHbDataTuples = self.getUnsatisfiedHydrogenBondData( unsHbData, object, model )

			#	# tuple has: (group_type, resi_name, rosetta_resi_num)
			#	# since "COO", "ASP" has two O's that should be marked, use the group2atoms array
			#	# to figure out which atoms to make sphere objects on
			#	for tuple in self.unsHbDataTuples:
			#		# I'm assuming that unsatisfied IGNOR groups should be ignored...
			#		if ( tuple[0]=="IGNOR" ):
			#			continue
			#		if ( tuple[0]=="BBO" ):
			#			atoms = ["O"]
			#		elif ( tuple[0]=="BBH" ):
			#			atoms = ["H"]
			#		else:
			#			try:
			#				atoms = group2atoms[(tuple[0],tuple[1])]
			#			except:
			#				print "Rosetta Tools Plugin: No atom translation available for group: '" + tuple[0] + "-" + tuple[1] + "'. Omitting this group."
			#				continue
			#
			#		for atom in atoms:
			#			sele_exp = '/' + '/'.join([ object, '', '*', tuple[2], atom])
			#			unsatisfied.extend( self.makeSphereObject( sele_exp ) )

			#	# now that we have all the spheres created, load them all in
			#	cmd.load_cgo( unsatisfied, object + "-unsat" )


			print "Rosetta Tools Plugin: Object " + object + " finished."

		# restore the saved view
		cmd.set_view(my_view)

	def displayLargestHydrophobicPatch(self):
		""" Method which gets executed when the user selected "display largest hydrophobic patch" in the Rosetta Tools GUI"""

		# save the view so we can restore it later
		my_view = cmd.get_view()

		# we have to pass the path to the QUILT executable to this method somehow
		# the easiest way is to have a config file in the users home directory
		self.path_to_quilt_executable = 0
		
		init_file_name = os.environ['HOME'] + '/.rosettatoolsplugin'
		if ( os.path.exists(init_file_name) ):
			init_file = open( init_file_name, "r" )
			for line in init_file.xreadlines():
				if ( line.find("quilt_executable") != -1 ):
					self.path_to_quilt_executable = line.split()[1]
			init_file.close()

		if ( self.path_to_quilt_executable == 0 ):
			showerror('Rosetta Tools Plugin: ERROR', "You need to specify a path to a QUILT executable in a file named '.rosettatoolsplugin' " \
						"in your home directory.\nquilt_executable path/to/executable", parent=self.parent )
			return

		print "Rosetta Tools Plugin: Running QUILT to find hydrophobic patches..."

		# if no objects are loaded, do nothing except print to the log
		if ( len(cmd.get_object_list()) == 0 ):   # cmd.get_names could also be used
			print "Rosetta Tools Plugin: ERROR 001: No structure found. Please load a structure before using this plugin."
			return

		# first check for files in the current working directory
		# if not there, skip that structure
		# at some point in the future, we should pop up a small window that asks the user to input the path to the file
		path_to_pdb_files = os.getcwd()

		# if multiple objects have been loaded, draw the patches for all of them
		object_list = cmd.get_names("objects", enabled_only=1)
		selection_list = cmd.get_names("all")
		for object in object_list:

			# don't try to find patches for things that are not structures - this is a consequence of calling get_names() for
			# all enabled objects
			if ( object.endswith("-hb") ):
				continue
			
			# don't recreate the patches if they're already there. this allows one to load a structure,
			# visualize the patches in that structures, load another structure and visualize again
			# but the quilt run is not performed again.
			if ( ( object + "-hpatch" ) in selection_list ):
				print "Rosetta Tools Plugin: Largest hydrophobic patch for object '" + object + "' already exists. Delete CGO object to redetermine patch."
				continue

			# concatenate the path to the filename and check if the file exists
			filename = object + ".pdb"
			absolute_file_path = os.path.join( path_to_pdb_files, filename )

			if ( os.path.exists(path_to_pdb_files) != True ):
				print "Rosetta Tools Plugin: ERROR 002: Path doesn't exist. A valid path to the PDB file must be specified. Skipping this structure."
				continue
			if ( os.path.isfile(absolute_file_path) != True ):
				print "Rosetta Tools Plugin: ERROR 003: Structure file " + filename + " not found. Issue fetch PDBID and retry."
				continue

			# Need to start a new thread/process to run rosetta with the right flags and create the hb output.
			# There are two ways to parallelize in python: forks and threads.  Threads are more lightweight
			# but created threads can't run on their own - the main thread must persist or the child threads
			# are killed. Also, PyMOL doesn't appear to run plugin code in a separate thread/process. When
			# a plugins code is executed, the main PyMOL GUI window stop responding until the plugin code 
			# finishes.  What this means is that it would do no good to make the Rosetta call a thread because 
			# the main plugin code would still have to return before control is returned to the main GUI window.
			# So, I believe a fork is the way to go. Can a fork exist on its own?  Presumably, the parent could
			# exit and the child process will be orphaned. Since orphaned processes are reaped somehow by the 
			# kernel, this should be ok. The child process is what will go on to display the CGO objects since
			# the hbond information will be created by it.  The parent process will not be able to do anything
			# really.

			command = "%s -n 252 -ep 1.0 -P CS \"%s\"" % (self.path_to_quilt_executable, absolute_file_path)
			print "Rosetta Tools Plugin: Running command: '" + command + "'"
			try:
				p = Popen(command, stdout=PIPE, stderr=PIPE, shell=True)
				output,error = p.communicate()
			except:
				showerror('Rosetta Tools Plugin: ERROR', 'An error occurred while trying to run QUILT to generate hydrophobic patch information.', parent=self.parent)
				return

			#rj not sure how to check the return status and capture the output with this new Popen/subprocess python module
			#rj so I'm just going to assume that if there's output the command finished successfully.

			# parse the output of the command
			patch_data_start = 0
			individual_patch_data = []
			all_patch_data = []
			for line in output.split("\n"):
				if ( string.find( line,"# 0") != -1 ): # start reading after we see the 1st largest patch line
					patch_data_start = 1
					continue
				if ( string.find( line, "# " ) != -1 ): # stop if we get to the 2nd largest patch
					all_patch_data.append( individual_patch_data )
					individual_patch_data = []
					continue

				if ( patch_data_start ):
					if ( string.find( line, "% " ) != -1 ): # ignore the lines that start with '%'
						continue
					else:
						individual_patch_data.append(line)

			if ( len(all_patch_data) == 0 ):
				print "Rosetta Tools Plugin: Hydrophobic patch info not found in output from quilt command."
			else:
				#print allPatchData
				print "Rosetta Tools Plugin: The following atoms are in hydrophobic patches:"
				patch_num = 1
				for patch_data in all_patch_data:
				
					print "Patch " + str(patch_num) + ": "
					# parse up the lines into 2-tuples containing the residue number and atom name
					# don_chain, don_resi, don_resn, don_res_atom, acc_chain, acc_resi, acc_resn, acc_atom, energy, distance) )
					patch_data_tuples = self.getPatchData( patch_data )

					sele_expr = ""
					
					# PyMOL can't handle long (~100-200 character) selection expressions. have to break down the really
					# long expressions into multiple selections
					line_count = 0
					
					# now go through the tuples and add them to the selection expression
					for tuple in patch_data_tuples:
					
						# we need the object name && residue and atom name for all atoms in the patch
						resnum = tuple[0]
						chain = tuple[1]
						atomname = tuple[2]
						
						if ( line_count != 0 ):
							sele_expr = object + "-hpatch" + str(patch_num) + " /" + '/'.join([ object, '', chain, resnum, atomname ])
						else:
							sele_expr = '/' + '/'.join([ object, '', chain, resnum, atomname ])
							
						#print sele_expr
						patch_sele = object + "-hpatch" + str(patch_num)

						try:
							cmd.select( patch_sele, sele_expr )
						except:
							print sele_expr
							print "Rosetta Tools Plugin: ERROR 007: Error in coloring patch."
							continue
						sele_expr = ""
						
						line_count += 1

					# only want to display the three largest hydrophobic patches - NOT all of them
					patch_num += 1
					if ( patch_num > 3 ):
						break

				cmd.show("spheres", object)
				cmd.color("sand", object)
				cmd.util.cnc("all")
				cmd.color("green", object + "-hpatch1")
				cmd.color("palegreen", object + "-hpatch2")
				cmd.color("gray80", object + "-hpatch3")
				cmd.hide("(all and hydro)")
				cmd.deselect()
			
			print "Rosetta Tools Plugin: Object " + object + " finished."

		# restore the saved view
		cmd.set_view(my_view)


#rj
#rj test main. __name__ gets set to __main__ when the module is run as a script (and not as a library).
#rj Thus, the code below can be used as a C-style test main, though with limited functionality.
#rj Since pymol is not importable as a module (as of 8/1/2006 and as far as the author knows), anything
#rj that uses pymol commands will not work.  Instead, whatever is defined for those commands in the
#rj except block at the top will be run.  Nonetheless, the code below can be used for rapid GUI debugging.
#rj Finally, this main has to be at the bottom of the file.  Otherwise, classes are undefined, and (as
#rj far as the author knows) there's no way to prototype functions/classes in Python.
#rj
if __name__ == '__main__':
	import Pmw

	from Tkinter import *
	# something has to serve as the parent application which the __init__ method of
	# HBondViewer will belong to. that's why we have to have a little Tk window with an
	# exit button pop up during testing.
	class App:
		pass

	app = App()
	app.root = Tkinter.Tk()
	Pmw.initialise(app.root)
	app.root.title('A PyMol application substitute')
	widget = HBondViewer(app)

	# the exitButton is not really necessary
	exitButton = Tkinter.Button(app.root, text = 'Exit', command = app.root.destroy)
	exitButton.pack()
	app.root.mainloop()


