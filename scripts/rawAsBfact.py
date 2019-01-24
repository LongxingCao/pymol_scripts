#! /usr/bin/env python

from pymol import cmd
from numpy import *

 
def zero_residues(sel1,offset=0,chains=0):
	offset = int(offset)

	# variable to store the offset
	stored.first = None
	# get the names of the proteins in the selection

	names = ['(model %s and (%s))' % (p, sel1)
					for p in cmd.get_object_list('(' + sel1 + ')')]

	if int (chains):
			names = ['(%s and chain %s)' % (p, chain)
							for p in names
							for chain in cmd.get_chains(p)]

	# for each name shown
	for p in names:
			# get this offset
			ok = cmd.iterate("first %s and polymer and n. CA" % p,"stored.first=resv")
			# don't waste time if we don't have to
			if not ok or stored.first == offset:
					continue;
			# reassign the residue numbers
			#cmd.alter("%s" % p, "resi=str(int(resi)-%s)" % str(int(stored.first)-offset))
			# update pymol

	cmd.rebuild()
 

def drab( nlim=-2.0, plim=2.0 ):
#zero the indexes
	#zero_residues("all", 1, 1)
#Read the values in
	fileName=cmd.get_names()[0];
	print fileName
	fileName=fileName+".dat"
	cmd.alter("all", "b=0.0")
	#file = open(fileName, 'r')
	#table = [row.strip().split('\t') for row in file]
	bVals = []
	#print table
	for line in file(fileName):
		line=line.strip().split()
		if ( line[0][0] != "#" ):
			selection="resi %s" % line[0]
			#Normalize the values
			bval=(float(line[1])-nlim)/(plim-nlim)
			#print bval
			cmd.alter( selection, "b=%f"%bval )
			print ( selection, "b=%f"%bval )
			bVals.append([selection, bval])
	cmd.hide("all")
	cmd.show("cartoon")
	cmd.spectrum("b", "blue_white_red", "all", "0.0", "1.0" , "1")
	cmd.ramp_new("rawFac_ramp", cmd.get_names()[0], [nlim, plim], color="[blue, white, red ]")
	cmd.recolor()
