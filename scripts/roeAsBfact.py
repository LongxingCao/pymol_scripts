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
                cmd.alter("%s" % p, "resi=str(int(resi)-%s)" % str(int(stored.first)-offset))
                # update pymol
 
        cmd.rebuild()
 

def rener(inName=""):
#zero the indexes
	zero_residues("all", 1, 1)
#Read the values in
	infileName=""
	if (inName==""):
		infileName=cmd.get_names()[0];
	else:
		infileName=inName
	fileName=infileName+".pdb"
	cmd.alter("all", "b=0.0")
	file = open(fileName, 'r')
	table = [row.strip().split() for row in file]
	bVals = []
	did_start=False
	rowcount=0
	for i in table:
		if ( i[0] == "#BEGIN_POSE_ENERGIES_TABLE" ):
			did_start=True
			rowcount=0
			continue
		elif did_start:
			rowcount+=1
			if (rowcount>3):
				if (i[0] == "#END_POSE_ENERGIES_TABLE"):
                                	did_start=False
					continue
				resn=i[0].split('_')[-1]	
				selection="resi %s" % resn
				val=i[-1]
				#if (val > 1.67):
				#	val=1.67
				#elif (val < 0.08):
				#	val=0.08
				bfac=val
				bVals.append([selection, bfac])
	numres=len(bVals)
	bValsAvg = zeros(numres)
	for i in range( 0, numres):
		bValsAvg[i]=float(bVals[i][1])
#		for j in range( 0, frag_size):
#			tmpval=float(bVals[i][1])
#			if(bValsAvg[i+j] < tmpval):
#				bValsAvg[i+j] = tmpval
#			bValsAvg[i]+= (float(bVals[i+j][1])/float(frag_size))
#	bValsAvg=bValsAvg+bValsAvg[i].min()
	for i in range( 0, numres):
		print bVals[i][0], bValsAvg[i]
		bfac = ("b=%f" % bValsAvg[i])
		cmd.alter(bVals[i][0], bfac)
	#cmd.hide("all")
	#cmd.show("cartoon")
#	cmd.spectrum("b", "blue_white_red", infileName, bValsAvg.min(), bValsAvg.max(), "1")
	cmd.spectrum("b", "blue_white_red", infileName, -2.0, 0.0, "1")
