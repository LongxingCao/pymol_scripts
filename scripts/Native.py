import pymol
from pymol import cmd

import sys,os


def nat(pdbid,chainid):
	
	cmd.view('v', 'store');
	
	object = cmd.get_names()[0]
		
	cmd.fetch(pdbid)
	cmd.select("design_nat","%s and chain %s and not hydro"%(object,"B"))
	cmd.select("native_nat","%s and chain %s and not hydro"%(pdbid,chainid))
	
	cmd.select("other_nat","%s and not chain %s"%(pdbid,chainid))
	cmd.hide("everything","other")
	cmd.hide('(resn HOH)')
	
	cmd.super("design_nat","native_nat")
	
	cmd.select("none")
	cmd.orient(object)
	
	cmd.hide("lines","all");
	cmd.show("sticks","native_nat");
	cmd.show("cartoon","native_nat");

	cmd.system("rm %s.pdb"%(pdbid));
	
	cmd.view('v', 'recall')

cmd.extend("nat",nat)
	
def na():
	
	cmd.view('v', 'store');
	
	
	object = cmd.get_names()[0]

	if object[0] == 'd':
		pdbid = object[1:5]
		chainid = object[5:6]
	else:
		pdbid = object[0:4]
		chainid = object[4:5]
	
	cmd.fetch(pdbid)
	cmd.select("design_na","%s and chain %s and not hydro"%(object,"B"))
	cmd.select("native_na","%s and chain %s and not hydro"%(pdbid,chainid))
	
	cmd.select("other_na","%s and not chain %s"%(pdbid,chainid))
	cmd.hide("everything","other")
	cmd.hide('(resn HOH)')
	
	cmd.super("native_na","design_na")
	
	cmd.hide("lines","all");
	cmd.show("sticks","native_na");
	cmd.show("cartoon","native_na");
	
	cmd.select("none")
	cmd.orient(object)
	
	pdbid = pdbid.lower()
	cmd.system("rm %s.pdb"%(pdbid));
	
	cmd.view('v', 'recall')

cmd.extend("na",na)

def info():

	object = cmd.get_names()[0]

	if object[0] == 'd':
		pdbid = object[1:5]
		chainid = object[5:6]
	else:
		pdbid = object[0:4]
		chainid = object[4:5]
	
	pdbid = pdbid.lower()
	cmd.system("open http://pdb.org/pdb/explore/explore.do?structureId=%s"%(pdbid))

cmd.extend("info",info)
