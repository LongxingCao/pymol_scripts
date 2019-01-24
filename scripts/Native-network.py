import pymol
from pymol import cmd

import sys,os


def nat(pdbid,chainid):
	
	cmd.view('v', 'store');
	
	object = cmd.get_names()[0]
		
	cmd.fetch(pdbid)
	cmd.select("design_nat","%s and chain %s"%(object,"B"))
	cmd.select("native_nat","%s and chain %s"%(pdbid,chainid))
	
	cmd.select("other_nat","%s and not chain %s"%(pdbid,chainid))
	cmd.hide("everything","other")
	cmd.hide('(resn HOH)')
	
	cmd.super("design_nat","native_nat")
	
	cmd.select("none")
	cmd.orient(object)

	cmd.system("rm %s.pdb"%(pdbid));
	
	cmd.view('v', 'recall')

cmd.extend("nat",nat)
	
def na():
	
	cmd.view('v', 'store');
	
	object = cmd.get_names()[0]
	
	pdbid = object[0:4]
	chainid = object[4:5]
	
	cmd.fetch(pdbid)
	cmd.select("design_na","%s and chain %s"%(object,"B"))
	cmd.select("native_na","%s and chain %s"%(pdbid,chainid))
	
	cmd.select("other_na","%s and not chain %s"%(pdbid,chainid))
	cmd.hide("everything","other")
	cmd.hide('(resn HOH)')
	
	cmd.super("design_na","native_na")
	
	cmd.select("none")
	cmd.orient(object)
	
	cmd.system("rm %s.pdb"%(pdbid));
	
	cmd.view('v', 'recall')

cmd.extend("na",na)
