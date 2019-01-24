
import pymol
from pymol import cmd

import sys,os

from MoleculeUtils import selectPolarProtons
from MoleculeUtils import selectApolarProtons
from MoleculeUtils import colorCPK
from PackingMeasureUtils import loadPackingPDB

def loadInterfacePDB(file,name=None,native=None,wt=None):

	print " Loading interface PDB %s"%(file)

	if name is None:
		name = name = os.path.basename(file)
	if name.endswith('.pdb'):
		name = name[:-4] 

	# Call Will's packing PDB loading function
	# Note: proteins will be cartoons, cavities will be spheres
	# Rosetta radii will be enabled
	name = loadPackingPDB(file,name,native)

	cavselname  = name+"cavities"
	protselname = name+"protein"
	cmd.hide("lines",protselname)

	backbone_colorlist = ['green', 'yellow', 'violet', 'cyan',    \
		   'salmon', 'lime', 'slate', 'magenta', 'orange', 'marine', \
		   'olive', 'forest', 'firebrick', 'chocolate' ]
	curr_bb_color = 0

 	carbon_colorlist = ['teal', 'wheat', 'grey', 'pink' ]
	curr_carbon_color = 0

	# Derive selections for the interface, color by chain
	cmd.select("interface", "none")
	cmd.select("heavy_interface", "none")
	selectPolarProtons("polar_protons")
	selectApolarProtons("apolar_protons")
	alphabet = list(('abcdefghijklmnopqrstuvwxyz').upper())
	for letter in alphabet:
		chainname = "chain"+letter
		cmd.select( chainname, "chain %s and not hetatm and not symbol w"%(letter) )

		# Check whether any protein atoms exist with this chain ID
		# JK Later, put in a special "non-interface" case for L/H antibody chains
		if cmd.count_atoms("chain%s"%(letter))>0:
			interfacename = "interface"+letter
			cmd.select("not_this_chain", "not hetatm and not symbol w and not %s"%(chainname) )
			cmd.select(interfacename, "byres %s and (not_this_chain around 4.0)"%(chainname) )
			cmd.select("heavy_%s"%(interfacename), "%s and not apolar_protons"%(interfacename))
			cmd.select("interface", "interface or %s"%(interfacename) )
			cmd.select("heavy_interface", "heavy_interface or heavy_%s"%(interfacename) )
			cmd.delete("not_this_chain")

			cmd.color(backbone_colorlist[curr_bb_color], chainname)
			curr_bb_color = curr_bb_color+1
			if(curr_bb_color == len(backbone_colorlist)):
				curr_bb_color = 0

			colorCPK(interfacename,carbon_colorlist[curr_carbon_color])
			curr_carbon_color = curr_carbon_color+1
			if(curr_carbon_color == len(carbon_colorlist)):
				curr_carbon_color = 0
			cmd.color("white", "%s and polar_protons"%(interfacename))

		else:
			cmd.delete(chainname)

	cmd.delete("apolar_protons")
	cmd.delete("polar_protons")

	# Show the interface in sticks, colored cpk
	#cmd.hide( "cartoon", "interface" )
	cmd.show( "sticks", "heavy_interface" )
	cmd.zoom("interface")

	cmd.show( "cartoon", "(not interface) or byres(neighbor(interface)) or byres(neighbor(byres(neighbor(interface))))" )

	# Show interface waters as small purple spheres
	cmd.select( "interface_water", "(symbol w or resn HOH) and (interface around 8.0)")
	if cmd.count_atoms("interface_water")>0:
		# Put the waters in a separate object, so that we can scale their radii
		newwatername = name+"waters"
		cmd.create(newwatername, "interface_water")
		cmd.remove("interface_water")
		cmd.color("purple", newwatername)
		cmd.show( "spheres", newwatername )
		#cmd.set( "sphere_scale", 0.5, newwatername )
		cmd.set( "sphere_scale", 0.1, newwatername )
	else:
		cmd.delete("interface_water")

	# Show interface ligands as pink sticks
	cmd.select( "interface_hetero", "(not symbol w and not resn HOH) and (hetatm and not symbol w and not resn WSS) and (interface around 4.5)")
	if cmd.count_atoms("interface_hetero")>0:
		cmd.color("pink", "interface_hetero")
		cmd.show( "sticks", "interface_hetero" )
	else:
		cmd.delete("interface_hetero")

	cmd.select("none")
   	cmd.load( wt, "wt" ) 
   	cmd.create( "wt_A", "wt and chain A" ) 
   	cmd.create( "wt_B", "wt and chain B" ) 
   	cmd.select("none") 
   	cmd.create( "des_A", "not wt and chain A" ) 
   	cmd.show( "surface", "des_A" ) 
   	cmd.create( "des_B", "not wt and chain B") 
   	cmd.show( "surface", "des_B" ) 
   	cmd.align( "wt_A", "des_A" ) 
   	cmd.align( "wt_B", "des_B" ) 
#   cmd.show( "lines", "wt_A" ) 
#   cmd.show( "lines", "wt_B" ) 
  	cmd.set( "transparency", 1 ) 
	cmd.zoom("interface") 
	return name

cmd.extend("loadInterfacePDB",loadInterfacePDB)		   




