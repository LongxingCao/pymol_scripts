
import pymol
from pymol import cmd
from pymol import stored

import sys,os

from MoleculeUtils import selectPolarProtons
from MoleculeUtils import selectApolarProtons
from MoleculeUtils import colorCPK
from PackingMeasureUtils import loadPackingPDB

def loadSurfaceInterfacePDB( file,name=None,native=None,wt=None):

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

        backbone_colorlist = [ 'forest','gold', 'violet', 'cyan',    \
                   'salmon', 'lime', 'slate', 'magenta', 'orange', 'marine', \
                   'olive', 'forest', 'firebrick', 'chocolate' ]
        curr_bb_color = 0

        carbon_colorlist = ['titanium', 'wheat', 'grey', 'pink' ]
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
                        #changing....
			colorCPK(interfacename,backbone_colorlist[curr_bb_color]) 
			curr_bb_color = curr_bb_color+1
                        if(curr_bb_color == len(backbone_colorlist)):
                                curr_bb_color = 0

                        #colorCPK(interfacename,carbon_colorlist[curr_carbon_color])
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
        #cmd.zoom("interface")

        cmd.create("design", "chain B")
        cmd.create("target", "chain A")
        cmd.show( "surface", "target" )
        cmd.show( "surface", "design" )
        cmd.set( "transparency", 0 )

        #if loaded together with wt 
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
#       cmd.show( "lines", "wt_A" ) 
#       cmd.show( "lines", "wt_B" ) 

        #cmd.show( "cartoon", "(interface) or byres(neighbor(interface)) or byres(neighbor(byres(neighbor(interface))))" )
	cmd.show( "cartoon")

        # Show interface waters as small purple spheres
        cmd.select( "interface_water", "(symbol w or resn HOH) and (interface around 8.0)")
        if cmd.count_atoms("interface_water")>0:
                # Put the waters in a separate object, so that we can scale their radii
                newwatername = name+"waters"
                cmd.create(newwatername, "interface_water")
                cmd.remove("interface_water")
                cmd.color("purple", newwatername)
                cmd.show( "spheres", newwatername )
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
        return name

cmd.extend("loadSurfaceInterfacePDB",loadSurfaceInterfacePDB)

def surfaceInterfacePDB():
        cmd.hide("lines")

        backbone_colorlist = [ 'forest','gold', 'violet', 'cyan',    \
                   'salmon', 'lime', 'slate', 'magenta', 'orange', 'marine', \
                   'olive', 'forest', 'firebrick', 'chocolate' ]
        curr_bb_color = 0

        carbon_colorlist = ['titanium', 'wheat', 'grey', 'pink' ]
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
                        #changing....
			colorCPK(interfacename,backbone_colorlist[curr_bb_color]) 
			curr_bb_color = curr_bb_color+1
                        if(curr_bb_color == len(backbone_colorlist)):
                                curr_bb_color = 0

                        #colorCPK(interfacename,carbon_colorlist[curr_carbon_color])
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

        cmd.create("design", "chain B")
        cmd.create("target", "chain A")
        cmd.show( "surface", "target" )
        cmd.show( "surface", "design" )

        cmd.set( "transparency", 0 )

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
        
        # Show polar contacts
        cmd.distance("hbonds","interfaceA","interfaceB",3.2,mode=2)
        cmd.hide("labels")
        
        # Show vacuum electrostatics
        cmd.util.protein_vacuum_esp("design", mode=2, quiet=0)
        cmd.util.protein_vacuum_esp("target", mode=2, quiet=0)

        cmd.disable("design_e_chg")
        cmd.disable("design_e_map")
        cmd.disable("design_e_pot")

        cmd.disable("target_e_chg")
        cmd.disable("target_e_map")
        cmd.disable("target_e_pot")

        cmd.disable("design")
        cmd.disable("target")

        #cmd.zoom("interface")
        cmd.remove("sele")
        cmd.select("none")


cmd.extend("surfaceInterfacePDB",surfaceInterfacePDB)


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

	backbone_colorlist = ['plutonium','wheat','green', 'yellow', 'violet', 'cyan',    \
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
			colorCPK(interfacename,backbone_colorlist[curr_bb_color])
			curr_bb_color = curr_bb_color+1
			if(curr_bb_color == len(backbone_colorlist)):
				curr_bb_color = 0

			#colorCPK(interfacename,carbon_colorlist[curr_carbon_color])
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
	#cmd.zoom("interface")

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
	return name

cmd.extend("loadInterfacePDB",loadInterfacePDB)		   

def interfacePDB():

	backbone_colorlist = ['plutonium','wheat','green', 'yellow', 'violet', 'cyan',    \
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
			colorCPK(interfacename,backbone_colorlist[curr_bb_color])
			curr_bb_color = curr_bb_color+1
			if(curr_bb_color == len(backbone_colorlist)):
				curr_bb_color = 0

			#colorCPK(interfacename,carbon_colorlist[curr_carbon_color])
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
	cmd.show( "sticks", "heavy_interface and not hydro" )
	#cmd.zoom("interface")

	cmd.show( "cartoon", "(not interface) or byres(neighbor(interface)) or byres(neighbor(byres(neighbor(interface))))" )

	# Show interface waters as small purple spheres
	cmd.select( "interface_water", "(symbol w or resn HOH) and (interface around 8.0)")
	if cmd.count_atoms("interface_water")>0:
		# Put the waters in a separate object, so that we can scale their radii
		newwatername = "interfacewater"
		cmd.create(newwatername, "interface_water")
		cmd.remove("interface_water")
		cmd.color("purple", newwatername)
		cmd.show( "spheres", newwatername )
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
	
	# Show polar contacts
	#cmd.distance("hbonds","interfaceA","interfaceB",3.2,mode=1)
	cmd.distance("hbonds","interfaceA","interfaceB",3.5,mode=2)
	#cmd.distance("hbonds","interfaceA","interfaceB",3.2,mode=3)
	
	cmd.hide("labels")
	cmd.hide("lines")
	cmd.select("none")

cmd.extend("interfacePDB",interfacePDB)



def interfaceResidues(cmpx, cA='c. A', cB='c. B', cutoff=1.0, selName="interface"):
	"""
	interfaceResidues -- finds 'interface' residues between two chains in a complex.
	
	PARAMS
		cmpx
			The complex containing cA and cB
		
		cA
			The first chain in which we search for residues at an interface
			with cB
		
		cB
			The second chain in which we search for residues at an interface
			with cA
		
		cutoff
			The difference in area OVER which residues are considered
			interface residues.  Residues whose dASA from the complex to
			a single chain is greater than this cutoff are kept.  Zero
			keeps all residues.
			
		selName
			The name of the selection to return.
			
	RETURNS
		* A selection of interface residues is created and named
			depending on what you passed into selName
		* An array of values is returned where each value is:
			( modelName, residueNumber, dASA )
			
	NOTES
		If you have two chains that are not from the same PDB that you want
		to complex together, use the create command like:
			create myComplex, pdb1WithChainA or pdb2withChainX
		then pass myComplex to this script like:
			interfaceResidues myComlpex, c. A, c. X
			
		This script calculates the area of the complex as a whole.  Then,
		it separates the two chains that you pass in through the arguments
		cA and cB, alone.  Once it has this, it calculates the difference
		and any residues ABOVE the cutoff are called interface residues.
			
	AUTHOR:
		Jason Vertrees, 2009.		
	"""
	# Save user's settings, before setting dot_solvent
	oldDS = cmd.get("dot_solvent")
	cmd.set("dot_solvent", 1)
	
	# set some string names for temporary objects/selections
	tempC, selName1 = "tempComplex", selName+"1"
	chA, chB = "chA", "chB"
	
	# operate on a new object & turn off the original
	cmd.create(tempC, cmpx)
	cmd.disable(cmpx)
	
	# remove cruft and inrrelevant chains
	cmd.remove(tempC + " and not (polymer and (%s or %s))" % (cA, cB))
	
	# get the area of the complete complex
	cmd.get_area(tempC, load_b=1)
	# copy the areas from the loaded b to the q, field.
	cmd.alter(tempC, 'q=b')
	
	# extract the two chains and calc. the new area
	# note: the q fields are copied to the new objects
	# chA and chB
	cmd.extract(chA, tempC + " and (" + cA + ")")
	cmd.extract(chB, tempC + " and (" + cB + ")")
	cmd.get_area(chA, load_b=1)
	cmd.get_area(chB, load_b=1)
	
	# update the chain-only objects w/the difference
	cmd.alter( "%s or %s" % (chA,chB), "b=b-q" )
	
	# The calculations are done.  Now, all we need to
	# do is to determine which residues are over the cutoff
	# and save them.
	stored.r, rVal, seen = [], [], []
	cmd.iterate('%s or %s' % (chA, chB), 'stored.r.append((model,resi,b))')

	cmd.enable(cmpx)
	cmd.select(selName1, 'none')
	for (model,resi,diff) in stored.r:
		key=resi+"-"+model
		if abs(diff)>=float(cutoff):
			if key in seen: continue
			else: seen.append(key)
			rVal.append( (model,resi,diff) )
			# expand the selection here; I chose to iterate over stored.r instead of
			# creating one large selection b/c if there are too many residues PyMOL
			# might crash on a very large selection.  This is pretty much guaranteed
			# not to kill PyMOL; but, it might take a little longer to run.
			cmd.select( selName1, selName1 + " or (%s and i. %s)" % (model,resi))

	# this is how you transfer a selection to another object.
	cmd.select(selName, cmpx + " in " + selName1)
	# clean up after ourselves
	cmd.delete(selName1)
	cmd.delete(chA)
	cmd.delete(chB)
	cmd.delete(tempC)
	# show the selection
	cmd.enable(selName)
	
	# reset users settings
	cmd.set("dot_solvent", oldDS)
	
	return rVal

cmd.extend("interfaceResidues", interfaceResidues)

def stix():

	backbone_colorlist = ['plutonium','wheat','green', 'yellow', 'violet', 'cyan',    \
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
			cmd.select(interfacename, "byres %s and (not_this_chain around 6.0)"%(chainname) )
			cmd.select("heavy_%s"%(interfacename), "%s and not apolar_protons"%(interfacename))
			cmd.select("interface", "interface or %s"%(interfacename) )
			cmd.select("heavy_interface", "heavy_interface or heavy_%s"%(interfacename) )
			cmd.delete("not_this_chain")

			cmd.color(backbone_colorlist[curr_bb_color], chainname)
			colorCPK(chainname,backbone_colorlist[curr_bb_color])
			curr_bb_color = curr_bb_color+1
			if(curr_bb_color == len(backbone_colorlist)):
				curr_bb_color = 0

			#colorCPK(interfacename,carbon_colorlist[curr_carbon_color])
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
	cmd.show( "sticks", "not hydro" )
	#cmd.zoom("interface")

	cmd.show( "cartoon", "(not interface) or byres(neighbor(interface)) or byres(neighbor(byres(neighbor(interface))))" )
	cmd.show( "cartoon" )
	

	# Show interface waters as small purple spheres
	cmd.select( "interface_water", "(symbol w or resn HOH) and (interface around 8.0)")
	if cmd.count_atoms("interface_water")>0:
		# Put the waters in a separate object, so that we can scale their radii
		newwatername = name+"waters"
		cmd.create(newwatername, "interface_water")
		cmd.remove("interface_water")
		cmd.color("purple", newwatername)
		cmd.show( "spheres", newwatername )
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
	
	# Show polar contacts
        #cmd.distance("hbonds","interfaceA","interfaceB",3.2,mode=1)
        cmd.distance("hbonds","interfaceA","interfaceB",4.5,mode=2)
        #cmd.distance("hbonds","interfaceA","interfaceB",3.2,mode=3)


        cmd.hide("labels")
        cmd.hide("lines")
	cmd.select("none")
	

cmd.extend("stix",stix)


