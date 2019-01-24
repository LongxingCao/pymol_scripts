
import pymol
from pymol import cmd

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
        cmd.zoom("interface")

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

        cmd.zoom("interface")
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
			cmd.select(interfacename, "byres %s and (not_this_chain around 4.0)"%(chainname) )
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


