# rosetta radii for sphere representation in PyMOL
# values taken from etable.cc

from pymol import cmd

def useRosettaRadii():
	cmd.alter("element C", "vdw=2.00")
	cmd.alter("element N", "vdw=1.75")
	cmd.alter("element O", "vdw=1.55")
#	cmd.alter("element H", "vdw=1.00") # rosetta polar H
	cmd.alter("element H", "vdw=1.20") #ja rosetta nonpolar H, probably more appropriate for packing
	cmd.alter("element P", "vdw=1.90")
	cmd.alter("element S", "vdw=1.90")
	cmd.alter("element O and r. HOH","1.4") # water
	cmd.alter("element W", "vdw=0.5") # Rosetta rotameric waters (small to avoid dominating structure)

	cmd.set("sphere_scale", 1.0) # forces spheres to update (I think)

cmd.extend('useRosettaRadii', useRosettaRadii)
