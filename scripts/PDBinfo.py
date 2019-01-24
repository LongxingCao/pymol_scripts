# Derrick
# This script should create pymol objects containing all residues belonging to a PODBinfo label


import re, string, gzip
from pymol import cmd, cgo

def create_PDBinfo_objs( lines, name ):

        model = cmd.get_model(name)

        PDBinfo_objects = []
        PDBinfo_dict = {}
	for line in lines:
                split_line = line.split()
		split_line = split_line[2:]
		for info_name in split_line[1:]:
			PDBinfo_dict.setdefault(info_name, []).append(split_line[0])
	for key in PDBinfo_dict:
		print key, "residues = ",PDBinfo_dict[key]
	for key in PDBinfo_dict:
		print key
		residues = ""
		for residue in PDBinfo_dict[key]:
			residues = residues + "," + str(residue)
			print residues
		cmd.select(key, "resi %s" % (residues[1:]))
		cmd.create("%s_obj" %key, "%s" %key)
		cmd.delete("%s" %key)

#		match = hb_re.search(line)
#                if match == None: continue
#                (d_resi, d_chain, d_atom, a_resi, a_chain, a_atom, energy) = match.groups()
#                energy = float(energy)
#
#                if energy < -0.05: # ingores very weak "hydrogen bonds"
#
#                        d_addr = '/%s//%s/%s/%s' % ( name, d_chain, d_resi, d_atom )
#                        a_addr = '/%s//%s/%s/%s' % ( name, a_chain, a_resi, a_atom )
#
#                        d_atm = model.atom[ cmd.index(d_addr)[0][1] - 1 ]
#                        a_atm = model.atom[ cmd.index(a_addr)[0][1] - 1 ]
#
#                        if energy <= -0.9: colorscale = 1.0
#                        else: colorscale = -1 * energy + 0.1 # ratio to strong, with offset
#                        hbonds.extend( hbond( d_atm, a_atm, colorscale, colorscale, 0.0 ) )
#        cmd.load_cgo( hbonds, 'hb_%s' % name )

###########################
### begin main function ###
###########################

def get_PDBinfo():

        # look in object list for rosetta pdb's with source files in pwd
        for name in cmd.get_names():

                source = None
                if os.path.exists( '%s.pdb' % name ): source = file( name + '.pdb', 'r' )
                elif os.path.exists( '%s.pdb.gz' % name ): source = gzip.open( '%s.pdb.gz' % name, 'r' )
                else: print 'cannot find source pdb file for', name; continue

                PDBinfo_lines = [];

                for line in source:
                        if line.startswith('REMARK PDBinfo-LABEL'): PDBinfo_lines.append(line)

                if PDBinfo_lines == []: print 'no PDBinfo lines found for %s' % name
                else:
                        print 'Showing Rosetta PDBinfo for %s...' % name
                        create_PDBinfo_objs( PDBinfo_lines, name )

cmd.extend('PDBinfo',get_PDBinfo)			
