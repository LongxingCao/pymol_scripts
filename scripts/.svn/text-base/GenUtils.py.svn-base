import os
import pymol
from pymol import cmd

import sys
import zlib
import glob
import urllib
import gzip
import re


# JK DEBUG - Load a BAFF PDB
#def baff_load(fname,name=None):
#		
#	cmd.load(fname)
#	cmd.hide("lines")
#	cmd.show( "cartoon" )
#	cmd.show( "sticks", "not element h and chain c and resi 65" )
#	cmd.show( "sticks", "not element h and chain i and byres ( resi 65 around 6.0 )" )
#	cmd.center("chain c and resi 65")
#
#cmd.extend("baff_load",baff_load)


# Load a PDB, with support for gzipped PDBs
def zload(fname,name=None):

	if name is None:
		name = name = os.path.basename(fname)
	if name.endswith('.gz'):
		name = name[:-3] 
	if name.endswith('.pdb'):
		name = name[:-4] 
	if name.endswith('.'):
		name = name[:-1] 

	# Check if the file exists. If not 4 chars, give a warning and do nothing else
	if not os.path.isfile(fname):
		# If fname is not 4 chars and the file doesn't exist, print an error message
		if len(name) != 4:
			print "Error, file does not exist"
			return

		# Same error if the first char is not a number
		if not re.match('[0-9]', name):
			print "Error, file does not exist"
			return

		# If the length is 4 and the first char is a number, this is presumably a standard PDB ID
		print "Using PDB file "+name+" from repository..."

		# If the filename is 4 chars, goto the PDB homedir repository looking for it
		# Assemble the PDB name (with .gz) in the homedir repository
		pdbdir = os.path.expanduser("PDBs/")
		if not os.path.isdir(pdbdir): os.mkdir(pdbdir)
		fname = pdbdir + name.lower() + ".pdb.gz"

		# If it's not found there, download it and add it to the PDB homedir repository so we can continue
		if not os.path.isfile(fname):
			pdburl = "http://www.rcsb.org/pdb/files/"+name+".pdb"
			response = urllib.urlopen(pdburl)
			f = gzip.open(fname, 'wb')
			f.write( response.read() )
			f.close()
			print "Downloaded PDB "+name

	pdbstr=''

	if fname.endswith('.gz'):
		FTEXT, FHCRC, FEXTRA, FNAME, FCOMMENT = 1, 2, 4, 8, 16

		input=open(fname)
		magic=input.read(2)
		if magic!='\037\213':
			print 'Not a gzipped file'
		if ord(input.read(1))!=8:
			print 'Unknown compression method' ; sys.exit(0)
		flag=ord(input.read(1))
		input.read(4+1+1) # Discard modification time, extra flags, and OS byte.

		if flag & FEXTRA:
			# Read & discard the extra field, if present
			xlen=ord(input.read(1))
			xlen=xlen+256*ord(input.read(1))
			input.read(xlen)
		if flag & FNAME:
			# Read and discard a null-terminated string containing the filename
			while (1):
				s=input.read(1)
				if s=='\000': break
		if flag & FCOMMENT:
			# Read and discard a null-terminated string containing a comment
			while (1):
				s=input.read(1)
				if s=='\000': break
		if flag & FHCRC:
			input.read(2)                   # Read & discard the 16-bit header CRC

		decompobj=zlib.decompressobj(-zlib.MAX_WBITS)
		crcval=zlib.crc32("")
		length=0
		while (1):
			data=input.read(1024)
			if data=="": break
			pdbstr+=decompobj.decompress(data)
			crcval=zlib.crc32(pdbstr, crcval)
		pdbstr+=decompobj.flush()
		input.close()
		crcval=zlib.crc32(pdbstr, crcval)

	else:

		f=open(fname, 'r')
		pdbstr=f.read()
		f.close()
		
	cmd.read_pdbstr(pdbstr, name)

# call useRosettaRadii to reset
def hilightPolar(sel='all'):
	N = "elem N"# and (not name N)"
	O = "elem O"# and (not name O)" # want to inlude look BB atoms...
	S = "elem S"
	sel = sel + " and " + N+" or "+O+" or "+S
	cmd.alter(sel,"vdw=0.3")
	cmd.rebuild()
	cmd.show('spheres',sel)

def unhilightPolar(sel='all'):
	N = "elem N"# and (not name N)"
	O = "elem O"# and (not name O)" # want to inlude look BB atoms...
	S = "elem S"
	sel = sel + " and " + N+" or "+O+" or "+S
	cmd.alter(sel,"vdw=1.0")
	cmd.rebuild()
	cmd.hide('spheres',sel)

# Load all PDBs in the current directory
def load_all_pdbs():

	for fname in glob.glob('*.pdb'):
		cmd.load(fname)

cmd.extend("load_all_pdbs",load_all_pdbs)


cmd.extend("zload",zload)
cmd.extend('hilightPolar',hilightPolar)
cmd.extend('unhilightPolar',unhilightPolar)

