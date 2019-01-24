# ashworth
# useful selection groups for protein-DNA interfaces

from pymol import cmd,util

# shortcut command for pymol's "color by chains (e. c)"
def color_by_chains():
	for obj in cmd.get_names('objects'):
		util.color_chains('%s and e. c' %obj)
cmd.extend('cbce',color_by_chains)

#class DNA_selections:
#	def __init__(self,display=True):
def DNA_selections(display='all'):
	bbatoms = 'name C2\*+C3\*+C4\*+C5\*+P+O3\*+O4\*+O5\*+O1P+O2P+H1\*+1H2\*+2H2\*+H3\*+H4\*+1H5\*+2H5\*+c2\'+c3\'+c4\'+c5\'+o3\'+o4\'+o5\'+op2+op1+h1\'+1h2\'+2h2\'+h3\'+h4\'+1h5\'+2h5\''
	waters = 'n. wo6+wn7+wn6+wn4+wo4 or r. hoh'
	cmd.select('DNA', 'r. g+a+c+t+gua+ade+cyt+thy+da+dc+dg+dt+5mc',enable=0)
	cmd.select('notDNA','not DNA',enable=0)
	cmd.select('DNAbases','DNA and not %s' % bbatoms ,enable=0)
	cmd.select('DNAbb','DNA and %s' % bbatoms ,enable=0)
	cmd.select('sc_base','byres notDNA w. 7 of DNAbases',enable=0)
	cmd.select('sc_base','sc_base and not n. c+n+o',enable=0)
	cmd.select('dna_h2o','%s w. 3.6 of DNAbases' %waters ,enable=0)
	cmd.set('sphere_transparency','0.5'); cmd.color('marine','dna_h2o')
	cmd.do('selectPolarProtons')
#	color_by_chains()
	cmd.color('gray','e. c')

	cmd.select('pbb','notDNA and n. c+n+ca',enable=0)

	if display != 'none':
		cmd.label('n. c1\*+c1\' and DNA','\'%s%s(%s)\' % (chain,resi,resn)')
		cmd.set('label_color','white')

	if display == 'all':
		# display things
		cmd.show('sticks','DNAbases or sc_base')
		cmd.show('ribbon','DNAbb')
		cmd.show('cartoon','notDNA')
		cmd.show('spheres','dna_h2o')
		cmd.hide('everything','e. h and not polar_protons')

cmd.extend('DNAselections', DNA_selections )
cmd.extend('DNAselections_nodisplay', lambda: DNA_selections('none') )
cmd.extend('DNAselections_labelsonly', lambda: DNA_selections('labels') )
