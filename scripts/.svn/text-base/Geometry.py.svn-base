from pymol import cmd
from math import asin,acos,atan,cos,sin,tan,sqrt



def getAtomProjectionOnPlane(sel,a,b,c,d):# specify d as distance along norm from origin rather than -distance (ax+by+cz=d)
	"""
	returns a list of all projected balls ordered by signed distance along
	the normal to the plane ax + by + cz + d = 0
	"""
	p = Plane(float(a),float(b),float(c),-float(d)) 
	slice = [s for s in p.getProjections(sel) if s is not None]
	slice.sort()
	return slice

def getAtomProjectionAsR(sel='all',fname='points.R',a=1,b=0,c=0,d=0,varname='slice'):	
	print "write projections to file",fname
	o = open(fname,'w')
	o.write( projectionAsRString(sel,a,b,c,d,varname) )
	o.close()

def projectionAsRString(sel='all',a=0,b=0,c=1,d=0,varname='slice'):	
	slice = getAtomProjectionOnPlane(sel,a,b,c,d)
	Npoints = len(slice)
	s  = ''
	s += "%(varname)s = array(NA,c(%(Npoints)i,9)) \n "%vars() 
	s += "%(varname)s = as.data.frame(%(varname)s) \n "%vars()
	s += "names(%(varname)s) = c('x','y','z','r','xorig','yorig','zorig','rorig','signdis' ) \n"%vars()
	s += "%(varname)s$atype = rep('H',%(Npoints)i) \n "%vars() 

	for ii in range(1,len(slice)+1):
		pa = slice[ii-1]
		#print pa
		x,y,z = pa.coordproj
		r = pa.radiusproj
		ox,oy,oz = pa.atom.coord
		origr = pa.atom.vdw
		signdis = pa.signeddis
		atype = pa.atom.symbol
		occ = pa.atom.q
		bfac = pa.atom.b
		s += "%(varname)s[%(ii)i,1] = "%vars() + `x` + '\n' 
		s += "%(varname)s[%(ii)i,2] = "%vars() + `y` + '\n' 
		s += "%(varname)s[%(ii)i,3] = "%vars() + `z` + '\n' 
		s += "%(varname)s[%(ii)i,4] = "%vars() + `r` + '\n' 
		s += "%(varname)s[%(ii)i,5] = "%vars() + `ox` + '\n' 
		s += "%(varname)s[%(ii)i,6] = "%vars() + `oy` + '\n' 
		s += "%(varname)s[%(ii)i,7] = "%vars() + `oz` + '\n' 
		s += "%(varname)s[%(ii)i,8] = "%vars() + `origr` + '\n' 
		s += "%(varname)s[%(ii)i,9] = "%vars() + `signdis` + '\n' 
		s += "%(varname)s[%(ii)i,9] = "%vars() + `signdis` + '\n' 
		s += "%(varname)s$atype[%(ii)i] = "%vars() + `atype` + '\n' 
		s += "%(varname)s$occ[%(ii)i] = "%vars() + `occ` + '\n' 
		s += "%(varname)s$bfac[%(ii)i] = "%vars() + `bfac` + '\n' 
	return s

def multiSliceAsR(sel='all',fname='/tmp/points.R',a=0,b=0,c=1,interval=1,varname='slice'):
	origin = Plane(float(a),float(b),float(c),0)
	mn =  99999999
	mx = -99999999
	for atom in cmd.get_model(sel).atom:
		d = origin.signedDisTo(atom)
		if mn > d:
			mn = d
		if mx < d:
			mx = d
	
	o = open(fname,'w')
	o.write( "%(varname)s = list() \n"%vars() )
	for ii,d in enumerate(range(mn,mx,interval)):
		print 'slice',ii,d
		variter = varname+"[[%i]]"%(ii+1)
		o.write( projectionAsRString(sel,a,b,c,d,varname=variter) )
	o.close()

#cmd.extend('getAtomProjectionOnPlane',getAtomProjectionOnPlane)
cmd.extend('getAtomProjectionAsR',getAtomProjectionAsR)
cmd.extend('multiSliceAsR',multiSliceAsR)



# argh! why can't I get this to work?
def scaleCoords(amount=1.0,sel='all'):
	m = cmd.get_model(sel)	
	for a in m.atom:
		if a.index == 80:
			delta = [float(k)*(amount-1) for k in a.coord]
			print amount, a.coord, delta
			cmd.translate( [ delta[0], delta[1], delta[2] ] , "index "+`a.index`)

class ProjectedAtom(object):
	"""The ProjectedAtom class."""
	def __init__( self, atom, coordproj, radiusproj, signeddis ):
		self.atom       = atom
		self.coordproj  = coordproj
		self.radiusproj = radiusproj
		self.signeddis  = signeddis
	def __cmp__(self,other):
		return cmp(self.signeddis,other.signeddis)
	def __repr__(self):
		return "TODO!!!"

# lets try something more successful...
class Plane(object):
	"""The Plane class."""

	def __init__(self, a=1,b=0,c=0,d=0):  # (a,b,c).(x,y,z) + d == ax + by + cz + d == 0
		l = sqrt(a**2+b**2+c**2)
		self.a = a/l
		self.b = b/l
		self.c = c/l
		self.d = d/l

	def signedDisTo(self,atom): # just (x0,y0,z0).(a,b,c) + d (assumes a,b,c,d are normalized)
 		return self.a*atom.coord[0] + self.b*atom.coord[1] + self.c*atom.coord[2] + self.d

	def disTo(self,atom):
		return abs( self.signedDisTo(atom) )

	def getProjection(self,atom):
		d = self.signedDisTo(atom)
		r = atom.vdw
		if abs(d) >= r: # atom doesn't intersect plain
			None
		else:
			rproj = r * cos( asin( d/r ) )
			xproj = atom.coord[0] - d * self.a   # XYZ - dist_to_plain * Normal == projected point on plain
			yproj = atom.coord[1] - d * self.b	
			zproj = atom.coord[2] - d * self.c
			pa = ProjectedAtom( atom, (xproj,yproj,zproj), rproj, d )
			return pa

	def getProjections(self,sel='all'):
		return [ self.getProjection(a) for a in cmd.get_model(sel).atom ]


