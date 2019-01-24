from pymol import cmd
from pymol.cgo import *
from random import randrange
import glob
import sets
# from vecmat import *

numcom = 0


def getres(sele):
   try:
      r = list(sets.Set([(x.chain,int(x.resi)) for x in cmd.get_model(sele).atom]))
      r.sort()
      return r
   except:
      return []


class Vec(object):
   def __init__(self,x,y=None,z=None):
      if type(x) is type(self):
         self.x,self.y,self.z = x.x,x.y,x.z
      elif type(x) in (type([]),type((1,))):
         self.x,self.y,self.z = x[0],x[1],x[2]
      elif y is None:
         if type(x) in (type(0),type(0.0)):
            self.x,self.y,self.z = x,x,x
         elif type(x) is type(Vec):
            self.x,self.y,self.z = x.x,x.y,x.z
      else:
         self.x,self.y,self.z = float(x),float(y),float(z)
   def dot(u,v):
      return u.x*v.x+u.y*v.y+u.z*v.z
   def length(u):
      return math.sqrt(u.dot(u))
   def cross(u,v):
      return Vec(u.y*v.z-u.z*v.y,u.z*v.x-u.x*v.z,u.x*v.y-u.y*v.x)
   def __mul__(u,a):
      if type(a) is type(0) or type(a) is type(0.0):
         return Vec(u.x*a,u.y*a,u.z*a)
      elif type(a) is Vec: 
         return u.dot(a)
      else:
         # print type(a)
         assert False
   def __rmul__(u,a):
      return u*a
   def __add__(u,v):
      if type(v) is type(u):
         return Vec(u.x+v.x,u.y+v.y,u.z+v.z)
      else:
         assert type(v) in (type(0),type(0.0))
         # print type(v),type(u.x),type(u.y),type(u.z)
         return Vec(u.x+v,u.y+v,u.z+v)
   def __radd__(u,v):
      return u+v
   def __sub__(u,v):
      return u+(-v)
   def __rsub__(u,v):
      return u+(-v)
   def __neg__(u):
      return Vec(-u.x,-u.y,-u.z)
   def __div__(u,a):
      return u*(1.0/a)
   def __str__(self):
      return "%f, %f, %f"%(self.x,self.y,self.z)
   def __repr__(self):
      return "Vec( %f, %f, %f )"%(self.x,self.y,self.z)
   def normalize(u):
      l = u.length()
      u.x /= l
      u.y /= l
      u.z /= l
   def normalized(u):
      v = Vec(u)
      u.normalize()
      return u
   def cgofrompoint(a,c):
      return [
            COLOR, 1.0, 1.0, 1.0,     
            SPHERE,  c.x, c.y, c.z, 1.0,
            CYLINDER,c.x    ,c.y    ,c.z    ,
                     c.x+a.x,c.y+a.y,c.z+a.z,0.5,
                     1,1,1,1,1,1,
            ]
      
      


class Mat(object):
   """docstring for Mat"""
   def __init__(self, xx=None, xy=None, xz=None, yx=None, yy=None, yz=None, zx=None, zy=None, zz=None):
      super(Mat, self).__init__()
      if xx is None: # identity default
         self.xx = 1.0
         self.yy = 1.0
         self.zz = 1.0
      if type(xx) in (type(0),type(0.0)):
         self.xx = float(xx)
         self.xy = float(xy)
         self.xz = float(xz)
         self.yx = float(yx)
         self.yy = float(yy)
         self.yz = float(yz)
         self.zx = float(zx)
         self.zy = float(zy)
         self.zz = float(zz)
      elif type(xx) is type(Vec(0)): # cols specified as vectors
         self.xx = xx.x
         self.xy = xy.x
         self.xz = xz.x
         self.yx = xx.y
         self.yy = xy.y
         self.yz = xz.y
         self.zx = xx.z
         self.zy = xy.z
         self.zz = xz.z     
      else:
         assert false
   def row(m,i):
      assert type(i) is type(1)
      if   i is 0: return Vec(m.xx,m.xy,m.xz)
      elif i is 1: return Vec(m.yx,m.yy,m.yz)
      elif i is 2: return Vec(m.zx,m.zy,m.zz)
      else: assert 0 <= i and i <= 2
   def col(m,i):
      assert type(i) is type(1)
      if   i is 0: return Vec(m.xx,m.yx,m.zx)
      elif i is 1: return Vec(m.xy,m.yy,m.zy)
      elif i is 2: return Vec(m.xz,m.yz,m.zz)
      else: assert 0 <= i and i <= 2
   def rowx(m): return m.row(0)
   def rowy(m): return m.row(1)
   def rowz(m): return m.row(2)      
   def colx(m): return m.col(0)
   def coly(m): return m.col(1)
   def colz(m): return m.col(2)      
   def __mul__(m,v):
      if type(v) in(type(0),type(0.0)):
         return Mat( v*m.xx, v*m.xy, v*m.xz, v*m.yx, v*m.yy, v*m.yz, v*m.zx, v*m.zy, v*m.zz )
      elif type(v) is Vec:
         return Vec( m.rowx()*v, m.rowy()*v, m.rowz()*v )
      elif type(v) is Mat:
         return Mat( m.rowx()*v.colx(), m.rowy()*v.colx(), m.rowz()*v.colx(),
                     m.rowx()*v.coly(), m.rowy()*v.coly(), m.rowz()*v.coly(),
                     m.rowx()*v.colz(), m.rowy()*v.colz(), m.rowz()*v.colz() )
      else:
         try:
            return v.__rmul__(m)
         except:
            print type(v)
            raise NotImplementedError
   def __div__(m,v):
      return m*(1/v)
   # def __add__(m,v):
   #    if type(v) in(type(0),type(0.0)):
   #       return Mat( v+m.xx, v+m.xy, v+m.xz, v+m.yx, v+m.yy, v+m.yz, v+m.zx, v+m.zy, v+m.zz )
   #    elif type(v) is Mat:
   #       return Mat( v.xx+m.xx, v.xy+m.xy, v.xz+m.xz, v.yx+m.yx, v.yy+m.yy, v.yz+m.yz, v.zx+m.zx, v.zy+m.zy, v.zz+m.zz )
   #    else:
   #       try:
   #          return v.__rmul__(m)
   #       except:
   #          print type(v)
   #          raise NotImplementedError
   def __str__(m):
      return "Mat( "+str(m.rowx())+"\n     "+str(m.rowy())+"\n     "+str(m.rowz()) + "  )"
   def transpose(m):
      return Mat( m.xx, m.yx, m.zx, m.xy, m.yy, m.zy, m.xz, m.yz, m.zz )
   def trace(m):
      return m.xx+m.yy+m.zz
   def axisofrot(R):
      cos_theta = (R.trace() - 1.0 ) / 2.0
      tolerance = 0.0000000000001
      if cos_theta > -1.0 + tolerance and cos_theta < 1.0 - tolerance:
         x = math.sqrt( ( R.xx - cos_theta ) / ( 1.0 - cos_theta ) )
         y = math.sqrt( ( R.yy - cos_theta ) / ( 1.0 - cos_theta ) )
         z = math.sqrt( ( R.zz - cos_theta ) / ( 1.0 - cos_theta ) )
         if R.zy < R.yz: x = -x
         if R.xz < R.zx: y = -y
         if R.yx < R.xy: z = -z
         theta = math.acos( cos_theta )
         return Vec( x, y, z ),theta
      elif cos_theta >= 1.0 - tolerance:
         theta = 0.0;
         return Vec( 1.0, 0.0, 0.0 ), theta

def projection_matrix(v):
   m = Mat( v.x * v.x, v.x * v.y, v.x * v.z, v.y * v.x, v.y * v.y, v.y * v.z, v.z * v.x, v.z * v.y, v.z * v.z )
   return m / v.dot(v)

def rotation_matrix_radians(axis,angle):
   n = axis.normalized()
   sin_theta = math.sin( angle )
   cos_theta = math.cos( angle )
   R = projection_matrix(n)
   R *= 1.0 - cos_theta
   R.xx += cos_theta;       R.xy -= sin_theta * n.z; R.xz += sin_theta * n.y
   R.yx += sin_theta * n.z; R.yy += cos_theta;       R.yz -= sin_theta * n.x
   R.zx -= sin_theta * n.y; R.zy += sin_theta * n.x; R.zz += cos_theta
   return R;

def rotation_matrix(axis,angle):
   return rotation_matrix_radians(axis,angle*math.pi/180.0)

def com(sel="all"):
   ## assumes equal weights (best called with "and name ca" suffix)
   model = cmd.get_model(sel)
   c = Vec(0)
   for a in model.atom:
      c += Vec(a.coord)
   c = c/len(model.atom)
   return c


def showcom(sel="all"):
   global numcom
   c=com(sel)
   print "Center of mass: ",c
   cgo = [pymol.cgo.COLOR, 1.0, 1.0, 1.0, SPHERE, c.x, c.y, c.z, 1.0] ## white sphere with 3A radius
   cmd.load_cgo(cgo, "com%i"%numcom)
   numcom += 1


class Stub(object):
   """docstring for Stub"""
   def __init__(self, frame, cen):
      super(Stub, self).__init__()
      self.frame = frame
      self.cen = cen
   def from_four_points(s,cen,a,b,c):
      s.cen = cen
      e1 = Vec(a-b).normalized()
      e3 = Vec(e1.cross(c-b)).normalized()
      e2 = Vec(e1.cross(e3)).normalized()
      s.frame = Mat(e1.x,e2.x,e3.x,e1.y,e2.y,e3.y,e1.z,e2.z,e3.z)
      return s
   def to_frame(s,x):
      return s.frame * (x - s.cen)
   def from_frame(s,x):
      return (s.frame.transpose() * x) + s.cen
      

class Jump(object):
   """docstring for Jump"""
   def __init__(self, rot, trans):
      super(Jump, self).__init__()
      assert type(rot) is Mat and type(trans) is Vec
      self.rot = rot
      self.trans = trans


class ResBB(object):
   """docstring for ResBB"""
   def __init__(self, n, ca=None, c=None, ss=None):
      super(ResBB, self).__init__()
      if type(n) is ResBB:
         self.n  = n.n
         self.ca = n.ca
         self.c  = n.c
         self.ss = n.ss
         return
      elif type(n) is type(""):
         assert len(getres(n)) is 1
         m = cmd.get_model(n)
         self.n  = Vec( m.atom[0].coord )
         self.ca = Vec( m.atom[1].coord )
         self.c  = Vec( m.atom[2].coord )
         self.ss = m.atom[1].ss
         return 
      assert type(n) is Vec
      assert type(ca) is Vec
      assert type(c) is Vec
      self.n  = n
      self.ca = ca
      self.c  = c
      self.ss = ss
   def rms(s,o):
      assert type(o) is ResBB
      return math.sqrt( (s.n-o.n)*(s.n-o.n) + (s.ca-o.ca)*(s.ca-o.ca) + (s.c-o.c)*(s.c-o.c) )
   def __rmul__(r,m):
      return ResBB( m*r.n, m*r.ca, m*r.c, r.ss )
   def __add__(r,v):
      return ResBB( r.n+v, r.ca+v, r.c+v, r.ss )
   def __sub__(r,v):
      return r + (-v)
   def __str__(r):
      return "n "+str(r.n)+", ca "+str(r.ca)+", c "+str(r.c)+", ss "+r.ss
   def stub(r):
      return Stub.from_four_points(Stub(None,None),r.ca,r.n,r.ca,r.c)


class DisulfLib(object):
   """docstring for DisulfLib"""
   def __init__(self, fn):
      super(DisulfLib, self).__init__()
      self.fn = fn
      self.lib = {"EE":[],"EH":[],"EL":[],"HE":[],"HH":[],"HL":[],"LE":[],"LH":[],"LL":[]}
      for l in open(fn).readlines():
         ss1,ss2,sep,rt,xx,xy,xz,yx,yy,yz,zx,zy,zz,x,y,z = l.split()
         self.lib[ss1+ss2].append(Jump(Mat(xx,xy,xz,yx,yy,yz,zx,zy,zz),Vec(x,y,z)))
   
   def disulf_rms(self,r1,r2):
      assert type(r1) is ResBB
      assert type(r2) is ResBB
      minrms = 9e9
      if (r1.ca-r2.ca).length() > 10: return minrms
      ss1,ss2 = r1.ss,r2.ss
      if r1.ss not in "HEL" or len(r1.ss) != 1: ss1 = "L"
      if r2.ss not in "HEL" or len(r2.ss) != 1: ss2 = "L"
      for j in self.lib[ss1+ss2][:1]:
         s = r1.stub()
         r1_in_s = s.to_frame(r1)
         r3_in_s = r1_in_s + Vec(1,1,1) #( j.rot * r1_in_s ) + j.trans
         r3      = s.from_frame(r3_in_s)
         rms = r2.rms(r3)
         if rms < minrms: minrms = rms
      return minrms


def chirality(fe1,fe2,fe3,fe4):
   a = fe2-fe1
   b = fe3-fe1
   c = a.cross(b)
   return c.dot(fe4-fe1)


def dihedral(p1,p2,p3,p4):
   a = ( p2 - p1 ).normalized()
   b = ( p3 - p2 ).normalized()
   c = ( p4 - p3 ).normalized()
   x = -a.dot(c) + a.dot(b) * b.dot(c)
   y =  a.dot( b.cross(c) );
   return abs(math.atan2(y,x)) * 180.0 / 3.14159


def angle(p1,p2,p3=None):
	if p3 is None:
   		return math.acos( p1.normalized().dot(p2.normalized()) ) * 180.0 / 3.14159
	else:
   		a = ( p2 - p1 ).normalized()
   		b = ( p2 - p3 ).normalized()
   		return math.acos( a.dot(b) ) * 180.0 / 3.14159











def dimeraxis1(m1,m2):
   if len(m1.atom) != len(m2.atom):
      print "selections must contain the same number of atoms!"
      return
   i = j = randrange(len(m1.atom))
   while j == i: j = randrange(len(m1.atom))
   a1,a2 = Vec(m1.atom[i].coord),Vec(m2.atom[i].coord)
   u = (a1+a2)/3
   a1,a2 = Vec(m1.atom[j].coord),Vec(m2.atom[j].coord)
   v = (a1+a2)/3
   if v.x-u.x < 0:
      u,v = v,u
   return (v-u).normalized()
   
   


def dimeraxis(sel1,sel2,nsamp=100):
   m1 = cmd.get_model(sel1)
   m2 = cmd.get_model(sel2)
   a,wtot = Vec(0),0
   for i in range(int(nsamp)):
      u = dimeraxis1(m1,m2)
      # print u
      a  += u
      wtot += u.length()
   a /= wtot
   return a
  



def axisofrot1(m1,m2,m3):
   if len(m1.atom) != len(m2.atom) or len(m2.atom) != len(m3.atom) or len(m3.atom) != len(m1.atom):
      print "selections must contain the same number of atoms!"
      return
   i = randrange(len(m1.atom))
   a1,a2,a3 = Vec(m1.atom[i].coord),Vec(m2.atom[i].coord),Vec(m3.atom[i].coord)
   u = (a1+a2+a3)/3
   i = randrange(len(m1.atom))
   a1,a2,a3 = Vec(m1.atom[i].coord),Vec(m2.atom[i].coord),Vec(m3.atom[i].coord)
   v = (a1+a2+a3)/3
   if u.length() > v.length():
      u,v = v,u
   return v-u
   
   


def axisofrot(sel1,sel2,sel3,nsamp=100):
   m1 = cmd.get_model(sel1).atom[0]
   m2 = cmd.get_model(sel2).atom[0]
   m3 = cmd.get_model(sel3).atom[0]
   m1 = Vec(m1.coord)
   m2 = Vec(m2.coord)
   m3 = Vec(m3.coord)
   return ((m2-m1).cross(m2-m3)).normalized()
   # a,wtot = Vec(0),0
   # for i in range(int(nsamp)):
   #    u = axisofrot1(m1,m2,m3)
   #    a  += u
   #    wtot += u.length()
   # a /= wtot
   return a
   
#   0.373875069479 0.733798169028 0.567236881341 24.7756576007
   
   
cmd.extend("com",com)
cmd.extend("showcom",showcom)
cmd.extend("axisofrot",axisofrot)


def bond_zns(sel):
   cmd.unbond(sel+" and resn ZNS",sel+" and not resn ZNS")
   # for c,i in getres("resn ZNS"):
   #    cmd.bond(sel+" and chain %s and resi %i and name ZN1"%(c,i), 
   #             sel+" and chain %s and resi %i and name S*"%(c,i))
   for c,i in getres("name ZN1"):
      cmd.bond(sel+" and chain %s and resi %i and name ZN1"%(c,i),
               sel+" and chain %s and resi %i and name S*,N2,N4"%(c,i))
   allsg = cmd.get_model(sel+" and (resn CYS and (name CB))").atom
   allsf = cmd.get_model(sel+" and (resn ZHC and (name S*,C1,C4))").atom
   while allsg:
      closest = [9e9,None,None]
      for sga in allsg:
         for sfa in allsf:
            sg = Vec(sga.coord)
            sf = Vec(sfa.coord)
            if (sg-sf).length() < closest[0]:
               closest = [(sg-sf).length(),sga,sfa]
      sga,sfa = closest[1:]
      allsg.remove(sga)
      allsf.remove(sfa)
      if closest[0] > 10.0: break
      cmd.bond(sel+" and resi %s and chain %s and name %s"%(sga.resi,sga.chain,sga.name),
               sel+" and resi %s and chain %s and name %s"%(sfa.resi,sfa.chain,sfa.name) )

cmd.extend("bond_zns",bond_zns)
 
   
def trans(sel,x,y=None,z=None):
   if type(x) is type(Vec(0)):
      x,y,z = x.x,x.y,x.z
   cmd.translate([x,y,z],sel,0,0)
   # cmd.do("alter_state 1,%s,x=x+ %f"%(sel,float(x)))
   # cmd.do("alter_state 1,%s,y=y+ %f"%(sel,float(y)))
   # cmd.do("alter_state 1,%s,z=z+ %f"%(sel,float(z)))
cmd.extend("trans",trans)

def rot(sel,axis,ang,cen=Vec(0,0,0)):
   cmd.rotate([axis.x,axis.y,axis.z],ang,sel,0,0,None,[cen.x,cen.y,cen.z])
   # R = rotation_matrix(axis,ang)
   # m = cmd.get_model(sel)
   # for i in range(len(m.atom)):
   #     tmp = ( R * (Vec(m.atom[i].coord)-cen) ) + cen
   #     m.atom[i].coord = [tmp.x,tmp.y,tmp.z]
   # cmd.load_model(m,sel,1,discrete=1)

def rotrad(sel,axis,ang,cen=None):
   return rot(sel,axis,ang*180.0/math.pi,cen)

def pointaxis(sel):
   u = Vec(1,1,1)
   m = cmd.get_model(sel)
   c = com(sel)
   axis = Vec(0)
   for a in m.atom:
      v = Vec(a.coord)-c
      if v.dot(u) < 0: v = -v
      axis += v
   return axis.normalized()
   

def alignaxis(sel,newaxis,oldaxis=None,cen=None):
   if oldaxis is None: oldaxis = pointaxis(sel)
   if cen is None: cen = com(sel)
   newaxis.normalize()
   oldaxis.normalize()
   axis = newaxis.cross(oldaxis).normalized()
   ang  = -math.acos(max(-1.0,min(1.0,newaxis.dot(oldaxis))))*180/math.pi
   # print "rot around",axis,"by",ang
   if str(ang) == "nan": return
   rot(sel,axis,ang,cen)

def showvec(x,y,z,l=100,c=None):   
   if c is None: c = com("all")
   v = Vec(x,y,z).normalized()
   a1 = c + v*l/2
   a2 = c - v*l/2
   obj = [
         BEGIN, LINES,
         COLOR, 1.0, 1.0, 1.0, 
         VERTEX, a1.x, a1.y, a1.z,
         VERTEX, a2.x, a2.y, a2.z,
         END
   ]                                                                                            
   cmd.load_cgo(obj,'vec')
   

def proj(u,v):
   return projection_matrix(u)*v

def projperp(u,v):
   return v - proj(u,v)



def mysetview(look,up,pos=None,cen=None,ncp=None,fcp=None):
   Xaxis = -look.cross(up).normalized()
   Yaxis = projperp(look,up).normalized()   
   Zaxis = -look.normalized()
   oldv = cmd.get_view()
   v = list(oldv)
   r = Mat(Xaxis,Yaxis,Zaxis)
   v[0] = r.xx; v[1] = r.xy; v[2] = r.xz
   v[3] = r.yx; v[4] = r.yy; v[5] = r.yz
   v[6] = r.zx; v[7] = r.zy; v[8] = r.zz
   if pos is not None: v[ 9:12] = pos.x,pos.y,pos.z
   if cen is not None: v[12:15] = cen.x,cen.y,cen.z
   if ncp is not None: v[15]    = ncp 
   if fcp is not None: v[16]    = fcp 
   cmd.set_view(v)
   return Yaxis

def meancoords(sel1,sel2,n="mix",w=0.5):
	cmd.delete(n)
	a = [Vec(x.coord) for x in cmd.get_model(sel1).atom]
	b = [Vec(x.coord) for x in cmd.get_model(sel2).atom]
	d = [(a[i]-b[i]).length() for i in range(len(a))]
	#print max(d), argmax(d), cmd.get_model("179L").atom[argmax(d)].resi
	#print min(d)
	cmd.create(n,sel1)
	m = cmd.get_model(sel1)
	for j in range(len(m.atom)):
		m.atom[j].coord[0] = w*b[j].x + (1.0-w)*a[j].x
		m.atom[j].coord[1] = w*b[j].y + (1.0-w)*a[j].y
		m.atom[j].coord[2] = w*b[j].z + (1.0-w)*a[j].z
	cmd.load_model(m,n,1)

cmd.extend("meancoords",meancoords)

	

def mygetview():
	v = cmd.get_view()
   	return Vec(v[2],v[5],v[8])

def swell():
	a = [Vec(x.coord) for x in cmd.get_model("177L").atom]
	b = [Vec(x.coord) for x in cmd.get_model("179L").atom]
	d = [(a[i]-b[i]).length() for i in range(len(a))]
	#print max(d), argmax(d), cmd.get_model("179L").atom[argmax(d)].resi
	#print min(d)
	for i in range(0,101):
		cmd.create("mix%03i"%i,"177L")
		m = cmd.get_model("177L")
		r = float(i)/100.0
		print "mixing",r
		for j in range(len(m.atom)):
			# for k in range(len(m.atom)):
			# 	d = (a[j]-a[k]).length() - (b[j]-b[k]).length()
			# print j,m.atom[j].coord
			m.atom[j].coord[0] = r*b[j].x + (1.0-r)*a[j].x
			m.atom[j].coord[1] = r*b[j].y + (1.0-r)*a[j].y
			m.atom[j].coord[2] = r*b[j].z + (1.0-r)*a[j].z
			# print j,m.atom[j].coord
		cmd.load_model(m,"mix%03i"%i,1)
		cmd.save("177L_179L_raw_mix%03i.pdb"%i,"mix%03i"%i)
	
cmd.extend("swell",swell)

def mkhelix(sel,t,r,n):
	v = cmd.get_view()
	for i in range(1,n):
		cmd.delete("h%i"%i)
		cmd.create("h%i"%i,sel)
		rot  ("h%i"%i,Vec(0,0,1),i*r,Vec(0,0,0))
		trans("h%i"%i,Vec(0,0,i*t))
	cmd.set_view(v)

def mkhelix4(sel,t,r,n):
	for i in range(n):
		cmd.delete("h%i"%i)
		cmd.create("h%i"%i,sel)
		rt = 90.0 * (i%4)
		print i,rt
		if i%4==0 and i > 0:
			cmd.create("hem%i"%i,"basehem")
			trans("hem%i"%i,Vec(0,0,t*(int(i/4))))
			rot  ("hem%i"%i,Vec(0,0,1),r*(int(i/4)),Vec(0,0,0))
			if i in (4,12,20,28,36,44,52):
				rot  ("hem%i"%i,Vec(0,0,1),90,Vec(0,0,0))
		rot  ("h%i"%i,Vec(0,0,1),rt,Vec(0,0,0))
		trans("h%i"%i,Vec(0,0,t*(int(i/4))))
		rot  ("h%i"%i,Vec(0,0,1),r*(int(i/4)),Vec(0,0,0))
	cmd.hide("ev","not (name N,CA,C,O,CB) and (h1,h3,h4,h6,h9,h11,h12,h14,h17,h19,h20,h22,h25,h27)")
	cmd.hide("ev","name SG and not (h0,h5,h8,h13,h16,h21,h24,h29,h32,h37,h40,h45,h48,base)")
#PyMOL>run /Users/sheffler/pymol_scripts/util.py; mkhelix4("4hel* and vis",10.2,-23.9,12)	


def mirror(sel,nname):
	cmd.delete(nname)
	a = [Vec(x.coord) for x in cmd.get_model(sel).atom]
	#print max(d), argmax(d), cmd.get_model("179L").atom[argmax(d)].resi
	#print min(d)
	cmd.create(nname,sel)
	m = cmd.get_model(sel)
	for j in range(len(m.atom)):
		m.atom[j].coord[0] *= -1
	cmd.load_model(m,nname,1)
	
	
def mkc4(sel,a=Vec(1,0,0),c=Vec(0,0,0)):
	cmd.delete("c1")
	cmd.delete("c2")
	cmd.delete("c3")
	cmd.delete("c4")			
	cmd.create("c1",sel)
	cmd.create("c2",sel)
	cmd.create("c3",sel)
	cmd.create("c4",sel)			
	rot("c2",a, 90,c)
	rot("c3",a,180,c)
	rot("c4",a,270,c)		
	
def alignall(sel="",obj=None):
	l = cmd.get_object_list()
	if not obj: obj = l[0]
	if obj not in l: 
		print "ERROR object",obj,"not found!!!"
		return
	for o in l:
		if o==obj: continue
		cmd.do("align "+o+" and ("+sel+"),"+obj+" and ("+sel+")")
	return

	
def showcst(fname):
	for l in open(fname).readlines():
		kind,a1,r1,a2,r2 = l.split()[:5]
		a1 = "(resi %s and name %s)"%(r1,a1)
		a2 = "(resi %s and name %s)"%(r2,a2)
		cmd.bond(a1,a2)
		cmd.show("lines",a1)
		cmd.show("lines",a2)		
	
	
def bondzn():
	for o in cmd.get_object_list():
		for c,r in getres(o+" and elem ZN"):
			print "add zn bonds for",o,"resi",r,"chain",c
			cmd.bond("(%s and resi %s and chain %s)"%(o,r,c), "(%s and elem N) within 2.5 of (%s and resi %s and chain %s)"%(o,o,r,c) )

def mktetrafrom3(sel):
	cmd.delete("tetra*")
	a2 = Vec( 0.942809043336065, 0                ,1-0.333333328372267).normalized()
	a3 = Vec(-0.471404521668032, 0.816496582446736,1-0.333333328372267).normalized()
	a4 = Vec(-0.471404521668032,-0.816496582446736,1-0.333333328372267).normalized()
	
	# cmd.create('tetra0',sel)
	cmd.create("tetra2",sel)
	cmd.create("tetra3",sel)
	cmd.create("tetra4",sel)
	rot("tetra2",a2,180,Vec(0,0,0))
	rot("tetra3",a3,180,Vec(0,0,0))
	rot("tetra4",a4,180,Vec(0,0,0))
	cmd.zoom()
	
	
def printtetra3():
	vecs = [
		Vec(0,0,+1                    ), Vec(1,0,0 ),
		Vec(0,0,+1                    ), Vec(1,0,0 ),
		Vec( 1,0,0                    ), Vec(0,0,+1),
		Vec( 1,0,0                    ), Vec(0,0,+1),
		Vec( 1,0,0                    ), Vec(0,0,+1),
		Vec( 1,0,0                    ), Vec(0,0,-1),
		Vec(-0.5,+0.866025403784439,0 ), Vec(0,0,+1),
		Vec(-0.5,+0.866025403784439,0 ), Vec(0,0,+1),
		Vec(-0.5,+0.866025403784439,0 ), Vec(0,0,+1),
		Vec(-0.5,+0.866025403784439,0 ), Vec(0,0,-1),
		Vec(-0.5,-0.866025403784438,0 ), Vec(0,0,+1),
		Vec(-0.5,-0.866025403784438,0 ), Vec(0,0,+1),
		Vec(-0.5,-0.866025403784438,0 ), Vec(0,0,+1),
		Vec(-0.5,-0.866025403784438,0 ), Vec(0,0,-1),
		
		
		
	]
	rots = [rotation_matrix(Vec( 0.942809043336065, 0                ,1-0.333333328372267).normalized(),  0.0),
	        rotation_matrix(Vec( 0.942809043336065, 0                ,1-0.333333328372267).normalized(),180.0),
			rotation_matrix(Vec(-0.471404521668032, 0.816496582446736,1-0.333333328372267).normalized(),180.0),
	 		rotation_matrix(Vec(-0.471404521668032,-0.816496582446736,1-0.333333328372267).normalized(),180.0)]	
	for i,r in enumerate(rots):
		for j,v in enumerate(vecs):
			t = r*v
			if j%2==1:
				print "%+.10F,%+.10F,%+.10F 0,0,0"%(t.x,t.y,t.z)			
			else:
				print "xyz %i %+.10F,%+.10F,%+.10F "%(i+1,t.x,t.y,t.z),
			
	print





def inversion(sel,nname="inv"):
	cmd.delete(nname)
	a = [Vec(x.coord) for x in cmd.get_model(sel).atom]
	#print max(d), argmax(d), cmd.get_model("179L").atom[argmax(d)].resi
	#print min(d)
	cmd.create(nname,sel)
	m = cmd.get_model(sel)
	for j in range(len(m.atom)):
		m.atom[j].coord[0] *= -1
		m.atom[j].coord[1] *= -1
		m.atom[j].coord[2] *= -1
	cmd.load_model(m,nname,1)
	

def setchains(sele="all"):
	aa="ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789"
	i = 0
	for o in enumerate(cmd.get_object_list()):
		for c in cmd.get_chains():
			cmd.alter("chain %s and %s"%(c,o),"chain='%s'"%aa[i%len(aa)])
			i += 1
	return
	
def avg(sel1,sel2,sel3,sel4):
	cmd.delete("avg")
	cmd.create("avg",sel1)
	m1 = cmd.get_model(sel1).atom
	m2 = cmd.get_model(sel2).atom
	m3 = cmd.get_model(sel3).atom
	m4 = cmd.get_model(sel4).atom			
	m  = cmd.get_model("avg")
	for i in range(len(m1)):
		m.atom[i].coord[0] = (m1[i].coord[0]+m2[i].coord[0]+m3[i].coord[0]+m4[i].coord[0])/4.0
		m.atom[i].coord[1] = (m1[i].coord[1]+m2[i].coord[1]+m3[i].coord[1]+m4[i].coord[1])/4.0
		m.atom[i].coord[2] = (m1[i].coord[2]+m2[i].coord[2]+m3[i].coord[2]+m4[i].coord[2])/4.0
	cmd.load_model(m,"avg",1)
	

def randvecsphere():
	v = 2*Vec(random.random(),random.random(),random.random())-1
	while v.length() > 1: v = 2*Vec(random.random(),random.random(),random.random())-1
	return v.normalized()

def gen180(sele,n):
	cmd.delete("oneighty*")
	for i in range(n):
		cmd.create("oneighty%i"%i,sele)
		rot("oneighty%i"%i,randvecsphere(),180,Vec(0,0,0))


def alignvis():
	l = [k for k,v in cmd.get_vis().items() if v[0]]
	l.sort()
	b = l[0]
	toprint = ""
	for i in l[1:]:
		if cmd.align(i,b)[0] > 2:
			if toprint=="": 
				for jj,j in enumerate(cmd.get_object_list()):
					if j.startswith(i):
						toprint = str(jj+1)
			cmd.delete(i)
	print "NEXT:",toprint
			

cmd.extend("alignvis",alignvis)

def gensymc(sele,n,axis=Vec(0,0,1)):
	cmd.delete("symc*")
	d = 360.0/n	
	for i in range(n):
		cmd.create("symc%i"%i,sele)
		rot("symc%i"%i,axis,d)


def aligndimerz(sele,alignsele=None):
	if alignsele is None: alignsele = sele
	cmd.remove(sele+" and resn HOH")
	trans(sele,-com(alignsele))
	a = cmd.get_model(alignsele+" and chain A and name CA").atom
	b = cmd.get_model(alignsele+" and chain B and name CA").atom
	if len(a) != len(b) or len(a) == 0:
		print "ERROR on %s: subunits are not the same size!"%alignsele
		return False
	axis = Vec(0,0,0)
	for i in range(len(a)):
		axis1 = ( Vec(a[i].coord)+Vec(b[i].coord) ).normalized()
		if axis.length() > 0.0001 and axis.dot(axis1) < 0:
			axis1 *= -1
		axis += axis1		
		# print axis1
	axis.normalize()
	# print "axis of rotation:",axis
	alignaxis(sele,Vec(0,0,1),axis,Vec(0,0,0))
	return True
	

def aligntrimerz(sele,alignsele=None):
	if alignsele is None: alignsele = sele
	cmd.remove(sele+" and resn HOH")
	trans(sele,-com(alignsele))
	a = cmd.get_model(alignsele+" and chain A and name CA").atom
	b = cmd.get_model(alignsele+" and chain B and name CA").atom
	c = cmd.get_model(alignsele+" and chain C and name CA").atom
	# print "subunit lengths:",len(a),len(b),len(c)
	if len(a) != len(b) or len(a) != len(c) or len(a) == 0:
		print "ERROR on %s: subunits are not the same size!"%alignsele
		return False
	axis = Vec(0,0,0)
	for i in range(len(a)):
		axis1 = ( Vec(a[i].coord)+Vec(b[i].coord)+Vec(c[i].coord) ).normalized()
		if axis.length() > 0.0001 and axis.dot(axis1) < 0:
			axis1 *= -1
		axis += axis1		
		# print axis1
	axis.normalize()
	# print "axis of rotation:",axis
	alignaxis(sele,Vec(0,0,1),axis,Vec(0,0,0))
	return True
	

def aligntrimersindir(dir):
	for fn in glob.glob(dir+"/*.pdb1"):
		cmd.delete("all")
		cmd.load(fn)
		if aligntrimerz(fn[len(fn)-9:len(fn)-5]):
			cmd.save(fn+"_alignz.pdb","chain A")
			print "SUCCESS on",fn,fn[len(fn)-9:len(fn)-5]
		else:
			print "ERROR on",fn



def loadmov(fn):
	cmd.delete("mov")
	for i in glob.glob(fn):
		cmd.load(i,"mov")
	
def maked2(sele):
	cmd.delete("S11")
	cmd.delete("S12")
	cmd.delete("S21")
	cmd.delete("S22")			
	cmd.create("S11",sele)
	cmd.create("S12",sele)
	cmd.create("S21",sele)
	cmd.create("S22",sele)
	rot("S12",Vec(1,0,0),180)
	rot("S21",Vec(0,1,0),180)
	rot("S22",Vec(1,0,0),180)
	rot("S22",Vec(0,1,0),180)
	cmd.alter("S11","chain='A'")
	cmd.alter("S12","chain='B'")
	cmd.alter("S21","chain='C'")
	cmd.alter("S22","chain='D'")
		
		
def makec2(sele):
	cmd.delete("S1")
	cmd.delete("S2")
	cmd.create("S1",sele)
	cmd.create("S2",sele)
	rot("S1",Vec(0,0,1),0)
	rot("S2",Vec(0,0,1),180)
	cmd.alter("S1","chain='A'")
	cmd.alter("S2","chain='B'")
		
		
def makec3(sele):
	cmd.delete("S1")
	cmd.delete("S2")
	cmd.delete("S3")
	cmd.create("S1",sele)
	cmd.create("S2",sele)
	cmd.create("S3",sele)
	rot("S1",Vec(0,0,1),0)
	rot("S2",Vec(0,0,1),120)
	rot("S3",Vec(0,0,1),240)
	cmd.alter("S1","chain='A'")
	cmd.alter("S2","chain='B'")
	cmd.alter("S3","chain='C'")
		
def getxform(sel1,sel2):
	cmd.create("TEMPORARY",sel1)
	sel1 = "TEMPORARY"
	c1a = com(sel1+" and resi 1 and name C")
	c1b = com(sel2+" and resi 1 and name C")
	v = c1b-c1a
	trans(sel1,v)
	c2a = (com(sel1+" and resi 2 and name C") - c1b).normalized()
	c2b = (com(sel2+" and resi 2 and name C") - c1b).normalized()
	axis1 = c2a.cross(c2b)
	ang1 = math.acos(c2a.dot(c2b))
	rot(sel1,axis1,180.0/math.pi*ang1,c1b)
	axis2 = (com(sel2+" and resi 2 and name C")-c1b).normalized()
	P = projection_matrix(axis2)
	c3a = (com(sel1+" and resi 3 and name C") - c1b)
	c3b = (com(sel2+" and resi 3 and name C") - c1b)
	# print c1b
	# print com(sel1+" and resi 3 and name C")
	# print com(sel2+" and resi 3 and name C")
	# print c3a,"|||",c3b,"|||",  c3a.dot(axis2), c3b.dot(axis2)
	c3a -= P*c3a
	c3b -= P*c3b
	axis2 = c3a.cross(c3b).normalized()
	# print c3a,"|||",c3b,"|||",  c3a.dot(axis2), c3b.dot(axis2)
	ang2 =   math.acos(c3a.normalized().dot(c3b.normalized()))
	# print ang2*180.0/math.pi
	rot(sel1,axis2,180.0/math.pi*ang2,c1b)
	M = rotation_matrix_radians(axis2,ang2)*rotation_matrix_radians(axis1,ang1)
	axis,ang = M.axisofrot()
	# trans(sel1,-v)
	# rot(sel1,axis,-ang,sel1+" and resi 1 and name C")
	# cmd.delete("TEMPORARY")
	return v,axis,ang*180.0/math.pi



def propagatexform(sel1,sel2,n=10):
	cmd.delete("PX*")
	for i in range(n): cmd.create("PX%i"%i,sel1)
	v,axis,ang = getxform(sel1,sel2)
	for i in range(1,n):
		for j in range(i,n):
			trans("PX%i"%j,v)
			rot("PX%i"%j,axis,ang,com("PX%i and name C and resi 1"%j))
	
	
	
	
	
	
	
	
	
		