# 
# class Vec(object):
#    def __init__(self,x,y=None,z=None):
#       if type(x) is type(self):
#          self.x,self.y,self.z = x.x,x.y,x.z
#       elif type(x) in (type([]),type((1,))):
#          self.x,self.y,self.z = x[0],x[1],x[2]
#       elif y is None:
#          assert type(x) in (type(0),type(0.0))
#          self.x,self.y,self.z = x,x,x
#       else:
#          self.x,self.y,self.z = float(x),float(y),float(z)
#    def dot(u,v):
#       return u.x*v.x+u.y*v.y+u.z*v.z
#    def length(u):
#       return math.sqrt(u.dot(u))
#    def cross(u,v):
#       return Vec(u.y*v.z-u.z*v.y,u.z*v.x-u.x*v.z,u.x*v.y-u.y*v.x)
#    def __mul__(u,a):
#       if type(a) is type(0) or type(a) is type(0.0):
#          return Vec(u.x*a,u.y*a,u.z*a)
#       elif type(a) is Vec: 
#          return u.dot(a)
#       else:
#          # print type(a)
#          assert False
#    def __rmul__(u,a):
#       return u*a
#    def __add__(u,v):
#       if type(v) is type(u):
#          return Vec(u.x+v.x,u.y+v.y,u.z+v.z)
#       else:
#          assert type(v) in (type(0),type(0.0))
#          # print type(v),type(u.x),type(u.y),type(u.z)
#          return Vec(u.x+v,u.y+v,u.z+v)
#    def __radd__(u,v):
#       return u+v
#    def __sub__(u,v):
#       return u+(-v)
#    def __rsub__(u,v):
#       return u+(-v)
#    def __neg__(u):
#       return Vec(-u.x,-u.y,-u.z)
#    def __div__(u,a):
#       return u*(1.0/a)
#    def __str__(self):
#       return "%f, %f, %f"%(self.x,self.y,self.z)
#    def __repr__(self):
#       return "Vec( %f, %f, %f )"%(self.x,self.y,self.z)
#    def normalize(u):
#       l = u.length()
#       u.x /= l
#       u.y /= l
#       u.z /= l
#       return u
#    def cgofrompoint(a,c):
#       return [
#             COLOR, 1.0, 1.0, 1.0,     
#             SPHERE,  c.x, c.y, c.z, 1.0,
#             CYLINDER,c.x    ,c.y    ,c.z    ,
#                      c.x+a.x,c.y+a.y,c.z+a.z,0.5,
#                      1,1,1,1,1,1,
#             ]
#       
#       
# 
# 
# class Mat(object):
#    """docstring for Mat"""
#    def __init__(self, xx, xy, xz, yx, yy, yz, zx, zy, zz):
#       super(Mat, self).__init__()
#       self.xx = float(xx)
#       self.xy = float(xy)
#       self.xz = float(xz)
#       self.yx = float(yx)
#       self.yy = float(yy)
#       self.yz = float(yz)
#       self.zx = float(zx)
#       self.zy = float(zy)
#       self.zz = float(zz)
#    def row(m,i):
#       assert type(i) is type(1)
#       if   i is 0: return Vec(m.xx,m.xy,m.xz)
#       elif i is 1: return Vec(m.yx,m.yy,m.yz)
#       elif i is 2: return Vec(m.zx,m.zy,m.zz)
#       else: assert 0 <= i and i <= 2
#    def col(m,i):
#       assert type(i) is type(1)
#       if   i is 0: return Vec(m.xx,m.yx,m.zx)
#       elif i is 1: return Vec(m.xy,m.yy,m.zy)
#       elif i is 2: return Vec(m.xz,m.yz,m.zz)
#       else: assert 0 <= i and i <= 2
#    def rowx(m): return m.row(0)
#    def rowy(m): return m.row(1)
#    def rowz(m): return m.row(2)      
#    def colx(m): return m.col(0)
#    def coly(m): return m.col(1)
#    def colz(m): return m.col(2)      
#    def __mul__(m,v):
#       if type(v) in(type(0),type(0.0)):
#          return Mat( v*m.xx, v*m.xy, v*m.xz, v*m.yx, v*m.yy, v*m.yz, v*m.zx, v*m.zy, v*m.zz )
#       elif type(v) is Vec:
#          return Vec( m.rowx()*v, m.rowy()*v, m.rowz()*v )
#       elif type(v) is Mat:
#          return Mat( m.rowx()*v.colx(), m.rowy()*v.colx(), m.rowz()*v.colx(),
#                      m.rowx()*v.coly(), m.rowy()*v.coly(), m.rowz()*v.coly(),
#                      m.rowx()*v.colz(), m.rowy()*v.colz(), m.rowz()*v.colz() )
#       else:
#          try:
#             return v.__rmul__(m)
#          except:
#             print type(v)
#             raise NotImplementedError
#    def __str__(m):
#       return "Mat( "+str(m.rowx())+"\n     "+str(m.rowy())+"\n     "+str(m.rowz()) + "  )"
#    def transpose(m):
#       return Mat( m.xx, m.yx, m.zx, m.xy, m.yy, m.zy, m.xz, m.yz, m.zz )
#       
