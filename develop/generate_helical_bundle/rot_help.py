def from_four_points(s,cen,a,b,c):
    import numpy as np
    from numpy import linalg as LA
    s.t = cen
    print "Rotation and translation",s
    ##e1 = (a-b).normalized()
    ##e3 = LA.norm(e1.cross(c-b))
    e1 = (a-b) / LA.norm(a-b)
    e3 = LA.norm(np.cross(e1,c-b))
    e2 = LA.norm(np.cross(e1,e3) )
    print "Rotation and translation",s
    # e2 = LA.norm(e1.cross(e3) )
    # e3 = e1.cross(c-b).normalized()
    # e2 = e1.cross(e3).normalized()
    # print "from_four_points"
    # print e1
    # print e2
    # print e3
    s.R =  Mat(e1.x,e2.x,e3.x,e1.y,e2.y,e3.y,e1.z,e2.z,e3.z)



def __init__(self, xx=None, xy=None, xz=None, yx=None, yy=None, yz=None, zx=None, zy=None, zz=None):
    super(Mat, self).__init__()
    if xx is None: # identity default
        self.xx, self.xy, self.xz = 1.0,0.0,0.0
        self.yx, self.yy, self.yz = 0.0,1.0,0.0
        self.zx, self.zy, self.zz = 0.0,0.0,1.0
    
    elif xy is None and ismat(xx):
        self.xx, self.xy, self.xz = xx.xx, xx.xy,xx.xz
        self.yx, self.yy, self.yz = xx.yx, xx.yy,xx.yz
        self.zx, self.zy, self.zz = xx.zx, xx.zy,xx.zz
    elif yx is None and isvec(xx) and isvec(xy) and isvec(xz):
        self.xx, self.xy, self.xz = xx.x, xy.x, xz.x
        self.yx, self.yy, self.yz = xx.y, xy.y, xz.y
        self.zx, self.zy, self.zz = xx.z, xy.z, xz.z
    elif isnum(xx):
        self.xx, self.xy, self.xz = float(xx), float(xy), float(xz)
        self.yx, self.yy, self.yz = float(yx), float(yy), float(yz)
        self.zx, self.zy, self.zz = float(zx), float(zy), float(zz)
    else:
        raise NotImplementedError
    assert isfloat(self.xx) and isfloat(self.xy) and isfloat(self.xz)
    assert isfloat(self.yx) and isfloat(self.yy) and isfloat(self.yz)
    assert isfloat(self.zx) and isfloat(self.zy) and isfloat(self.zz)
