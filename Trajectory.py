import math
from numpy import arctan2, hypot, sin, cos, tan
class Trajectory:
    """Tracjectory"""
    def minabs(x,y):
        if abs(x) <= abs(y):
            return x
        return y

    def __init__(self, name="", weight=1., startx=0., starty=0., endx=1., endy=1.  ):
        self.name = name
        self.weight = weight
        self.startx = startx
        self.starty = starty
        self.endx = endx
        self.endy = endy
        self.angle = arctan2(endy-starty, endx-startx) * 180 / math.pi
        #print name, weight, startx, starty, endx, endy, self.angle
        self.length = hypot(endy-starty, endx-startx)
        
        self.slope = 0
        
        if endx - startx == 0:
            self.slope = float('inf')
        else:
            self.slope = (endy - starty) / (endx - startx)
        self.segments = []
        self.xstep = []
        self.ystep = []

    def make_segments(self, segment_length=100., ):
        """return list of lines"""
        
        self.xstep = segment_length*cos(self.angle * math.pi / 180.0 )
        self.ystep = segment_length*sin(self.angle * math.pi / 180.0)
        nsegs = int(math.ceil(self.length/segment_length + 0.001))
        #print "step", self.xstep, self.ystep, self.length, self.angle
        seg_startx = self.startx
        seg_starty = self.starty
        
        for i in range(0,nsegs):
            self.segments.append(self.TrajectorySegment(self, seg_startx, seg_starty, seg_startx+self.xstep, seg_starty+self.ystep))
            seg_startx += self.xstep
            seg_starty += self.ystep
    #    print self.startx, self.starty, self.endx, self.endy, " segments: ", self.segments
  
    def get_segment_at(self, pointx, pointy):
#        print "segment_at ", pointx, pointy, "start", self.startx, self.starty,  "end", self.endx, self.endy, "step", self.xstep, self.ystep, "angle", self.angle, "\n"
        if self.xstep == 0:
           return self.segments[int((pointy-self.starty)/self.ystep)] 
        return self.segments[int((pointx-self.startx)/self.xstep)]



    class TrajectorySegment:
        """Trajectory Segment"""
        def __init__(self, parent, startx, starty, endx, endy):
            self.parent = parent
            self.startx = startx
            self.starty = starty
            self.endx = endx 
            self.endy = endy
            self.id = parent.name + ":" + str(startx) + ":" +  str(starty)
        def __str__(self):
            return self.id
        def __repr__(self):
            return self.__str__()
