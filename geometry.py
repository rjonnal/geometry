from matplotlib import pyplot as plt
from matplotlib import patches
import numpy as np
import os,sys
import numpy as np
import matplotlib.colors as colors
import matplotlib.cm as cmx

DEFAULT_MARKER = 'y-'

# convenience classes for geometrical calculations and visualization later:
class Point:
    def __init__(self,x,y):
        self.x = x
        self.y = y

    def distance(self,other_point):
        return np.sqrt((self.x-other_point.x)**2+(self.y-other_point.y)**2)

    def lies_on(self,line):
        line_m,line_b = line.m,line.b
        return self.y==line_m*self.x+line_b

    def get_line(self,other_point):
        m = (self.y-other_point.y)/(self.x-other_pont.x)
        b = self.y - m*self.x
        return Line(m,b)

    def __str__(self):
        return 'x=%0.2f,y=%0.2f'%(self.x,self.y)

    def __repr__(self):
        return 'x=%0.2f,y=%0.2f'%(self.x,self.y)

    def __add__(self,point):
        return Point(self.x+point.x,self.y+point.y)

    def __sub__(self,point):
        return Point(self.x-point.x,self.y-point.y)

    def __div__(self,scalar):
        return Point(self.x/scalar,self.y/scalar)

    def plot(self):
        ph = plt.plot(self.x,self.y,'bo',mec='b',mfc='b',ms=6)

class Line:
    def __init__(self,m=0.0,b=0.0):
        self.m = float(m)
        self.b = float(b)

    def from_points(self,p1,p2):
        # factory method, as follows:
        # p1 = Point(10,3)
        # p2 = Point(4,7)
        # l = Line().from_points(p1,p2)
        if p1.x==p2.x:
            m = np.inf
        else:
            m = (p1.y-p2.y)/(p1.x-p2.x)
        b = p1.y - m*p1.x
        return Line(m,b)

    def __str__(self):
        return 'y = %0.2f x + %0.2f'%(self.m,self.b)

    def __repr__(self):
        return 'y = %0.2f x + %0.2f'%(self.m,self.b)
    
class LineSegment:
    def __init__(self,p1,p2):
        self.p1 = p1
        self.p2 = p2
        self.line = Line().from_points(p1,p2)
        self.length = p1.distance(p2)
        self.midpoint = self.get_midpoint()
        
    def __lt__(self,ls):
        return self.length<ls.length

    def __gt__(self,ls):
        return self.length>ls.length

    def __div__(self,ls):
        return self.length/ls.length

    def get_midpoint(self):
        return (self.p1 + self.p2)/2.0

    def plot(self,marker=DEFAULT_MARKER):
        plt.plot([self.p1.x,self.p2.x],[self.p1.y,self.p2.y],marker)


class Polygon:
    def __init__(self,point_list):
        # point list is assumed to be 'open', i.e. w/o final repetition of first point
        self.vertices = point_list
        self.area = self.get_area()
        self.n = len(point_list)
        self.edges = []
        xmin = np.inf
        xmax = -np.inf
        ymin = np.inf
        ymax = -np.inf
        for v1,v2 in zip(self.vertices[:-1],self.vertices[1:]):
            self.edges.append(LineSegment(v1,v2))
            xmin = min(xmin,min(v1.x,v2.x))
            xmax = max(xmax,max(v1.x,v2.x))
            ymin = min(ymin,min(v1.y,v2.y))
            ymax = max(ymax,max(v1.y,v2.y))
        
        self.edges.append(LineSegment(self.vertices[-1],self.vertices[0]))
        self.xmin = xmin
        self.xmax = xmax
        self.ymin = ymin
        self.ymax = ymax
        self.width = self.xmax - self.xmin
        self.height = self.ymax - self.ymin
        
    def get_circular_vertices(self):
        temp = self.vertices
        temp.append(self.vertices[0])
        return temp
    
    def get_area(self):
        corners = self.vertices
        n = len(corners) # of corners
        area = 0.0
        for i in range(n):
            j = (i + 1) % n
            area += corners[i].x * corners[j].y
            area -= corners[j].x * corners[i].y
        area = abs(area) / 2.0
        return area


    def get_nn(self):
        return self.n

    def centroid(self):
        # don't use self.get_area, which returns only positive areas;
        # this algorithm requires that the area be negative sometimes.
        A = 0.0 
        
        cvertices = self.get_circular_vertices()
        tempx = 0.0
        tempy = 0.0
        for v1,v2 in zip(cvertices[:-1],cvertices[1:]):
            t2 = (v1.x*v2.y - v2.x*v1.y)
            tempx = tempx + (v1.x + v2.x)*t2
            tempy = tempy + (v1.y + v2.y)*t2
            A = A + t2/2.0
        return Point(float(tempx)/6.0/A,float(tempy)/6.0/A)

    def patch(self,color='r'):
        vtx = np.array([[p.x,p.y] for p in self.vertices])
        return patches.Polygon(vtx,fill=True,facecolor=color)


class Rectangle(Polygon):
    def __init__(self,Point1,Point2):
        self.upper_left = Point1
        self.lower_right = Point2
        self.lower_left = Point(Point1.x,Point2.y)
        self.upper_right = Point(Point2.x,Point1.y)
        self.width = Point2.x - Point1.x
        self.height = Point2.y - Point1.y
        Polygon.__init__(self,[self.upper_left,self.upper_right,self.lower_right,self.lower_left])
        
    

class Square(Rectangle):
    def __init__(self,Point_center,edge_length):
        hel = edge_length/2.0
        upper_left = Point(Point_center.x-hel,Point_center.y-hel)
        lower_right = Point(Point_center.x+hel,Point_center.y+hel)
        Rectangle.__init__(self,upper_left,lower_right)



class Region(Polygon,Point):
    def __init__(self,center,points_list):
        Polygon.__init__(self,points_list)
        Point.__init__(self,center.x,center.y)
        self.center = center
        self.x = center.x
        self.y = center.y
        self.short_radii = []
        for e in self.edges:
            self.short_radii.append(LineSegment(self.center,e.midpoint))

        self.long_radii = []
        for v in self.vertices:
            self.long_radii.append(LineSegment(self.center,v))

        self.area = self.get_area()
            
    def __str__(self):
        return 'Region with area %0.2f'%(self.area)


    def get_nn(self):
        return float(len(self.edges))
    
    def get_cc(self):
        return self.center.distance(self.centroid())

    def get_roundness(self):
        return np.min(self.short_radii)/np.max(self.long_radii)


    def plot(self,marker=DEFAULT_MARKER,plot_type='edge'):

        if plot_type is None:
            plot_type = 'edge'

        if plot_type.lower().find('edge')>-1:
            for e in self.edges:
                e.plot(marker)

        if plot_type.lower().find('short')>-1:
            for sr in self.short_radii:
                sr.plot(marker)

        if plot_type.lower().find('long')>-1:
            for lr in self.long_radii:
                lr.plot(marker)

        if plot_type.lower().find('center')>-1:
            for lr in self.long_radii:
                self.center.plot()

    def get_patch(self):
        points = []
        for v in self.vertices:
            points.append([v.x,v.y])
        points = np.array(points)
        return patches.Polygon(points)

    def computeArea(self):
        return self.get_area()
    
        
