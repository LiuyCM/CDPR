#Created by Liu yifan 2023-05-26
#Connect me by email: liu.y.cm@titech.ac.jp
#Last modified time: 2023-07-20
#This file contains some geometric elements
import math
import numpy

class Point3D:
    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z

    def __add__(self,p):
        return Point3D(self.x + p.x, self.y + p.y, self.z + p.z)

    def __sub__(self,p):
        return Point3D(self.x - p.x, self.y - p.y, self.z - p.z)

    def __truediv__(self, p):
        magnitude_self = math.sqrt(self.x ** 2 + self.y ** 2 + self.z ** 2)
        magnitude_p = math.sqrt(p.x ** 2 + p.y ** 2 + p.z ** 2)

        if magnitude_p != 0:
            return magnitude_self / magnitude_p
        else:
            raise ZeroDivisionError("Division by zero is not allowed.")

    def __getitem__(self, index):
        if index == 0:
            return self.x
        elif index == 1:
            return self.y
        elif index == 2:
            return self.z
        else:
            raise IndexError("Point3D object index out of range")

    def norm(self,ignoreX = None):
        if ignoreX is None:
            return math.sqrt(self.x ** 2 + self.y ** 2 + self.z ** 2)
        else :
            return math.sqrt(self.y ** 2 + self.z ** 2)

    #Beware! Our plenary is y-z plane
    def div2D(self,p):
        magnitude_self = math.sqrt(self.y ** 2 + self.z ** 2)
        magnitude_p = math.sqrt(p.y ** 2 + p.z ** 2)
        if magnitude_p != 0:
            return magnitude_self / magnitude_p
        else:
            raise ZeroDivisionError("Division by zero is not allowed.")

    def copy(self):
        return Point3D(self.x, self.y, self.z)

    def debugcheck(self):
        print("x: ", self.x,end = ",")
        print("y: ", self.y,end = ",")
        print("z: ", self.z)

    def tolist(self):
        return [self.x, self.y, self.z]

    def np(self):
        return numpy.array([self.x, self.y, self.z])

    def __str__(self):
        return "x: " + str(self.x) + " y: " + str(self.y) + " z: " + str(self.z)


#This class is used to store the ticks of the axis
class Ticks:
    def __init__(self):
        self.x = []
        self.y = []
        self.z = []
    def set_xticks(self, x):
        self.x = x
    def set_yticks(self, y):
        self.y = y
    def set_zticks(self, z):
        self.z = z
    def set_ticks(self,x,y,z):
        self.x = x
        self.y = y
        self.z = z
    def get_xticks(self):
        return self.x
    def get_yticks(self):
        return self.y
    def get_zticks(self):
        return self.z
