#Created by Liu yifan 2023-07-19
#Connect me by email: liu.y.cm@m.titech.ac.jp
#Last modified time: 2023-08-05
import geometric_elements as ge
import numpy as np
import simulation_utils as su
class move_platform:
    def __init__(self):
        self.x = 0;
        self.y = 0;
        self.z = 0;
        self.phi = 0;
        self.theta = 0;
        self.psi = 0;
        self.fin_length = 0;
        self.base_points = {};
        self.local_base_points = {};
        self.fin_points = {};
        self.g = 9.8
        self.m = 3.5
        self.weight = self.m*self.g

    def initial(self,width,height,length):
        self.set_local_base_points(width,height)
        self.set_fin_length(length)
        self.set_position(0,0,0)
        self.set_orientation(0,0,0)
        #for i in range(1,5):
        #    print(self.local_base_points[i])

    def set_local_base_points(self, width,height):
        self.local_base_points[1] = ge.Point3D(0, -width/2, height/2)
        self.local_base_points[2] = ge.Point3D(0, -width/2, -height/2)
        self.local_base_points[3] = ge.Point3D(0, width/2, -height/2)
        self.local_base_points[4] = ge.Point3D(0, width/2, height/2)
    
    def set_fin_length(self,length):
        self.fin_length = length

    def set_fin_points(self,finpoints):
        self.fin_points = finpoints

    def set_position(self,x,y,z):
        self.x = x
        self.y = y
        self.z = z

    def set_orientation(self,phi,theta,psi):
        self.phi = phi
        self.theta = theta
        self.psi = psi

    def set_pose(self,pose):
        self.set_position(pose[0],pose[1],pose[2])
        self.set_orientation(pose[3],pose[4],pose[5])
        self.update_base_points()

    #convert the base point to the global coordinate
    #using the rotation matrix
    def update_base_points(self,get=None):
        #rotation matrix
        R = su.get_rotation_matrix([self.phi,self.theta,self.psi])
        for i in range(1,5):
            x = self.x + R[0][0]*self.local_base_points[i].x + R[0][1]*self.local_base_points[i].y + R[0][2]*self.local_base_points[i].z
            y = self.y + R[1][0]*self.local_base_points[i].x + R[1][1]*self.local_base_points[i].y + R[1][2]*self.local_base_points[i].z
            z = self.z + R[2][0]*self.local_base_points[i].x + R[2][1]*self.local_base_points[i].y + R[2][2]*self.local_base_points[i].z
            self.base_points[i] = ge.Point3D(x,y,z)
        if(get == True):
            return self.base_pointsget_fin_length

    def get_rotation_matrix(self):
        R = su.get_rotation_matrix([self.phi,self.theta,self.psi])
        return R

    def get_base_points(self):
        return self.base_points

    def get_fin_length(self):
        return self.fin_length

    def get_fin_points(self):
        return self.fin_points

    def get_fin_point(self,id_code):
        return self.fin_points[id_code]

    def get_pose(self):
        return [self.x, self.y, self.z, self.phi, self.theta, self.psi]

    def get_local_base_points(self):
        return self.local_base_points

    def get_position(self):
        return [self.x,self.y,self.z]

    def get_mass(self):
        return self.m
    




