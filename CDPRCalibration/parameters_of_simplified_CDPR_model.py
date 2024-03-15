#Created by Liu yifan 2023-05-26
#Connect me by email: liu.y.cm@m.titech.ac.jp
#Last modified time: 2023-08-06
#This file contains some parameters of the simplified CDPR model
#后续工作：把现在的class CDPR拆成两个class，一个是CDPR，保留机构的参数，另一个继承CDPR类但是添加关于CDPR的方法
#这样可以把CDPR的参数和方法分开，方便后续的开发

import geometric_elements as ge
import simulation_utils as su
import draw_utils as du
import numpy as np
import sensor
import CDPR_fuction as CF
#import move_platform as mp
#designed parameters
#长度单位mm
#重量单位kg
DESINGED_HEIGHT = 3440 #3440
DESINGED_WIDTH = 1430  #1465.2
DESIGNED_PULLEY_R = 17.6
DESIGNED_EE_HEIGHT = 49
DESIGNED_EE_WIDTH = 86
DESIGNED_FIN_LENGTH = 287
#很可疑的重量
DESIGNED_WEIGHT = 3.5


#DESINGED_HEIGHT += 2 * DESIGNED_PULLEY_R
DESINGED_WIDTH += 2 * DESIGNED_PULLEY_R

INITIAL_MP_X = 0
INITIAL_MP_Y = DESINGED_WIDTH/2
INITIAL_MP_Z = DESINGED_HEIGHT/2


class Pulley:
    def __init__(self,positon=None,orientation = None,tangent = None,clock = None):
        self.radius = DESIGNED_PULLEY_R
        self.diameter = 2*self.radius
        #self.weight = 0
        #self.tension = 0
        self.position = ge.Point3D(0, 0, 0)
        self.orientation = ge.Point3D(0, 0, 0)
        self.tangent = 0
        self.center = self.position + ge.Point3D(0, self.radius, 0)
        #用于表示该滑轮是否是顺时针旋转时为收线
        self.clock = True
        if(positon != None):
            self.position = positon
        if(orientation != None):
            self.orientation = orientation
        if(tangent != None):
            self.tangent = tangent
        if(clock != None):
            self.clock = clock
        self.update()
        self.cable_point = ge.Point3D(0,0,0)
        self.attach_point = np.array([0, 0, 0])
        self.attach_angle = np.pi
        
        
    def update(self,attachpoint = None,attachangle = None):
        #update the center of the pulley
        rotation_matrix = su.get_rotation_matrix(self.get_orientation())
        new_relative_center = np.dot(rotation_matrix, np.array([self.radius, 0, 0]))
        self.center = self.position + ge.Point3D(new_relative_center[0], new_relative_center[1], new_relative_center[2])
        if attachpoint is not None:
            self.attach_point = attachpoint
        if attachangle is not None:
            self.attach_angle = attachangle

    def set(self,postion = None,orientation = None,tangent = None):
        if(postion):
            self.position = postion
        if(orientation):
            self.orientation = orientation
        if(tangent):
            self.tangent = tangent
        self.update()

    def set_rotation(self,angle):
        self.tangent = angle
        self.update()

    def set_tangent(self,angle):
        self.tangent = angle

    #注意，这里返回的是经过tangent变换以后的orientation，不是原始方向，原始方向请直接调用成员
    def get_orientation(self):
        return self.orientation + ge.Point3D(0, 0, self.tangent)

    def get_rotation_matrix(self):
        R = su.get_rotation_matrix(self.orientation)
        return R

    def get_rotated_matrix(self):
        R = su.get_rotation_matrix(self.get_orientation())
        return R

    def get_unit_vector(self):
        R = self.get_rotated_matrix()
        #将x，y，z三个单位向量经过旋转矩阵变换为u，w，k三个新的单位向量
        u = np.dot(R, np.array([1, 0, 0]))
        w = np.dot(R, np.array([0, 1, 0]))
        k = np.dot(R, np.array([0, 0, 1]))
        return u,w,k

    def set_attach_point(self, point):
        self.attach_point = point
        
#- - - - - - - - - - - - - - - - - - - - - - - - - -
#--------------------------CDPR---------------------------
#- - - - - - - - - - - - - - - - - - - - - - - - - -
class CDPR:
    def __init__(self,Moveplatform):
        self.width = DESINGED_WIDTH
        self.height = DESINGED_HEIGHT
        self.ticks = ge.Ticks()
        self.ticks.set_xticks([-100, -50, 0, 50, 100])
        self.ticks.set_yticks([-1000, 0, 1000, 2000, 3000])
        self.ticks.set_zticks([0, 1000, 2000, 3000, 4000])
        self.MP = Moveplatform
        #CDPR的相对中心
        self.move_platfrom_initial_position = ge.Point3D(INITIAL_MP_X, INITIAL_MP_Y, INITIAL_MP_Z)
        #the base points' position of the cables
        self.cable_base_points_designed = {}
        self.cable_base_points = {}
        self.cable_base_points_estimated = {}
        self.cable_base_point2Ds_estimated = {}
        self.cable_attach_points = {}
        self.cable_base_points_designed[1] = ge.Point3D(0, 0, 0)
        self.cable_base_points_designed[2] = ge.Point3D(0, 0, DESINGED_HEIGHT)
        self.cable_base_points_designed[3] = ge.Point3D(0, DESINGED_WIDTH, DESINGED_HEIGHT)
        self.cable_base_points_designed[4] = ge.Point3D(0, DESINGED_WIDTH, 0)
        self.manuError = {}
        self.manuError[1] = ge.Point3D(10.223,16.173,-1.874)
        self.manuError[2] = ge.Point3D(-17.916,-10.282,27.52)
        self.manuError[3] = ge.Point3D(19.32,-18.628,54.429)
        self.manuError[4] = ge.Point3D(0.472,-70.355,58.181)
        for i in range(1, 5):
            self.cable_base_points_designed[i] += self.manuError[i]
        self.cable_length_designed = np.sqrt((DESINGED_HEIGHT*0.5) ** 2 + (DESINGED_WIDTH*0.5-DESIGNED_PULLEY_R) ** 2) - self.MP.fin_length + np.pi * (80/180) * DESIGNED_PULLEY_R
        print(self.cable_length_designed)
        self.cable_length_designed = 1605

        #the pulleys
        self.pulleys = {}
        self.pulleys[1] = Pulley(positon=self.cable_base_points_designed[1],orientation=ge.Point3D(0,0,0),tangent=-90)
        self.pulleys[2] = Pulley(positon=self.cable_base_points_designed[2],orientation=ge.Point3D(0,0,0),tangent=-90,clock=False)
        self.pulleys[3] = Pulley(positon=self.cable_base_points_designed[3],orientation=ge.Point3D(0,0,0),tangent=90,clock=False)
        self.pulleys[4] = Pulley(positon=self.cable_base_points_designed[4],orientation=ge.Point3D(0,0,0),tangent=90)
        #the encoders
        self.encoders = {}
        for i in range(1,5):
            self.encoders[i] = sensor.Encoder(pulley=self.pulleys[i],id_code = i)
        #the force sensors
        self.force_sensors = {}
        for i in range(1,5):
            self.force_sensors[i] = sensor.ForceSensor(pulley=self.pulleys[i],id_code = i)
        
        #这些内容都被集成进move platform类里面了，以下仅作为存档保留
        '''
        #the base points' position of the moving platform
        self.move_platfrom_base_points = {}
        self.move_platfrom_base_points[1] = ge.Point3D(0, -DESIGNED_EE_WIDTH/2, DESIGNED_EE_HEIGHT/2)
        self.move_platfrom_base_points[2] = ge.Point3D(0, -DESIGNED_EE_WIDTH/2, -DESIGNED_EE_HEIGHT/2)
        self.move_platfrom_base_points[3] = ge.Point3D(0, DESIGNED_EE_WIDTH/2, -DESIGNED_EE_HEIGHT/2)
        self.move_platfrom_base_points[4] = ge.Point3D(0, DESIGNED_EE_WIDTH/2, DESIGNED_EE_HEIGHT/2)

        #the position and gesture of the moving platform(or end effector)
        self.move_platfrom_position_initial = ge.Point3D(INITIAL_MP_X, INITIAL_MP_Y, INITIAL_MP_Z)
        

        #move_platfrom_gesture = ge.Point3D(0, 0, 0)
        self.cable_length_init = []
        self.cable_length_true = []
        self.cable_length_estimated = []
        for i in range(4):
            self.cable_length_init.append(i)
            self.cable_length_estimated.append(i)
            self.cable_length_true.append(i)
        #the length of the fins
        self.fins_length_designed = {}
        self.fins_length_designed[1] = DESIGNED_FIN_LENGTH
        self.fins_length_designed[2] = DESIGNED_FIN_LENGTH
        self.fins_length_designed[3] = DESIGNED_FIN_LENGTH
        self.fins_length_designed[4] = DESIGNED_FIN_LENGTH
        '''


    #generate the geometric errors
    #the errors are generated based on the designed parameters
    #and the results with errors are stored in the self.cable_base_points
    #当前版本认为仅有存在pulley的位置存在误差
    def init_error(self):
        #avg = ge.Point3D(0,0,0)
        #-------------------------------DEBUG-------------------
        print("\nWarning!Debug mode started, the cable base point is pre-setted\n")
#        self.cable_base_points[1] = ge.Point3D(27.34, -14.44, 16.78)
#        self.cable_base_points[2] = ge.Point3D(-15.76, 15.44, 3497.95)
#        self.cable_base_points[3] = ge.Point3D(-21.29, 1481.96, 3460.16)
#        self.cable_base_points[4] = ge.Point3D(9.25, 1446.53, 11.49)
        self.cable_base_points[1] = ge.Point3D(0, 0, 0)
        self.cable_base_points[2] = ge.Point3D(0, 7.6, 3441.3)
        self.cable_base_points[3] = ge.Point3D(0, 1475.7, 3439.3)
        self.cable_base_points[4] = ge.Point3D(0, 1471.6, 2)
        #-------------------------------DEBUG-------------------
        '''
        for label,point in self.cable_base_points_designed.items():
            #self.cable_base_points[label] = su.generate_point_assembling_error(point)
            #self.cable_base_points[label] = su.generate_point_assembling_error(point,debug='x')
            avg.x += self.cable_base_points[label].x
            avg.y += self.cable_base_points[label].y
            avg.z += self.cable_base_points[label].z
        avg.x /= 4
        avg.y /= 4
        avg.z /= 4
        self.move_platfrom_initial_position = avg
        '''
        
    
    def init_without_error(self):
        for label,point in self.cable_base_points_designed.items():
            self.cable_base_points[label] = ge.Point3D(point.x,point.y,point.z)

    def init_with_input(self,set_errors):
        for label,point in self.cable_base_points_designed.items():
            self.cable_base_points[label] = ge.Point3D(point.x + set_errors[label][0],point.y + set_errors[label][1],point.z + set_errors[label][2])

    def initial(self, error = True, set_errors = None):
        if set_errors is None:
            if error:
                self.init_error()
            else:
                self.init_without_error()
        else:
            self.init_with_input(set_errors)
        self.initial_pulleys()

    def initial_pulleys(self):
        for label,point in self.cable_base_points.items():
            self.pulleys[label].position = point
            self.pulleys[label].update()

    #更新滑轮的信息，在解算时别更新，浪费资源
    def update_pulleys(self,MP):
        angles = CF.IK(self,MP,pulley=True)
        for i in range(1,5):
            self.pulleys[i].update(attachpoint = self.cable_attach_points[i])
            #print("pulley " + str(i) + " angle: " + str(np.degrees(np.pi-abs(angles[i]))))
            self.pulleys[i].attach_angle = abs(np.pi-abs(angles[i]))
    
    #更新编码器的初始长度信息
    def initial_encoders(self,init_length = None):
        if init_length is not None:
            for i in range(1,5):
                self.encoders[i].update_init(init_length[i])
        else:
            for i in range(1,5):
                self.encoders[i].update_init(self.read_encoder(id_code=i,isRela=False))
            

    def debug_check_position(self):
        #print("move_platfrom_initial_position: ")
        #print(str(self.move_platfrom_initial_position))
        for label,point in self.cable_base_points.items():
            print(str(label) + ": ",end="")
            point.debugcheck()

    def debug_check_designed_position(self):
        print("move_platfrom_initial_position: ")
        print(str(self.move_platfrom_initial_position))
        for label,point in self.cable_base_points_designed.items():
            print(label + ": ")
            point.debugcheck()

    def get_cable_base_points(self):
        return self.cable_base_points

    def set_cable_attach_points(self,cable_attach_points):
        self.cable_attach_points = cable_attach_points

    def get_cable_attach_points(self):
        return self.cable_attach_points

    def read_encoder(self,id_code,isRela = None,error_level = None):
        if isRela is None or isRela != False:
            isRela = True
        if error_level == None:
            return self.encoders[id_code].read(self.MP,isRela=isRela)
        else:
            return self.encoders[id_code].read(self.MP,isRela=isRela,error_level=error_level)

    def read_encoders(self,isRela = None,error_level = None):
        encoders = {}
        for i in range(1,5):
            if error_level == None:
                encoders[i] = self.read_encoder(i,isRela=isRela)
            else:
                encoders[i] = self.read_encoder(i,isRela=isRela,error_level=error_level)
        return encoders

    def read_force_sensor(self,id_code,error_level = None):
        if error_level == None:
            return self.force_sensors[id_code].read()
        else:
            return self.force_sensors[id_code].read(error_level=error_level)

    def read_force_sensors(self,error_level = None):
        forces = {}
        for i in range(1,5):
            if error_level == None:
                forces[i] = self.read_force_sensor(i)
            else:
                forces[i] = self.read_force_sensor(i,error_level=error_level)
        return forces

    def write_force_sensor(self,forces):
        for i in range(1,5):
            self.force_sensors[i].write(forces[i-1])

    def update_force_sensors(self,MP=None,check=None,forces = None):
        #！！！！！！！！！！！！！！！！！！！！！！！！！
        #如果要优化的话，别用通用的IS方法，那个方法存在大量在此处没必要的残差计算
        if forces is not None:
            self.write_force_sensor(forces)
        else:
            res = CF.IS(self,MP)
            if check is not None:
                residuals = res[1]
                residual = 0
                residuals[0] = (residuals[0]*8)**2
                residuals[1] = (residuals[1]/125)**2
                residuals[2] = (residuals[2]/125)**2
                for r in residuals:
                    residual += r
                if residual > 0.1:
                    raise Exception("the residual of IS is too large, beforce update the F-sensor, pls make true the EE is in the Workingspace!")
                else:
                    self.write_force_sensor(res[0])
            else:
                self.write_force_sensor(res[0])

    def check_sensors(self):
        for i in range(1,5):
            length = self.read_encoder(i)
            force = self.read_force_sensor(i)
            print("Encoder",i,"length:",length," mm")
            print("F-Sensor",i,"force:",force, " N")

    def set_cable_base_points(self,cable_base_points):
        self.cable_base_points = cable_base_points

    def set_initial_length(self,cable_init_length):
        for i in range(1,5):
            self.encoders[i].init_distance = cable_init_length[i]

    #------------------------------------------------------------------------------------------------------------
    #--------------------------the functions below are out of date! Do not use them!-----------------------------
    #------------------------------------------------------------------------------------------------------------
    '''
    #Function name: inverse_kinematic   
    #the IKmode is used to choose the inverse kinematic method
    #method 1 usage: calculate the all length of the cables based on the position of the moving platform
    #method 2 usage: calculate the length of the ith cable based on the position of the moving platform
    def inverse_kinematic(self,IKmode,MP_position = None, MP_pose = None, base_points = None, cable_index = None):
        if MP_position == None:
            MP_position = self.move_platfrom_position
        if base_points == None:
            base_points = self.cable_base_points
        #1和2代表简化的逆运动学，只把platform当做点计算
        if IKmode == 1:
            cable_length = []
            for point in base_points.values():
                length = su.distance_between_two_points(point, MP_position)
                cable_length.append(length)
            return cable_length
        elif IKmode == 2:
            if cable_index == None:
                print("The cable index is not defined")
            point = base_points[cable_index]
            length = su.distance_between_two_points(point, MP_position)
            return length
        #3代表复杂的逆运动学，把platform当做模型计算,不仅返回长度，还返回了平台上连接点的位置
        elif IKmode == 3:
            cable_length = []
            for point in base_points.values():
                length = su.distance_between_two_points(point, MP_position)
                cable_length.append(length)
            return cable_length
        else:
            print("The IKmode is not defined")

    def forward_kinematic (self,FKmode,cable_length, base_points = None):
        if base_points == None:
            base_points = self.cable_base_points
        if FKmode != 0:
            p = su.find_min_grav_2point_intersection(base_points[2],base_points[3],cable_length[1],cable_length[2])
            disA1 = su.distance_between_two_points(base_points[1],p)
            disA4 = su.distance_between_two_points(base_points[4],p)
            if disA1 - cable_length[0] > 0 or disA4 - cable_length[3] > 0:
                p = su.find_min_tension_intersection(base_points,cable_length)
            return p
        else:
            print("The FKmode is Wrong")
            return None
    
    def random_sampling(self,num = None):
        if num is None:
            num = 50
        samling_list = su.generate_points_in_workspace(self.cable_base_points.values(),num)
        return samling_list

    def after_calibration(self,post_processing):
        for i in range(4):
            self.cable_base_points_estimated[i+1] = post_processing.cable_base_points_result[i]
            self.cable_base_point2Ds_estimated[i+1] = post_processing.cable_base_point2Ds_result[i]
            self.cable_length_true[i] = post_processing.cable_init_length_true[i]
            self.cable_length_estimated[i] = post_processing.cable_init_length_result[i]
            self.cable_length_init[i] = post_processing.cable_init_length[i]

    #simulate the relative length from encoder
    #the length_init is the TRUE initial length of the cable
    #the MP_position is the simpling position of the moving platform
    def simu_encoder(self,length_init,MP_position):
        encoder_read = {}
        length = {}
        for index in self.cable_base_points.keys():
            cable_length = self.inverse_kinematic(IKmode = 2,cable_index = index,base_points = self.cable_base_points,MP_position = MP_position)
            length[index] = cable_length
        if len(length.values()) != len(length_init.values()):
            print("The length of the cable is not equal to the initial length")
            return encoder_read
        for index in length.keys():
            encoder = length[index] - length_init[index]
            #encoder = su.generate_encoder_random_error(encoder)
            encoder_read[index]=encoder
        return encoder_read

    def draw(self,ax,mode = None):
        if mode == None:
            mode = 'simple'
        #绘制简化模型（末端执行器当作点处理）
        if mode == 'simple':
            du.draw_points(self.cable_base_points, ax)
            du.draw_point(self.move_platfrom_position, ax,'b','o')
            du.draw_labels(self.cable_base_points, ax)
            du.draw_lines(self.cable_base_points,self.move_platfrom_position, ax)
        #绘制复杂模型（末端执行器当作模型处理）
        elif mode == 'complex':
            du.draw_points(self.cable_base_points, ax)
            du.draw_point(self.move_platfrom_position, ax,'b','o')
            for i in range(1,5):
                du.draw_point(self.move_platfrom_base_points[i]+self.move_platfrom_position,ax,'g','o')
            for i in range(1,5):
                du.draw_line(self.cable_base_points[i],self.move_platfrom_base_points[i]+self.move_platfrom_position,ax)
                if(i<4):
                    du.draw_line(self.move_platfrom_base_points[i]+self.move_platfrom_position,self.move_platfrom_base_points[i+1]+self.move_platfrom_position,ax)
                else:
                    du.draw_line(self.move_platfrom_base_points[i]+self.move_platfrom_position,self.move_platfrom_base_points[1]+self.move_platfrom_position,ax)
    '''


