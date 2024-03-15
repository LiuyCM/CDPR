#Created by Liu yifan 2023-09-21
#Connect me by email: liu.y.cm@m.titech.ac.jp
#Last modified time: 2023-11-12
#模拟传感器

import simulation_utils as su

#编码器
class Encoder():
    def __init__(self,pulley,id_code,init_distance = None):
        #注意这里的所有误差均为置信区间，使用时需要配合置信度
        #以防你未来忘了，这里举个例子说明一下，比如置信区间为0.06，置信度为0.9
        #就代表着如果真值为x，则测量值有90%的概率处于x-0.03到x+0.03之间
        #默认为正态分布
        #默认随机误差，单位均为毫米mm
        self.default_random_error = 0.06
        self.pulley = pulley
        self.id = id_code
        if init_distance is None:
            self.init_distance = 0
        #默认置信度，这里使用函数自带的默认置信度90%
        #DEFAULT_CONFIDENCE_LEVEL = 0.90x

    def update_init(self,init_distance):
        self.init_distance = init_distance
        
    def read(self,MP,isRela = None,error_level = None):
        #懒得写就直接忽视没滑轮的情况了
        #根据自己的id查到对应的pulley和fin
        fin_point = MP.get_fin_point(self.id)
        pulley_attach_point = self.pulley.attach_point
        #计算测量值
        dis = su.distance_between_two_points(fin_point,pulley_attach_point)
        #dis += self.pulley.radius * self.pulley.attach_angle
        #和热海的算法保持一致
        #但这个其实有问题吧，这里少算四分之一个圆弧，那从winch到pulley的长度就该补上这个四分之一圆
        dis += self.pulley.radius * self.pulley.attach_angle - (self.pulley.radius * 3.1415926/2)
        #加入随机误差
        if error_level != None:
            dis += su.generate_norm_error(error_level)
        if isRela == None:
            isRela = True
        if isRela:
            return dis - self.init_distance
        else:
            return dis

#力传感器
class ForceSensor():
    def __init__(self,pulley,id_code):
        #注意这里的所有误差均为置信区间，使用时需要配合置信度
        #单位均为N
        self.default_random_error = 1
        self.pulley = pulley
        self.id = id_code
        self.force = 0.0

    def read(self,error_level = None):
        readforce = self.force
        if error_level != None:
            readforce += su.generate_norm_error(error_level)
        return readforce
    pass

    def write(self,force):
        self.force = force




