import numpy as np
log_file_path = './ProcessedLog/samplepoint.txt'
class ExResult():
    def __init__(self):
        self.log_file_path = log_file_path
        self.point_position = []
        self.cable_tension = []
        self.cable_length = []
        self.samnum = 0
        self.dataflag = False
        self.init_position = []
        self.init_length = {} 
        self.init_tension = {}
        self.pulley_r = 17.6
        self.Kw = 59600 #N/m
        self.length_sp = {}
        self.length_sp[1] = 6309
        self.length_sp[2] = 6503
        self.length_sp[3] = 6503
        self.length_sp[4] = 6309
        with open(log_file_path,'r') as f:
            for line in f:
                data = line.strip().split()
                if len(data) >= 11:
                    self.point_position.append([0.0,float(data[1])+self.pulley_r,float(data[2])])
                    self.cable_tension.append({})
                    self.cable_length.append({})
                    for i in range(4):
                        self.cable_tension[self.samnum][i+1] = float(data[3+i])
                        self.cable_length[self.samnum][i+1] = float(data[7+i])
                    self.samnum += 1
                    self.dataflag = True
                elif len(data) == 10:
                    self.init_position = [float(0.0),float(data[0])+self.pulley_r,float(data[1])]
                    for i in range(4):
                        self.init_tension[i+1] = float(data[2+i])
                        self.init_length[i+1] = float(data[6+i])
        if self.dataflag == False:
            print('No data in the file!')
            exit(0)

    def check(self):
        print('The number of sampling points is:',self.samnum)
        for i in range(self.samnum):
            print("Point ",i)
            print(self.point_position[i])
            print(self.cable_tension[i])
            print(self.cable_length[i])

    def read_relative_lengths(self,point_index):
        if point_index < 0 or point_index >= self.samnum:
            print('The index is out of range!')
            return None
        else:
            relative = {}
            for i in range(1,5):
                relative[i] = self.cable_length[point_index][i] - self.init_length[i]
            return relative

    def elongation_compensation(self):
        compensation = {}
        for i in range(1,5):
            compensation[i] = self.init_tension[i] * (self.length_sp[i]+self.init_length[i]) / self.Kw
        return compensation
        


