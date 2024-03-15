#Created by Liu yifan 2023-05-26
#Connect me by email: liu.y.cm@m.titech.ac.jp
#Last modified time: 2023-08-02
#This file contains some functions about time

import time

def time_diff(timestamp):
    current_time = time.perf_counter()
    return (current_time - timestamp)*1000

def time_diff_str(timestamp,second = None):
    if(second):
        formatted_str = "{:.3f}s".format(time_diff(timestamp)/1000)
    else:
        formatted_str = "{:.3f}ms".format(time_diff(timestamp))
    return formatted_str

def time_stamp():
    return time.perf_counter()

#This function is wrong
#Because timestamp is a local variable, doesn't have pass by reference
#pls use time_stamp() instead
def time_update(timestamp):
    current_time = time.perf_counter()
    timestamp = current_time

def wait(x):
    time.sleep(x)