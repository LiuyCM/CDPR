#Created by Liu yifan 2023-05-26
#Connect me by email: liu.y.cm@titech.ac.jp
#Last modified time: 2023-12-19
#This file contains some functions that can be used to draw specific geometric elements
import numpy as np
import simulation_utils as su

def draw_points(points, ax, color=None, mark=None,label = None):
    if color is None:
        color = 'r'

    if mark is None:
        mark = 'o'
    #debug info
    #print('color: ', color)
    #print('marker: ', mark)
    x = []
    y = []
    z = []
    if isinstance(points, list):
        for point in points:
            x.append(point[0])
            y.append(point[1])
            z.append(point[2])
    elif isinstance(points, dict):
        for point in points.values():
            x.append(point[0])
            y.append(point[1])
            z.append(point[2])
    else:
        raise('Warring:points must be list or dict')
    if label is None:
        ax.scatter(x, y, z, c=color, marker = mark)
    else:
        ax.scatter(x, y, z, c=color, marker = mark,label = label)

def draw_point(point, ax, color=None, mark=None):
    draw_points([point], ax, color, mark)

def draw_label(point,ax,label,color=None):
    if color is None:
        color = 'black'
    ax.text(point.x, point.y, point.z, label, color=color)

def draw_labels(points, ax, color=None):
    for label,point in points.items():
        draw_label(point,ax,label,color)
    

def draw_line(point1, point2, ax, color=None, mark=None):
    if color is None:
        color = 'black'
    if mark is None:
        mark = '-'
    #debug info
    #print('color: ', color)
    #print('marker: ', mark)
    line_x = [point1[0], point2[0]]
    line_y = [point1[1], point2[1]]
    line_z = [point1[2], point2[2]]
    ax.plot(line_x,line_y,line_z,color = color,linestyle = mark)

def draw_lines(points, point0, ax, color=None, mark=None):
    for point in points:
        draw_line(point, point0, ax, color, mark)

def draw_lines2(points0,points1, ax, color=None, mark=None):
    for i in range(1,5):
        draw_line(points0[i], points1[i], ax, color, mark)

def draw_circle(center,radius,ax,orientation = None):
    #print("start drawing circle,centre location is" + str(center))
    numpoints = 100
    theta = np.linspace(0,2*np.pi,numpoints)
    x = radius*np.cos(theta)
    y = np.zeros(numpoints)
    z = radius*np.sin(theta)
    if(orientation):
        RotationMatrix = su.get_rotation_matrix(orientation)
        rotated_pointes = np.dot(RotationMatrix,np.array([x,y,z]))
        x,y,z = rotated_pointes
    x += center[0]
    y += center[1]
    z += center[2]
    ax.plot(x,y,z,color = 'm')

def draw_CDPR(CDPR,MP,ax,Pulley=None, Fin = None,showlegend = None):
    p_color1 = 'r'
    p_color2 = 'g'
    p_color3 = 'b'
    l_color1 = 'k'
    l_color2 = 'g'
    l_color3 = 'b'
    c_color1 = 'm'
    if Pulley is None and Fin is None:
        draw_points(CDPR.get_cable_base_points(),ax,color = 'r',mark='o')
        draw_lines2(CDPR.get_cable_base_points(),MP.get_fin_points(),ax,mark='-')
        draw_points(MP.get_fin_points(),ax,color = 'g',mark='o')
        draw_lines2(MP.get_fin_points(),MP.get_base_points(),ax,color = 'g',mark='-')
        draw_point(MP.get_position(),ax,color = 'b',mark='o')
        for i in range(1,5):
            if(i<4):
                draw_line(MP.base_points[i],MP.base_points[i+1],ax)
            else:
                draw_line(MP.base_points[i],MP.base_points[1],ax)
    elif Pulley is None:
        draw_points(CDPR.get_cable_base_points(),ax,color = 'r',mark='o')
        draw_lines2(CDPR.get_cable_base_points(),MP.get_fin_points(),ax,mark='-')
        draw_points(MP.get_fin_points(),ax,color = 'g',mark='o')
        draw_lines2(MP.get_fin_points(),MP.get_base_points(),ax,color = 'g',mark='-')
        draw_point(MP.get_position(),ax,color = 'b',mark='o')
        for i in range(1,5):
            if(i<4):
                draw_line(MP.base_points[i],MP.base_points[i+1],ax)
            else:
                draw_line(MP.base_points[i],MP.base_points[1],ax)
        for i in range(1,5):
            draw_circle(center = CDPR.pulleys[i].center, radius = CDPR.pulleys[i].radius, ax = ax,orientation=CDPR.pulleys[i].get_orientation())
    elif Fin is None:
        draw_points(CDPR.get_cable_base_points(),ax,color = 'b',mark='o')
        draw_points(CDPR.get_cable_attach_points(),ax,color = 'r',mark='o')
        draw_lines2(CDPR.get_cable_attach_points(),MP.get_base_points(),ax,mark='-')
        draw_points(MP.get_base_points(),ax,color = 'g',mark='o')
        draw_point(MP.get_position(),ax,color = 'b',mark='o')
        for i in range(1,5):
            if(i<4):
                draw_line(MP.base_points[i],MP.base_points[i+1],ax)
            else:
                draw_line(MP.base_points[i],MP.base_points[1],ax)
        for i in range(1,5):
            draw_circle(center = CDPR.pulleys[i].center, radius = CDPR.pulleys[i].radius, ax = ax,orientation=CDPR.pulleys[i].get_orientation())
    else:
        draw_points(CDPR.get_cable_base_points(),ax,color = 'b',mark='o')
        draw_points(CDPR.get_cable_attach_points(),ax,color = 'r',mark='o')
        draw_lines2(CDPR.get_cable_attach_points(),MP.get_fin_points(),ax,mark='-')
        draw_points(MP.get_fin_points(),ax,color = 'g',mark='o')
        draw_lines2(MP.get_fin_points(),MP.get_base_points(),ax,color = 'g',mark='-')
        draw_point(MP.get_position(),ax,color = 'b',mark='o')
        for i in range(1,5):
            if(i<4):
                draw_line(MP.base_points[i],MP.base_points[i+1],ax)
            else:
                draw_line(MP.base_points[i],MP.base_points[1],ax)
        for i in range(1,5):
            draw_circle(center = CDPR.pulleys[i].center, radius = CDPR.pulleys[i].radius, ax = ax,orientation=CDPR.pulleys[i].get_orientation())
    if showlegend is not None:
        ax.legend(loc='upper right')

def draw_cube(points,ax, color=None, mark=None):
    if color is None:
        color = 'black'
    if mark is None:
        mark = '-'
    
    #draw points
    draw_points(points,ax,color = 'r',mark='o')


def draw_point_cloud(points,values,ax,fig):
    x_coords = [point[0] for point in points]
    y_coords = [point[1] for point in points]
    z_coords = [point[2] for point in points]
    #print("start drawing point cloud")
    img = ax.scatter(x_coords,y_coords,z_coords,c=values,cmap='jet')
    #print("finish drawing point cloud")
    fig.colorbar(img)
    #plt.show()