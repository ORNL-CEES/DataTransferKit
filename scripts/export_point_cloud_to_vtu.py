import pyevtk
import numpy as np

# Input and output file are hardcoded at this time but that can be changed if
# necessary
input_filename = 'point_cloud.txt'
# If output_path is ./point_cloud, the script will output in the current
# directory a file named point_cloud.vtu
output_path = './point_cloud'

point_cloud = np.genfromtxt(input_filename, delimiter=' ')
x = np.ascontiguousarray(point_cloud[:,0])
y = np.ascontiguousarray(point_cloud[:,1])
z = np.ascontiguousarray(point_cloud[:,2])
pyevtk.hl.pointsToVTK(output_path, x, y, z, data={"val" :
    np.ones([len(point_cloud)])})
