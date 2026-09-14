"""
Authors: Zane Thornburg

Functions to interact with and manipulate triangulated surfaces from FreeDTS and TS2CG
"""

import numpy as np
import os


#########################################################################################
def rewriteQfile(in_q_file, out_q_file, rescale_factor=1):
    """
    Inputs:
    Returns:
    Called by:
    Description:
    """
    
    vertices, triangles = readQ(in_q_file)
    
    vertices_RS = rotateScaleVertices(vertices, rescale_factor=rescale_factor)
    
    writeQ(out_q_file, vertices_RS, triangles)
    
    return None
#########################################################################################


#########################################################################################
def readQ(in_q_file):
    """
    Inputs:
    Returns:
    Called by:
    Description:
    """
    
    with open(in_q_file) as f:
        data_lines = f.readlines()
        
    data = []
    
    for line in data_lines:
        line_data = line.split('\n')[0].split(' ')
        line_data = [dat for dat in line_data if dat != '']
        data.append(line_data)
        
#     print(data)
    box_center = [int(float(data[0][0])/2), int(float(data[0][1])/2), int(float(data[0][2])/2)]
#     print(box_center)
    vertex_count = int(data[1][0])
#     print(vertex_count)
    triangle_count = int(data[vertex_count+2][0])
#     print(triangle_count)
    
    vertices = []
    triangles = []
    
    for i in range(vertex_count):
        
        vert_data = data[i+2]

        vert = [float(vert_data[1])-box_center[0], float(vert_data[2])-box_center[1], float(vert_data[3])-box_center[2]]
#         print(vert)
        vertices.append(vert)
        
    vertices = np.array(vertices)
    
    for i in range(triangle_count):
        
        tri_data = data[i+vertex_count+3]
#         print(tri_data)
        tri = [int(float(tri_data[1])), int(float(tri_data[2])), int(float(tri_data[3]))]
#         print(vert)
        triangles.append(tri)
        
    triangles = np.array(triangles)
    
#     print(len(vertices), len(triangles))
    
    return vertices, triangles
#########################################################################################


#########################################################################################
def createRotationMatrix(vertices, x=0, y=0, z=1):
    """
    Inputs:
    Returns:
    Called by:
    Description:
    """
    
    max_dist = 0
    
    for vertex1 in vertices:

        for vertex2 in vertices:

            vector = vertex2 - vertex1

            dist = np.linalg.norm(vector)

            if dist > max_dist:

                max_dist = dist
                n = vector
    
    n_z = np.array([x,y,z])

#     print(n_z)

    n_z = n_z/np.linalg.norm(n_z)

#     print(n_z)

    n_norm = n/np.linalg.norm(n)

    axis = np.cross(n_z,n_norm)
    
    axis = axis/np.linalg.norm(axis)

    theta = np.arccos(np.dot(n_z,n_norm))

    a = np.cos(theta/2)
    b, c, d = axis * np.sin(theta/2)
    aa, bb, cc, dd = a * a, b * b, c * c, d * d
    bc, ad, ac, ab, bd, cd = b * c, a * d, a * c, a * b, b * d, c * d

    rotation_matrix = np.array([[aa + bb - cc - dd, 2 * (bc + ad), 2 * (bd - ac)],
                     [2 * (bc - ad), aa + cc - bb - dd, 2 * (cd + ab)],
                     [2 * (bd + ac), 2 * (cd - ab), aa + dd - bb - cc]])

    return rotation_matrix
#########################################################################################


#########################################################################################
def rotateScaleVertices(vertices, rescale_factor=1, new_center = 750):
    """
    Inputs:
    Returns:
    Called by:
    Description:
    """

    vertices_S = vertices*rescale_factor

    rotation_matrix = createRotationMatrix(vertices_S)
    
    vertices_R = []
    
    for vertex in vertices_S:

        vertex_T = np.transpose(vertex)

        rotated_vertex = np.dot(rotation_matrix,vertex_T)

        vertices_R.append(rotated_vertex)

    vertices_RS = np.array(vertices_R) + new_center
    
    return vertices_RS
#########################################################################################


#########################################################################################
def writeQ(qfilename, vertices, triangles, box_size=[1500, 1500, 1500]):
    """
    Inputs:
    Returns:
    Called by:
    Description:
    """
    
    NV = len(vertices)
    NT = len(triangles)
    
    vertices_array_t = vertices.T
    triangle_array_t = triangles.T

    with open(qfilename, 'w') as f:
        f.write(str(box_size[0])+".00000   " + str(box_size[1]) + ".00000   "  + str(box_size[2]) + ".00000\n")
        f.write(str(NV)+"\n")

        for i in range(NV):
            f.write(str(i)+"  ") #carbon atoms
            np.savetxt(f, [vertices_array_t[:,i].T], fmt='%1.3f', newline=' ', delimiter=' ') #coordinates
            f.write("0\n")

        f.write(str(NT)+"\n")

        for j in range(NT):
            f.write(str(j)+"  ") #carbon atoms
            np.savetxt(f, [triangle_array_t[:,j].T], fmt='%d', newline=' ', delimiter=' ') #coordinates
            f.write("1\n")

    return None
########################################################################################


#########################################################################################
def TS2CG_PLM(inFile, outFile, rescaleFactor=[1, 1, 1], Mashno=4):
    """
    Inputs:
    Returns:
    Called by:
    Description:
    """
    
#     '/home/zane/MinCell/4DWCM/external_programs/TS2CG/PLM -TSfile '
    
    os.system('/home/zane/Software/TS2CG/TS2CG1.1/PLM -TSfile ' + inFile + ' -r PLM -Mashno ' 
              + str(Mashno) + ' -rescalefactor ' + str(rescaleFactor[0]) + ' '
              + str(rescaleFactor[1]) + ' ' + str(rescaleFactor[2]) + ' -o ' + outFile)
    
    return None
#########################################################################################


#########################################################################################
def readTS2CG(fname, rescaleFactor=1):
    """
    Inputs:
    Returns:
    Called by:
    Description:
    """
    
    data_file = fname # sim_properties['working_directory'] + 'membrane/OuterBM.dat'
    
    with open(data_file) as f:
        data_lines = f.readlines()

    data_lines = data_lines[4:]

    mem_coords = []

    for line in data_lines:

        all_data = line.split(' ')
        all_data = [dat for dat in all_data if dat != '']

        x, y, z = all_data[3:6]

        mem_coords.append([float(x),float(y),float(z)])

    memCoordArray = np.array(mem_coords)*rescaleFactor

    vert_x, vert_y, vert_z = memCoordArray.T

    x_CoM = np.average(vert_x)
    y_CoM = np.average(vert_y)
    z_CoM = np.average(vert_z)

    vert_x = vert_x - x_CoM
    vert_y = vert_y - y_CoM
    vert_z = vert_z - z_CoM

    memCoordArray = np.array([vert_x,vert_y,vert_z])

    memCoordArray = memCoordArray.T
    
    memBoundaryCoords = np.reshape(memCoordArray,(int(len(memCoordArray)),3),order='F')
    
    return memBoundaryCoords
#########################################################################################


#########################################################################################
def getCoordConversion(volFrac, small_radius=200.0):
    """
    Inputs:
    Returns:
    Called by:
    Description:
    """

#     small_radius = 200.0

    SA1 = 4*np.pi*small_radius**2
    SA3 = 2*SA1

    smallVol = 4/3*np.pi*small_radius**3
    cellVol = 2*smallVol
    radius = (cellVol*3/4/np.pi)**(1/3)
#     print(radius)
    # radius = 255.0
    
    SA2 = 4*np.pi*radius**2
    
    reducedVol = volFrac*cellVol
    
    SA_V_reducedV = SA2/reducedVol
    
    SA_real = SA_V_reducedV*cellVol
    
    conversionRadius = (SA_real/4/np.pi)**(1/2)

    return conversionRadius
#########################################################################################


