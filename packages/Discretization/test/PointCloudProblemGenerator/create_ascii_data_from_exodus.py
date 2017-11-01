############################################################################
# Copyright (c) 2012-2017 by the DataTransferKit authors                   #
# All rights reserved.                                                     #
#                                                                          #
# This file is part of the DataTransferKit library. DataTransferKit is     #
# distributed under a BSD 3-clause license. For the licensing terms see    #
# the LICENSE file in the top-level directory.                             #
############################################################################

 # This script reads an exodus file and dumps its node coordinates and
# connectivities of all elements in all element blocks into 2 seperate files
# to be used by driver test code

import netCDF4
import numpy
import sys

# exodus file
exodus_file = sys.argv[1]

# prefix of output files
output_prefix = sys.argv[2]

# load the exodus file
nc = netCDF4.Dataset( exodus_file )

# extract the node coordinates
x_coords = nc.variables['coordx']
y_coords = nc.variables['coordy']
z_coords = nc.variables['coordz']

# get the number of nodes
num_node = nc.dimensions['num_nodes'].size

# write the coordinates to a file. first value is a global id (starting at 1)
# next values are the coordinates
coord_file = open( output_prefix + '_node_coords.dat', 'w' )
value = str(num_node) + '\n'
coord_file.write( value )
for n in xrange(num_node):
    value = str(n+1) + ' ' + str(x_coords[n]) + ' ' + str(y_coords[n]) + ' ' + str(z_coords[n]) + '\n'
    coord_file.write( value )
coord_file.close()

# get the number of element blocks
num_block = nc.dimensions['num_el_blk'].size

# get the total number of cell
total_num_cell = nc.dimensions['num_elem'].size

# write the connectivity to a file
connect_file = open( output_prefix + '_connectivity.dat', 'w' )
value = str(total_num_cell) + '\n'
connect_file.write( value )
for b in xrange(num_block):

    # get the element connectivity in the block. these appear to start at 1
    # instead of 0
    connect_name = 'connect' + str(b+1)
    connectivity = nc.variables[connect_name]

    # get the number of cells in the block
    num_cell_name = 'num_el_in_blk' + str(b+1)
    num_cell = nc.dimensions[num_cell_name].size

    # get the number of nodes per cell in the block
    node_per_cell_name = 'num_nod_per_el' + str(b+1)
    nodes_per_cell = nc.dimensions[node_per_cell_name].size

    # write the connectivity of the cells in the block to file. First entry is
    # the global id (starting at 1) second entry is the topology, third entry
    # is number of nodes making the element, next entries are the node ids
    for c in xrange(num_cell):
        value = str(c+1) + ' ' + connectivity.elem_type + ' ' + str(nodes_per_cell) + ' '
        for n in xrange(nodes_per_cell):
            value += str(connectivity[c][n]) + ' '
        value += '\n'
        connect_file.write( value )

connect_file.close()
