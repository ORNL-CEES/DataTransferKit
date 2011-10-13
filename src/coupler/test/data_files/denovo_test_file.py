###############################################################################
## denovo_test_file.py
## Stuart R. Slattery
## Fri Jun 10 11:14:21 2011
## $Id: template.py,v 1.3 2008/01/02 17:18:47 9te Exp $
##---------------------------------------------------------------------------##
# Problem
#
# - infinite medium of U235
# - 10x10x10 mesh
# - all boundary conditions are reflecting
#
###############################################################################
## Copyright (C) 2008 Oak Ridge National Laboratory, UT-Battelle, LLC.
###############################################################################

import sys, os, math

from sc import *

##---------------------------------------------------------------------------##

initialize(sys.argv)

if node() == 0:
    print "Denovo - pykba Python Front-End"
    print "-------------------------------"
    print "Release      : %16s" % (release())
    print "Release Date : %16s" % (release_date())
    print "Build Date   : %16s" % (build_date())
    print

##---------------------------------------------------------------------------##
## DATABASE
##---------------------------------------------------------------------------##

db = DB("coupler")

db.insert("problem_name", "coupler")
db.insert("problem_type", "EIGENVALUE")
db.insert("downscatter", 1, 1)
db.insert("Pn_order", 0)
db.insert("num_groups", 4)
db.insert("tolerance", 1.0e-6)
db.insert("partition_upscatter", 0, 1)

reflect = [1, 1, 1, 1, 1, 1]
db.insert("boundary", "reflect")
#db.insert("boundary", "vacuum")

db.add_db("boundary_db", "bnd_conditions")
db.insert("boundary_db", "reflect", reflect, 1)

db.insert("tolerance", 1.0e-6)
db.insert("max_itr", 1000)
db.insert("aztec_diag", 0)
db.insert("aztec_output", 0)

db.insert("mg_solver", "krylov")

## MESH

db.insert("num_cells_i", 4)
db.insert("num_cells_j", 4)
db.insert("num_cells_k", 4)

db.insert("delta_x", 2.5)
db.insert("delta_y", 2.5)
db.insert("delta_z", 2.5)

## DECOMPOSITION

db.insert("num_z_blocks", 1)

db.insert("num_blocks_i", 1)
db.insert("num_blocks_j", 1)

db.insert("num_sets", 1)

# Angular options
db.add_db("quadrature_db", "quad_options")
db.insert("quadrature_db", "Sn_order", 4)

# Make Eigenvalue database
db.insert("eigen_solver", "power_iteration")
db.add_db("eigenvalue_db", "eigenvalue")
db.insert("eigenvalue_db", "L2_tolerance", 1.0e-4)
db.insert("eigenvalue_db", "arnoldi_kspace", 20)

db.output()

##---------------------------------------------------------------------------##

manager = Manager()
mat     = Mat()
angles  = Angles()

manager.partition(db, mat, angles)
indexer = manager.get_indexer()
mesh    = manager.get_mesh()

Gx = indexer.num_global(X)
Gy = indexer.num_global(Y)
Gz = mesh.num_cells_dim(Z)

if node() == 0:
    print ">>> Partitioned global mesh with %i x %i x %i cells" \
        % (Gx, Gy, Gz)

##---------------------------------------------------------------------------##
## Material setup

mat.set_num(1)

S_1 = [[1.0]]
S_2 = [[0.0], [1.0]]
S_3 = [[0.0], [0.0], [1.0]]
S_4 = [[0.0], [0.0], [0.0], [1.0]]

mat.assign_xs(0, 0, 2.0, S_1)
mat.assign_xs(0, 1, 2.0, S_2)
mat.assign_xs(0, 2, 2.0, S_3)
mat.assign_xs(0, 3, 2.0, S_4)

nuf = [1.0, 1.0, 1.0, 1.0]
chi = [0.25, 0.25, 0.25, 0.25]
kappa_sigma = [205.0, 205.0, 205.0, 205.0]
mat.assign_fission_kappa(0, nuf, chi, kappa_sigma)

cells = []
mat.assign_id(0, cells)

manager.partition_energy(mat, angles)

##---------------------------------------------------------------------------##
## Verify setup
source = Zero_Source()
manager.setup(source)
manager.verify()

##---------------------------------------------------------------------------##
## HPC OUTPUT
##---------------------------------------------------------------------------##

# Make an HPC output file
hpc_output = HPC_Problem_Output(1, Gx, Gy, Gz)
hpc_output.open("NeutronicsTest")

# Write the database
hpc_output.write_db(db)

# Write the material cross sections
hpc_output.write_xs(mat)

# Write the mixed cross sections
# material ids
matids = Vec_Int(hpc_output.chunk())
hpc_output.start_field_loop()
while not hpc_output.finished_field_loop():
    k = hpc_output.current_chunk()
    if k < Gz:
        for j in xrange(Gy):
            for i in xrange(Gx):
                cell  = indexer.g2g(i, j, k)
                index = indexer.g2g(i, j, 0)
                matids[index] = 0
    hpc_output.write_matids(matids)
    hpc_output.advance_loop()
del matids

# Output the source
shapes = Vec_Dbl()
hpc_output.write_src_info(ZERO_SOURCE, 0, shapes)

# Close
hpc_output.close()

manager.close()
finalize()

###############################################################################
##                            end of denovo_test_file.py
###############################################################################

