from dolfin import *
import os

src_idx = 436
filename = '/home/fenics/shared/results/anis/def_7_4.h5'

print 'Loading meshes'
mesh =  Mesh(os.path.join("mesh", "mesh_setup.xml"))
V = FunctionSpace(mesh, "CG", 2)
phi = Function(V)

x = phi.vector()
dump_file = HDF5File(mesh.mpi_comm(), filename, 'r')
dump_file.read(x, str(src_idx), True)

dump_file = File("anis_"+str(src_idx)+".pvd")
dump_file << phi
