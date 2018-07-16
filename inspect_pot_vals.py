from dolfin import *
import os

print 'Loading meshes'
mesh =  Mesh(os.path.join("mesh", "mesh_setup.xml"))
V = FunctionSpace(mesh, "CG", 2)
phi = Function(V)

x = phi.vector()
dump_file = HDF5File(mesh.mpi_comm(), 'traub_hom_1.h5', 'r')
dump_file.read(x, '642', True)

dump_file = File("hom_642.pvd")
dump_file << phi
