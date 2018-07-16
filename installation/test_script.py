import dolfin as d
import vtk

print(vtk.VTK_VERSION)
dump_file = d.HDF5File(d.mpi_comm_world(), 'test.h5', 'w')
dump_file_2 = d.File('test.xml.gz')
solver = d.KrylovSolver("cg", "hypre_amg")
mesh = d.UnitCubeMesh(10, 10, 10)
V = d.FunctionSpace(mesh, "CG", 2)
phi = d.Function(V)

dump_file.write(phi, 'phi_func')
dump_file.close()
        
dump_file_2 << phi

