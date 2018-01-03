import numpy as np
import dolfin as d
import parameters as params


def extract_pots(phi, positions):
    compt_values = np.zeros(positions.shape[0])
    for ii in range(positions.shape[0]):
        compt_values[ii] = phi(positions[ii, :])
    return compt_values


def load_vector(filename, label, mesh):
    V = d.FunctionSpace(mesh, "CG", 2)
    phi = d.Function(V)
    x = phi.vector()
    dump_file = d.HDF5File(mesh.mpi_comm(), filename, 'r')
    dump_file.read(x, label, True)
    return phi

# pos_1 = [[5.196, 22.913, -4.9957]]
# pos_2 = [[5.5945, 22.699, -5.6637]]
# pos_2 = [[8.4, 31.4, -6.15]]


mesh = params.load_just_mesh()
pos_list = [[5.196, 22.913, -4.9957], [8.4, 31.4, -6.151],
            [5.5945, 22.699, -5.6637]]
phi_1 = load_vector('test_del_ih.h5', '0', mesh)
phi_2 = load_vector('test_del_ih.h5', '1', mesh)
phi_3 = load_vector('test_del_ih.h5', '2', mesh)


print(extract_pots(phi_1, np.array([pos_list[1], pos_list[2]])))
print(extract_pots(phi_2, np.array([pos_list[0], pos_list[2]])))
print(extract_pots(phi_3, np.array([pos_list[0], pos_list[1]])))
