
import os
import numpy as np
import dolfin as d
import parameters as params
from glob import glob


def extract_pots(phi, positions):
    compt_values = np.zeros(positions.shape[0])
    for ii in range(positions.shape[0]):
        compt_values[ii] = phi(positions[ii, :])
    return compt_values


def load_file(filename, mesh):
    print(filename)
    dump_file = d.HDF5File(mesh.mpi_comm(), filename, 'r')
    return dump_file


def load_vector(dump_file, mesh, label):
    V = d.FunctionSpace(mesh, "CG", 2)
    phi = d.Function(V)
    x = phi.vector()
    dump_file.read(x, label, True)
    return phi


def test_few_points():
    mesh = params.load_just_mesh()
    pos_list = [[5.196, 22.913, -4.9957], [8.4, 31.4, -6.151],
                [5.5945, 22.699, -5.6637]]
    dump_file = load_file(os.path.join(params.results_path, 'test_del_ih.h5'),
                          mesh)
    phi_1 = load_vector(dump_file, mesh, '0')
    phi_2 = load_vector(dump_file, mesh, '1')
    phi_3 = load_vector(dump_file, mesh, '2')
    print(extract_pots(phi_1, np.array([pos_list[1], pos_list[2]])))
    print(extract_pots(phi_2, np.array([pos_list[0], pos_list[2]])))
    print(extract_pots(phi_3, np.array([pos_list[0], pos_list[1]])))


def how_many_procs(path, sbspt):
    op_files = glob(os.path.join(path, sbspt+'*.h5'))
    prefix = op_files[0].split(sbspt)[1]
    return int(prefix.split('_')[0])


def obtain_unsorted_srcVele(conductivity, save=False):
    mesh = params.load_just_mesh()
    pos_list, conductivity, path, sbspt = params.default_run(conductivity)
    num_pts = len(pos_list)
    phi_mat = np.zeros((num_pts, num_pts))
    num_proc = how_many_procs(path, sbspt)
    proc_vals = np.linspace(0, num_pts, num_proc + 1).astype(int)
    for kk in range(num_proc):
        proc_idx = kk + 1  # as proc ids start from 1
        pt_idxs = range(proc_vals[proc_idx - 1], proc_vals[proc_idx])
        f_string = sbspt + str(num_proc) + '_' + str(proc_idx) + '.h5'
        dump_file = load_file(os.path.join(path, f_string), mesh)
        for pt in pt_idxs:
            this_phi = load_vector(dump_file, mesh, str(pt))
            phi_mat[pt, :] = extract_pots(this_phi, np.array(pos_list))
        dump_file.close()
        print('Finished processing proc number: ', proc_idx)
    if save:
        phi_mat_fname = sbspt + conductivity + '_phi_mat.npy'
        np.save(os.path.join(path, phi_mat_fname), phi_mat)
    return phi_mat


obtain_unsorted_srcVele('anisotropic', save=True)
