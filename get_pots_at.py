
import os
import meshes

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
    mesh = meshes.load_just_mesh()
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


def find_nearest_avail_pt(pos, xx, yy, zz):
    '''Fetches the closes available point in pos'''
    k = np.array((xx, yy, zz))
    dists = np.linalg.norm(pos-k, axis=1)
    if np.min(dists) != 0.0:
        print('Exact point is unavailable')
        print('Picking the closest point to the selection')
    xx, yy, zz = pos[np.argmin(dists), :]
    return xx, yy, zz


def how_many_procs(path, sbspt):
    '''returns how many procs were used in for that particular run'''
    op_files = glob(os.path.join(path, sbspt+'*.h5'))
    prefix = op_files[0].split(sbspt)[1]
    return int(prefix.split('_')[0])


def which_file(path, sbspt, pos, src_idx):
    'returns the file name string where the corresp src matrix is stored'''
    num_proc = how_many_procs(path, sbspt)
    num_pts = len(pos)
    proc_vals = np.linspace(0, num_pts, num_proc + 1).astype(int)
    proc_idx = np.argmax(proc_vals > src_idx)-1
    if proc_idx == -1:
        proc_idx = len(proc_vals) - 1
    f_string = os.path.join(path, sbspt + str(num_proc)
                            + '_' + str(proc_idx) + '.h5')
    return f_string
    
def obtain_orth_planes(src_loc, xx, yy, zz, conductivity, save=False):
    ''' Sources placed at src_loc, with orth planes meeting at xx, yy, zz 
    The pots observed at the corresponding xx, yy, zz planes is reported
    - dumps three npy arrays corresp to each plane'''
    mesh = meshes.load_just_mesh()
    pos_list, conductivity, path, sbspt = params.default_run(conductivity)
    src_x, src_y, src_z = src_loc
    src_loc = find_nearest_avail_pt(pos_list, src_x, src_y, src_z)
    k = np.array((src_loc))
    dists = np.linalg.norm(pos_list-k, axis=1)
    src_idx = np.argmin(dists)    
    f_string = which_file(path, sbspt, pos_list, src_idx)
    dump_file = load_file(os.path.join(path, f_string), mesh)
    this_phi = load_vector(dump_file, mesh, str(src_idx))
    planes = fetch_cut_planes(xx, yy, zz)
    ### instead get this from params as a special planes
    ### Planes get dumped from probe_points file instead - so premtively obtain these planes.
    phi_planes = []
    #  plane pots 
    for ii, label in enumerate('x', 'y', 'z'):
        ele_list = planes[ii]
        num_ele = len(ele_list)
        phi_mat = np.zeros((1, num_ele))
        phi_mat[0, :] = extract_pots(this_phi, np.array(ele_list))
        phi_planes.append(phi_mat)
    dump_file.close()
    print(phi_planes)
    if save:
        phi_mat_fname = sbspt + conductivity + '_plane_phi_mat.npy'
        np.save(os.path.join(path, phi_mat_fname), phi_mat)
    return phi_mat
    

def obtain_unsorted_srcVele(conductivity, save=False):
    '''Fetch the potentials for the default run - SEEG section to test reciprocity'''
    mesh = meshes.load_just_mesh()
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

def obtain_traubs_srcVele(conductivity, save=False):
    '''Fetch the potentials for the default run - SEEG section For Traubs column'''
    mesh = meshes.load_just_mesh()
    pos_list, conductivity, path, sbspt = params.default_run(conductivity)
    num_pts = len(pos_list)
    ele_list = params.load_traubs_points() # points from traubs colum
    num_ele = len(ele_list)
    phi_mat = np.zeros((num_pts, num_ele))
    num_proc = how_many_procs(path, sbspt)
    proc_vals = np.linspace(0, num_pts, num_proc + 1).astype(int)
    for kk in range(num_proc):
        proc_idx = kk + 1  # as proc ids start from 1
        pt_idxs = range(proc_vals[proc_idx - 1], proc_vals[proc_idx])
        f_string = sbspt + str(num_proc) + '_' + str(proc_idx) + '.h5'
        dump_file = load_file(os.path.join(path, f_string), mesh)
        for pt in pt_idxs:
            this_phi = load_vector(dump_file, mesh, str(pt))
            phi_mat[pt, :] = extract_pots(this_phi, np.array(ele_list))
        dump_file.close()
        print('Finished processing proc number: ', proc_idx)
    if save:
        phi_mat_fname = sbspt + conductivity + '_traub_phi_mat.npy'
        np.save(os.path.join(path, phi_mat_fname), phi_mat)
    return phi_mat


# obtain_unsorted_srcVele('anisotropic', save=True)
# obtain_unsorted_srcVele('homogeneous', save=True)
# obtain_unsorted_srcVele('inhomogeneous', save=True)
# obtain_traubs_srcVele('anisotropic', save=True)
obtain_orth_planes(np.array((5., 23., -5.)), 5, 23., -6., 'anisotropic', save=False)
