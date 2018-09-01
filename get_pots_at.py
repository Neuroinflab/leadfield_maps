
import os
import meshes

import numpy as np
import dolfin as d
import parameters as params
from glob import glob


def fetch_cut_planes(xx, yy, zz, res=1.):
    '''gets you planes (as np.array points ) for sag, ax and trav planes'''
    xmin, xmax, ymin, ymax, zmin, zmax = params.limits_brain(ipsilateral=False)
    plane_pt_idx = []
    # X plane
    pl_y, pl_z = np.mgrid[ymin:ymax:res,
                          zmin:zmax:res]
    pl_y = pl_y.flatten()
    pl_z = pl_z.flatten()
    pl_x = np.zeros_like(pl_y) + xx
    x_plane = np.vstack((pl_x, pl_y, pl_z)).T
    #x_plane = inside_brain(x_plane)
    # Y plane
    pl_x, pl_z = np.mgrid[xmin:xmax:res,
                          zmin:zmax:res]
    pl_x = pl_x.flatten()
    pl_z = pl_z.flatten()
    pl_y = np.zeros_like(pl_x) + yy
    y_plane = np.vstack((pl_x, pl_y, pl_z)).T
    #y_plane = inside_brain(y_plane)
    # Z plane
    pl_x, pl_y = np.mgrid[xmin:xmax:res,
                          ymin:ymax:res]
    pl_x = pl_x.flatten()
    pl_y = pl_y.flatten()
    pl_z = np.zeros_like(pl_y) + zz
    z_plane = np.vstack((pl_x, pl_y, pl_z)).T
    #z_plane = inside_brain(z_plane)
    return x_plane, y_plane, z_plane


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
    '''Fetches the closes available point in pos - returns index'''
    print('Checking for :', xx, yy, zz)
    k = np.array((xx, yy, zz))
    dists = np.linalg.norm(pos-k, axis=1)
    if np.min(dists) != 0.0:
        print('Exact point is unavailable')
        print('Picking the closest point to the selection')
        print('Closest point index is:', np.argmin(dists))
        print('Closest point is: ', pos[np.argmin(dists)])
    #xx, yy, zz = pos[np.argmin(dists)]
    return np.argmin(dists)


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
    proc_idx = np.argmax(proc_vals > src_idx) - 1
    if proc_idx == -1:  # perhaps the last entry
        proc_idx = len(proc_vals) - 1
    f_string = os.path.join(path, sbspt + str(num_proc)
                            + '_' + str(proc_idx + 1) + '.h5')   # + 1 because numbering for proc idx starts from 1
    return f_string


def obtain_orth_planes(src_label, orig_label, conductivity, save=False):
    ''' Sources placed at src_label, with orth planes meeting at orig_label 
    The pots observed at the corresponding xx, yy, zz planes is reported
    - dumps an npz with arrays corres to orig_label corresp to each plane'''
    mesh = meshes.load_just_mesh()
    pos_list, conductivity, path, sbspt = params.default_run(conductivity)
    sp_pts = params.load_special_points()
    src_x, src_y, src_z = sp_pts[src_label]
    xx, yy, zz = sp_pts[orig_label]
    print('Testing if src location has exact probe point')
    src_idx = find_nearest_avail_pt(pos_list, src_x, src_y, src_z)
    k = np.array((pos_list[src_idx]))
    dists = np.linalg.norm(pos_list-k, axis=1)
    src_idx = np.argmin(dists)    
    f_string = which_file(path, sbspt, pos_list, src_idx)
    dump_file = load_file(os.path.join(path, f_string), mesh)
    this_phi = load_vector(dump_file, mesh, str(src_idx))
    planes = fetch_cut_planes(xx, yy, zz, 0.1)
    ### TODO - FUTURE
    ### instead get this from params as a special planes
    ### Planes get dumped from probe_points file instead - so premtively obtain these planes.
    phi_planes = []
    #  plane pots
    for ii, label in enumerate(['x', 'y', 'z']):
        ele_list = planes[ii]
        num_ele = len(ele_list)
        phi_mat = np.zeros((1, num_ele))
        phi_mat[0, :] = extract_pots(this_phi, np.array(ele_list))
        phi_planes.append(phi_mat)
    dump_file.close()
    if save:
        phi_mat_fname = sbspt + conductivity + '_' + src_label + '_' + orig_label + '_phi_mat.npz'
        np.savez(os.path.join(path, phi_mat_fname),
                 xx=phi_planes[0], yy=phi_planes[1], zz=phi_planes[2],
                 x_plane=planes[0], y_plane=planes[1], z_plane=planes[2])
    return
    

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
    return

def obtain_traubs_srcVele(conductivity, run_type='seeg', save=False):
    '''Fetch the potentials for the default run - SEEG section For Traubs column'''
    mesh = meshes.load_just_mesh()
    if run_type == 'seeg':
        pos_list, conductivity, path, sbspt = params.default_run(conductivity)
    elif run_type == 'ecog':
        pos_list, conductivity, path, sbspt = params.ecog_run(size=0.98)
    elif run_type == 'heeg':
        pos_list, conductivity, path, sbspt = params.hippo_eeg_run()
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
        phi_mat_fname = sbspt + conductivity + '_' + run_type + '_traub_phi_mat.npy'
        np.save(os.path.join(path, phi_mat_fname), phi_mat)
    return


# obtain_unsorted_srcVele('anisotropic', save=True)
# obtain_unsorted_srcVele('homogeneous', save=True)
# obtain_unsorted_srcVele('inhomogeneous', save=True)
# obtain_traubs_srcVele('anisotropic', run_type='seeg', save=True)
# obtain_traubs_srcVele('anisotropic', run_type='ecog', save=True)
# obtain_traubs_srcVele('anisotropic', run_type='heeg', save=True)
# obtain_traubs_srcVele('anisotropic', run_type='heeg', save=True)
# obtain_orth_planes('sp1', 'sp2', 'anisotropic', save=False)


# # Testing for inspect_pot_vals
pos_list, conductivity, path, sbspt = params.default_run('anisotropic')
sp_pts = params.load_special_points()
src_x, src_y, src_z = sp_pts['sp3']
src_idx = find_nearest_avail_pt(pos_list, src_x, src_y, src_z)
print(which_file(path, sbspt, pos_list, src_idx))
