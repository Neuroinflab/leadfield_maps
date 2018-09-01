import os
import csv
import numpy as np


curr_dir = os.getcwd()
dti_path = os.path.join(curr_dir, os.pardir, 'dti_package')
points_path = os.path.join(curr_dir, 'points')
mesh_path = os.path.join(curr_dir, 'mesh')
sigmas_path = os.path.join(curr_dir, 'sigmas')

anis_path = os.path.join(sigmas_path, 'anis')
inhom_path = os.path.join(sigmas_path, 'inhom')
hom_path = os.path.join(sigmas_path, 'hom')

results_path = os.path.join(curr_dir, 'results')
results_ani_path = os.path.join(results_path, 'anis')
results_inhom_path = os.path.join(results_path, 'inhom')
results_hom_path = os.path.join(results_path, 'hom')


def limits_brain(ipsilateral=True):
    if not ipsilateral:
        xmin, xmax = [1.8, 18.0]   # LEFT - RIGHT
        print('Points do not inc. cerebellum.')
        print('Both hemispheres')
    else:
        xmin, xmax = [1.8, 9.85]   # LEFT - RIGHT
        print('Points do not inc. cerebellum.')
        print('ipsilateral hemisphere (L side) only')
    ymin, ymax = [13.25, 34.35]  # posterior - anterior
    zmin, zmax = [-13.2, -2.3]   # inferior - superior
    # all inclusive.
    # x : 0 to 20
    # y : 4 to 36
    # z : -18.0 to 0
    # obtaining a fine grid first
    return xmin, xmax, ymin, ymax, zmin, zmax


def load_probe_points(filename):
    probe_points = np.load(filename)
    return probe_points


def load_barycenters():
    all_pts = np.load(os.path.join(points_path, 'mesh_midpts.npz'))
    points_to_sample = np.vstack((all_pts['x'],
                                  all_pts['y'],
                                  all_pts['z'])).T
    return points_to_sample


def load_barycenters_ids():
    ids = np.load(os.path.join(points_path, 'special_mesh_midpts.npz'))
    idx_out = ids['idx_out']
    idx_grnd = ids['idx_grnd']
    idx_csf = ids['idx_csf']
    np_fa = ids['np_fa']
    return idx_out, idx_grnd, idx_csf, np_fa


def load_eigen_vectors():
    eig_vecs = np.load(os.path.join(points_path, 'eigen_vecs.npz'))
    np_v1 = eig_vecs['np_v1']
    np_v2 = eig_vecs['np_v2']
    np_v3 = eig_vecs['np_v3']
    return np_v1, np_v2, np_v3


def default_run(conductivity='anisotropic'):
    pos = np.load(os.path.join(points_path, 'probe_points_ipsi_1.0.npy'))
    pos_list = pos.tolist()
    if conductivity == 'anisotropic':
        path = results_ani_path
    elif conductivity == 'inhomogeneous':
        path = results_inhom_path
    else:
        path = results_hom_path
    sbspt = 'def_'   # meaning default_run
    return pos_list, conductivity, path, sbspt


def load_special_points():
    # DO NOT EDIT THESE POINTS
    sp_pts = {}
    sp_pts['sp1'] = np.array((5., 23., -5.))   
    sp_pts['sp2'] = np.array((5., 23., -6))
    sp_pts['sp3'] = np.array((6.7, 22.6, -4.55))  # point in cortex_L
    return sp_pts


def load_traub_morph_props():
    # Fetch traub points
    num_cmpts = [74, 74, 59, 59, 59, 59, 61, 61, 50, 59, 59, 59]
    cell_range = [0, 1000, 1050, 1140, 1230, 1320, 1560,
                  2360, 2560, 3060, 3160, 3260, 3360]
    pop_names = ['pyrRS23', 'pyrFRB23', 'bask23', 'axax23', 'LTS23',
                 'spinstel4', 'tuftIB5', 'tuftRS5', 'nontuftRS6',
                 'bask56', 'axax56', 'LTS56']
    num_cells = np.diff(cell_range) / 10  # 10% MODEL
    total_cmpts = list(num_cmpts * num_cells)
    return num_cmpts, cell_range, pop_names, num_cells, total_cmpts


def load_traubs_points():
    pos_list = []
    with open(os.path.join(points_path, 'traub_post_transform.csv'), 'rb') as csvfile:
        next(csvfile, None)
        spamreader = csv.reader(csvfile, delimiter=',')
        for row in spamreader:
            pos_list.append([float(row[1]), float(row[2]), float(row[3])])
    return pos_list


def load_hippocampus_points():
    pos_list = []
    with open(os.path.join(points_path, 'hippocampus_L_side.csv'), 'r') as csvfile:
        spamreader = csv.reader(csvfile, delimiter=' ')
        for row in spamreader:
            pos_list.append([float(row[0]), float(row[1]), float(row[2])])
    refined_list = pos_list[::100]
    #print(refined_list, len(refined_list))
    return refined_list

def load_cortex_points(size=0.95):
    pos_list = []
    filename = 'cortex_L_'+str(size)+'.csv'
    with open(os.path.join(points_path, filename), 'r') as csvfile:
        spamreader = csv.reader(csvfile, delimiter=' ')
        for row in spamreader:
            pos_list.append([float(row[0]), float(row[1]), float(row[2])])
    refined_list = pos_list[::10]
    return refined_list

def ecog_run(size=0.95):
    conductivity = 'anisotropic'
    path = results_ani_path
    pos_list = load_cortex_points(size)
    sbspt = 'ecog' + str(size) + '_'
    return pos_list, conductivity, path, sbspt


def hippo_eeg_run():
    conductivity = 'anisotropic'
    pos_list = load_hippocampus_points()
    path = results_ani_path
    sbspt = 'hippo_'
    return pos_list, conductivity, path, sbspt
