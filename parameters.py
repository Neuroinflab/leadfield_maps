import os
import dolfin as d
import numpy as np


curr_dir = os.getcwd()
dti_path = os.path.join(curr_dir, os.pardir, 'dti_package')
points_path = os.path.join(curr_dir, 'points')
mesh_path = os.path.join(curr_dir, 'mesh')
sigmas_path = os.path.join(curr_dir, 'sigmas')

anis_path = os.path.join(sigmas_path, 'anis')
inhom_path = os.path.join(sigmas_path, 'inhom')
hom_path = os.path.join(sigmas_path, 'hom')


def load_meshes():
    mesh = d.Mesh(os.path.join(mesh_path, "mesh_setup.xml"))
    subdomain = d.MeshFunction("size_t", mesh,
                               os.path.join(mesh_path,
                                            "mesh_setup_physical_region.xml"))
    boundaries = d.MeshFunction("size_t", mesh,
                                os.path.join(mesh_path,
                                             "mesh_setup_facet_region.xml"))
    return mesh, subdomain, boundaries


def load_just_mesh():
    mesh = d.Mesh(os.path.join(mesh_path, "mesh_setup.xml"))
    return mesh


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
