import vtk_functions as vtk_utils
from vtk.util.numpy_support import vtk_to_numpy
import parameters as params
import numpy as np
import meshes
import os
import sys
import csv


def probe_points_assign(res, ipsilateral, save=False):
    if not ipsilateral:
        xmin, xmax = [1.8, 18.0]
        print('Points do not inc. cerebellum.')
        print('Both hemispheres')
    else:
        xmin, xmax = [1.8, 9.85]
        print('Points do not inc. cerebellum.')
        print('ipsilateral hemisphere (L side) only')
    ymin, ymax = [13.25, 34.35]
    zmin, zmax = [-13.2, -2.3]
    # all inclusive.
    # x : 0 to 20
    # y : 4 to 36
    # z : -18.0 to 0
    # obtaining a fine grid first
    xx, yy, zz = np.mgrid[xmin:xmax:res,
                          ymin:ymax:res,
                          zmin:zmax:res]
    xx = xx.flatten()
    yy = yy.flatten()
    zz = zz.flatten()
    points_to_sample = np.vstack((xx, yy, zz)).T
    tpoints = points_to_sample.shape[0]
    print('Total points inside cubic meshgrid: ', tpoints)
    print('Now probing for points inside brain')
    # probe points where its non zero to give points inside brain.
    seg = img_wrapper('brain_mask_refined.nii.gz', 'mask')
    seg.probe_points(points_to_sample)
    print('Done probing for points inside brain mask refined')
    np_seg = vtk_to_numpy(seg.get_array())
    idx = np.where(np_seg > 0)
    print('Total points inside brain volume: ', len(idx[0]))
    selected_points = points_to_sample[idx]
    if save:
        if ipsilateral:
            save_as = 'probe_points_ipsi_'
        else:
            save_as = 'probe_points_'
        np.save(os.path.join(params.points_path, save_as + str(res) + '.npy'),
                selected_points)
        print('Dumped points inside brain in ' + save_as + str(res) + '.npy')
    return selected_points


def probe_points_labels(points, save=False):
    brain_atlas_path = os.path.join(params.dti_path, 'brain_atlas.nii.gz')
    seg = vtk_utils.image_wrapper(brain_atlas_path, 'mask')
    seg.probe_points(points)
    np_labels = vtk_to_numpy(seg.get_array())
    if save:
        np.save(os.path.join(params.points_path,
                             'probe_pts_labels.npy'), np_labels)
    return np_labels


def compute_point_wts(points, save=False):
    labels_idx = probe_points_labels(points)
    unique_idx, counts = np.unique(labels_idx, return_counts=True)
    # Idea is to give a dummy weight to the point location based on freq occ.
    # The extension of this to to enable a sorting based on amplitudes
    # as factions
    unique_ordered = [x for _, x in sorted(zip(counts, unique_idx),
                                           key=lambda pair: pair[0],
                                           reverse=True)]
    first_wt = [(ii+1) * jj for ii, jj in enumerate(np.sort(counts))][::-1]
    assert(len(np.unique(first_wt)) == len(counts))
    points_wt = [first_wt[unique_ordered.index(kk)] for kk in labels_idx]
    # print([labels_txt[ii] for ii in unique_ordered])
    # print(np.sort(counts)[::-1])
    # print(counts_wt, unique_ordered, labels_idx, points_wt)
    if save:
        np.save(os.path.join(params.points_path,
                             'probe_pts_wts.npy'), points_wt)
    return points_wt


def compute_barycenters(save=False):
    ''' Loads mesh, gets cells of the mesh, and evals barycenters
    dumps this in mesh_midpts.npz'''
    mesh, subdomains, boundaries = meshes.load_meshes()
    cells = params.d.cells(mesh)
    count = 0
    pts_x = []
    pts_y = []
    pts_z = []
    for cell in cells:
        count += 1
        pts_x.append(cell.midpoint().x())
        pts_y.append(cell.midpoint().y())
        pts_z.append(cell.midpoint().z())
    if save:
        np.savez(os.path.join(params.points_path, 'mesh_midpts.npz'),
                 x=np.array(pts_x),
                 y=np.array(pts_y),
                 z=np.array(pts_z))
        print('Obtained bary centers of the meshes, dumped in mesh_midpts.npz')
    print('Total Tets: ', count)
    points_to_sample = np.vstack((np.array(pts_x),
                                  np.array(pts_y),
                                  np.array(pts_z))).T
    return points_to_sample


def img_wrapper(name, id):
    '''Helper function to load vtk nifti files'''
    ii = vtk_utils.image_wrapper(os.path.join(params.dti_path, name), id)
    return ii


def id_barycenters(points_to_sample, save=False):
    '''finds the indices of the barycenters corresp to csf, air, ele, brain
    dumps these values as special_mesh_points.npz'''
    # Electrode - ID-ing ground electrode
    grnd = img_wrapper('ground_ele.nii.gz', 'grnd')
    grnd.probe_points(points_to_sample)
    np_grnd = vtk_to_numpy(grnd.get_array())
    idx_grnd = np.where(np_grnd == 2)[0]
    # Mask points - Here ID-ing air
    seg = img_wrapper('brain_mask_refined.nii.gz', 'segs')
    seg.probe_points(points_to_sample)
    np_seg = vtk_to_numpy(seg.get_array())
    idx_out = np.where(np_seg == 0)[0]
    # Fractional Anisotropy points - 3vox smoothed - ID-ing wm-ness
    fa_map = img_wrapper('fa_smooth_3vox_masked.nii.gz', 'FA')
    # fa_map = fetch_img_wrapper('fa_masked.nii.gz', 'FA')
    fa_map.probe_points(points_to_sample)
    np_fa = vtk_to_numpy(fa_map.get_array())
    idx_wm = np.where(np_fa > 0.4)[0]
    # Mean Diffusivity - ID-ing csf
    md_map = img_wrapper('md_masked.nii.gz', 'MD')
    md_map.probe_points(points_to_sample)
    np_md = vtk_to_numpy(md_map.get_array())
    idx_csf = np.where(np_md >= 0.0006)[0]
    if save:
        np.savez(os.path.join(params.points_path, 'special_mesh_midpts.npz'),
                 idx_out=idx_out,
                 idx_grnd=idx_grnd,
                 idx_csf=idx_csf, md_csf=np_md,
                 idx_wm=idx_wm, np_fa=np_fa)
        print('Dumped in special_mesh_midpts.npz')
    del seg, np_seg, md_map, np_md
    print('Number of Tets outside brain vol: ', len(idx_out))
    print('Number of Tets in ground ele: ', len(idx_grnd))
    print('Number of Tets in CSF vol: ', len(idx_csf))
    print('Number of Tets in WMvol-not hard cutoff: ', len(idx_wm))
    return [idx_out, idx_csf, idx_grnd]


def eigen_vectors(points_to_sample, save=False):
    '''useful for computing the anisotropic conductivity tensor'''
    v1 = img_wrapper('v1_masked.nii.gz', 'v1')
    v2 = img_wrapper('v2_masked.nii.gz', 'v2')
    v3 = img_wrapper('v3_masked.nii.gz', 'v3')
    v1.probe_points(points_to_sample)
    v2.probe_points(points_to_sample)
    v3.probe_points(points_to_sample)
    np_v1 = vtk_to_numpy(v1.get_array())
    np_v2 = vtk_to_numpy(v2.get_array())
    np_v3 = vtk_to_numpy(v3.get_array())
    if save:
        np.savez(os.path.join(params.points_path, 'eigen_vecs.npz'),
                 np_v1=np_v1,
                 np_v2=np_v2,
                 np_v3=np_v3)
        print('Dumped Eigen vectors in eigen_vecs.npz')
    print('Done fetching eigen_vecs from the imaging data')
    return np_v1, np_v2, np_v3


if __name__ == '__main__':
    if sys.argv[-1] == 'redo_all':
        print('Preproduction phase')
        print('Previous files will be overwritten')
        points_to_sample = compute_barycenters(save=True)
        idx_out, idx_csf, idx_grnd = id_barycenters(points_to_sample,
                                                    save=True)
        np_v1, np_v2, np_v3 = eigen_vectors(points_to_sample, save=True)
        grid_sel = probe_points_assign(res=1.0, ipsilateral=True, save=True)
    elif sys.argv[-1] == 'default':
        print('Computing just probe_labels and points_wts')
        print('Previous files will be overwritten')
        grid_sel = probe_points_assign(res=1.0, ipsilateral=True, save=True)
        probe_labels = probe_points_labels(grid_sel, save=True)
        probe_wts = compute_point_wts(grid_sel, save=True)
    else:
        print('You feeling lucky punk?')
        pass
