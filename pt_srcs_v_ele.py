import os
import csv
import parameters as params
import numpy as np
import vtk_functions as vtk_utils
from vtk.util.numpy_support import vtk_to_numpy


def probe_points_labels(points, save=False):
    brain_atlas_path = os.path.join(params.dti_path, 'brain_atlas.nii.gz')
    seg = vtk_utils.image_wrapper(brain_atlas_path, 'mask')
    seg.probe_points(points)
    np_labels = vtk_to_numpy(seg.get_array())
    if save:
        np.save(os.path.join(params.points_path,
                             'probe_pts_labels.npy'), np_labels)
    return np_labels


def fetch_label_dicts():
    labels_txt = {}
    labels_clr = {}
    with open(os.path.join(params.points_path, 'labels.csv'), 'rb') as f:
        csv_file = csv.reader(f, delimiter=',')
        for row in csv_file:
            labels_txt[int(row[0])] = row[4]
            labels_clr[int(row[0])] = (int(row[1]), int(row[2]), int(row[3]))
    return labels_txt, labels_clr


points = np.load(os.path.join(params.points_path, 'probe_points_ipsi_1.0.npy'))
labels_idx = probe_points_labels(points)
unique_idx, counts = np.unique(labels_idx, return_counts=True)
labels_txt, labels_clr = fetch_label_dicts()  # This is already unique

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

