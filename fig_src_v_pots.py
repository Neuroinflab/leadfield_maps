# from scipy.spatial.distance import pdist, squareform
# from sklearn import datasets
# from fastcluster import linkage

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import parameters as params
import csv
import os


def fetch_label_dicts():
    labels_txt = {}
    labels_clr = {}
    with open(os.path.join(params.points_path, 'labels.csv'), 'rb') as f:
        csv_file = csv.reader(f, delimiter=',')
        for row in csv_file:
            labels_txt[int(row[0])] = row[4]
            labels_clr[int(row[0])] = (int(row[1]), int(row[2]), int(row[3]))
    return labels_txt, labels_clr


conds = ['homogeneous', 'inhomogeneous', 'anisotropic']
phi_mats = []
for case in conds:
    pos_list, conductivity, path, sbspt = params.default_run(case)
    fname = sbspt + conductivity + '_phi_mat.npy'
    data = np.load(os.path.join(path, fname))
    np.fill_diagonal(data, 0)
    phi_mats.append(data)

labels_idx = np.load(os.path.join(params.points_path, 'probe_pts_labels.npy'))
point_wts = np.load(os.path.join(params.points_path, 'probe_pts_wts.npy'))
unique_idx, counts = np.unique(point_wts, return_counts=True)
unique = [x for _, x in sorted(zip(counts, unique_idx),
                               key=lambda pair: pair[0],
                               reverse=True)]
phi_seg_list = []
unique_list_lens = []
for u_val in unique:
    idx = np.where(point_wts == u_val)
    phi_seg_list.extend(idx[0].tolist())
    unique_list_lens.append(len(idx[0]))


def draw_color_line(ax, orientation='vertical'):
    points = [[0, 0]]
    colors = []
    labels = []
    labels_txt, labels_clr = fetch_label_dicts()
    ii = 0
    for kk, val in enumerate(unique_list_lens):
        if orientation == 'vertical':
            points.append([0, val + ii])
        else:
            points.append([val + ii, 0])
        qq = np.where(point_wts == unique[kk])
        tt = labels_idx[qq[0][0]]
        colors.append([cc / 255. for cc in labels_clr[tt]])
        labels.append(labels_txt[tt])
        ii += val
    points = np.array(points)
    print(points[:, 1])
    for pp in range(len(colors)):
        print(points[pp:pp+1, 1])
        ax.plot([points[pp, 0], points[pp+1, 0]],
                [points[pp, 1], points[pp+1, 1]],
                color=colors[pp], linewidth=5.)
    return ax
    
# fig = plt.figure()
# im = plt.imshow(phi_sorted)
# #                norm=LogNorm(vmin=0.01, vmax=np.max(4.)))
# plt.colorbar()

gs = gridspec.GridSpec(2, 3, height_ratios=[1, 0.05])
for jj, case in enumerate(conds):
    ax = plt.subplot(gs[0, jj])
    data = phi_mats[jj][phi_seg_list, ][:, phi_seg_list, ]
    np.fill_diagonal(data, 0)
    im = plt.imshow(data[::-1], cmap=plt.cm.Greens, vmin=0.0, vmax=4.0)
    cax = plt.subplot(gs[1, jj])
    plt.colorbar(im, cax=cax, orientation='horizontal')
    draw_color_line(ax)
    
plt.show()
