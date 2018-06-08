# from scipy.spatial.distance import pdist, squareform
# from sklearn import datasets
# from fastcluster import linkage

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import colors as clr
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


#conds = ['homogeneous', 'inhomogeneous', 'anisotropic']
conds = ['anisotropic']
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
    points = [0]
    colors = []
    labels = []
    labels_txt, labels_clr = fetch_label_dicts()
    ii = 0
    for kk, val in enumerate(unique_list_lens):
        points.append(val + ii)
        qq = np.where(point_wts == unique[kk])
        tt = labels_idx[qq[0][0]]
        colors.append([cc / 255. for cc in labels_clr[tt]])
        labels.append(labels_txt[tt])
        ii += val
    points = np.array(points)
    for pp in range(len(colors)):
        ax.plot([1000, 1000], [points[pp], points[pp+1]],
                color=colors[pp], lw=5, clip_on=True)
        if pp < 10:
            ax.text(1050, (points[pp] + points[pp+1])/2., labels[pp])
    return ax
    
# fig = plt.figure()
# im = plt.imshow(phi_sorted)
# #                norm=LogNorm(vmin=0.01, vmax=np.max(4.)))
# plt.colorbar()

#gs = gridspec.GridSpec(2, 1, height_ratios=[1, 0.05])
fig = plt.figure(figsize=(12, 10))
for jj, case in enumerate(conds):
    ax = plt.subplot(111)
    data = phi_mats[jj][phi_seg_list, ][:, phi_seg_list, ]
    print(np.min(data), np.max(data))
    np.fill_diagonal(data, 0)
    #im = ax.pcolormesh(data, norm=clr.LogNorm(vmin=np.min(data), vmax=np.max(data)),  cmap=plt.cm.gray)#, vmin=0.0, vmax=3.0)
    im = ax.pcolormesh(data, norm=clr.PowerNorm(gamma=1./2.),  cmap=plt.cm.gray_r)#, vmin=0.0, vmax=3.0)
    draw_color_line(ax)
    ax.set_aspect('equal')
    ax.set_xlabel('Electrode number')
    ax.set_ylabel('Source number')
    ax.spines['top'].set_visible(False)                                                                
    ax.spines['right'].set_visible(False)
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("bottom", size="5%", pad=.7)
    #cax = plt.subplot(gs[1, jj])
    plt.colorbar(im, cax=cax, orientation='horizontal')
    

plt.savefig('points_v_points.png', dpi=300)
#plt.show()
