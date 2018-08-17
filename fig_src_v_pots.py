# from scipy.spatial.distance import pdist, squareform
# from sklearn import datasets
# from fastcluster import linkage

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import colors as clr
from figure_properties import *
import numpy as np
import parameters as params
import scipy.spatial
import scipy
import csv
import os


def inv_distance(src_pos, ele_pos):
    '''computes the inverse distance between src_pos and ele_pos'''
    dist_matrix = np.zeros((src_pos.shape[0], ele_pos.shape[0]))
    for ii, electrode in enumerate(ele_pos):
        dist_matrix[:, ii] = scipy.spatial.distance.cdist(src_pos,
                                                          electrode.reshape(1, 3)).flatten()
    dist_matrix = 1 / dist_matrix  # inverse distance matrix
    return dist_matrix


def reciprocity_norm_measure(data):
    np.fill_diagonal(data, 0)
    mat = np.matrix(mat)
    symm = 0.5 * (mat  + mat.T)
    asymm = 0.5 * (mat - mat.T)
    return np.linalg.norm(asymm) / np.linalg.norm(symm)


def rdm_vals(anis, data):
    rdm = np.zeros(anis.shape[0])  # number of source locations
    for ii in range(len(rdm)):
        arr = anis[ii, :]
        cmpr = data[ii, :]
        this_anis = np.hstack((arr[:ii], arr[ii+1:]))  # Skip self electrode
        this_cmpr = np.hstack((cmpr[:ii], cmpr[ii+1:]))
        ani_denom = np.linalg.norm(this_anis)
        cmpr_denom = np.linalg.norm(this_cmpr)
        sum_diff = 0
        for jj in range(len(this_anis)):
            f = this_anis[jj] / ani_denom
            s = this_cmpr[jj] / cmpr_denom
            sum_diff += np.abs((f-s))
        rdm[ii] = np.sqrt(sum_diff)
    return rdm
        
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
#conds = ['anisotropic']
phi_mats = {}
for case in conds:
    pos_list, conductivity, path, sbspt = params.default_run(case)
    fname = sbspt + conductivity + '_phi_mat.npy'
    data = np.load(os.path.join(path, fname))
    np.fill_diagonal(data, 0)
    phi_mats[case] = data

# infinite homogeneous case
phi_mats_hom = inv_distance(np.array(pos_list), np.array(pos_list)).T
phi_mats_hom *= (1 / (4 * np.pi * 0.3))
np.fill_diagonal(phi_mats_hom, 0)
print('infinite', np.min(phi_mats_hom), np.max(phi_mats_hom))
conds.append('infinite')
phi_mats['infinite'] =phi_mats_hom

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


def draw_color_line(ax, orientation='horizontal'):
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
        if orientation == 'vertical':
            ax.plot([1000, 1000], [points[pp], points[pp+1]],
                    color=colors[pp], lw=5, clip_on=True)
            if pp < 10:
                ax.text(1050, (points[pp] + points[pp+1])/2., labels[pp])
        elif orientation == 'horizontal':
            ax.plot([points[pp], points[pp+1]], [-20, -20],
                    color=colors[pp], lw=5, clip_on=True)
            if pp < 10:
                if len(labels[pp]) > 15:
                    size_txt = 8
                else:
                    size_txt = 10
                if len(labels[pp].split()) > 3:
                    labels[pp] = ' '.join(labels[pp].split()[-3:])
                ax.text((points[pp] + points[pp+1])/2., -70,  labels[pp],
                        horizontalalignment='right', rotation=45, size=size_txt)
    return ax
    
# fig = plt.figure()
# im = plt.imshow(phi_sorted)
# #                norm=LogNorm(vmin=0.01, vmax=np.max(4.)))
# plt.colorbar()

gs = gridspec.GridSpec(1, 3, width_ratios=[1, 0.05, 1], height_ratios=[1])
fig = plt.figure(figsize=(20, 10))
case = 'anisotropic'
jj = conds.index(case)
#ax = plt.subplot(111)
ax = plt.subplot(gs[0, 0])
data = phi_mats[case][phi_seg_list, ][:, phi_seg_list, ]
print(case, np.min(data), np.max(data))
np.fill_diagonal(data, 0)
im = ax.pcolormesh(data, vmax=3.,  cmap=plt.cm.Greens, vmin=0.0)
#im = ax.pcolormesh(data, norm=clr.PowerNorm(gamma=1./2.),  cmap=plt.cm.gray_r)#, vmin=0.0, vmax=3.0)
draw_color_line(ax)
ax.set_aspect('equal')
#ax.set_xlabel('')
ax.set_xticks([])
ax.set_ylabel('Electrode number')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_visible(False)
ax.spines['bottom'].set_visible(False)    
plt.title('Leadfield matrix - anisotropic conductivity')
divider = make_axes_locatable(ax)
#cax = divider.append_axes("right", size="5%", pad=1.)
cax = plt.subplot(gs[0, 1])
plt.colorbar(im, cax=cax, orientation='vertical', extend='max')
ax.text(1200, 400, 'A/S')
ax = plt.subplot(gs[0, 2])
hom = rdm_vals(phi_mats['anisotropic'], phi_mats['homogeneous'])
inhom = rdm_vals(phi_mats['anisotropic'], phi_mats['inhomogeneous'])
#infi = rdm_vals(phi_mats['anisotropic'], phi_mats['infinite'])
plt.plot(hom, label='Anis-Hom')
plt.plot(inhom, label='Anis-Inhom')
#plt.plot(infi, label='Anis-Inf')
plt.legend()
plt.show()

    
#plt.savefig('points_v_points.png', dpi=300)
#plt.show()
