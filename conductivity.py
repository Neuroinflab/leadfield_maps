from math import erf
import dolfin as d
import numpy as np
import parameters as params
import os
import sys

sig_ele = 1e+07
sig_air = 1e-14
sig_zero = 0.00

sig_csf = 1.5  # See M.Rullmann, NeuroImage 44 (2009) Table 1
sig_brain = 0.33
sig_wm = 0.142
sig_wm_long = 0.65  # See C.H.Wolters, NeuroImage 30 (2006) Table 1
sig_wm_trans = 0.065


def scale_gm_wm(fa, s_gm, s_wm):
    # (x) -1. --> 1. :: (fa) 0.3 --> 0.5
    # (y) -1. --> 1. :: (scale) s_gm --> s_wm
    s_diff = s_gm - s_wm
    k = (np.sqrt(np.pi) * 10) / 2.
    scaled_val = (s_diff * ((-1 * erf(k * (fa - 0.4)) + 1) / 2.))
    return s_wm + scaled_val


def tensor_components(mesh):
    ''' c00 c01 c02
            c11 c12
                c22 '''
    c00 = d.MeshFunction("double", mesh, 3)
    c01 = d.MeshFunction("double", mesh, 3)
    c02 = d.MeshFunction("double", mesh, 3)
    c11 = d.MeshFunction("double", mesh, 3)
    c12 = d.MeshFunction("double", mesh, 3)
    c22 = d.MeshFunction("double", mesh, 3)
    return c00, c01, c02, c11, c12, c22


def compute_conductivity(conductivity, save=False):
    print('Computing sigma tensor for case :', conductivity)
    mesh = params.load_just_mesh()
    cells_all = d.cells(mesh)
    idx_out, idx_grnd, idx_csf, np_fa = params.load_barycenters_ids()
    points_to_sample = params.load_barycenters()
    idx_all = np.arange(len(points_to_sample))
    if conductivity == 'anisotropic':
        np_v1, np_v2, np_v3 = params.load_eigen_vectors()
        longi_sigma = []
        trans_sigma = []
        for jj in idx_all:
            longi_sigma.append(scale_gm_wm(np_fa[jj], sig_brain, sig_wm_long))
            trans_sigma.append(scale_gm_wm(np_fa[jj], sig_brain, sig_wm_trans))
        print 'Did the scaling, now assigning vals'
        c00, c01, c02, c11, c12, c22 = tensor_components(mesh)
        for jj in idx_all:
            cell = cells_all.next()
            if jj in idx_grnd:  # Electrode
                c00[cell] = sig_ele
                c11[cell] = sig_ele
                c22[cell] = sig_ele
                c01[cell] = sig_zero
                c02[cell] = sig_zero
                c12[cell] = sig_zero
            elif jj in idx_out:  # Pt outside brain volume
                c00[cell] = sig_air
                c01[cell] = sig_zero
                c02[cell] = sig_zero
                c11[cell] = sig_air
                c12[cell] = sig_zero
                c22[cell] = sig_air
            elif jj in idx_csf:
                c00[cell] = sig_csf
                c01[cell] = sig_zero
                c02[cell] = sig_zero
                c11[cell] = sig_csf
                c12[cell] = sig_zero
                c22[cell] = sig_csf
            else:
                sig_wm_matrix = np.zeros((3, 3))
                np.fill_diagonal(sig_wm_matrix,
                                 (longi_sigma[jj],
                                  trans_sigma[jj],
                                  trans_sigma[jj]))
                S = np.matrix((np_v1[jj], np_v2[jj], np_v3[jj]))
                # This is correct. S is eigen vectors are row vectors
                sigma_mat = np.around(S.T * sig_wm_matrix * S, decimals=4)
                c00[cell] = sigma_mat[0, 0]
                c11[cell] = sigma_mat[1, 1]
                c22[cell] = sigma_mat[2, 2]
                c01[cell] = sigma_mat[0, 1]
                c02[cell] = sigma_mat[0, 2]
                c12[cell] = sigma_mat[1, 2]

        if save:
            c00_file = d.File(os.path.join(params.anis_path,
                                           "sigma_anis_d0.xml.gz"))
            c01_file = d.File(os.path.join(params.anis_path,
                                           "sigma_anis_d1.xml.gz"))
            c02_file = d.File(os.path.join(params.anis_path,
                                           "sigma_anis_d2.xml.gz"))
            c11_file = d.File(os.path.join(params.anis_path,
                                           "sigma_anis_d3.xml.gz"))
            c12_file = d.File(os.path.join(params.anis_path,
                                           "sigma_anis_d4.xml.gz"))
            c22_file = d.File(os.path.join(params.anis_path,
                                           "sigma_anis_d5.xml.gz"))
            c00_file << c00
            c01_file << c01
            c02_file << c02
            c11_file << c11
            c12_file << c12
            c22_file << c22

    elif conductivity == 'inhomogeneous':
        c = d.MeshFunction("double", mesh, 3)
        fa_scaled_sigma = []
        for jj in idx_all:
            fa_scaled_sigma.append(scale_gm_wm(np_fa[jj], sig_brain, sig_wm))
        print 'Did the scaling, now assigning vals'
        for jj in idx_all:
            cell = cells_all.next()
            if jj in idx_grnd:
                c[cell] = sig_ele
            elif jj in idx_out:
                c[cell] = sig_air
            elif jj in idx_csf:
                c[cell] = sig_csf
            else:
                c[cell] = fa_scaled_sigma[jj]
        if save:
            c00_file = d.File(os.path.join(params.inhom_path,
                                           "sigma_inhom.xml.gz"))
            c00_file << c

    elif conductivity == 'homogeneous':
        c = d.MeshFunction("double", mesh, 3)
        for jj in idx_all:
            cell = cells_all.next()
            if jj in idx_grnd:
                c[cell] = sig_ele
            elif jj in idx_out:
                c[cell] = sig_air
            elif jj in idx_csf:
                c[cell] = sig_csf
            else:
                c[cell] = sig_brain
        if save:
            c00_file = d.File(os.path.join(params.hom_path,
                                           "sigma_hom.xml.gz"))
            c00_file << c


if __name__ == '__main__':
    if sys.argv[-1] == 'redo_all':
        print('Preproduction phase')
        print('Previous files will be overwritten')
        compute_conductivity('anisotropic', save=True)
        compute_conductivity('inhomogeneous', save=True)
        compute_conductivity('homogeneous', save=True)
    else:
        print('What else is there to do?')
        pass
