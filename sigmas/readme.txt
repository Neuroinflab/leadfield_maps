Folder has sigma values as tensors.

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

S = np.matrix((np_v1[jj], np_v2[jj], np_v3[jj]))
 # This is correct. S is eigen vectors are row vectors
sigma_mat = np.around(S.T * sig_wm_matrix * S, decimals=4)

hom is constant of sigma_gm
inhom is scaled according to FA(P_i) between s_wm and s_gm
Anis is scaled according to sigma_wm_long and simga_wm_trans along the principle eigen vectors.
