import os
import h5py
import numpy as np
import parameters as params

# Fetch tansmembrane current values.
def find_map_idx(h, pop_name, field_name):
    '''Find the corresponding, locations of morphology in field'''
    mor = h['/data/static/morphology/' + pop_name]
    i_data = h['/data/uniform/' + pop_name + '/' + field_name]
    mor_names = mor.dims[0].values()[0]  # entry in /map
    i_names = i_data.dims[0].values()[0]  # entry in /map
    if np.array_equal(i_names, mor_names):
        idx = range(i_names.shape[0])
    else:
        idx = [np.where(i_names.value == entry)[0][0] for entry in mor_names]
    return idx

def fetch_pop_pots(pop_name):
    pop_idx = pop_names.index(pop_name)
    if pop_idx > 0 :
        start = np.sum(total_cmpts[:pop_idx])
        end = np.sum(total_cmpts[:pop_idx + 1])
    else:
        start = 0
        end = total_cmpts[0]
    pots_pop = pots[:, start:end]
    return pots_pop

def write_hdf5(path, ele_pots, ele_coords, start_t, end_t):
    with h5py.File(os.path.join(path, 'traub_pots_v_t.h5'), 'w') as hf:
        hf.create_dataset("ele_pos",  data=ele_coords)
        for ii in range(start_t, end_t):
            tt = ii-start_t
            hf.create_dataset(str(ii), data=ele_pots[:, tt])
        ele_coords = np.array(ele_coords)
        hf.create_dataset('x', data=ele_coords[:, 0])
        hf.create_dataset('y', data=ele_coords[:, 1])
        hf.create_dataset('z', data=ele_coords[:, 2])
    return


num_cmpts, cell_range, pop_names, num_cells, total_cmpts = params.load_traub_morph_props()
# Fetch electrode points
ele_coords, conductivity, path, sbspt = params.default_run('anisotropic')
# fetch pot
phi_mat_fname = sbspt + conductivity + '_traub_phi_mat.npy'
pots = np.load(os.path.join(path, phi_mat_fname))  # shape = ele x traub_pos_list
start_t, end_t  = 2750, 3750
field_names = ['i']

tot_ele_pots = np.zeros((len(ele_coords), end_t - start_t))
h = h5py.File('/home/chaitanya/ratty/KESI/traub_syn.h5', 'r')
for pop_name in pop_names:   # ['pyrFRB23']
    pop_pots = fetch_pop_pots(pop_name)
    src_time = np.zeros((pop_pots.shape[1], end_t - start_t))
    for field_name in field_names:
        idx = find_map_idx(h, pop_name, field_name)
        this_field = h['/data/uniform/' + pop_name + '/' + field_name].value
        src_time += this_field[idx][:, start_t:end_t]  # Order according to correct indices
    tot_ele_pots += np.dot(pop_pots, src_time)
    print('Done for population :', pop_name)

write_hdf5(path, tot_ele_pots, ele_coords, start_t, end_t)   # Either or
#np.save(os.path.join(path, 'traubs_pots_v_t.npy'), tot_ele_pots)
