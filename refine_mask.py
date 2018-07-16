import nifti
#import nibabel as nib
import os
import numpy as np

exclude_labels = [76, 44, 42, 41]
exclude_names =  ['spinal trigeminal tract',
                  'inner ear', 'optic tract and optic chiasm',
                  'optic nerve']
old_mask = nifti.NiftiImage(os.path.join("../dti_package/brain_mask.nii.gz"))
d = nifti.NiftiImage(os.path.join("../dti_package/brain_atlas.nii.gz"))

mask = old_mask.data
r = d.data

# for ii, label in enumerate(exclude_labels):
#     print 'Excluding: ', exclude_names[ii]
#     mask[r == label] = 0

# # mask[r > 0] = 1  # Create new refined mask

# save_as = os.path.join("../dti_package/brain_mask_refined.nii.gz")
# nifti.NiftiImage(mask, header=d.header).save(save_as)



include_hippo_labels = [95, 96, 97, 98, 99, 100]
hippo_names = ['cornu ammonis 1', 'dentate gyrus',
               'cornu ammonis 2', 'cornu ammonis 3',
               'fasciola cinereum', 'subiculum']

mask[:] = 0
for ii, label in enumerate(include_hippo_labels):
    print('Hippocampus : ', hippo_names[ii])
    mask[r == label] = 1

save_as = os.path.join("../dti_package/hippocampus.nii.gz")
nifti.NiftiImage(mask, header=d.header).save(save_as)
