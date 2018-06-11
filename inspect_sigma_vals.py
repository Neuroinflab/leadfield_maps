import dolfin as d
import parameters as params
import conductivity_c as cc
import meshes
import os

def sigma_principle_axis(mesh, conductivity):
    def assign(path_to_xml):
        return d.MeshFunction("double", mesh, path_to_xml)
    if conductivity == 'anisotropic':
        c00 = assign(os.path.join(params.anis_path, "sigma_anis_d0.xml.gz"))
        print('Anisotropic case')
    else:
        if conductivity == 'homogeneous':
            c00 = assign(os.path.join(params.hom_path, "sigma_hom.xml.gz"))
            print('Homogeneous case')
        elif conductivity == 'inhomogeneous':
            c00 = assign(os.path.join(params.inhom_path, "sigma_inhom.xml.gz"))
            print('Inhomogeneous case')
        # C - code assignment remains the same for inhom and hom.
    return c00



print('Loading meshes')
conductivity = 'anisotropic'
mesh =  d.Mesh(os.path.join(params.mesh_path, "mesh_setup.xml"))
c00 = sigma_principle_axis(mesh, conductivity)

dump_file = d.File(os.path.join(params.anis_path, "sigmas_princ.pvd"))
dump_file << c00
