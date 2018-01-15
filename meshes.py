
import os
import dolfin as d
import parameters as params


def load_meshes():
    mesh = d.Mesh(os.path.join(params.mesh_path, "mesh_setup.xml"))
    subdomain = d.MeshFunction("size_t", mesh,
                               os.path.join(params.mesh_path,
                                            "mesh_setup_physical_region.xml"))
    boundaries = d.MeshFunction("size_t", mesh,
                                os.path.join(params.mesh_path,
                                             "mesh_setup_facet_region.xml"))
    return mesh, subdomain, boundaries


def load_just_mesh():
    mesh = d.Mesh(os.path.join(params.mesh_path, "mesh_setup.xml"))
    return mesh
