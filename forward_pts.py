import numpy as np
import dolfin as d
import parameters as params
import conductivity_c as cc
import os


def sigma_tensor(mesh, conductivity):
    def assign(path_to_xml):
        return d.MeshFunction("double", mesh, path_to_xml)
    if conductivity == 'anisotropic':
        c00 = assign(os.path.join(params.anis_path, "sigma_anis_d0.xml.gz"))
        c01 = assign(os.path.join(params.anis_path, "sigma_anis_d1.xml.gz"))
        c02 = assign(os.path.join(params.anis_path, "sigma_anis_d2.xml.gz"))
        c11 = assign(os.path.join(params.anis_path, "sigma_anis_d3.xml.gz"))
        c12 = assign(os.path.join(params.anis_path, "sigma_anis_d4.xml.gz"))
        c22 = assign(os.path.join(params.anis_path, "sigma_anis_d5.xml.gz"))
        print('Anisotropic case')
        c = d.Expression(cppcode=cc.anisotropy_code, degree=0)
        c.c00 = c00
        c.c01 = c01
        c.c02 = c02
        c.c11 = c11
        c.c12 = c12
        c.c22 = c22
        C = d.as_matrix(((c[0], c[1], c[2]),
                         (c[1], c[3], c[4]),
                         (c[2], c[4], c[5])))
    else:
        if conductivity == 'homogeneous':
            c00 = assign(os.path.join(params.hom_path, "sigma_hom.xml.gz"))
            print('Homogeneous case')
        elif conductivity == 'inhomogeneous':
            c00 = assign(os.path.join(params.inhom_path, "sigma_inhom.xml.gz"))
            print('Inhomogeneous case')
        # C - code assignment remains the same for inhom and hom.
        c = d.Expression(cppcode=cc.homogeneous_code, degree=0)
        c.c00 = c00
        C = d.as_matrix(((c[0], d.Constant(0.0), d.Constant(0.0)),
                         (d.Constant(0.0), c[0], d.Constant(0.0)),
                         (d.Constant(0.0), d.Constant(0.0), c[0])))
    return C


def extract_pots(phi, positions):
    compt_values = np.zeros(positions.shape[0])
    for ii in range(positions.shape[0]):
        compt_values[ii] = phi(positions[ii, :])
    return compt_values


def fem_pts(conductivity, pos_list, save_as, ele_list=None):
    if save_as.find('.h5') > -1:
        HDF5Save = True
        dump_file = d.HDF5File(d.mpi_comm_world(), save_as, 'w')
    elif save_as.find('.xml.gz'):
        HDF5Save = False
        dump_file = d.File(save_as)
    else:
        print('Only .h5 or .xml.gz supported')
    mesh, subdomain, boundaries = params.load_meshes()
    sigma = sigma_tensor(mesh, conductivity=conductivity)
    print('Done loading meshes and sigma')
    V = d.FunctionSpace(mesh, "CG", 2)
    v = d.TestFunction(V)
    u = d.TrialFunction(V)

    phi = d.Function(V)
    dx = d.Measure("dx")(subdomain_data=subdomain)
    ds = d.Measure("ds")(subdomain_data=boundaries)

    a = d.inner(sigma * d.grad(u), d.grad(v))*dx
    L = d.Constant(0)*v*dx()
    A = d.assemble(a)
    # "hypre_amg") #"hypre_euclid") # "hypre_amg") # "petsc_amg" "petsc_amg"
    solver = d.KrylovSolver("cg", "hypre_amg")
    solver.parameters["maximum_iterations"] = 1000
    solver.parameters["absolute_tolerance"] = 1E-8
    solver.parameters["error_on_nonconvergence"] = True
    solver.parameters["monitor_convergence"] = True
    # solver.parameters["divergence_limit"] = 1E+6
    # solver.parameters["nonzero_initial_guess"] = True
    x = phi.vector()
    d.info(solver.parameters, verbose=True)
    d.set_log_level(d.TRACE)
    for pos_idx, position in enumerate(pos_list):
        print('Started computing for,at: ', pos_idx, position)
        b = d.assemble(L)
        bc = d.DirichletBC(V, d.Constant(0), boundaries, 1030)
        # Surface of the grnd ele = 1030
        bc.apply(A, b)
        xx, yy, zz = position
        point = d.Point(xx, yy, zz)
        delta = d.PointSource(V, point, 1.)
        delta.apply(b)
        solver.solve(A, x, b)
        # file = File("pots_anis.pvd")
        # file << phi
        if HDF5Save:
            dump_file.write(x.array(), str(pos_idx))
            dump_file.flush()
        else:
            dump_file << x
        print('Finished computing for :', pos_idx)
    # if HDF5Save:
    #     dump_file.write(phi, 'phi_func')
    #     dump_file.close()
    # else:
    #     pass
    return


if __name__ == '__main__':
    # import sys
    # import csv
    # pos_list = []
    # with open('traub_post_transform.csv', 'rb') as csvfile:
    #     next(csvfile, None)
    #     spamreader = csv.reader(csvfile, delimiter=',')
    #     for row in spamreader:
    #         pos_list.append([float(row[1]), float(row[2]), float(row[3])])
    pos_list = [[5.196, 22.913, -4.9957]]
    fem_pts('anisotropic', pos_list, 'test_ani_i.h5')
    # if len(sys.argv) == 3:
    #     print('Running process ', sys.argv[-1], 'of ', sys.argv[-2])
    #     big_loop(pos_list, int(sys.argv[-1]), int(sys.argv[-2]))
    # else:
    #     print('Running default, 1 process, first two point locs')
    #     big_loop(pos_list[0:2], 1, 1)
