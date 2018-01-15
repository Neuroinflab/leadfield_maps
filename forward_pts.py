
import numpy as np
import dolfin as d
import parameters as params
import conductivity_c as cc
import meshes
import os
import sys


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


def set_solver():
    # "hypre_amg") #"hypre_euclid") # "hypre_amg") # "petsc_amg" "petsc_amg"
    solver = d.KrylovSolver("cg", "hypre_amg")
    solver.parameters["maximum_iterations"] = 1000
    solver.parameters["absolute_tolerance"] = 1E-8
    solver.parameters["error_on_nonconvergence"] = True
    solver.parameters["monitor_convergence"] = True
    # solver.parameters["divergence_limit"] = 1E+6
    # solver.parameters["nonzero_initial_guess"] = True
    d.info(solver.parameters, verbose=True)
    d.set_log_level(d.PROGRESS)
    return solver


def fem_pts(conductivity, pos_list, save_dest, ele_list=None, sel_idx=None):
    NPYSave = False
    HDF5Save = False
    if save_dest.find('.h5') > -1:
        HDF5Save = True
        dump_file = d.HDF5File(d.mpi_comm_world(), save_dest, 'w')
    elif save_dest.find('.npy') > -1:
        if not ele_list:
            print('Expecting ele_list argument')
        else:
            NPYSave = True
    else:
        print('Only .h5 (entire mesh space), .npy (known ele_pos) supported')
        print('Will not save anything this time')
    if not sel_idx:
        sel_idx = range(len(pos_list))
    print('On this process no. pt. srcs = ', len(sel_idx))
    mesh, subdomain, boundaries = meshes.load_meshes()
    sigma = sigma_tensor(mesh, conductivity=conductivity)
    print('Done loading meshes and conductivity')
    V = d.FunctionSpace(mesh, "CG", 2)
    v = d.TestFunction(V)
    u = d.TrialFunction(V)
    dx = d.Measure("dx")(subdomain_data=subdomain)
    # ds = d.Measure("ds")(subdomain_data=boundaries)
    a = d.inner(sigma * d.grad(u), d.grad(v))*dx
    L = d.Constant(0)*v*dx()
    A = d.assemble(a)
    # Surface of the grnd ele = 1030
    bc = d.DirichletBC(V, d.Constant(0), boundaries, 1030)

    for curr_idx in sel_idx:
        solver = set_solver()
        phi = d.Function(V)
        x = phi.vector()
        print('Started computing for,at: ', curr_idx, pos_list[curr_idx])
        b = d.assemble(L)
        bc.apply(A, b)
        xx, yy, zz = pos_list[curr_idx]
        point = d.Point(xx, yy, zz)
        delta = d.PointSource(V, point, 1.)
        delta.apply(b)
        solver.solve(A, x, b)
        # file = d.File("pots_anis.pvd")
        # file << phi
        if HDF5Save:
            dump_file.write(x.array(), str(curr_idx))
            dump_file.flush()
        if NPYSave:
            vals = extract_pots(phi, np.array(ele_list))
            np.save(save_dest, vals)
        print('Finished computing for :', curr_idx)
    return


if __name__ == '__main__':
    if len(sys.argv) == 3:
        print('Running process ', sys.argv[-1], 'of ', sys.argv[-2])
        pos_list, conductivity, path, sbspt = params.default_run('homogeneous')
        num_proc = int(sys.argv[-2])
        proc_idx = int(sys.argv[-1])
        proc_vals = np.linspace(0, len(pos_list), num_proc + 1).astype(int)
        # diff_proc = proc_vals[proc_idx] - proc_vals[proc_idx-1]
        pt_idxs = range(proc_vals[proc_idx - 1], proc_vals[proc_idx])
        save_dest = os.path.join(path, sbspt + str(num_proc)
                                 + '_' + str(proc_idx) + '.h5')
        fem_pts(conductivity, pos_list, save_dest, sel_idx=pt_idxs)
    else:
        print('Running default, 1 process, three point locs')
        pos_list = [[5.196, 22.913, -4.9957], [8.4, 31.4, -6.151],
                    [5.5945, 22.699, -5.6637]]
        conductivity = 'anisotropic'
        path = params.results_path
        save_dest = os.path.join(params.results_path, 'test_del_a.h5')
        fem_pts(conductivity, pos_list, save_dest)
