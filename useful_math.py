from mpl_toolkits.mplot3d import Axes3D
import mpl_toolkits.mplot3d as m3d
from matplotlib import cm
import matplotlib.tri as mtri
import matplotlib.pyplot as plt
from math import sin, cos, acos, sqrt
import numpy as np
import csv


def get_plane(pt2, pt3):
    ''' Get a plane passing through two points and origin '''
    A = pt2[1] * pt3[2] - pt3[1] * pt2[2]
    B = pt2[2] * pt3[0] - pt3[2] * pt2[0]
    C = pt2[0] * pt3[1] - pt3[0] * pt2[1]
    D = 0
    return A, B, C, D


def is_on_plane(A, B, C, D, pt):
    ''' check if point on a plane error is 1e-5 '''
    return np.abs((A * pt[0]) + (B * pt[1]) + (C * pt[2]) - D) < 1e-5


def get_plane_thr_pt(A, B, C, pt):
    ''' get the factor D for a given plane normal and pt passing'''
    return (A * pt[0]) + (B * pt[1]) + (C * pt[2])


def get_irreg_pts_plane(A, B, C, D, grid_x, grid_y):
    ''' Obtain an irregular (obtuse) grid of pts on a given plane'''
    num_pts = grid_x.shape[0]**2
    grid_z = np.zeros((num_pts)) - D
    grid_z += A * grid_x.reshape(num_pts,)
    grid_z += B * grid_y.reshape(num_pts,)
    grid_z /= -1. * C
    return grid_z.reshape(grid_x.shape[0], grid_x.shape[0])


def get_perp_plane(A, B, C, pt2, pt3):
    '''get a plane perperdicular to a given plane and passing thr 2pts'''
    # Px+Qy+Rz = S
    P = C * (pt3[1] - pt2[1]) - B * (pt3[2] - pt2[2])
    Q = A * (pt3[2] - pt2[2]) - C * (pt3[0] - pt2[0])
    R = B * (pt3[0] - pt2[0]) - A * (pt3[1] - pt2[1])
    S = P * pt2[0] + Q * pt2[1] + R * pt2[2]
    return P, Q, R, S


def obtain_pt_on_line(pt1, pt2, t=0):
    ''' Obtain the point on a line passing through two points'''
    diff_line = [ii - jj for ii, jj in zip(pt1, pt2)]
    return [ii + (t * jj) for ii, jj in zip(pt1, diff_line)]


def rotate_abt_normal(vec, normal):
    ''' assumes the rotation is about z axis - points in z=0 plane
    https://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula'''
    k = np.array((0, 0, 1))
    normal = np.array(normal)
    vec = np.array(vec)
    mag = np.linalg.norm(k) * np.linalg.norm(normal)
    theta = np.arccos(np.dot(k, normal) / mag)

    vec_rot = (vec * np.cos(theta)) + (np.cross(k, vec) * np.sin(theta)) + \
              (k * np.dot(k, vec) * (1 - np.cos(theta)))
    return vec_rot


def normalize(v, tolerance=0.00001):
    '''Taken from
    https://stackoverflow.com/questions/4870393/rotating-coordinate-system-via-a-quaternion'''
    mag2 = sum(n * n for n in v)
    if abs(mag2 - 1.0) > tolerance:
        mag = sqrt(mag2)
        v = tuple(n / mag for n in v)
    return v


def q_mult(q1, q2):
    '''quaterion multiplication'''
    w1, x1, y1, z1 = q1
    w2, x2, y2, z2 = q2
    w = w1 * w2 - x1 * x2 - y1 * y2 - z1 * z2
    x = w1 * x2 + x1 * w2 + y1 * z2 - z1 * y2
    y = w1 * y2 + y1 * w2 + z1 * x2 - x1 * z2
    z = w1 * z2 + z1 * w2 + x1 * y2 - y1 * x2
    return w, x, y, z


def q_conjugate(q):
    '''quaterion conjugate'''
    w, x, y, z = q
    return (w, -x, -y, -z)


def qv_mult(q1, v1):
    '''quaterion vector multiplcation'''
    q2 = (0.0,) + v1
    return q_mult(q_mult(q1, q2), q_conjugate(q1))[1:]


def axisangle_to_q(v, theta):
    '''Angle in rad and vector in quaterion'''
    v = normalize(v)
    x, y, z = v
    theta /= 2
    w = cos(theta)
    x = x * sin(theta)
    y = y * sin(theta)
    z = z * sin(theta)
    return w, x, y, z


def q_to_axisangle(q):
    '''Quaterion to axis angle'''
    w, v = q[0], q[1:]
    theta = acos(w) * 2.0
    return normalize(v), theta


def obtain_angle_btwn_vecs(v1, v2):
    '''Angle between 2 vectors'''
    mag = np.linalg.norm(v1) * np.linalg.norm(v2)
    return np.arccos(np.dot(v1, v2) / mag)


def obtain_n_pts_btwn_pts(start_pt, end_pt, dist_btwn_gridpts):
    '''Equidisant points between 2 points'''
    dist = np.linalg.norm(np.array(start_pt) - np.array(end_pt))
    num_planes = np.int(np.ceil(dist / dist_btwn_gridpts)) - 1
    pts_on_line = []
    pts_on_line.append(start_pt)
    for ii in range(1, num_planes):
        m = float(ii)
        n = float(num_planes - ii)
        x_coord = ((m * end_pt[0]) + (n * start_pt[0])) / num_planes
        y_coord = ((m * end_pt[1]) + (n * start_pt[1])) / num_planes
        z_coord = ((m * end_pt[2]) + (n * start_pt[2])) / num_planes
        pts_on_line.append((x_coord, y_coord, z_coord))
    pts_on_line.append(end_pt)
    return pts_on_line


def circular_grid_on_z(rad, num_pts):
    gt = np.complex(0, num_pts)  # density  of the points
    xx, yy = np.mgrid[-rad:rad:gt, -rad:rad:gt]
    idx_in = xx**2 + yy**2 <= rad**2  # Cylindrical
    xx = xx[idx_in]
    yy = yy[idx_in]
    zz = np.zeros_like(xx)
    return xx, yy, zz


def obtain_cylindrical_stack_btwn_pts(start_pt, end_pt, rad=1., num_pts=5):
    '''Fill a cylinder between two points,radius in with mgrid style points'''
    gt = np.complex(0, num_pts)  # density  of the points
    xx, yy = np.mgrid[-rad:rad:gt, -rad:rad:gt]
    idx_in = xx**2 + yy**2 <= rad**2  # Cylindrical
    xx = xx[idx_in]
    yy = yy[idx_in]
    zz = np.zeros_like(xx)
    all_vecs = np.vstack((xx.flatten(), yy.flatten(), zz.flatten())).T
    all_vecs = all_vecs.astype(float) # shape = num of pts x 3
    rot_vecs = np.zeros_like(all_vecs)
    k = (0, 0, 1)  # z axis - normal to the plane xy
    orientation = tuple(jj - ii for ii, jj in zip(start_pt, end_pt))
    theta = obtain_angle_btwn_vecs(k, orientation)
    r1 = axisangle_to_q(np.cross(k, orientation), theta)
    for ii in range(xx.shape[0]):
        rot_vecs[ii] = qv_mult(r1, tuple((all_vecs[ii])))
    rot_xx, rot_yy, rot_zz = rot_vecs.T
    dist_btwn_gridpts = 2 * rad / num_pts
    pts_on_line = obtain_n_pts_btwn_pts(start_pt, end_pt, dist_btwn_gridpts)
    return (rot_xx, rot_yy, rot_zz), pts_on_line


def area_triangle_three_pts(p1, p2, p3):
    '''Area between 3 points'''
    v1 = p1 - p2
    v2 = p2 - p3
    area = np.linalg.norm(np.cross(v1, v2)) / 2.
    return area


def radius_btw_line_pts(pt1, pt2, pos_list):
    '''Distance mean+3sigma between line and a list of points'''
    dist = []
    dist_base = np.linalg.norm(pt1 - pt2)
    for pos in pos_list:
        area = area_triangle_three_pts(pt1, pt2, np.array(pos))
        dist.append(2 * area / dist_base)
    dists = np.array(dist)
    return dists.mean() + (3 * np.std(dists))


def closest_pt_vector(vec, vec_pt, pt):
    '''vec passes thru vec_pt and the closest pt on this vec'''
    fact = np.dot((pt - vec_pt), vec)
    return vec_pt + (fact * vec)


def fit_a_line(pos_list):
    '''http://stackoverflow.com/a/2333251/603292'''
    pos_array = np.array(pos_list)
    pos_mean = pos_array.mean(axis=0)
    uu, dd, vv = np.linalg.svd(pos_array - pos_mean)
    all_dists = np.linalg.norm(pos_array, axis=1)
    mx_idx = np.where(all_dists == all_dists.max())[0][0]
    mn_idx = np.where(all_dists == all_dists.min())[0][0]
    start_pt = closest_pt_vector(vv[0], pos_mean, pos_array[mn_idx])
    end_pt = closest_pt_vector(vv[0], pos_mean, pos_array[mx_idx])
    linepts = np.vstack((start_pt, end_pt))
    return linepts


pos_list = []
with open('traub_post_transform.csv', 'rb') as csvfile:
    next(csvfile, None)
    spamreader = csv.reader(csvfile, delimiter=',')
    for row in spamreader:
        pos_list.append([float(row[1]), float(row[2]), float(row[3])])
linepts = fit_a_line(pos_list)
radius = radius_btw_line_pts(linepts[0], linepts[1], pos_list)

start_pt = linepts[0]
end_pt = linepts[1]
rot_vecs, pts_on_line = obtain_cylindrical_stack_btwn_pts(start_pt, end_pt, rad=radius)
rot_xx, rot_yy, rot_zz = rot_vecs

num_pts = np.complex(0, 5)
xx, yy = np.mgrid[-radius:radius:num_pts, -radius:radius:num_pts]
idx_in = xx**2 + yy**2 <= radius**2
xx = xx[idx_in]
yy = yy[idx_in]
zz = np.zeros_like(xx)

plt3d = plt.figure().gca(projection='3d')
tri = mtri.Triangulation(rot_xx, rot_yy)
for ii in pts_on_line:
    plt3d.plot_trisurf(rot_xx + ii[0],
                       rot_yy + ii[1],
                       rot_zz + ii[2], triangles=tri.triangles, color='blue')
    # plt3d.scatter3D(rot_xx + ii[0],
    #                 rot_yy + ii[1],
    #                 rot_zz + ii[2], c='blue')

plt3d.plot(xs=(start_pt[0], end_pt[0]),
           ys=(start_pt[1], end_pt[1]),
           zs=(start_pt[2], end_pt[2]), color='green')
#plt3d.scatter3D(*np.array(pos_list).T, c='red')

plt.show()

# Obtain the planes passing through a point and origin and a perp plan to that
# Used in obtaining the points on planes of interest when computing CSD
# Uses Ax + By + Cz = D equation of the plane. Where (A, B, C) is the normal of the plane

# grid_x_plane, grid_y_plane = np.mgrid[-0.9:0.9:42j, -0.9:0.9:42j]
# pl_A, pl_B, pl_C, pl_D = get_plane(src_pos, snk_pos)
# print is_on_plane(pl_A, pl_B, pl_C, pl_D, snk_pos)
# grid_z_plane = get_pts_plane(pl_A, pl_B, pl_C, pl_D, grid_x_plane, grid_y_plane)
# idx_out = grid_x_plane**2 + grid_y_plane**2 + grid_z_plane**2 < 0.81 # Ensures that the points are within the sphere.
# grid_x_plane, grid_y_plane, grid_z_plane = grid_x_plane[idx_out], grid_y_plane[idx_out], grid_z_plane[idx_out]

# pl_P, pl_Q, pl_R, pl_S = get_perp_plane(pl_A, pl_B, pl_C, dpole.src_pos, dpole.snk_pos)
# grid_x_perp, grid_y_perp = np.mgrid[-0.9:0.9:42j, -0.9:0.9:42j]
# grid_z_perp = get_pts_plane(pl_P, pl_Q, pl_R, pl_S, grid_x_perp, grid_y_perp)
# idx_out = grid_x_perp**2 + grid_y_perp**2 + grid_z_perp**2 < 0.81
# grid_x_perp, grid_y_perp, grid_z_perp = grid_x_perp[idx_out], grid_y_perp[idx_out], grid_z_perp[idx_out]

# grid_x = np.hstack((grid_x_plane, grid_x_perp))
# grid_y = np.hstack((grid_y_plane, grid_y_perp))
# grid_z = np.hstack((grid_z_plane, grid_z_perp))
