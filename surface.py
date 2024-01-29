"""
Copyright (c) 2020 by Martin Kern

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
USA or visit their web site at www.gnu.org .
"""

#import pdb_reader
#import networkx as nx
import numpy as np
#import time
#from matplotlib import pyplot
#from scipy.spatial import distance_matrix, Delaunay
#from scipy.interpolate import BSpline
#from mpl_toolkits.mplot3d import Axes3D
from sklearn.preprocessing import PolynomialFeatures
from sklearn.linear_model import LinearRegression
#from math import sqrt
#import sys
#from geomdl import fitting, operations
#from geomdl.visualization import VisMPL as vis
#import nurbspy as nrb


#def sample_average(pdb_file, key_atoms, sample_number, radius):
#    # fig = pyplot.figure()
#    # ax = Axes3D(fig)
#    np.set_printoptions(threshold=sys.maxsize)

#    box = pdb_reader.box(pdb_file)
#    positions = pdb_reader.atom_positions(pdb_file)
#    np_pos = np.zeros([len(positions), 3])
#    sample_dist_x = box[0] / sample_number
#    sample_dist_y = box[1] / sample_number
#    x_samples = sample_number + 1  # int(box[0] / sample_number)
#    y_samples = sample_number + 1  # int(box[1] / sample_number)
#    i = 0
#    for pos in positions:
#        np_pos[i, 0] = pos[0]
#        np_pos[i, 1] = pos[1]
#        np_pos[i, 2] = pos[2]
#        i += 1

#    np_grid = np.zeros([x_samples * y_samples, 2])
#    for j in range(x_samples):
#        for k in range(y_samples):
#            np_grid[j * y_samples + k][0] = sample_dist_x * j
#            np_grid[j * y_samples + k][1] = sample_dist_y * k

#    dm = distance_matrix(np_pos[:, :2], np_grid)
#    bool_matrix = dm <= radius
#    normalizer = np.sum(bool_matrix, 0)
#    normalizer[normalizer == 0] = 1
#    grid_z = np.dot(np_pos[:, 2], bool_matrix) / normalizer  # np.sum(bool_matrix, 0)

#    surf_data = np.c_[np_grid, grid_z]

#    key_positions = pdb_reader.atom_positions(key_atoms)
#    np_keys = np.zeros([len(key_positions), 3])
#    i = 0
#    for pos in key_positions:
#        np_keys[i, 0] = pos[0]
#        np_keys[i, 1] = pos[1]
#        np_keys[i, 2] = pos[2]
#        i += 1

#    points_tuple = tuple(map(tuple, surf_data))

#    # print(x_samples, y_samples)
#    surf = fitting.interpolate_surface(points_tuple, x_samples, y_samples, 3, 3)
#    print(surf.knotvector_u, surf.knotvector_v)
#    n_dim = 3
#    n = surf.ctrlpts_size_u
#    m = surf.ctrlpts_size_v
#    p = np.zeros((n_dim, n, m))
#    ctl_points = np.array(surf.ctrlpts)

#    for i in range(n):
#        for j in range(m):
#            p[:, j, i] = ctl_points[i * m + j]

#    w = np.array(surf.weights)
#    nrb_surf = nrb.NurbsSurface(control_points=p,
#                                u_degree=3, v_degree=3,
#                                u_knots=np.array(surf.knotvector_u), v_knots=np.array(surf.knotvector_v))

#    # print(nrb_surf.get_curvature(0.5, 0.5))

#    # nrb_surf.plot_curvature()
#    # nrb_surf.plot(control_points=True)
#    # pyplot.show()

#    # d = surf.derivatives(0.5, 0.5, 2)
#    # for i in range(3):
#    #     for j in range(3):
#    #         print(i, j, d[i][j])

#    pic_res = 62
#    res_size = 1 / pic_res

#    x_coord = res_size
#    y_coord = res_size
#    curvature_pic = np.zeros([pic_res - 2, pic_res - 2])
#    for i in range(pic_res-2):
#        for j in range(pic_res-2):
#            y_coord += res_size
#            if y_coord > 1:
#                y_coord = 1
#            d = surf.derivatives(x_coord, y_coord, 2)
#            s_u = np.array(d[1][0])
#            s_v = np.array(d[0][1])
#            s_uu = np.array(d[2][0])
#            s_uv = np.array(d[1][1])
#            s_vv = np.array(d[0][2])
#            normal = np.array(operations.normal(surf, [x_coord, y_coord])[1])

#            e_1 = np.sum(s_u * s_u, axis=0)
#            f_1 = np.sum(s_u * s_v, axis=0)
#            g_1 = np.sum(s_v * s_v, axis=0)

#            e_2 = np.sum(s_uu * normal, axis=0)
#            f_2 = np.sum(s_uv * normal, axis=0)
#            g_2 = np.sum(s_vv * normal, axis=0)

#            # v = np.linalg.norm(np.cross(df1, df2)) / (np.linalg.norm(df1)**3)
#            # print(v)
#            # if v > 1:
#            #     v = 0

#            s_u = d[1][0][2] / d[1][0][0]
#            s_uu = d[2][0][2] / d[2][0][0]
#            c_u = s_uu / ((1 + s_u**2)**1.5)
#            s_v = d[0][1][2] / d[0][1][1]
#            s_vv = d[0][2][2] / d[0][2][1]
#            c_v = s_vv / ((1 + s_v**2)**1.5)
#            c_m = 0.5 * (c_u + c_v)
#            c_g = c_u * c_v
#            curvature_pic[i, j] = nrb_surf.get_curvature(x_coord, y_coord)[0]

#            # curvature_pic[j, i] = - 0.5 * (g_2 * e_1 - 2 * f_2 * f_1 + e_2 * g_1) / (e_1 * g_1 - f_1**2)
#            # curvature_pic[i, j] = (e_2 * g_2 - f_2**2) / (e_1 * g_1 - f_1**2)
#        y_coord = 0
#        x_coord += res_size
#        if x_coord > 1:
#            x_coord = 1

#    # # Sample points
#    # fig = pyplot.figure()
#    # ax = Axes3D(fig)
#    # ax.scatter(np_keys[:, 0], np_keys[:, 1], np_keys[:, 2])
#    # ax.scatter(surf_data[:, 0], surf_data[:, 1], surf_data[:, 2])
#    # pyplot.show()

#    # # curvature
#    pyplot.imshow(curvature_pic)
#    pyplot.colorbar()
#    pyplot.show()

#    # # surface
#    surf.delta = 0.05
#    vis.VisConfig(ctrlpts=False, bbox=False)
#    surf.vis = vis.VisSurface(vis.VisConfig(ctrlpts=False))
#    surf.render()


#def regression_test():
#    p = np.array([[0.1, -0.05, 0.9], [1.3, 0.054, 2.3], [-0.9, 0.2, 1.94], [0.01, 2.4, 3.6], [-0.03, -1.76, -3.2],
#                  [4.2, 0.54, 6.32], [-3.76, 0.043, 5.76], [2.54, 3.24, 5.34]])
#    z = p[:, 2]
#    x = np.array([np.ones([8]), p[:, 0], p[:, 1], p[:, 0]**2, p[:, 0] * p[:, 1], p[:, 1]**2])
#    b = np.matmul(np.matmul(np.linalg.inv(np.matmul(x, np.transpose(x))), x), np.transpose(z))
#    print(b)

#    poly = PolynomialFeatures(degree=2)
#    x_t = poly.fit_transform(p[:, :2])
#    clf = LinearRegression()
#    clf.fit(x_t, p[:, 2])
#    coefficients = clf.coef_
#    coefficients[0] = clf.intercept_
#    print(coefficients)


def curve_fit(positions, box, order, pbc):
    #start = time.process_time()
    atom_data = positions # see pdb_reader
    data_np = np.zeros([len(atom_data), 3])
    box = box # see pdb_reader
    i = 0
    for atom in atom_data:
        data_np[i, 0] = atom[0]
        data_np[i, 1] = atom[1]
        data_np[i, 2] = atom[2]
        i += 1

    pbx = box[0] * pbc / 100
    pby = box[1] * pbc / 100

    pbc_1 = data_np[data_np[:, 0] < pbx]
    pbc_1[:, 0] = pbc_1[:, 0] + box[0]
    pbc_2 = data_np[data_np[:, 0] > box[0] - pbx]
    pbc_2[:, 0] = pbc_2[:, 0] - box[0]
    pbc_3 = data_np[data_np[:, 1] < pby]
    pbc_3[:, 1] = pbc_3[:, 1] + box[1]
    pbc_4 = data_np[data_np[:, 1] > box[1] - pby]
    pbc_4[:, 1] = pbc_4[:, 1] - box[1]
    pbc_13 = data_np[(data_np[:, 0] < pbx) & (data_np[:, 1] < pby)]
    pbc_13[:, 0] = pbc_13[:, 0] + box[0]
    pbc_13[:, 1] = pbc_13[:, 1] + box[1]
    pbc_14 = data_np[(data_np[:, 0] < pbx) & (data_np[:, 1] > box[1] - pby)]
    pbc_14[:, 0] = pbc_14[:, 0] + box[0]
    pbc_14[:, 1] = pbc_14[:, 1] - box[1]
    pbc_23 = data_np[(data_np[:, 0] > box[0] - pbx) & (data_np[:, 1] < pby)]
    pbc_23[:, 0] = pbc_23[:, 0] - box[0]
    pbc_23[:, 1] = pbc_23[:, 1] + box[1]
    pbc_24 = data_np[(data_np[:, 0] > box[0] - pbx) & (data_np[:, 1] > box[1] - pby)]
    pbc_24[:, 0] = pbc_24[:, 0] - box[0]
    pbc_24[:, 1] = pbc_24[:, 1] - box[1]
    pbc_data = np.concatenate((data_np, pbc_1, pbc_2, pbc_3, pbc_4, pbc_13, pbc_14, pbc_23, pbc_24))
    # pbc_data[:, 0] = pbc_data[:, 0] - 100

    x, y = np.meshgrid(np.arange(0, box[0], 2), np.arange(0, box[1], 2))
    xx = x.flatten()
    yy = y.flatten()

    poly = PolynomialFeatures(degree=order)
    x_t = poly.fit_transform(pbc_data[:, :2])
    clf = LinearRegression()
    clf.fit(x_t, pbc_data[:, 2])

    poly_matrix = np.ones(xx.shape)
    for n in range(order):
        for m in range(n + 2):
            poly_matrix = np.c_[poly_matrix, xx**(n + 1 - m) * yy**m]

    order2 = 7
    poly_matrix2 = np.ones(xx.shape)
    for n in range(order2):
        for m in range(n + 2):
            poly_matrix2 = np.c_[poly_matrix2, xx ** (n + 1 - m) * yy ** m]

    coefficients = clf.coef_
    coefficients[0] = clf.intercept_
    # print(coefficients)
    z = np.dot(poly_matrix, coefficients).reshape(x.shape)

    # c, r, rank, s = numpy_regression(pbc_data, order2)
    # print(c)
    # z2 = np.dot(poly_matrix2, c).reshape(x.shape)

    # min_pen = 10000

    ##pd_x = np.ones(xx.shape)
    ##pd_y = np.ones(xx.shape)
    ##pd_xx = np.ones(xx.shape)
    ##pd_yy = np.ones(xx.shape)
    ##pd_xy = np.ones(xx.shape)
    ##for n in range(order - 1):
    ##    for m in range(n + 2):
    ##        pd_x = np.c_[pd_x, xx**(n + 1 - m) * yy**m]
    ##        pd_y = np.c_[pd_y, xx**(n + 1 - m) * yy**m]
    ##for n in range(order - 2):
    ##    for m in range(n + 2):
    ##        pd_xx = np.c_[pd_xx, xx**(n + 1 - m) * yy**m]
    ##        pd_yy = np.c_[pd_yy, xx**(n + 1 - m) * yy**m]
    ##        pd_xy = np.c_[pd_xy, xx**(n + 1 - m) * yy**m]
    ##z_x = np.dot(pd_x, partial_derivative_x(coefficients, order)).reshape(x.shape)
    ##z_y = np.dot(pd_y, partial_derivative_y(coefficients, order)).reshape(x.shape)
    ##z_xx = np.dot(pd_xx, partial_derivative_xx(coefficients, order)).reshape(x.shape)
    ##z_yy = np.dot(pd_yy, partial_derivative_yy(coefficients, order)).reshape(x.shape)
    ##z_xy = np.dot(pd_xy, partial_derivative_xy(coefficients, order)).reshape(x.shape)

    ##curvature = ((1 + z_y**2) * z_xx - 2 * z_x * z_y * z_xy + (1 + z_x**2) * z_yy)\
    ##    / (2 * np.sqrt(1 + z_x**2 + z_y**2)**3)

    ##curvature_x = z_xx / np.sqrt((1 + z_x**2)**3)
    ##curvature_y = z_yy / np.sqrt((1 + z_y**2)**3)
    ##mean_curvature = (curvature_x + curvature_y) / 2

    ##curvature = (z_xx * (1 + z_y**2) - 2 * z_x * z_y * z_xy + z_yy * (1 + z_x**2)) / (2 * np.sqrt(1 + z_x**2 + z_y**2))

    #end = time.process_time()
    #print('time: ', end-start)

    # fig, ax = pyplot.subplots()
    # pyplot.imshow(-curvature)
    # pyplot.colorbar()

    ##fig = pyplot.figure()
    ##ax = Axes3D(fig)
    ##ax.set_zlim(55, 118)
    ##ax.plot_surface(x, y, z, color='purple', antialiased=True, alpha=0.2)
    ##ax.plot_surface(x, y, z)
    # ax.plot_surface(x, y, z_y)
    # ax.plot_surface(x, y, z2)
    #
    # key_pos = pdb_reader.atom_positions(key_atoms)
    # x_key = []
    # y_key = []
    # z_key = []
    # for key in key_pos:
    #     x_key.append(key[0])
    #     y_key.append(key[1])
    #     z_key.append(key[2])
    #
    # ax.scatter(x_key, y_key, z_key)
    # print(min_pen)
    ##pyplot.show()

    return ([x,y,z])


#def numpy_regression(data, degree: int):
#    x = data[:, 0].flatten()
#    y = data[:, 1].flatten()
#    z = data[:, 2].flatten()
#    poly_matrix = np.ones(x.shape)
#    for n in range(degree):
#        for m in range(n + 2):
#            poly_matrix = np.c_[poly_matrix, x ** (n + 1 - m) * y ** m]
#    return np.linalg.lstsq(poly_matrix, z, rcond=None)


#def partial_derivative_x(coefficients, degree):
#    derivative = []
#    i = 0
#    for n in range(degree + 1):
#        for m in range(n + 1):
#            if (n - m) != 0:
#                derivative.append(coefficients[i] * (n - m))
#            i += 1
#    return np.array(derivative)


#def partial_derivative_y(coefficients, degree):
#    derivative = []
#    i = 0
#    for n in range(degree + 1):
#        for m in range(n + 1):
#            if m != 0:
#                derivative.append(coefficients[i] * m)
#            i += 1
#    return np.array(derivative)


#def partial_derivative_xx(coefficients, degree):
#    d1 = partial_derivative_x(coefficients, degree)
#    return partial_derivative_x(d1, degree - 1)


#def partial_derivative_yy(coefficients, degree):
#    d1 = partial_derivative_y(coefficients, degree)
#    return partial_derivative_y(d1, degree - 1)


#def partial_derivative_xy(coefficients, degree):
#    d1 = partial_derivative_x(coefficients, degree)
#    return partial_derivative_y(d1, degree - 1)


#def graph_surface():
#    y = np.arange(0, 10, 0.05)
#    sin = np.sin(y)
#    np.random.seed(0)
#    x = np.sin(y) + 0.5 - np.random.random(y.shape)
#    m = np.column_stack((x, y))
#    del_tri = Delaunay(m)
#    tri_points = del_tri.points
#    indices = del_tri.vertex_neighbor_vertices[1]
#    indptr = del_tri.vertex_neighbor_vertices[0]
#    spanning_tree = nx.Graph()
#    spanning_tree.add_node(indices[0])  # Insert first node into spanning tree
#    visited = [indices[0]]
#    print(indices[0])

#    for i in range(len(m) - 1):
#        closest_node1 = -1
#        closest_node2 = -1
#        shortest_dist = 100000000000000000000000
#        for node in spanning_tree.nodes:
#            neighbor_nodes = list(indices[indptr[node]:indptr[node + 1]])
#            neighbor_nodes = [n for n in neighbor_nodes if n not in visited]
#            if neighbor_nodes:
#                pos1 = tri_points[node]
#                for index in neighbor_nodes:
#                    pos2 = del_tri.points[index]
#                    new_dist = sqrt((pos1[0] - pos2[0]) ** 2 + (pos1[1] - pos2[1]) ** 2)
#                    if new_dist < shortest_dist:
#                        shortest_dist = new_dist
#                        closest_node1 = index
#                        closest_node2 = node
#        visited.append(closest_node1)
#        spanning_tree.add_node(closest_node1)
#        spanning_tree.add_edge(closest_node2, closest_node1)
#    nx.draw(spanning_tree, m)
#    # pyplot.triplot(m[:, 0], m[:, 1], del_tri.simplices)
#    # pyplot.plot(m[:, 0], m[:, 1], 'o')
#    pyplot.show()


#def benchmark(size, x_box, y_box, z_box, iterations, order, sample_distance, radius):
#    data_list = []
#    for i in range(iterations):
#        data = np.random.random([size, 3])
#        data[:, 0] *= x_box
#        data[:, 1] *= y_box
#        data[:, 2] *= z_box
#        data_list.append(data)

#    print('Surface benchmark.\nData size:  ' + str(size) + '\nIterations: ' + str(iterations) + '\n')
#    # Parameterized regression
#    p_start = time.time()
#    for data_set in data_list:
#        x, y = np.meshgrid(np.arange(0, x_box, 2), np.arange(0, y_box, 5))
#        xx = x.flatten()
#        yy = y.flatten()
#        poly = PolynomialFeatures(degree=order)
#        x_t = poly.fit_transform(data_set[:, :2])
#        clf = LinearRegression()
#        clf.fit(x_t, data_set[:, 2])
#        poly_matrix = np.ones(xx.shape)
#        for n in range(order):
#            for m in range(n + 2):
#                poly_matrix = np.c_[poly_matrix, xx ** (n + 1 - m) * yy ** m]

#        coefficients = clf.coef_
#        coefficients[0] = clf.intercept_
#        z = np.dot(poly_matrix, coefficients).reshape(x.shape)
#    p_end = time.time()
#    print('Parameterized - Time:   ' + str(p_end - p_start))

#    # Unparameterized regression
#    u_start = time.time()
#    for data_set in data_list:
#        x_samples = int(x_box / sample_distance)
#        y_samples = int(y_box / sample_distance)
#        np_grid = np.zeros([x_samples * y_samples, 2])
#        for j in range(x_samples):
#            for k in range(y_samples):
#                np_grid[j * x_samples + k][0] = sample_distance * j
#                np_grid[j * x_samples + k][1] = sample_distance * k

#        dm = distance_matrix(data_set[:, :2], np_grid)
#        bool_matrix = dm <= radius
#        grid_z = np.dot(data_set[:, 2], bool_matrix) / np.sum(bool_matrix, 0)
#    u_end = time.time()
#    print('Unparameterized - Time: ' + str(u_end - u_start))


#def test_spline():
#    x = np.arange(0.0, 3.0, 0.1)
#    tau = np.array([0, 0, 0, 1, 2, 3, 4, 4, 4])
#    control_points = np.array([[0, 0], [1, 0], [2, 1], [1, 1], [0.5, 1], [0, 0]])
#    order = 3
#    spline = BSpline(tau, control_points, order)
#    y = spline(x)
#    pyplot.plot(y[:, 0], y[:, 1])
#    pyplot.show()


#def poly_val(coefficients, degree, x, y):
#    """
#    Calculate the value of a polynomial at point x, y.
#    :param coefficients: Coefficients of the polynomial.
#    :param degree: Degree of the polynomial.
#    :param x: x-coordinate
#    :param y: y-coordinate
#    :return: Value at point x, y.
#    """
#    i = 0
#    value = 0
#    for n in range(degree + 1):
#        for m in range(n + 1):
#            value += coefficients[i] * x**(n - m) * y**m
#            i += 1
#    return value
