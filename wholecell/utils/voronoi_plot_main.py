"""
Plot voronoi diagram

Algorithm Reference: A. Nocaj and U. Brandes, "Computing voronoi Treemaps:
Faster, Simpler, and Resolution-independent, Computer Graphics Forum,
vol. 31, no. 31, pp. 855-864, 2012. Available: 10.1111/j.1467-8659.2012.03078.x.

This program generate a voronoi diagram, in which the canvas of a plot is
divided in to n polygons. Each data point (or "Points" in this program)
corresponds to one polygon, and the centroid of each polygon is called a "site".
The area of each polygon is proportional to the value each point wants to
represent. For example, it can be the mass compositions of a cell, the relative
abundance of  different proteins, or the energy cost of different cellular
processes.

Dividing a canvas (in this program we use np.array([[0,0],[4,0],[4,4],[0,4]]))
into multiple polygons is a complicated computational geometry problem. The main
working horse of this task is the function "_voronoi_treemap".
In "_voronoi_treemap", we did the following things:
(1) initialize random sites within the canvas (using "random_points_in_canvas")
and assign an initial weight to each site. The initial weights of each site are
set equal so that it will be easier to compute a draft.
(2) compute the draft of our voronoi diagram, which is called "Power Diagram"
following the notation of our algorithm reference.
(using "_compute_power_diagram")
(3) adjust the positions of each site and the weights of each site until the
program reaches the maximum iterations (i_max = 75) or until error < err_thres.
(4) return the final error of the whole voronoi diagram.

In voronoi diagram, certain amount of error in area representation is expected.
According to the analysis in our algorithm reference, an error between 5~23% is
expected. However, an error up to 20% can severely distort the diagram. We
therefore set a stricter limit (15%) on our error in the function
"_voronoi_main_function". If the newly computed voronoi diagram has a lower
error compared to the previous one, it will continuously increase the number of
iterations until the error falls below 15%. However, if the error of a newly
computed diagram is greater than the previous one, the program will plot the
diagram all over again. Generally, a final error <=10% will give you a nice
representation.

Using this algorithm, the program is capable of dividing 32 polygons in 10
seconds, and the program scales with O(n*log(n)). As a result, it would be
preferable to introduce layered structure if the number of polygons are more
than 50.

If you want to implement the functions that are created here and make a new
voronoi plot with your own data, please read the following instructions:
(1)Place this in the heading, or other ways that can import every function in
this file.
from wholecell.utils.voronoiPlotMain import VoronoiMaster

(2)Copy the following code and modify it appropriately:
vm = VoronoiMaster()
vm.plot("YOUR DICTIONARY", title = "YOUR TITLE")
exportFigure(plt, plotOutDir, plotOutFileName, metadata)
plt.close("all")

"""
from __future__ import absolute_import, division, print_function

from typing import cast, Iterable, List, Optional, Union

import numpy as np
import math as math
from scipy.spatial import distance_matrix
from scipy.spatial import ConvexHull
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
import matplotlib.pyplot as plt
from six.moves import range, zip


COLORS_256 = [
    [55, 126, 184],
    [255, 127, 0],
    [152, 78, 163],
    [77, 175, 74],
    [255, 255, 51],
    [166, 86, 40],
    [247, 129, 191],
    [228, 26, 28],
    ] # From colorbrewer2.org, qualitative 8-class set 1
COLORS = [[colorValue/255. for colorValue in color] for color in COLORS_256]


def angle(v1, v2):
    def dot_product(v1, v2):
        return sum((a * b) for a, b in zip(v1, v2))

    cosine = dot_product(v1, v2) / (
            math.hypot(v1[0], v1[1]) * math.hypot(v2[0], v2[1]))
    if (cosine > 1) or (cosine < -1):
        return 100
    else:
        return math.acos(cosine)

def is_on_segment(p, edge):
    """
    Determine whether a point is on a segment by checking if Ax+By-C == 0
    and falls between the two corners which define the edge.
    """
    [[x1, y1], [x2, y2]] = edge
    [x, y] = p
    # convert to ax + by = c
    a = (y2 - y1)
    b = - (x2 - x1)
    c = x1*(y2 - y1) - y1*(x2 - x1)
    if (a**2 + b**2) == 0:
        result = (x == x1) and (y == y1)

    else:
        test = (a*x + b*y - c)
        if int(test*(10**9))/(10.**9) == 0:
            x = int(x*(10**9) + 0.5)/(10.**9)
            x1 = int(x1*(10**9) + 0.5)/(10.**9)
            x2 = int(x2*(10**9) + 0.5)/(10.**9)
            y = int(y*(10**9) + 0.5)/(10.**9)
            y1 = int(y1*(10**9) + 0.5)/(10.**9)
            y2 = int(y2*(10**9) + 0.5)/(10.**9)
            result = ((x >= min(x1, x2)) and (x <= max(x1, x2))
                and (y >= min(y1, y2)) and (y <= max(y1, y2)))
        else:
            result = False

    return result

class PolygonClass(object):
    def __init__(self, xy):
        '''
        A polygon object stores the coordinates of the corners, the edges, and
        the area of the polygon.
        '''
        def _reorder_points(corners):
            """
            This function reorders the corners of a polygon in a
            counterclockwise manner. The input should be a numpy array of nx2.
            """
            ordered_points = corners
            com = ordered_points.mean(axis = 0) # find center of mass
            ordered_points = ordered_points[
                np.argsort(np.arctan2((ordered_points - com)[:, 1],
                (ordered_points - com)[:, 0]))]
            return ordered_points

        self.xy = _reorder_points(xy)
        self.n_corners = len(xy)
        self.edges = []
        for i in range(self.n_corners):
            self.edges.append(self.xy[[i-1, i], :])
        self.area = self._polygon_area(self.xy)

    def _polygon_area(self, corners):
        """
        Calculate polygon area using shoelace formula.
        Please make sure that the corners are reordered before calling
        _polygon_area function!
        """
        n = len(corners)
        area = 0.0
        for i in range(n):
            j = (i + 1) % n
            area += corners[i][0] * corners[j][1]
            area -= corners[j][0] * corners[i][1]
        area = abs(area) / 2.0
        return area

    def random_points_in_canvas(self, n_points):
        """
        This function assigns n random points within the canvas. The algorithm
        is as followed:
        (1) divide the canvas into multiple triangels
        (2) weight each triangle according to the area
        (3) randomly placed points within the triangles
        Please make sure that the points have been reordered properly!
        """
        vec_all = self.xy[1:, :] - self.xy[0, :]
        n_triangle = self.n_corners - 2
        area_triangle = np.zeros(n_triangle)

        for i in range(n_triangle):
            area_triangle[i] = self._polygon_area(
                np.vstack([[0, 0], vec_all[i], vec_all[i+1]]))
        rand_scale = np.sum(np.tril(area_triangle), axis = 1)/sum(area_triangle)
        rand_num = np.hstack(
            [np.random.rand(n_points, 3), np.zeros([n_points, 1])])

        sites = np.zeros([n_points, 2])
        for i in range(n_triangle):
            mask = (rand_num[:, 0] <= rand_scale[i]) * (rand_num[:, 3] == 0)
            vec1_tile_xy = np.tile(vec_all[i, :], (sum(mask), 1))
            vec2_tile_xy = np.tile(vec_all[i+1, :], (sum(mask), 1))
            rand_num_masked = rand_num[mask, 1].reshape((-1, 1))
            rand_len_masked = np.sqrt(rand_num[mask, 2].reshape((-1, 1)))
            sites[mask, :] = (rand_num_masked*vec1_tile_xy*rand_len_masked 
                + (1 - rand_num_masked)*vec2_tile_xy*rand_len_masked)
            rand_num[mask, 3] = 1

        sites += np.tile(self.xy[0, :], (n_points, 1))
        return sites

    def point_is_within_canvas(self, p):
        """
        Check if a point is within canvas by computing the sum of the angle
        formed by adjacent corners and the point. If a point is within canvas,
        the sum of the angle should be 2*pi.
        """
        def _point_is_on_corner_of_canvas(p):
            """
            Check if a point is one of the corner of the canvas.
            """
            a = (self.xy == p)
            result = len(np.nonzero(a[:, 0] * a[:, 1])[0]) == 1
            return result

        if _point_is_on_corner_of_canvas(p):
            result = True

        else:
            sum_angle = 0
            for edge_canvas in self.edges:
                a = edge_canvas - p
                sum_angle += angle(a[0, :], a[1, :])

            if not math.isnan(sum_angle):
                result = int((sum_angle - 2*math.pi)*(10**9))/(10.**9) == 0
            else:
                result = False

        return result

    def find_min_dist_to_border(self, site):
        """
        find the minimum distance of the site to the borders of its polygon.
        """
        def _find_min_dist_to_edge(site, edge):
            """
            find the minimum distance of the site to one specific edge of a
            polygon.
            """
            [[x1, y1],[x2, y2]] = edge
            [xs, ys] = site
            # convert to ax + by + c = 0
            a = (y2 - y1)
            b = - (x2 - x1)
            c = y1*(x2 - x1) - x1*(y2 - y1)
            dist_perp = abs(a*xs + b*ys + c)/math.sqrt(a**2 + b**2)
            xp = (b*(b*xs - a*ys) - a*c)/(a**2 + b**2)
            yp = (a*(- b*xs + a*ys) - b*c)/(a**2 + b**2)
            if is_on_segment([xp, yp], edge):
                min_dist = dist_perp

            else:
                min_dist = min(math.sqrt((x1 - xs)**2 + (y1 - ys)**2), 
                    math.sqrt((x2 - xs)**2 + (y2 - ys)**2))

            return min_dist

        min_dist_list = []
        for edge in self.edges:
            min_dist_list.append(_find_min_dist_to_edge(site, edge))
        distance_border = min(min_dist_list)
        return distance_border

class VoronoiClass(object):
    def __init__(self, polygons, sites, weights, values, canvas_obj):
        '''
        A voronoi object stores all the polygons of a voronoi diagram, the
        coordinates/weights/values of each site, and the canvas of the whole
        plot.
        '''
        self.polygons = polygons
        self.n_sites = len(polygons)
        self.sites = sites
        self.weights = weights
        self.values = values
        self.canvas_obj = canvas_obj

    def find_sites_causing_displacement(self):
        """
        If a site need to expand its area greatly in the next round, the other
        sites surrounding it should be displaced a little bit in order to
        facilitate the process of optimizing the voronoi diagram. This is part
        of the speed-up heuristic.
        """
        self.sites_area_ratio = [-np.inf] * self.n_sites
        self.move_me = [True] * self.n_sites
        for i in range(self.n_sites):
            if not self.polygons[i]: #rescue the sites
                self.sites[i, :] = self.canvas_obj.random_points_in_canvas(1)
                self.weights[i] = self.canvas_obj.area/self.n_sites

            else:
                area_current = self.polygons[i].area
                area_target = self.canvas_obj.area*self.values[i]/sum(self.values)
                self.sites_area_ratio[i] = area_target/area_current
                if area_current/self.canvas_obj.area < 0.05:
                    self.move_me[i] = False

        try:
            self.site_causing_disp = cast(int, np.argmax(self.sites_area_ratio))
            self.move_me[self.site_causing_disp] = False

        except ValueError:
            self.site_causing_disp = 0

        return

    def adapt_positions_weights(self):
        """
        Adapt the positions and weights of each site within an existing voronoi
        diagram. This algorithm also considers the point with largest
        Area_{target}/Area_{current} ratio and displace its neighbor in order to
        speed up the optimization process.
        """

        # 1. move the location of each site to the centroid of each polygon.
        # If there is no closed polygon for that site, randomly reassign a
        # location for that site in order to rescue it.
        for i in range(self.n_sites):
            if not self.polygons[i]: #rescue the sites
                self.sites[i, :] = self.canvas_obj.random_points_in_canvas(1)
                self.weights[i] = self.canvas_obj.area/self.n_sites

            else: # adapt the position to the centroid
                self.sites[i, :] = self.polygons[i].xy.mean(axis = 0)

        # 2. speed up heuristic: if a site is expected to expand its area
        # greatly, its neighbor should be displaced in advance in order to speed
        # up the optimization process.
        for i in range(self.n_sites):
            if self.polygons[i]:
                if self.move_me[i]: #speed-up heuristic#
                    dist_vec = (self.sites[i, :]
                                - self.sites[self.site_causing_disp, :])
                    dist_abs = math.sqrt(dist_vec[0]**2 + dist_vec[1]**2)
                    if dist_abs > 0:
                        ds_abs = (math.sqrt(
                            self.polygons[i].area/math.pi))**2/(5*dist_abs)
                        site_new = self.sites[i, :] + dist_vec*ds_abs/dist_abs
                        if self.polygons[i].point_is_within_canvas(site_new):
                            self.sites[i, :] = site_new

        # 3. if the square root of weights exceed the minimum distance from the
        # site to the current border of the polygon, the weight should be
        # decreased.
                distance_border = self.polygons[i].find_min_dist_to_border(
                    self.sites[i, :])
                self.weights[i] = (
                    np.nanmin([math.sqrt(self.weights[i]), distance_border]))**2

        return

    def adapt_weights(self):
        """
        Adapt the weight of each site. increase the weight if current area is
        different from the targeted area of each site. However, if the weight is
        too big such that the site is encroaching its neighbor, the weight
        should be decreased.
        """

        def _nearest_neighbor(self):
            """
            find the nearest neighbor of each site and return the DISTANCE
            (rather than the index) between each site and its nearest neighbor
            """
            nn_dist_cal = np.nanmin(
                distance_matrix(self.sites, self.sites)
                + np.diag([np.nan]*self.n_sites), axis = 0)
            return nn_dist_cal

        nn_dist = _nearest_neighbor(self)
        epsilon = self.canvas_obj.area/1000 # minimum weight of each site
        for i in range(self.n_sites):
            if not self.polygons[i]: # rescue the sites
                self.sites[i, :] = self.canvas_obj.random_points_in_canvas(1)
                self.weights[i] = self.canvas_obj.area/self.n_sites

            else:
                area_current = self.polygons[i].area
                area_target = self.canvas_obj.area*self.values[i]/sum(self.values)
                f_adapt = area_target/area_current
                sqrt_w_new = math.sqrt(self.weights[i])*f_adapt
                sqrt_w_max = nn_dist[i]
                self.weights[i] = (np.nanmin([sqrt_w_new, sqrt_w_max]))**2
                self.weights[i] = np.nanmax([self.weights[i], epsilon])

        return

    def compute_error(self):
        """
        compute the error of the voronoi plot. The error is defined as 
        sum(abs(Area_current - Area_target))/(2*Area_canvas)
        """
        error = 0
        total_value = sum(self.values)
        for i in range(self.n_sites):
            if not self.polygons[i]: # rescue the sites if the polygon is empty
                self.sites[i, :] = self.canvas_obj.random_points_in_canvas(1)
                self.weights[i] = self.canvas_obj.area/self.n_sites
                error += np.inf
            else:
                error += abs(self.polygons[i].area - 
                    self.canvas_obj.area*self.values[i]/total_value)
        error = error/(2*self.canvas_obj.area)
        return error

class LineClass(object):
    def __init__(self, x_col, y_col, weights):
        '''
        A line object store the information of a line in this format of
        "dx + ey = f", using formula:
        2(x1-x2)x + 2(y1-y2)y = (x1^2+y1^2-w1) - (x2^2+y2^2-w2)
        '''
        [[x1], [x2]] = x_col
        [[y1], [y2]] = y_col
        [w1, w2] = weights
        self.d = 2*(x1 - x2)
        self.e = 2*(y1 - y2)
        self.f = (x1**2 + y1**2 - w1) - (x2**2 + y2**2 - w2)

    def find_intersect_with_canvas(self, canvas_obj):
        """
        This function finds the intersection points of a bisector with a canvas.
        """
        def _line_intersect_with_edge(self, edge_canvas):
            """
            This function finds the intersection points of a line with an edge.
            """
            # [P1, P2] = edge_canvas
            [[x1, y1],[x2, y2]] = edge_canvas
            # convert to ax + by = c, dx + ey = f
            a = (y2 - y1)
            b = - (x2 - x1)
            c = x1*(y2 - y1) - y1*(x2 - x1)
            det = a*self.e - b*self.d
            det_x = c*self.e - self.f*b
            det_y = a*self.f - c*self.d 

            # if their slopes are the same but they don't intersect,
            # they are parallel
            if det == 0:
                if (det_x == 0) and (det_y == 0):
                    result = True
                    intersect_point = [np.inf]
                else:
                    result = False
                    intersect_point = []
            else:
                x_int = det_x/det
                y_int = det_y/det
                if is_on_segment([x_int, y_int], edge_canvas):
                    result = True
                    intersect_point = np.array([x_int, y_int])
                else:
                    result = False
                    intersect_point = []

            return result, intersect_point

        intersect_TF_list = []
        intersect_coordinates = []
        for edge_canvas in canvas_obj.edges:
            result, intersect_point = _line_intersect_with_edge(self, edge_canvas)
            intersect_TF_list.append(result)
            if result:
                if len(intersect_point) == 2:
                    if len(intersect_coordinates) == 0:
                        intersect_coordinates = intersect_point
                    else:
                        intersect_coordinates = np.vstack(
                            [intersect_coordinates, intersect_point])

        edge = []
        if len(intersect_coordinates) > 0:
            intersect_coordinates = intersect_coordinates.reshape((-1, 2))
            intersect_coordinates = (
                            intersect_coordinates*(10**9)).astype(int)/(10.**9)
            unique_intersect_coordinate = np.unique(
                intersect_coordinates, axis = 0)
            if len(unique_intersect_coordinate) == 2:
                edge = unique_intersect_coordinate

        return edge

class RayClass(object):
    def __init__(self, origin, tangent, main_sites, adjunct_site):
        '''
        A ray object stores the origin and the tangent vector of a ray, and
        the site it represents. It is store in this format:
        ~ origin + a*tangent, a>=0
        ~ a_r*x + b_r*y + c_r = 0
        ~ t_y*x - t_x*y + (t_x*y_int - t_y*x_int) = 0
        '''
        self.origin = origin
        self.tangent = tangent
        self.a_r = tangent[1]
        self.b_r = - tangent[0]
        self.c_r = (tangent[0]*origin[1] - tangent[1]*origin[0])
        self.main_sites = main_sites
        self.adjunct_site = adjunct_site
        self.index_ip = -10000

    def is_on_ray(self, p):
        """
        Determine whether a point is on a ray by checking if Ax+By-C == 0 and
        its relative position with respect to the origin is on the same
        direction as the tangent vector of the ray.
        """
        t_x = self.tangent[0]
        t_y = self.tangent[1]
        test = (self.a_r*p[0] + self.b_r*p[1] + self.c_r)
        if math.isnan(test):
            result = False

        else:
            if int(test*(10**9))/(10.**9) == 0:
                if t_x != 0:
                    result = int((p[0] - self.origin[0])/t_x*(10**9))/(10.**9) >= 0

                else:
                    result = int((p[1] - self.origin[1])/t_y*(10**9))/(10.**9) >= 0

            else:
                result = False

        return result

    def ray_intersect_with_edge(self, edge_canvas):
        """
        This function determines if a ray is intersected with an edge.
        """
        [p1, p2] = edge_canvas
        [[x1, y1], [x2, y2]] = edge_canvas
        # convert to ax + by = c, dx+ey = f
        a = (y2 - y1)
        b = - (x2 - x1)
        c = x1*(y2 - y1) - y1*(x2 - x1)
        det = a*self.b_r - b*self.a_r
        det_x = c*self.b_r - b*(-self.c_r)
        det_y = a*(-self.c_r) - c*self.a_r
        if det == 0:
            if (det_x != 0) or (det_y != 0):
                result = False
                intersect_point = []
            else: # (det == 0) and (det_x == 0) and (det_y == 0)
                p1_is_on_ray = self.is_on_ray(p1)
                p2_is_on_ray = self.is_on_ray(p2)
                if p1_is_on_ray and p2_is_on_ray:
                    result = True
                    intersect_point = [np.inf]
                elif p1_is_on_ray != p2_is_on_ray:
                    if (np.array_equal(p1, self.origin)
                        or np.array_equal(p2, self.origin)):
                        result = True
                        intersect_point = (p1*np.array_equal(p1, self.origin) 
                            + p2*np.array_equal(p2, self.origin))
                    else:
                        result = True
                        intersect_point = [np.inf]
                else:
                    result = False
                    intersect_point = []            
        else:
            x_int = det_x/det
            y_int = det_y/det
            if (is_on_segment([x_int, y_int], edge_canvas)
                    and self.is_on_ray(np.array([x_int, y_int]))):
                result = True
                intersect_point = np.array([x_int, y_int])
            else:
                result = False
                intersect_point = []

        return result, intersect_point

    def keep_nearest_point_on_ray(self, intersect_coordinates):
        """
        If a ray is intersect with multiple ipoints, 
        only keep the one that is closest to the origin of the ray.
        """
        n_candidates = len(intersect_coordinates)
        # find candidates with shortest distance from the origin
        dist = np.sum((intersect_coordinates
                       - np.tile(self.origin,(n_candidates, 1)))**2, axis = 1)
        dist = (dist*(10**9)).astype(int)/(10.**9)
        # avoid the intersection point that is the origin itself
        dist[dist == 0] = np.inf
        intersect_coordinates = intersect_coordinates[np.argmin(dist), :]
        return intersect_coordinates

    def find_ray_intersect_with_canvas(self, canvas_obj):
        """
        This function finc the point where a ray object is intersected with
        canvas by checking the intersection with every edge of the canvas.
        """
        intersect_TF_list = []
        intersect_coordinates = []
        for edge_canvas in canvas_obj.edges:
            result, intersect_point = self.ray_intersect_with_edge(edge_canvas)
            intersect_TF_list.append(result)
            if result:
                if len(intersect_point) == 2:
                    if len(intersect_coordinates) == 0:
                        intersect_coordinates = intersect_point
                    else:
                        intersect_coordinates = np.vstack(
                            [intersect_coordinates, intersect_point])

        if len(intersect_coordinates) == 0:
            intersect_coordinates = self.origin
        else:
            intersect_coordinates = np.vstack(
                [intersect_coordinates, self.origin])

        edge = []
        if len(intersect_coordinates) > 0:
            intersect_coordinates = intersect_coordinates.reshape((-1, 2))
            intersect_coordinates = (
                            intersect_coordinates*(10**9)).astype(int)/(10.**9)
            unique_intersect_coordinate = np.unique(
                intersect_coordinates, axis = 0)
            if len(unique_intersect_coordinate) == 2:
                edge = unique_intersect_coordinate

        self.edge = edge
        return

    def find_ray_intersect_with_canvas_and_ipoints(self, canvas_obj, ipoints,
                                                   simplices_prune):
        """
        This function finds the terminal of each ray object.
        """
        # 0. get the label of each ray object
        intersect_coordinates = []
        simplex = simplices_prune[self.index_ip]
        a = ((simplices_prune == simplex[0]) 
            + (simplices_prune == simplex[1]) 
            + (simplices_prune == simplex[2]))
        arr = np.nonzero((a.sum(axis = 1) == 2))[0]

        # 1. check if the ray is intersected with any nearby ipoints
        for i in arr:
            ipoint = ipoints[i, 0:2]
            if self.is_on_ray(ipoint):
                if len(intersect_coordinates) == 0:
                    intersect_coordinates = ipoint
                else:
                    intersect_coordinates = np.vstack(
                        [intersect_coordinates, ipoint])

        # 2. if the ray is not intersected with any ipoints, find its
        # intersection with canvas
        if len(intersect_coordinates) == 0:
            intersect_TF_list = []
            for edge_canvas in canvas_obj.edges:
                result, intersect_point = self.ray_intersect_with_edge(
                    edge_canvas)
                intersect_TF_list.append(result)
                if result:
                    if len(intersect_point) == 2 :
                        if len(intersect_coordinates) == 0:
                            intersect_coordinates = intersect_point
                        else:
                            intersect_coordinates = np.vstack(
                                [intersect_coordinates, intersect_point])

        # 3. only keep the nearest point to the origin of that ray
        if len(intersect_coordinates) > 0:
            intersect_coordinates = intersect_coordinates.reshape((-1, 2))
        if len(intersect_coordinates) > 1:
            intersect_coordinates = self.keep_nearest_point_on_ray(
                intersect_coordinates)

        # 4. append origin and create an edge
        if len(intersect_coordinates) == 0:
            intersect_coordinates = self.origin
        else:
            intersect_coordinates = np.vstack(
                [intersect_coordinates, self.origin])

        # 5. round and leave only the unique point
        intersect_coordinates = intersect_coordinates.reshape((-1, 2))
        edge = []
        if len(intersect_coordinates) > 0:
            intersect_coordinates = (
                            intersect_coordinates*(10**9)).astype(int)/(10.**9)
            unique_intersect_coordinate = np.unique(
                intersect_coordinates, axis = 0)
            if len(unique_intersect_coordinate) == 2:
                edge = unique_intersect_coordinate

        self.edge = edge
        return

class VoronoiMaster(object):
    def __init__(self, i_max = 75, err_thres = 1E-6):
        self.i_max = i_max
        self.err_thres = err_thres

    def plot(self, dic, side_length=(4, 4), custom_shape_vertices=None,
            font_size=8, title=None, ax_shape=None, chained=None,
            verbose=False):
        '''
        Master function of generating layered or non-layered voronoi plot from a
        dictionary.
        gross_error = \sum_{i} |Area(i) - Expected Area(i)|
        error_all = \sum_{i} |Area(i) - Expected Area(i)|/(2 * total area)
        The factor 2 in the error formula is to correct the repeated calculation
        in area error.
        Args:
            dic: Layered dictionary which contains the labels and the values you
                intend to represent in a voronoi diagram. This can be a single
                dictionary or multiple dictionaries in a nested list. The shape
                of the list should be the same as the layout of the subplots.
            side_length: The side length of whole voronoi diagram. The whole
                diagram is defaulted to be rectangular.
            custom_shape_vertices: If you want a custom shaped voronoi diagram,
                enter the vertices of the voronoi diagram in a nested np array
                like this: np.array([[0, 0], [4, 0], [4, 1], [1.5, 3], [0, 2]])
                The shape of the whole voronoi diagram will then become a
                pentagon with vertices [0, 0], [4, 0], [4, 1], [1.5, 3], [0, 2].
            font_size: The font size of the labeling on the voronoi diagram.
            title: The title of the plot. This can be a single title or a nested
                list of titles.
            ax_shape: the shape of the subplot, in (nrows, ncols)
            chained: whether the new subplot should be based on the previous
                subplot or not.
            verbose: If true, print the magnitude of the error in the areas
                represented by the plot

        Returns:
            error_all: the error in total area representation.
        '''
        voronoi_list_old = []
        if ax_shape is not None:
            error_all = []
            nrows, ncols = ax_shape
            fig, axes = plt.subplots(
                nrows=nrows, ncols=ncols,
                figsize = (ncols*side_length[0], nrows*side_length[1]+1),
                dpi=200)
            axes = axes.reshape(ax_shape)
            plt.tight_layout()
            for i in range(nrows):
                for j in range(ncols):
                    dic_current = dic[i][j]

                    if (i == 0 and j == 0) or (chained is None):
                        voronoi_list_new, polygon_value_list_new, label_site_list_new = self._compute_boundaries(
                            dic_current, side_length, custom_shape_vertices)
                    else:
                        voronoi_list_new, polygon_value_list_new, label_site_list_new = self._compute_boundaries(
                            dic_current, side_length, custom_shape_vertices,
                            voronoi_list_old = voronoi_list_old)

                    axes[i][j].set_aspect('equal')
                    axes[i][j].axis('off')
                    axes[i][j].title.set_text(title[i][j])
                    self._generate_plot(voronoi_list_new, axes[i][j])
                    total_value = sum(voronoi_list_new[0].values)
                    total_area = voronoi_list_new[0].canvas_obj.area
                    self._add_labels(
                        label_site_list_new, font_size, axes[i][j])
                    gross_error = self._compute_error(
                        polygon_value_list_new, total_value, total_area)
                    error_new = gross_error / (
                            2 * voronoi_list_new[0].canvas_obj.area)
                    if verbose:
                        print('The error in the area representation of the whole '
                              'voronoi diagram (%i, %i): %.4f'
                              % (i, j, error_new,))
                    error_all.append(error_new)
                    voronoi_list_old = voronoi_list_new
        else:
            voronoi_list, polygon_value_list, \
            label_site_list = self._compute_boundaries(
                dic, side_length, custom_shape_vertices)

            fig = plt.figure(figsize=(side_length[0], side_length[1] + 1),
                             dpi=200)
            ax = fig.add_axes([0, 0, 1, 1])
            ax.set_aspect('equal')
            plt.axis('off')
            self._generate_plot(voronoi_list, ax)

            total_value = sum(voronoi_list[0].values)
            total_area = voronoi_list[0].canvas_obj.area

            self._add_labels(label_site_list, font_size, ax)
            gross_error = self._compute_error(
                polygon_value_list, total_value, total_area)
            error_all = gross_error / (2 * voronoi_list[0].canvas_obj.area)
            if verbose:
                print('The error in the area representation of the whole'
                      'voronoi diagram: %.4f' % (error_all,))
            plt.title(title)

        return error_all

    def _generate_plot(self, voronoi_list, ax, counter = 0):
        """
        Generate the layered voronoi diagram.
        """
        for voronoi in voronoi_list:
            if isinstance(voronoi, list):
                ax, counter = self._generate_plot(
                    voronoi, ax, counter = counter + 1)

            else:
                canvas_obj = voronoi.canvas_obj

                # plot canvas
                for edge_canvas in canvas_obj.edges:
                    ax.plot(edge_canvas[:, 0], edge_canvas[:, 1],
                            'black', lw = 1.5, solid_capstyle = 'round',
                            zorder = 2)

                # plot polygons
                patches = []
                colors_all = []
                for polygon in voronoi.polygons:
                    polygon_plot_obj = Polygon(polygon.xy, True)
                    ax.plot(polygon.xy[:, 0], polygon.xy[:, 1],
                            color = 'black', alpha = 1, linewidth = 1,
                            solid_capstyle = 'round', zorder = 2)
                    patches.append(polygon_plot_obj)
                    colors = COLORS[counter] + np.random.rand(3) / 5
                    if max(colors) >= 1:
                        colors = colors / max(colors)
                    colors_all.append(colors)
                p = PatchCollection(patches, facecolors = colors_all, alpha = 1)
                ax.add_collection(p)

        return ax, counter

    def _add_labels(self, label_site_list, font_size, ax):
        '''
        Create the label of the Voronoi plot.
        '''
        for element in label_site_list:
            if isinstance(element, list):
                self._add_labels(element, font_size, ax)
            else:
                ax.text(element[1][0], element[1][1], element[0],
                        fontsize = font_size, horizontalalignment = 'center',
                        verticalalignment = 'center')
        return

    def _compute_error(self, polygon_value_list, total_value, total_area):
        '''
        compute the gross error of the Voronoi plot.
        '''
        gross_error = 0
        for element in polygon_value_list:
            if isinstance(element, list):
                gross_error += self._compute_error(
                    element, total_value, total_area)
            else:
                gross_error += abs(element[0].area -
                                   total_area * element[1] / total_value)
        return gross_error

    def _find_total(self, val):
        # type: (Union[int, float, dict]) -> float
        '''
        Find the total value within a number or nestable dictionary.
        '''
        total = 0.0
        if isinstance(val, (float, int)):
            total += val
        else:
            for value in val.values():
                total += self._find_total(value)
        return total

    def _compute_boundaries(self, dic,
                            side_length = (4, 4),
                            custom_shape_vertices = None,
                            voronoi_list_old = None):
        """
        Main function used for computing layered voronoi diagram based on
        existing voronoi diagram to ensure the location of each block stays
        the same.
        """
        if custom_shape_vertices is None:
            canvas_vertices = np.array(
                [[0, 0], [side_length[0], 0],
                 [side_length[0], side_length[1]], [0, side_length[1]]])
            canvas_obj = PolygonClass(canvas_vertices)
        else:
            canvas_obj = PolygonClass(custom_shape_vertices)

        labels = list(dic.keys())
        values = [self._find_total(dic[key]) for key in dic]
        if voronoi_list_old is not None:
            voronoi_old = voronoi_list_old[0]
            voronoi_out, error_0 = self._voronoi_main_function(
                labels, values, canvas_obj, voronoi_old = voronoi_old)
        else:
            voronoi_out, error_0 = self._voronoi_main_function(
                labels, values, canvas_obj)
        voronoi_list = [voronoi_out]
        polygon_value_list = [[] for _ in dic]  # type: List[Iterable]
        label_site_list = [[] for _ in dic]  # type: List[Iterable]

        for i, value in enumerate(dic.values()):
            if isinstance(value, (float, int)):
                polygon_value_list[i] = (voronoi_out.polygons[i], values[i])
                label_site_list[i] = (labels[i], voronoi_out.sites[i])

            elif isinstance(value, dict):
                if voronoi_list_old is not None:
                    voronoi_old = voronoi_list_old[1]
                    voronoi_list_old.pop(0)
                    new_voronoi, polygon_value_list[i], label_site_list[i] = self._compute_boundaries(
                        value, voronoi_list_old = voronoi_old,
                        custom_shape_vertices = voronoi_out.polygons[i].xy)
                else:
                    new_voronoi, polygon_value_list[i], label_site_list[i] = self._compute_boundaries(
                        value, custom_shape_vertices = voronoi_out.polygons[i].xy)
                voronoi_list.append(new_voronoi)

        return voronoi_list, polygon_value_list, label_site_list

    def _voronoi_main_function(self, points, values, canvas_obj,
                               voronoi_old = None):
        """ 
        Call voronoi_treemap to compute the voronoi diagram, with baseline
        number of iterations = 75. If the computed voronoi diagram has error
        >= 0.15, increase the iterations to 150 and so on until the error
        < 0.15. If the error increases after increasing the number of iterations
        , recompute the whole voronoi diagram all over again.

        Please be careful when interpreting the error output. An error between
        0.5~0.15 is expected. If the error = 0, that means the whole voronoi
        diagram is not properly generated.
        """
        if voronoi_old is None:
            voronoi, error = self._voronoi_treemap(canvas_obj, points, values)
            counter = 0
            while error > 0.15:
                voronoi, error_new = self._voronoi_treemap_recal(voronoi)
                if error_new <= 0.15:
                    error = error_new
                    break
                else:
                    if error_new >= error:
                        voronoi, error_new = self._voronoi_treemap(
                            canvas_obj, points, values)
                error = error_new
                counter += 1
                if counter >= 50:
                    print('This random seed is not able to produce an error less '
                          'than 0.15. Please choose another random seed instead.')
                    break
        else:
            voronoi_merge = voronoi_old
            voronoi_merge.values = values
            voronoi_merge.canvas_obj = canvas_obj
            voronoi, error = self._voronoi_treemap_recal(voronoi_merge)
            counter = 0
            while error > 0.15:
                voronoi, error_new = self._voronoi_treemap_recal(voronoi)
                if error_new <= 0.15:
                    error = error_new
                    break
                else:
                    if error_new >= error:
                        voronoi, error_new = self._voronoi_treemap(
                            canvas_obj, points, values)
                error = error_new
                counter += 1
                if counter >= 50:
                    print('This random seed is not able to produce an error less '
                          'than 0.15. Please choose another random seed instead.')
                    break
        return voronoi, error

    def _voronoi_treemap(self, canvas_obj, points, values):
        """
        The main working horse of the whole function. This function compute the
        initial voronoi diagram, adapt the location of each site and the weight
        of each site, re-compute the voronoi diagram, and compute the error of
        the voronoi diagram. If the error is too high, repeat this process until
        the error falls below error threshold or the number of iteration exceeds
        i_max(= 75).

        Certain amount of error is intrinsic to this plot. Therefore, blindly
        increasing the number of maximum iteration does not guarantee a decrease
        in error.
        """
        n_points = len(points)

        if n_points == 1:
            # if there is only 1 point, we don't need to compute the diagram
            voronoi = VoronoiClass([canvas_obj],
                                   canvas_obj.xy.mean(axis = 0).reshape((-1, 2)),
                                   canvas_obj.area, values, canvas_obj)
            error = 0
            voronoi.points = points
            return voronoi, error

        else:
            # initialization
            sites = canvas_obj.random_points_in_canvas(n_points)
            weights = np.ones(n_points)*canvas_obj.area/n_points
            error = float("inf")

            # calculate initial voronoi plot
            voronoi = self._compute_power_diagram(
                canvas_obj, sites, weights, values)

            # adjust position and weight until error becomes acceptable
            i = 0
            while i < self.i_max:
                voronoi.adapt_positions_weights()
                voronoi = self._recompute_power_diagram(voronoi)
                voronoi.adapt_weights()
                voronoi = self._recompute_power_diagram(voronoi)
                error = voronoi.compute_error()
                i += 1
                if error < self.err_thres:
                    break
            voronoi.points = points

        return voronoi, error

    def _voronoi_treemap_recal(self, voronoi):
        """
        Basically the same function as _voronoi_treemap but skip the
        initialization part. This ensures the new plot is built on previous
        results and the previous iterations are not wasted.
        """
        points = voronoi.points
        i = 0
        error = np.inf
        while i < self.i_max:
            voronoi.adapt_positions_weights()
            voronoi = self._recompute_power_diagram(voronoi)
            voronoi.adapt_weights()
            voronoi = self._recompute_power_diagram(voronoi)
            error = voronoi.compute_error()
            i += 1
            if error < self.err_thres:
                break
        voronoi.points = points
        return voronoi, error

    def _compute_power_diagram(self, canvas_obj, sites, weights, values):
        """
        Compute the draft of the voronoi diagram, or the power diagram following
        the nomeclature of our code reference. This function seperate the cases
        for n = 2, 3, 4 and above. For n = 2, since there should always be a
        solution, the function will rescue itself by resetting the initial sites
        until a solution is found.
        """
        n_polygon = len(weights)
        if n_polygon < 4:
            # directly compute by intersecting the bisectors
            if n_polygon == 2:
                x_col = sites[:, 0].reshape((-1, 1))
                y_col = sites[:, 1].reshape((-1, 1))
                bisector = LineClass(x_col, y_col, weights)

                # find intersect point with canvas
                edge_new = bisector.find_intersect_with_canvas(canvas_obj)
                while len(edge_new) == 0: 
                    sites = canvas_obj.random_points_in_canvas(n_polygon)
                    x_col = sites[:, 0].reshape((-1, 1))
                    y_col = sites[:, 1].reshape((-1, 1))
                    bisector = LineClass(x_col, y_col, weights)
                    edge_new = bisector.find_intersect_with_canvas(canvas_obj)
                polygons_all = self._divide_polygons(
                    sites, canvas_obj, edge_new)

            else:
                # n_polygon = 3
                site_label = [0, 1, 2]
                sites, ray_obj_all = self._find_bisector_positive_ray(
                    sites, weights, site_label, canvas_obj)
                polygons_all = self._divide_polygons(
                    sites, canvas_obj, ray_obj_all)

        else:
            x_col = sites[:, 0].reshape((-1, 1))
            y_col = sites[:, 1].reshape((-1, 1))

            # 1. map sites into dual space
            dual_sites = np.hstack(
                (x_col, y_col, x_col**2 + y_col**2 - weights.reshape(-1, 1)))

            # 2. find the convex hull of the points in dual space
            faces = ConvexHull(dual_sites)
            simplices = faces.simplices

            # 3. select only the lower convex hull
            # and convert the triangle in dual space into the intersection
            # points(ipoints) in normal space
            ipoints, simplices_prune = self._prune_and_convert_to_ipoints(
                dual_sites, simplices, canvas_obj)

            # 4. compute the ray objects of each ipoints
            ray_obj_all = self._find_multiple_positive_ray(
                sites, ipoints, simplices_prune)

            # 5. compute the coordinate of each polygon
            polygons_all = self._divide_polygons(
                sites, canvas_obj, (ray_obj_all, simplices_prune, ipoints))  

        # convert to object and find the site with largest
        # area_target/area_current ratio, which will displace its neighbor
        voronoi = VoronoiClass(polygons_all, sites, weights, values, canvas_obj)
        voronoi.find_sites_causing_displacement()
        return voronoi

    def _recompute_power_diagram(self, voronoi):
        """
        the same function as _compute_power_diagram but accept an existing
        voronoi object as input.
        """
        sites = voronoi.sites
        weights = voronoi.weights
        if voronoi.n_sites < 4:
            # directly compute by intersecting the bisectors
            if voronoi.n_sites == 2:
                x_col = sites[:, 0].reshape((-1, 1))
                y_col = sites[:, 1].reshape((-1, 1))
                bisector = LineClass(x_col, y_col, weights)

                # find intersect point with canvas
                edge_new = bisector.find_intersect_with_canvas(
                    voronoi.canvas_obj)
                while len(edge_new) == 0: 
                    sites =  voronoi.canvas_obj.random_points_in_canvas(
                        voronoi.n_sites)
                    x_col = sites[:, 0].reshape((-1, 1))
                    y_col = sites[:, 1].reshape((-1, 1))
                    bisector = LineClass(x_col, y_col, weights)
                    edge_new = bisector.find_intersect_with_canvas(
                        voronoi.canvas_obj)

                polygons_all = self._divide_polygons(
                    sites, voronoi.canvas_obj, edge_new)

            else:
                # n_polygon = 3
                site_label = [0, 1, 2]
                sites, ray_obj_all = self._find_bisector_positive_ray(
                    sites, weights, site_label, voronoi.canvas_obj)
                polygons_all = self._divide_polygons(
                    sites, voronoi.canvas_obj, ray_obj_all)

        else:
            x_col = sites[:, 0].reshape((-1, 1))
            y_col = sites[:, 1].reshape((-1, 1))

            # 1. map sites into dual space
            dual_sites = np.hstack(
                (x_col, y_col, x_col**2 + y_col**2 - weights.reshape(-1, 1)))

            # 2. find the convex hull of the points in dual space
            faces = ConvexHull(dual_sites)
            simplices = faces.simplices

            # 3. select only the lower convex hull and convert the triangle in
            # dual space into the intersection points(ipoints) in normal space
            ipoints, simplices_prune = self._prune_and_convert_to_ipoints(
                dual_sites, simplices, voronoi.canvas_obj)

            # 4. compute the ray objects of each ipoints
            ray_obj_all = self._find_multiple_positive_ray(
                sites, ipoints, simplices_prune)

            # 5. compute the coordinate of each polygon
            polygons_all = self._divide_polygons(
                sites, voronoi.canvas_obj,
                (ray_obj_all, simplices_prune, ipoints))

        # convert to object and find the site with largest
        # area_target/area_current ratio, which will displace its neighbor
        voronoi = VoronoiClass(
            polygons_all, sites, weights, voronoi.values, voronoi.canvas_obj)
        voronoi.find_sites_causing_displacement()

        return voronoi

    def _divide_polygons(self, sites, canvas_obj, args):
        n_corners = canvas_obj.n_corners
        n_sites = len(sites)

        if n_sites == 2:
            edge = args
            # 1. convert edge to ax + by = c
            [[x1, y1], [x2, y2]] = edge
            a = (y2 - y1)
            b = - (x2 - x1)
            c = x1*(y2 - y1) - y1*(x2 - x1)

            # 2. calculate for each site and each corner of the canvas,
            # ax + by - c > 0 or < 0 to determine if they are on the same side
            # to the edge.
            above_below_sites = (np.dot(sites, np.array([a, b])) - c) >= 0
            above_below_corners = (np.dot(
                canvas_obj.xy, np.array([a, b])) - c) >= 0
            polygons_all = []

            # 3. assign the corners which is on the same side as each site, and
            # determine the final coordinates of the 2 polygons.
            if above_below_sites[0] != above_below_sites[1]:
                for i in range(n_sites):
                    corner_polygon = canvas_obj.xy[
                        above_below_corners == above_below_sites[i]]
                    corner_polygon = np.vstack([corner_polygon, edge])
                    corner_polygon = (
                        corner_polygon*(10**9)).astype(int)/(10.**9)
                    corner_polygon = np.unique(corner_polygon, axis = 0)
                    polygon = PolygonClass(corner_polygon)
                    polygons_all.append(polygon)

        elif n_sites == 3:
            ray_obj_all = args

            # 0. convert each ray to a_r1*x + b_r1*y + c_r1 = 0
            a_r1, a_r2, a_r3 = ray_obj_all[0].a_r, ray_obj_all[1].a_r, \
                               ray_obj_all[2].a_r
            b_r1, b_r2, b_r3 = ray_obj_all[0].b_r, ray_obj_all[1].b_r, \
                               ray_obj_all[2].b_r
            c_r1, c_r2, c_r3 = ray_obj_all[0].c_r, ray_obj_all[1].c_r, \
                               ray_obj_all[2].c_r
            coefficient_matrix = np.array(
                [[a_r1, a_r2, a_r3], [b_r1, b_r2, b_r3]])
            c_r_matrix = np.array([c_r1, c_r2, c_r3])

            # 1. calculate the edge created by each ray by finding the
            # intersection point of each ray with canvas
            for i in range(3):
                ray_obj_all[i].find_ray_intersect_with_canvas(canvas_obj)

            # 2. determine if each site and each corner is above or below each
            # ray by calculating Ax+By+C >0 or <0 (AB table = above/below table)
            above_below_corners = (np.dot(canvas_obj.xy, coefficient_matrix) 
                + np.tile(c_r_matrix, (n_corners, 1)) >= 0)
            above_below_sites = (np.dot(sites, coefficient_matrix) 
                + np.tile(c_r_matrix, (n_sites, 1)) >= 0)

            # 3. group the corners, edges, intersection points that are on the
            # same side as each site and form the coordinates of each polygon.
            polygons_all = [None for _ in range(3)]  # type: List[Optional[PolygonClass]]
            indices = np.arange(3)
            for i in range(3):
                ab_corners_temp = above_below_corners[:,indices != i]
                ab_site_temp = above_below_sites[i, indices != i]
                target_corner = (np.equal(
                    ab_corners_temp, ab_site_temp).sum(axis = 1) == 2)
                # append the included corners by 2 rays
                corner_polygon = canvas_obj.xy[target_corner, :]
                j_list = np.arange(3)
                j_list = j_list[j_list != i]
                for j in j_list:
                    if not len(ray_obj_all[j].edge) == 0:
                        if len(corner_polygon) == 0:
                            corner_polygon = ray_obj_all[j].edge

                        else:
                            corner_polygon = np.vstack(
                                [corner_polygon, ray_obj_all[j].edge])

                if len(corner_polygon) >= 3:
                    corner_polygon = (corner_polygon*(10**9)).astype(int)/(10.**9)
                    corner_polygon = np.unique(corner_polygon, axis = 0)
                    polygon = PolygonClass(corner_polygon)
                    polygons_all[i] = polygon

        else:
            ray_obj_all = args[0]
            simplices_prune = args[1]
            ipoints = args[2]
            n_ipoints = len(simplices_prune)

            # 1. for all the intersection points, we first determine whether
            # each site and each corners are above or below each ray.
            ab_corners_table = [np.arange(0) for _ in range(n_ipoints)]
            ab_sites_table = [np.arange(0) for _ in range(n_ipoints)]
            for i in range(n_ipoints):
                if ray_obj_all[i]:
                    # 1-1. convert each ray to a_r1*x + b_r1*y + c_r1 = 0
                    a_r1, a_r2, a_r3 = ray_obj_all[i][0].a_r, \
                                       ray_obj_all[i][1].a_r, \
                                       ray_obj_all[i][2].a_r
                    b_r1, b_r2, b_r3 = ray_obj_all[i][0].b_r, \
                                       ray_obj_all[i][1].b_r, \
                                       ray_obj_all[i][2].b_r
                    c_r1, c_r2, c_r3 = ray_obj_all[i][0].c_r, \
                                       ray_obj_all[i][1].c_r, \
                                       ray_obj_all[i][2].c_r
                    coefficient_matrix = np.array(
                        [[a_r1, a_r2, a_r3], [b_r1, b_r2, b_r3]])
                    c_r_matrix = np.array([c_r1, c_r2, c_r3])

                    # 1-2. calculate the edge created by each ray
                    for j in range(3):
                        ray_obj_all[i][j].find_ray_intersect_with_canvas_and_ipoints(
                            canvas_obj, ipoints, simplices_prune)

                    # 1-3. calculate Ax+By+C >0 or <0 (AB table = above/below table)
                    ab_corners_table[i] = (np.dot(canvas_obj.xy, coefficient_matrix) 
                        + np.tile(c_r_matrix, (n_corners, 1)) >= 0)
                    ab_sites_table[i] = (np.dot(sites, coefficient_matrix) 
                        + np.tile(c_r_matrix, (n_sites, 1)) >= 0)

            # 2. by referring to the ab_corners_table and ab_sites_table
            # find the coordinates of each polygon.
            polygons_all = [None for _ in range(n_sites)]
            corner_polygon_all = [np.arange(0) for _ in range(n_sites)]
            indices = np.arange(3)
            for k in range(n_sites):
                iP_index_required = np.where(
                    np.sum((simplices_prune == k), axis = 1) == 1)[0]
                target_corner_test = []
                for m in iP_index_required:
                    if ray_obj_all[m]:
                        simplex = simplices_prune[m]
                        [loc_012] = np.where(simplex == k)[0]
                        # 2-1. find the potential corners if they are on the
                        # same side as each site:
                        ab_corners_temp = ab_corners_table[m][:, indices != loc_012]
                        ab_site_temp = ab_sites_table[m][k, indices != loc_012]
                        target_corner_test.append(
                            set(np.where(
                                np.equal(
                                    ab_corners_temp, ab_site_temp).sum(axis = 1) == 2)[0]))

                        # 2-2. append the edges formed by nearby rays.
                        n_list = np.arange(3)
                        n_list = n_list[n_list != loc_012]
                        for n in n_list:
                            if not len(ray_obj_all[m][n].edge) == 0:
                                if len(corner_polygon_all[k]) == 0:
                                    corner_polygon_all[k] = ray_obj_all[m][n].edge
                                else:
                                    corner_polygon_all[k] = np.vstack(
                                        [corner_polygon_all[k], ray_obj_all[m][n].edge])

                # 2-3. reconcile the potential corners across different ipoints,
                # find target_corner (some corners may be on the same side as
                # the site with respect to one ipoint, but not other nearby
                # ipoints. We have to make sure that the final corners included
                # are on the same side as each site with respect to all ipoints
                # nearby.)
                if not len(target_corner_test) == 0:
                    for p in range(1, len(target_corner_test)):
                        target_corner_test[0] = target_corner_test[0].intersection(
                            target_corner_test[p])

                # 2-4. append the target corners
                if not len(target_corner_test) == 0:
                    if len(corner_polygon_all[k]) == 0:
                        corner_polygon_all[k] = canvas_obj.xy[
                                                list(target_corner_test[0]), :]
                    else:
                        corner_polygon_all[k] = np.vstack(
                            [corner_polygon_all[k], canvas_obj.xy[
                                                    list(target_corner_test[0]), :]])

                if not len(corner_polygon_all[k]) == 0:
                    # 2-5. remove redundant corners of each polygon
                    corner_polygon_all[k] = (
                        corner_polygon_all[k]*(10**9)).astype(int)/(10.**9)
                    corner_polygon_all[k] = np.unique(
                        corner_polygon_all[k], axis = 0)

                    # 2-6. create and store the polygon object for each site
                    if len(corner_polygon_all[k]) >= 3:
                        polygon = PolygonClass(corner_polygon_all[k])
                        polygons_all[k] = polygon

        return polygons_all

    def _find_bisector_positive_ray(self, sites, weights, site_label,
                                    canvas_obj):
        """
        This function is designed for n=3. This function compute the 3 bisectors
        among 3 sites, find the intersect point (ipoint) of these 3 bisectors,
        and define the positive direction of each ray.
        """
        def _compute_determinant(x_col, y_col, weights):
            [[x0], [x1], [x2]] = x_col
            [[y0], [y1], [y2]] = y_col
            [w0, w1, w2] = weights
            a1 = 2*x0
            a2 = 2*x1
            a3 = 2*x2
            b1 = 2*y0
            b2 = 2*y1
            b3 = 2*y2
            d1 = - (x0**2 + y0**2 - w0)
            d2 = - (x1**2 + y1**2 - w1)
            d3 = -(x2**2 + y2**2 - w2)
            det = - a1*b2 - a2*b3 - a3*b1 + a3*b2 + a1*b3 + a2*b1
            if det == 0:
                x_int = np.nan
                y_int = np.nan
            else:
                x_int = (- d1*b2 - d2*b3 - d3*b1 + d3*b2 + d1*b3 + d2*b1)/(- det)
                y_int = (- a1*d2 - a2*d3 - a3*d1 + a3*d2 + a1*d3 + a2*d1)/(- det)
            return det, x_int, y_int

        # 1. Find the 3 bisectors and the ipoint
        # If no solution is possible, reset the location of sites until there
        # is a solution
        x_col = sites[:, 0].reshape((-1, 1))
        y_col = sites[:, 1].reshape((-1, 1))
        det, x_int, y_int = _compute_determinant(x_col, y_col, weights)

        while det == 0:
            sites = canvas_obj.random_points_in_canvas(3)
            x_col = sites[:, 0].reshape((-1, 1))
            y_col = sites[:, 1].reshape((-1, 1))
            det, x_int, y_int = _compute_determinant(x_col, y_col, weights)

        ipoint = np.array([x_int, y_int])

        # If the intersection point is not within the canvas, 
        # reset the location of sites until it is within canvas
        while not canvas_obj.point_is_within_canvas(ipoint):
            det = 0
            while det == 0:
                sites = canvas_obj.random_points_in_canvas(3)
                x_col = sites[:, 0].reshape((-1, 1))
                y_col = sites[:, 1].reshape((-1, 1))
                det, x_int, y_int = _compute_determinant(x_col, y_col, weights)

            ipoint = np.array([x_int, y_int])

        # 2. find the positive direction of the 3 rays.
        ray_obj_all = []
        if self._point_is_within_triangle(ipoint, sites):
            # if the ipoint is within the triangle formed by 3 sites, the
            # positive direction of the ray is the one that moves away from the
            # opposite site.
            for k in range(3):
                site_remain = np.arange(3)
                site_remain = site_remain[site_remain != k]
                tag1 = site_remain[0]
                tag2 = site_remain[1]
                t_bisector = np.array(
                    [2*(y_col[tag2][0] - y_col[tag1][0]),
                     2*(x_col[tag1][0] - x_col[tag2][0])])
                if np.dot(t_bisector,(sites[k, :] - ipoint)) > 0:
                    t_bisector = - t_bisector
                ray_obj_all.append(
                    RayClass(ipoint, t_bisector,
                             [site_label[tag1], site_label[tag2]],
                             site_label[k]))

        else:
            # if the ipoint is not within the triangle formed by 3 sites, the
            # positive direction is defined in a way that can divide the canvas
            # into 3 polygons and each polygon contains 1 site.
            com_line = (sites[[1, 2, 0], :] + sites[[2, 0, 1], :])/2
            point_opposite = np.argmin(
                (com_line - ipoint)[:, 0]**2 + (com_line - ipoint)[:, 1]**2)
            for k in range(3):
                site_remain = np.arange(3)
                site_remain = site_remain[site_remain != k]
                tag1 = site_remain[0]
                tag2 = site_remain[1]
                t_bisector = np.array([2*(y_col[tag2][0] - y_col[tag1][0]), 
                    2*(x_col[tag1][0] - x_col[tag2][0])])
                if k == point_opposite:
                    if np.dot(t_bisector,(sites[k, :] - ipoint)) > 0:
                        t_bisector = - t_bisector
                    ray_obj_all.append(
                        RayClass(ipoint, t_bisector,
                                 [site_label[tag1], site_label[tag2]],
                                 site_label[k]))
                else:
                    if np.dot(t_bisector,(com_line[k, :] - ipoint)) < 0:
                        t_bisector = - t_bisector
                    ray_obj_all.append(
                        RayClass(ipoint, t_bisector,
                                 [site_label[tag1], site_label[tag2]],
                                 site_label[k]))

        return sites, ray_obj_all

    def _find_multiple_positive_ray(self, sites, ipoints, simplices_prune):
        """
        This function is designed for n>=4. This function define the positive
        direction of each ray in a way that is the same as
        _find_bisector_positive_ray.
        """    
        n_ipoints = len(simplices_prune)
        ray_obj_all = [[] for _ in range(n_ipoints)]

        for i in range(n_ipoints):
            ipoint = ipoints[i, 0:2]
            site_label = simplices_prune[i]
            sites_temp = sites[simplices_prune[i], :]
            x_col_temp = sites[simplices_prune[i],0].reshape((-1, 1)) 
            y_col_temp = sites[simplices_prune[i],1].reshape((-1, 1))

            ray_obj_ipoint = []
            if self._point_is_within_triangle(ipoint, sites_temp):
                for k in range(3):
                    site_remain = np.arange(3)
                    site_remain = site_remain[site_remain != k]
                    tag1 = site_remain[0]
                    tag2 = site_remain[1]
                    t_bisector = np.array(
                        [2*(y_col_temp[tag2][0] - y_col_temp[tag1][0]),
                        2*(x_col_temp[tag1][0] - x_col_temp[tag2][0])])
                    if np.dot(t_bisector, (sites_temp[k, :] - ipoint)) > 0:
                        t_bisector = - t_bisector
                    ray_obj = RayClass(
                        ipoint, t_bisector,
                        [site_label[tag1], site_label[tag2]], site_label[k])
                    ray_obj.index_ip = i
                    ray_obj.index_ray = k
                    ray_obj_ipoint.append(ray_obj)

            else:
                com_line = (sites_temp[[1, 2, 0], :] + sites_temp[[2, 0, 1], :])/2
                point_opposite = np.argmin(
                    (com_line - ipoint)[:, 0]**2 + (com_line - ipoint)[:, 1]**2)
                for k in range(3):
                    site_remain = np.arange(3)
                    site_remain = site_remain[site_remain != k]
                    tag1 = site_remain[0]
                    tag2 = site_remain[1]
                    t_bisector = np.array(
                        [2*(y_col_temp[tag2][0] - y_col_temp[tag1][0]),
                        2*(x_col_temp[tag1][0] - x_col_temp[tag2][0])])
                    if k == point_opposite:
                        if np.dot(t_bisector, (sites_temp[k, :] - ipoint)) > 0:
                            t_bisector = - t_bisector
                        ray_obj = RayClass(
                            ipoint, t_bisector,
                            [site_label[tag1], site_label[tag2]], site_label[k])
                        ray_obj.index_ip = i
                        ray_obj.index_ray = k
                        ray_obj_ipoint.append(ray_obj)
                    else:
                        if np.dot(t_bisector, (com_line[k, :] - ipoint)) < 0:
                            t_bisector = - t_bisector
                        ray_obj = RayClass(
                            ipoint, t_bisector,
                            [site_label[tag1], site_label[tag2]], site_label[k])
                        ray_obj.index_ip = i
                        ray_obj.index_ray = k
                        ray_obj_ipoint.append(ray_obj)

            ray_obj_all[i] = ray_obj_ipoint
        return ray_obj_all

    def _prune_and_convert_to_ipoints(self, dual_sites, simplices, canvas_obj):
        """
        This function determines which ipoints should be kept. 
        Only the ipoints that (1) belong to the lower convex hull and
        (2) locate within canvas will be kept.
        """
        com_dualsites = dual_sites.mean(axis = 0)
        simplices_prune = simplices
        ipoints = np.array([])
        prune_list = []
        for i, simplex in enumerate(simplices):
            x1, y1, z1 = dual_sites[simplex[0]]
            x2, y2, z2 = dual_sites[simplex[1]]
            x3, y3, z3 = dual_sites[simplex[2]]
            alpha = y1*(z2 - z3) + y2*(z3 - z1) + y3*(z1 - z2)
            beta = z1*(x2 - x3) + z2*(x3 - x1) + z3*(x1 - x2)
            gamma = x1*(y2 - y3) + x2*(y3 - y1) + x3*(y1 - y2)
            delta = x1*(y2*z3 - y3*z2) + x2*(y3*z1 - y1*z3) + x3*(y1*z2 - y2*z1)
            a = - alpha/gamma
            b = - beta/gamma
            c = delta/gamma # z = ax+by+c => nf = (a, b, -1)

            # (1) determine if the face belongs to the lower convex Hull
            com_face = (dual_sites[simplex]).mean(axis = 0)
            if np.dot([a, b, -1], (com_face - com_dualsites)) > 0:

                # (2) determine if the ipoint is within canvas
                if canvas_obj.point_is_within_canvas(np.array([a/2, b/2])):
                    ipoints = np.append(
                        ipoints, np.array([a/2, b/2, -c]), axis = 0)
                else:
                    prune_list.append(i)
            else:
                prune_list.append(i)

        ipoints = ipoints.reshape(-1, 3)
        simplices_prune = np.delete(simplices_prune, prune_list, axis = 0)
        return ipoints, simplices_prune

    def _point_is_within_triangle(self, p, triangle):
        """
        Check if a point is within a triangle. The algorithm used is the same as
        point_is_within_canvas. However, this function does not require the
        input to be a POLYGON object.
        """
        a = (triangle == p)
        if len(np.nonzero(a[:, 0] * a[:, 1])[0]) == 1:
            result = True

        else:
            sum_angle = 0
            for i in range(3):
                edge = triangle[[i-1, i], :]
                a = edge - p
                sum_angle += angle(a[0, :], a[1, :])

            if not math.isnan(sum_angle):
                result = int((sum_angle - 2*math.pi)*(10**9))/(10.**9) == 0

            else:
                result = False

        return result
