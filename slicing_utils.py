import math
import time
import Rhino
import rhinoscriptsyntax as rs
import extruder_turtle  
import turtle_utilities as tu
from extruder_turtle import *

import geometry_utils
from geometry_utils import *

import vertical_utils
from vertical_utils import *

import contour_utils
from contour_utils import *

import graph_utils
from graph_utils import *


def get_corner(t, outer_curve, inner_curve, points):
    offset = float(t.get_extrude_width())

    prec = 100
    pnts = rs.DivideCurve(inner_curve, prec)
    center = (sum([p.X for p in pnts])/prec, sum([p.Y for p in pnts])/prec, sum([p.Z for p in pnts])/prec)

    closest = {"point": 0, "distance": 1000000}

    for p in range(len(points)):
        pnt = points[p]
        prev_pnt = points[(p-1) % len(points)]
        next_pnt = points[(p+1) % len(points)]
        v1 = rs.VectorUnitize(rs.VectorSubtract(prev_pnt, pnt))
        v2 = rs.VectorUnitize(rs.VectorSubtract(next_pnt, pnt))
        dot = rs.VectorDotProduct(v1, v2)
        add = rs.VectorAdd(v1, v2)

        if not (add.X == 0 and add.Y == 0 and add.Z == 0):
            # acos has a domain of [-1, 1], sometimes rhinoscript
            # returns close to 1 due to a precision error
            angle = math.acos(min(max(dot, -1), 1)) * (180 / math.pi)

            if angle < 160 and angle > 20:
                dist = rs.Distance(pnt, center)
                # check that the corner isn't "inverted"
                inside = rs.PointInPlanarClosedCurve(rs.VectorAdd(pnt, rs.VectorScale(rs.VectorUnitize(add), offset/4)), outer_curve)
                if dist < closest["distance"] and inside:
                    closest["distance"] = dist
                    closest["point"] = p

    return closest["point"]


def spiral_contours(t, isocontours, start_index):
    if len(isocontours) == 0:
        print("Error: no isocontours passed in")
        return
    # a single spirallable region
    # connect each isocontour with the one succeeding it
    spiral = []
    offset = float(t.get_extrude_width())
    num_pnts = get_num_points(isocontours[0], offset)
    points = rs.DivideCurve(isocontours[0], num_pnts)

    # choose appropriate starting index on outermost contour
    # get center point of innermost contour to compare
    if not start_index:
        start_index = get_corner(t, isocontours[0], isocontours[-1], points)

    # spiral region
    start_point = None
    spiral_contour_indices = []
    for i in range(len(isocontours)):
        start_point = points[start_index]

        # get break point on isocontour for connecting spiral
        marching_order = range(start_index, -1, -1) + range(len(points)-1, start_index, -1)
        closest = {"point": None, "distance": 1000000}

        for j in marching_order:
            dist = offset - rs.Distance(start_point, points[j])
            if dist < -offset/5:
                break
            if abs(dist) < closest["distance"]:
                closest["distance"] = dist
                closest["point"] = j
        break_index = closest["point"]
        break_point = points[break_index]

        # append the points from the isocontour from start_point
        # to break_point to the spiral
        indices = []
        if start_index > break_index:
            indices = range(start_index, len(points))+range(0, break_index+1)
        elif start_index < break_index:
            indices =  range(start_index, break_index+1)
        else:
            print("Error: start and break index should not be the same")

        spiral = spiral + [points[j] for j in indices]
        spiral_contour_indices.append(len(spiral)-1)

        # find closest point in next contour to the break point
        # if we are not at the centermost contour
        if i < len(isocontours) - 1:
            next_points = rs.DivideCurve(isocontours[i+1], get_num_points(isocontours[i+1], offset))
            closest = {"point": None, "distance": 1000000}
            for j in range(len(next_points)):
                dist = rs.Distance(break_point, next_points[j])
                if dist < closest["distance"]:
                    closest["distance"] = dist
                    closest["point"] = j

            # set next start index and next points
            start_index = closest["point"]
            points = next_points

    # return spiral and corresponding contour index
    return spiral, spiral_contour_indices


def fermat_spiral(t, spiral, indices):
    offset = float(t.get_extrude_width())

    # spiral path inward
    in_spiral = []
    start_index = 0
    order = range(0, len(indices), 2)
    for i in order:
        connection = {"point": None, "distance": 1000000}
        for j in range(indices[i]-1, 0, -1):
            dist = offset - rs.Distance(spiral[indices[i]], spiral[j])
            if dist < -offset/5:
                break
            if abs(dist) < connection["distance"]:
                connection["distance"] = dist
                connection["point"] = j

        in_spiral = in_spiral + range(start_index, connection["point"])

        if i < len(indices) - 1:
            start_index = indices[i+1]
        else:
            start_index = len(spiral) - 1

    # if the number of isocontours is even, include a link
    # from the end point of the spiral to the second to last index
    out_spiral = []
    start_index = indices[0]
    order = range(1, len(indices), 2)
    for i in order:
        connection = {"point": None, "distance": 1000000}
        for j in range(indices[i]-1, 0, -1):
            dist = offset - rs.Distance(spiral[indices[i]], spiral[j])
            if dist < -offset/5:
                break
            if abs(dist) < connection["distance"]:
                connection["distance"] = dist
                connection["point"] = j

        out_spiral = out_spiral + range(start_index, connection["point"])

        if i < len(indices) - 1:
            start_index = indices[i+1]
        else:
            start_index = len(spiral) - 1

    new_spiral = in_spiral + list(reversed(out_spiral))
    new_spiral = [spiral[i] for i in new_spiral]

    return new_spiral


def get_leaves(node):
    if len(node["children"]) == 0:
        return [node]
    else:
        leaves = []
        for n in node["children"]:
            leaves = leaves + get_leaves(n)
        return leaves


def get_branches(node):
    if len(node["children"]) == 0:
        return []
    else:
        branches = [node]
        for n in node["children"]:
            branches = branches + get_branches(n)
        return branches


def get_all_nodes(node):
    if len(node["children"]) == 0:
        return [node]
    else:
        all_nodes = [node]
        for n in node["children"]:
            all_nodes = all_nodes + get_all_nodes(n)
        return all_nodes


def set_node_types(node):
    if len(node["children"]) == 0:
        node["type"] = 1
    if len(node["children"]) == 1:
        node["type"] = 1
        set_node_types(node["children"][0])
    if len(node["children"]) > 1:
        node["type"] = 2
        for child in node["children"]:
            set_node_types(child)


def segment_tree(root):
    region_root = {"guid":"0", "curves":[], "children":[]}
    fill_region(region_root, root, 0)
    return region_root


def fill_region(region_node, node, idx):
    # if the number of children is zero, append to curves
    if(len(node["children"]) == 0):
        # add curve to the region
        region_node["curves"].append(node["guid"])
    # if the number of children is one, proceed to
    # the next node in the tree
    elif len(node["children"]) == 1:
        # add curve to the region
        region_node["curves"].append(node["guid"])
        fill_region(region_node, node["children"][0], 0)
    # if the number of children is greater than one,
    # we have found a split and need to add a new node
    elif len(node["children"]) > 1:
        # add new region node and curve to the new region
        new_node = region_node
        if len(region_node['curves']) > 0:
            new_node = {"guid":region_node["guid"]+"_"+str(idx), "curves":[], "parent":region_node, "children":[]}
            region_node["children"].append(new_node)

        new_node["curves"].append(node["guid"])

        idx = 0
        for child in node["children"]:
            if len(child['children']) > 1:
                fill_region(new_node, child, idx)
                idx = idx + 1
            else:
                # add new region node
                new_new_node = {"guid":new_node["guid"]+"_"+str(idx), "curves":[], "parent":new_node, "children":[]}
                new_node["children"].append(new_new_node)

                fill_region(new_new_node, child, 0)
                idx = idx + 1


def connect_spiralled_nodes(t, root):
    find_connections(t, root)
    all_nodes = get_all_nodes(root)
    all_nodes = {node['guid']: node for node in all_nodes}

    path = connect_path(t, root, all_nodes, 0, [])
    spiral = []
    for p in range(len(path)):
        node = all_nodes[path[p][0]]
        indices = get_marching_indices(node, path[p][1], path[p][2], path[p][3])
        spiral = spiral + [node['fermat_spiral'][idx] for idx in indices]

    return spiral


def connect_path(t, node, all_nodes, start_idx, spiral):
    final_idx = len(node['fermat_spiral'])-1
    if node.get('parent'):
        final_idx = next(idx for idx in node['connection'][node['parent']['guid']] if idx != start_idx)

    marching_order, reverse = get_marching_order(node, start_idx, final_idx)
    if not node.get('reverse'):
        node['reverse'] = reverse

    if node['children']:
        everything_sorted = sorted([(k, n) for n in node['connection'] for k in node['connection'][n].keys()], key=lambda x: marching_order.index(x[0]))
        sorted_children = []
        for c in everything_sorted:
            include = True
            for child in sorted_children:
                if child[1] == c[1]:
                    include = False
                    break
            if include and c[0] != start_idx: sorted_children.append(c)

        for c in sorted_children:
            if len(c[1]) > len(node['guid']):
                child = all_nodes[c[1]]
                end_idx = c[0]

                # add path from node to child index to spiral
                spiral.append((node['guid'], start_idx, end_idx, reverse))

                # recursively call connect_path on child
                child_start_idx = node['connection'][child['guid']][end_idx]
                spiral = connect_path(t, child, all_nodes, child_start_idx, spiral)

                # find next start index for node
                other_idx = next(idx for idx in child['connection'][node['guid']] if idx != child_start_idx)
                start_idx = child['connection'][node['guid']][other_idx]

    # close off spiral with remaining portion of node's fermat spiral
    spiral.append((node['guid'], start_idx, final_idx, reverse))

    return spiral


def get_marching_order(node, start, end):
    points = node['fermat_spiral']

    other_idx = len(node['fermat_spiral'])-1
    if node.get('parent'):
        other_idx = next(idx for idx in node['connection'][node['parent']['guid']] if idx != start)

    marching_order = []
    reverse_marching_order = []
    if start > other_idx:
        marching_order = range(start, len(points)) + range(0, other_idx+1)
        reverse_marching_order = range(start, other_idx-1, -1)
    elif start < other_idx:
        marching_order = range(start, other_idx+1)
        reverse_marching_order = range(start, -1, -1) + range(len(points)-1, other_idx-1, -1)
    else:
        print("Error: start and end indices should not be the same")

    if len(marching_order) < len(reverse_marching_order):
        return get_marching_indices(node, start, end, True), True
    else:
        return get_marching_indices(node, start, end, False), False


def get_marching_indices(node, start, end, reverse):
    points = node['fermat_spiral']

    if reverse:
        if start > end:
            return range(start, end-1, -1)
        elif start < end:
            return range(start, -1, -1) + range(len(points)-1, end-1, -1)
        else:
            print("Error: start and end indices should not be the same")
    else:
        if start > end:
            return range(start, len(points)) + range(0, end+1)
        elif start < end:
            return range(start, end+1)
        else:
            print("Error: start and end indices should not be the same")   


def find_connections(t, node):
    for child in node["children"]:
        find_connections(t, child)

    parent = node.get('parent')
    if parent and parent.get('fermat_spiral'):
        connect_node_to_parent(t, node, parent)


def connect_node_to_parent(t, node, parent):
    offset = float(t.get_extrude_width())

    node_start, node_end = get_connection_indices(t, node)

    points = node['fermat_spiral']
    start_pnt = points[node_start]
    end_pnt = points[node_end]

    if node_start == node_end or rs.Distance(start_pnt, end_pnt) == 0:
        print("Error: could not find suitable connection indices", node_start, node_end)
        return

    direction = 90
    vec = rs.VectorSubtract(end_pnt, start_pnt)
    vec = rs.VectorScale(rs.VectorUnitize(rs.VectorRotate(vec, direction, [0, 0, 1])), offset)
    pnt1 = rs.VectorAdd(start_pnt, vec)
    pnt2 = rs.VectorAdd(end_pnt, vec)

    if rs.PointInPlanarClosedCurve(pnt1, node['curves'][0]) or rs.PointInPlanarClosedCurve(pnt2, node['curves'][0]):
        vec = rs.VectorSubtract(end_pnt, start_pnt)
        vec = rs.VectorScale(rs.VectorUnitize(rs.VectorRotate(vec, -direction, [0, 0, 1])), offset)
        pnt1 = rs.VectorAdd(start_pnt, vec)
        pnt2 = rs.VectorAdd(end_pnt, vec)

    # search parent node points for points
    # closest to computed intersection
    closest = {"start": {"point": None, "distance": 1000000}, "end": {"point": None, "distance": 1000000}}
    points = parent['fermat_spiral']
    for p in range(len(points)):
        dist = rs.Distance(points[p], pnt1)
        if dist < closest["start"]["distance"]:
            closest["start"]["distance"] = dist
            closest["start"]["point"] = p

    for p in range(len(points)):
        dist1 = rs.Distance(points[p], pnt2)
        dist2 = rs.Distance(points[p], pnt1) - offset
        if abs(dist2) < offset*1.2 and dist1 < closest["end"]["distance"]:
            closest["end"]["distance"] = dist1
            closest["end"]["point"] = p

    parent_start = closest['start']['point']
    parent_end = closest['end']['point']

    if not parent.get('connection'):
        parent['connection'] = {}
    parent['connection'][node['guid']] = { parent_start: node_start, parent_end: node_end }

    if not node.get('connection'):
        node['connection'] = {}
    node['connection'][parent['guid']] = { node_start: parent_start, node_end: parent_end }


def get_connection_indices(t, node):
    offset = float(t.get_extrude_width())

    start_index = 0
    end_index = len(node['fermat_spiral']) - 1

    if node['type'] == 2:
        points = node['fermat_spiral']
        # verify that we haven't already tried to connect to a node
        # at those indices, otherwise move along curve to a new spot
        available_indices = range(len(points))
        connections = node.get('connection')
        if connections:
            connection = [connections[n].keys() for n in connections]
            for c in connection:
                connect_indices = get_shortest_indices(c[0], c[1], points)
                available_indices = [x for x in available_indices if not x in connect_indices]

            # set start index
            start_index = available_indices[0]

        start_pnt = points[start_index]
        # find end index
        closest = {"point": None, "distance": 1000000}
        for i in available_indices:
            dist = rs.Distance(points[i], start_pnt) - offset
            if dist < -offset*2:
                break
            if abs(dist) < closest["distance"]:
                closest["distance"] = abs(dist)
                closest["point"] = i

        end_index = closest["point"]
    
    return start_index, end_index


def connect_curves(curves, offset):
    # order curves by surface area, outermost curve will
    # have the largest surface area
    curves = sorted(curves, key=lambda x: get_size(x), reverse=True)

    # find closest points between all curves
    closest = {curve: None for curve in curves}
    for c1 in range(len(curves)):
        curve1 = curves[c1]
        pnts1 = rs.DivideCurve(curve1, 100)

        closest[curve1] = {}
        for c2 in range(c1+1, len(curves)):
            curve2 = curves[c2]
            pnts2 = rs.DivideCurve(curve2, 100)

            minimum = 100000
            for pnt in pnts1:
                idx, dist = closest_point(pnt, pnts2)
                if dist < minimum:
                    minimum = dist
                    closest[curve1][curve2] = (pnt, pnts2[idx], dist)

    # construct a fully connected graph where each curve
    # represents a single node in the graph, and edges
    # are the minimum distance between each curve
    graph = Graph()
    for c in range(len(curves)):
        node = Graph_Node(curves[c])
        node.name = 'c' + str(c)
        graph.add_node(node)

    for curve1 in closest:
        for curve2 in closest[curve1]:
            weight = closest[curve1][curve2][2]
            node1 = graph.get_node(curve1)
            node2 = graph.get_node(curve2)
            graph.add_edge(Graph_Edge(node1, node2, weight))
            graph.add_edge(Graph_Edge(node2, node1, weight))

    #graph.print_graph_data()

    # trim graph by ordering edges from most weighted to least
    ordered_edges = get_edge_tuples(graph)

    # remove edges, starting with most weighted and determine
    # if nodes are still reachable from the start node
    start_node = graph.get_node(curves[0])
    for edge in ordered_edges:
        node1 = edge[0]
        node2 = edge[1]
        weight = edge[2]
        graph.edges[node1].pop(node2)
        graph.edges[node2].pop(node1)

        for node in graph.nodes:
            if node != start_node:
                if not graph.check_for_path(start_node, node)[0]:
                    # add edge back if unable to reach all nodes
                    # in graph from the start node (outermost curve)
                    graph.add_edge(Graph_Edge(node1, node2, weight))
                    graph.add_edge(Graph_Edge(node2, node1, weight))
                    break

    #graph.print_graph_data()
    #min_edges = get_edge_tuples(graph)
    #print([(tuple[0].name, tuple[1].name, tuple[2]) for tuple in min_edges])

    # split curves at shortest connection points to other curves
    new_curves = {}
    curve_ends = {}
    max_weight = 0
    for node1 in graph.edges:
        split_points = []
        for node2 in graph.edges[node1]:
            split_point, weight = get_connect_point(closest, node1, node2)
            if weight > max_weight: max_weight = weight
            split_points.append(split_point)

        split_curves, split_ends = split_curve_at(node1.data, split_points, tolerance=offset)
        new_curves[node1] = split_curves
        curve_ends[node1] = split_ends

    all_curves = [curve for n in new_curves for curve in new_curves[n]]

    # join new split curves together with lines between ends
    min_edges = get_edge_tuples(graph)
    for edge in min_edges:
        node1 = edge[0]
        node2 = edge[1]
        connect_point_1, weight = get_connect_point(closest, node1, node2)
        connect_point_2, weight = get_connect_point(closest, node2, node1)

        ends1 = sorted([(end, rs.Distance(end, connect_point_1)) for end in curve_ends[node1]], key=lambda x: x[1])
        pnt1_1 = ends1[0][0]
        pnt1_2 = ends1[1][0]
        ends2 = sorted([(end, rs.Distance(end, connect_point_2)) for end in curve_ends[node2]], key=lambda x: x[1])
        pnt2_1 = ends2[0][0]
        pnt2_2 = ends2[1][0]

        intersect = rs.PlanarClosedCurveContainment(rs.AddCurve([pnt1_1, pnt2_1]), rs.AddCurve([pnt1_2, pnt2_2]))

        if intersect == 0:
            all_curves.append(rs.AddCurve([pnt1_1, pnt2_1]))
            all_curves.append(rs.AddCurve([pnt1_2, pnt2_2]))
        else:
            all_curves.append(rs.AddCurve([pnt1_1, pnt2_2]))
            all_curves.append(rs.AddCurve([pnt1_2, pnt2_1]))

    return rs.JoinCurves(all_curves)


def get_connect_point(closest, node1, node2):
    split_point = None
    connection = closest[node1.data].get(node2.data)
    if connection:
        split_point = connection[0]
    else:
        connection = closest[node2.data].get(node1.data)
        split_point = connection[1]

    weight = connection[2]
    return split_point, weight


def get_edge_tuples(graph):
    ordered_edges = []
    for node1 in graph.edges:
        for node2 in graph.edges[node1]:
            weight = graph.edges[node1][node2]
            if (node1, node2, weight) not in ordered_edges and (node2, node1, weight) not in ordered_edges:
                ordered_edges.append((node1, node2, weight))

    return sorted(ordered_edges, key=lambda x: x[2], reverse=True)


def fill_curves_with_fermat_spiral(t, curves, start_pnt=None, wall_mode=False, walls=3, initial_offset=0.5):
    # connect curves if given more than one
    curve = curves[0]
    if len(curves) > 1:
        curve = connect_curves(curves, float(t.get_extrude_width())/8)

    # slice the shape
    # Generate isocontours
    root, isocontours = get_contours(t, curve, walls=walls, wall_mode=wall_mode, initial_offset=initial_offset)

    travel_paths = []

    region_tree = segment_tree(root)
    set_node_types(region_tree)
    all_nodes = get_all_nodes(region_tree)

    # Spiralling Regions
    for n in all_nodes:
        if len(n.get('curves')) > 0:
            num_pnts = get_num_points(n['curves'][0], float(t.get_extrude_width()))
            if n['type'] == 1:
                start_idx = None
                if start_pnt:
                    start_idx, d = closest_point(start_pnt, rs.DivideCurve(n['curves'][0], num_pnts))
                spiral, indices = spiral_contours(t, n["curves"], start_idx)
                n["fermat_spiral"] = fermat_spiral(t, spiral, indices)
            elif n['type'] == 2:
                n['fermat_spiral'] = rs.DivideCurve(n["curves"][0], num_pnts)
        else:
            print("Error: node with no curves in it at all", n)

    # Connect Spiralled Regions
    final_points = connect_spiralled_nodes(t, region_tree)
    final_curve = rs.AddCurve(final_points)
    final_spiral = rs.DivideCurve(final_curve, int(rs.CurveLength(final_curve)/t.get_resolution()))
   #print(len(final_points), len(final_spiral))
    t.pen_up()
    travel_paths.append(rs.AddCurve([t.get_position(), final_spiral[0]]))
    t.set_position(final_spiral[0].X, final_spiral[0].Y, t.get_position().Z)
    t.set_position(final_spiral[0].X, final_spiral[0].Y, final_spiral[0].Z)
    t.pen_down()
    for p in final_spiral:
        t.set_position(p.X, p.Y, p.Z)
    t.set_position(final_spiral[0].X, final_spiral[0].Y, final_spiral[0].Z)

    return travel_paths, final_spiral


def fill_curves_with_spiral(t, curves, start_pnt=None):
    # connect curves if given more than one
    curve = curves[0]
    if len(curves) > 1:
        curve = connect_curves(curves, float(t.get_extrude_width())/8)

    # slice the shape
    # Generate isocontours
    root, isocontours = get_contours(t, curve)

    region_tree = segment_tree(root)
    set_node_types(region_tree)
    all_nodes = get_all_nodes(region_tree)

    # Spiral Regions
    travel_paths = []
    for node in all_nodes:
        if 'root' in node['curves']: node['curves'].remove('root')
        if len(node['curves']) > 0:
            num_pnts = get_num_points(node['curves'][0], float(t.get_extrude_width()))
            if node['type'] == 1:
                start_idx = None
                if start_pnt:
                    start_idx, d = closest_point(start_pnt, rs.DivideCurve(node['curves'][0], num_pnts))
                spiral, indices = spiral_contours(t, node["curves"], start_idx)
                t.pen_up()
                travel_paths.append(rs.AddCurve([t.get_position(), spiral[0]]))
                t.set_position(spiral[0].X, spiral[0].Y, spiral[0].Z)
                t.pen_down()
                for p in spiral:
                    t.set_position(p.X, p.Y, p.Z)
            elif node['type'] == 2:
                points = rs.DivideCurve(node["curves"][0], num_pnts)
                t.pen_up()
                travel_paths.append(rs.AddCurve([t.get_position(), spiral[0]]))
                t.set_position(points[0].X, points[0].Y, t.get_position().Z)
                t.set_position(points[0].X, points[0].Y, points[0].Z)
                t.pen_down()
                for p in points:
                    t.set_position(p.X, p.Y, p.Z)

    return travel_paths


def fill_curves_with_contours(t, curves, start_pnt=None, wall_mode=False, walls=3, initial_offset=0.5):
    offset = float(t.get_extrude_width())

    # connect curves if given more than one
    curve = curves[0]
    if len(curves) > 1:
        curve = connect_curves(curves, offset/8)

    # slice the shape
    # Generate isocontours
    root, isocontours = get_contours(t, curve, walls=walls, wall_mode=wall_mode, initial_offset=initial_offset)

    isocontours = [rs.DivideCurve(i, int(rs.CurveLength(i)/t.get_resolution())) for i in isocontours]

    travel_paths = []
    start_idx = 0
    for i in range(len(isocontours)):
        points = isocontours[i]
        start_idx, d = closest_point(start_pnt, points)
        if start_idx == None: start_idx = 0

        travel_paths = travel_paths + [rs.AddCurve([t.get_position(), points[start_idx]])]
        t.pen_up()
        t.set_position(points[start_idx].X, points[start_idx].Y, t.get_position().Z)
        t.set_position(points[start_idx].X, points[start_idx].Y, points[start_idx].Z)
        t.pen_down()

        for p in (range(start_idx+1, len(points)) + range(0, start_idx)):
            t.set_position(points[p].X, points[p].Y, points[p].Z)
        t.set_position(points[start_idx].X, points[start_idx].Y, points[start_idx].Z)

    return travel_paths


def slice_fermat_fill(t, shape, start_pnt=None, start=0, end=None, wall_mode=False, walls=3, fill_bottom=False, bottom_layers=3, initial_offset=0.5):
    travel_paths = []
    layers = int(math.floor(get_shape_height(shape) / t.get_layer_height())) + 1

    print("Number of layers: "+str(layers-1))

    if end is None: end = layers

    for l in range(start, min(layers, end)):
        curve_groups = get_curves(shape, l*t.get_layer_height())

        for crvs in curve_groups:
            if start_pnt == None: start_pnt = t.get_position()
            try:
                if not wall_mode or (wall_mode and fill_bottom and l<bottom_layers):
                    travel_paths = travel_paths + fill_curves_with_fermat_spiral(t, crvs, start_pnt=start_pnt, initial_offset=initial_offset)[0]
                else:
                    travel_paths = travel_paths + fill_curves_with_fermat_spiral(t, crvs, start_pnt=start_pnt, wall_mode=wall_mode, walls=walls, initial_offset=initial_offset)[0]
            except:
                try:
                    print("Unable to fermat spiral layer, generating contours: "+str(l))
                    if not wall_mode or (wall_mode and fill_bottom and l<bottom_layers):
                        travel_paths = travel_paths + fill_curves_with_contours(t, crvs, start_pnt=start_pnt, initial_offset=initial_offset)
                    else:
                        travel_paths = travel_paths + fill_curves_with_contours(t, crvs, start_pnt=start_pnt, wall_mode=wall_mode, walls=walls, initial_offset=initial_offset)
                except:
                    print("Error: unable to slice layer "+str(l))

    return travel_paths


def slice_spiral_fill(t, shape, start=0, end=None):
    travel_paths = []
    layers = int(math.floor(get_shape_height(shape) / t.get_layer_height())) + 1

    print("Number of layers: "+str(layers-1))

    if end is None: end = layers

    for l in range(start, min(layers, end)):
        curve_groups = get_curves(shape, l*t.get_layer_height())

        for crvs in curve_groups:
            travel_paths = travel_paths + fill_curves_with_spiral(t, crvs, start_pnt=t.get_position())

    return travel_paths


def slice_contour_fill(t, shape, start=0, end=None, wall_mode=False, walls=3, fill_bottom=False, bottom_layers=3, initial_offset=0.5):
    travel_paths = []
    layers = int(math.floor(get_shape_height(shape) / t.get_layer_height())) + 1

    print("Number of layers: "+str(layers-1))

    if end is None: end = layers

    for l in range(start, min(layers, end)):
        curve_groups = get_curves(shape, l*t.get_layer_height())

        for crvs in curve_groups:
            if not wall_mode or (wall_mode and fill_bottom and l<bottom_layers):
                travel_paths = travel_paths + fill_curves_with_contours(t, crvs, start_pnt=t.get_position(), initial_offset=initial_offset)
            else:
                travel_paths = travel_paths + fill_curves_with_contours(t, crvs, start_pnt=t.get_position(), wall_mode=wall_mode, walls=walls, initial_offset=initial_offset)

    return travel_paths


def slice_vertical_and_fermat_fill(t, shape, wall_mode=False, walls=3, fill_bottom=False, bottom_layers=3, initial_offset=0.5):
    overall_start_time = time.time()

    travel_paths = []
    center_points = []
    #try:
    tree, path, center_points, bboxes, edges = best_vertical_path(t, shape)

    start_point = path[0].data.sub_nodes[0].start_point
    for sup_node in path:
        for node in sup_node.data.sub_nodes:
            #if node == sup_node.data.sub_nodes[0]: start_point = node.start_point
            #else: start_point = t.get_position()
            start_point = t.get_position()
            for curves in node.data:
                try:
                    if not wall_mode or (wall_mode and fill_bottom and node.height<bottom_layers):
                        travel_paths = travel_paths + fill_curves_with_fermat_spiral(t, curves, start_pnt=start_point, initial_offset=initial_offset)[0]
                    else:
                        travel_paths = travel_paths + fill_curves_with_fermat_spiral(t, curves, start_pnt=start_point, wall_mode=wall_mode, walls=walls, initial_offset=initial_offset)[0]
                except Exception as err:
                    print("Failed to fermat spiral layer "+str(node.height), err)
                    if not wall_mode or (wall_mode and fill_bottom and node.height<bottom_layers):
                        travel_paths = travel_paths + fill_curves_with_contours(t, curves, start_pnt=start_point, initial_offset=initial_offset)
                    else:
                        travel_paths = travel_paths + fill_curves_with_contours(t, curves, start_pnt=start_point, wall_mode=wall_mode, walls=walls, initial_offset=initial_offset)
    #except Exception as error:
        #print("Failed to find vertical travel path minimization, printing layer by layer. "+str(error))
        #travel_paths = travel_paths + slice_fermat_fill(t, shape, wall_mode=wall_mode, walls=walls, fill_bottom=fill_bottom, bottom_layers=bottom_layers)

    print("Full path generation: "+str(time.time()-overall_start_time)+" seconds")

    return travel_paths, tree, path, center_points, bboxes, edges


def slice_2_half_D_fermat(t, curves, layers=3, wall_mode=False, walls=3, fill_bottom=False, bottom_layers=3, initial_offset=0.5):
    groups = get_curve_groupings(curves)
    group = groups[0]

    curve = group[0]
    if len(group) > 1:
        curve = connect_curves(curves, float(t.get_extrude_width())/8)

    first_curve = sorted(get_isocontour(curve, float(t.get_extrude_width())*initial_offset), key=lambda x: get_area(x), reverse=True)[0]
    root = {"guid": first_curve, "depth": 0, "children":[]}
    isocontours = [] + [first_curve]
    new_curves = get_isocontours(t, first_curve, root, wall_mode=wall_mode, walls=walls)
    if new_curves:
        isocontours = isocontours + new_curves

    travel_paths = []

    region_tree = segment_tree(root)
    set_node_types(region_tree)
    all_nodes = get_all_nodes(region_tree)

    for l in range(layers):
        for n in all_nodes:
            if len(n.get('curves')) > 0:
                num_pnts = get_num_points(n['curves'][0], float(t.get_extrude_width()))
                if n['type'] == 1:
                    start_idx, d = closest_point(t.get_position(), rs.DivideCurve(n['curves'][0], num_pnts))
                    spiral, indices = spiral_contours(t, n["curves"], start_idx)
                    n["fermat_spiral"] = fermat_spiral(t, spiral, indices)
                elif n['type'] == 2:
                    n['fermat_spiral'] = rs.DivideCurve(n["curves"][0], num_pnts)
            else:
                print("Error: node with no curves in it at all", n)

        final_points = connect_spiralled_nodes(t, region_tree)
        final_curve = rs.AddCurve(final_points)
        final_spiral = rs.DivideCurve(final_curve, int(rs.CurveLength(final_curve)/t.get_resolution()))

        t.pen_up()
        #travel_paths.append(rs.AddCurve([t.get_position(), final_spiral[0]]))
        t.set_position(final_spiral[0].X, final_spiral[0].Y, t.get_layer_height()*l)
        t.pen_down()
        for p in final_spiral:
            t.set_position(p.X, p.Y, t.get_layer_height()*l)

    return travel_paths
