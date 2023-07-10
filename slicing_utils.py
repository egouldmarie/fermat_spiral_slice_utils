import math
import Rhino
import rhinoscriptsyntax as rs
import extruder_turtle  
import turtle_utilities as tu
from extruder_turtle import *


def get_size(shape):
    # bounding box of shape
    bb = rs.BoundingBox(shape)
    size = rs.Distance(bb[0], bb[2])
    return size


def get_shape_height(shape):
    # bounding box of shape
    bb = rs.BoundingBox(shape)
    height = rs.Distance(bb[0], bb[4])
    return height


def get_plane(z):
    return rs.PlaneFromNormal([0, 0, z], [0, 0, 1])


def get_surface(curve, z):
    size = get_size(curve)*2

    points = []
    points.append(rs.CreatePoint(-size,-size,z))
    points.append(rs.CreatePoint(-size,size,z))
    points.append(rs.CreatePoint(size,size,z))
    points.append(rs.CreatePoint(size,-size,z))
    surface = rs.AddSrfPt(points)
    return surface


def get_precision(t, curve):
    resolution = t.get_resolution()
    points = rs.DivideCurve (curve, 100)
    ll = tu.line_length(points)
    num_points = int(ll/resolution)+1
    return num_points


def closest_point(point, points):
    closest = {"point": None, "distance": 1000000}
    for p in range(len(points)):
        dist = rs.Distance(point, points[p])
        if dist < closest['distance']:
            closest['distance'] = dist
            closest['point'] = p
    
    return closest['point']


def get_corner(t, outer_curve, inner_curve, points):
    offset = t.get_extrude_width()
    center = rs.CurveAreaCentroid(inner_curve)[0]
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
                inside = rs.PointInPlanarClosedCurve(rs.VectorAdd(pnt, rs.VectorScale(rs.VectorUnitize(add), offset/10)), outer_curve)
                if dist < closest["distance"] and inside:
                    closest["distance"] = dist
                    closest["point"] = p

    return closest["point"]


def get_isocontours(t, curves, depth, precision, wall_mode=False, walls=3):
    new_curves = get_isocontour(t, curves, precision)
    if not new_curves:
        return []
    else:
        all_curves = [] + new_curves
        new_depth = depth+1
        if not wall_mode or (wall_mode and new_depth <= walls):
            new_curves = get_curve_groupings(new_curves)
            for cs in new_curves:
                new_new_curves = get_isocontours(t, cs, new_depth, precision, wall_mode=wall_mode, walls=walls)
                for nc in new_new_curves:
                    all_curves.append(nc)
        return all_curves


def get_isocontour(t, curves, precision):
    if not (curves and all([rs.IsCurveClosed(curve) for curve in curves])):
        print("Error: get_isocontours called with not a closed curve: ", curves)
        return None

    if precision <= 2:
        print("Error: precision too low")
        return None

    points = []
    for curve in curves:
        points = points + rs.DivideCurve(curve, precision)

    offset = t.get_extrude_width()

    winding_order, direction = get_winding_order(t, curves[0], rs.DivideCurve(curves[0], precision))
    if winding_order == None:
        print("Error: winding order could not return a direction")
        return None

    # determine each new p' at distance offset away from p
    new_points = [None]*len(points)
    discarded_points = [None]*len(points)
    for i in range(len(points)):
        prev_i = (i-1) % len(points)
        next_i = (i+1) % len(points)

        # get tangent vector
        tangent = rs.VectorSubtract(points[next_i], points[prev_i])
        # get vector orthogonal to tangent vector
        ortho = rs.VectorRotate(tangent, direction, [0, 0, 1])
        # normalize and scale orthogonal vector
        ortho = rs.VectorScale(ortho, offset/rs.VectorLength(ortho))
        # compute new point
        new_point = rs.VectorAdd(points[i], ortho)

        # check that distance from all points is >= offset
        include = True
        for j in range(len(points)):
            if not i == j and rs.Distance(points[j], new_point) < offset:
                include = False
                break
        if include:
            new_points[i] = new_point
        else:
            discarded_points[i] = new_point

    # if there are new points
    if not all(x is None for x in new_points):
        # if any points have been discarded
        if not all(x is None for x in discarded_points):
            # get list of lists of all sequential indices
            sequences = [[]]
            start_index = next((index for index, value in enumerate(new_points) if value != None and new_points[index-1] == None), -1)
            if start_index !=- 1:
                indices = range(start_index, len(new_points)) + range(0, start_index)
                for i in indices:
                    next_i = (i+1)%len(new_points)
                    if new_points[i] != None:
                        sequences[-1].append(i)
                    elif new_points[next_i] != None and next_i != start_index:
                        sequences.append([])

            # remove sequences that consist of a single point
            for seq in sequences:
                if len(seq) <= 1:
                    sequences.remove(seq)

            # get start and end points of all sequences
            start = [seq[0] for seq in sequences]
            end = [seq[-1] for seq in sequences]

            # get connections between sequences
            # prime connections dictionary with None
            # index and large initial minimum distance
            connections = {j: (None, 100000) for j in end}
            for i in start:
                for j in end:
                    dist = rs.Distance(new_points[i], new_points[j])
                    if dist < connections[j][1]:
                        connections[j] = (i, dist)

            # get rid of distances now that
            # connections have been made
            connections = {j: connections[j][0] for j in connections}

            # get sequence grouping indices of connections
            nodes = [[i] for i in range(len(sequences))]

            # check through groupings until no more
            # connections are found between sequences
            connections_found = True
            while connections_found:
                connections_found = False
                for node in nodes:
                    idx = nodes.index(node)
                    seq = sequences[node[-1]]
                    connect = connections[seq[-1]]
                    for node2 in nodes:
                        seq2 = sequences[node2[0]]
                        idx2 = nodes.index(node2)
                        if idx != idx2 and connect in seq2:
                            nodes[idx] = nodes[idx] + nodes[idx2]
                            nodes.pop(idx2)
                            connections_found = True

            # construct curves from groupings of indices
            curves = []
            for node in nodes:
                curves.append([])
                for idx in node:
                    curves[-1] = curves[-1] + [new_points[i] for i in sequences[idx]]

            # make sure curve is closed;
            # add start point to the end
            curves = [c+[c[0]] for c in curves]

            # Transform point lists into curves
            curves = [rs.AddCurve(c) for c in curves if len(c) > 5]

            return curves
        else:
            curves = []
            for c in range(0, len(new_points), precision):
                start = c
                end = c + precision
                curve_points = [new_points[i] for i in range(start, end)]+[new_points[start]]
                if not any(x is None for x in curve_points): curves.append(rs.AddCurve(curve_points))
                else: print("Error: All new points are none")
            return curves
    else:
        return None


def get_winding_order(t, curve, points):
    # get winding order, CW or CCW
    offset = t.get_extrude_width()
    winding_order = None

    for i in range(4):
        index = i*(len(points)/4)

        tangent = rs.VectorSubtract(points[index+1], points[index-1])
        pnt_cw = rs.VectorAdd(points[index], rs.VectorScale(rs.VectorUnitize(rs.VectorRotate(tangent, -90, [0, 0, 1])), offset/10))
        pnt_ccw = rs.VectorAdd(points[index], rs.VectorScale(rs.VectorUnitize(rs.VectorRotate(tangent, 90, [0, 0, 1])), offset/10))

        direction = None
        if rs.PointInPlanarClosedCurve(pnt_cw, curve):
            winding_order = "CW"
            direction = -90
            break
        elif rs.PointInPlanarClosedCurve(pnt_ccw, curve):
            winding_order = "CCW"
            direction = 90
            break

    return winding_order, direction


def spiral_contours(t, isocontours, precision, start_index):
    if len(isocontours) == 0:
        print("Error: no isocontours passed in")
        return
    # a single spirallable region
    # connect each isocontour with the one succeeding it
    spiral = []
    points = rs.DivideCurve(isocontours[0], precision)
    offset = t.get_extrude_width()

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
            next_points = rs.DivideCurve(isocontours[i+1], precision)
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
    offset = t.get_extrude_width()

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
        new_node = {"guid":region_node["guid"]+"_"+str(idx), "curves":[], "parent":region_node, "children":[]}
        new_node["curves"].append(node["guid"])

        region_node["children"].append(new_node)

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


def get_curve_groupings(curves):
    # find curve groupings from intersection of shape with plane
    # curves can represent the inside of a surface or potentially
    # a nested curve within another set of curves defining a surface
    all_points = [rs.DivideCurve(curve, 100) for curve in curves]
    inside = {c:{c2:all([rs.PointInPlanarClosedCurve(p, curves[c2]) for p in all_points[c]]) for c2 in range(len(curves)) if c2 != c} for c in range(len(curves))}
    outer_curves = [c for c in range(len(curves)) if not any([inside[c][k] for k in inside[c]])]
    inner_curves = [c for c in range(len(curves)) if c not in outer_curves]

    iteration = 0
    curve_groupings = {c:[] for c in outer_curves}
    while len(inner_curves) > 0:
        if iteration > 30: break
        iteration = iteration+1
        next_group = []
        for c in outer_curves:
            for c2 in inner_curves:
                if inside[c2][c] and all([inside[c2][c3] == inside[c][c3] for c3 in inside[c2] if c3 != c]):
                    curve_groupings[c].append(c2)
                    next_group.append(c2)
                    inner_curves.remove(c2)

        outer_curves = []
        for c in next_group:
            for c2 in inner_curves:
                if inside[c2][c] and all([inside[c2][c3] == inside[c][c3] for c3 in inside[c2] if c3 != c]):
                    outer_curves.append(c2)
                    curve_groupings[c2] = []
                    inner_curves.remove(c2)

    groups = [[c]+curve_groupings[c] for c in curve_groupings]
    curves = [[curves[c] for c in g] for g in groups]
    return curves


def tree_isocontours(t, root, isocontours):
    # find outermost contours
    points = [rs.DivideCurve(curve, 100) for curve in isocontours]
    inside = {c:{c2:all([rs.PointInPlanarClosedCurve(p, isocontours[c]) for p in points[c2]]) for c2 in range(len(isocontours)) if c2 != c} for c in range(len(isocontours))}
    #for i in inside: print(i, inside[i])

    node = {"guid":isocontours[0], "depth":root["depth"]+1, "children":[]}
    root["children"].append(node)

    has_parent = [False for i in range(len(isocontours))]
    has_parent[0] = True
    find_children(t, node, isocontours, points, inside, has_parent)
    print(has_parent)

    return root


def directly_inside(out_, in_, inside):
    if out_ != in_:
        if in_==1:
            print('')
            print(out_, in_)
            print(inside[out_])
            print(inside[in_][out_])
            print([inside[idx3][in_] for idx3 in inside[out_] if idx3 != in_])
        return inside[out_][in_] or (not inside[in_][out_] and not any([inside[idx3][in_] for idx3 in inside[out_] if idx3 != in_]))
    return False


def find_children(t, parent, isocontours, points, inside, has_parent):
    offset = t.get_extrude_width()
    pidx = isocontours.index(parent['guid'])
    for i in range(len(isocontours)):
        if i==1: print(i, pidx, directly_inside(pidx, i, inside))
        if not has_parent[i] and directly_inside(pidx, i, inside) and any([rs.Distance(p, p2) <= offset*1.5 for p in points[pidx] for p2 in points[i]]):
            node = {"guid":isocontours[i], "depth":parent["depth"]+1, "children":[]}
            parent['children'].append(node)
            has_parent[i] = True
            find_children(t, node, isocontours, points, inside, has_parent)


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
    if parent:
        connect_node_to_parent(t, node, parent)


def connect_node_to_parent(t, node, parent):
    offset = t.get_extrude_width()

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
    offset = t.get_extrude_width()

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
                connect_indices = []
                if c[0] > c[1]: connect_indices = range(c[0], c[1]-1, -1)
                else:  connect_indices = range(c[0], c[1]+1)
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


def fill_layer_with_fermat_spiral(t, shape, z, start_pnt=None, wall_mode=False, walls=3):
    # get outer curves of shape
    plane = get_plane(z)
    
    # if there is more than one curve when intersecting
    # with the plane, we will have travel paths between them
    curves = rs.AddSrfContourCrvs(shape, plane)
    curves = get_curve_groupings(curves)

    # slice the shape
    print("Generating Isocontours")
    precision = 300
    isocontours = []
    for curve_group in curves:
        isocontours.append(curve_group)
        new_curves = get_isocontours(t, curve_group, 1, precision, wall_mode, walls)
        if new_curves:
            isocontours[-1] = isocontours[-1] + new_curves

    print("Treeing Isocontours")
    root = {"guid": "root", "depth": 0, "children":[]}
    for i in isocontours:
        tree_isocontours(t, root, i)

    print(len(root['children']), len(get_all_nodes(root))-1)

    return root

    region_tree = segment_tree(root)
    set_node_types(region_tree)
    all_nodes = get_all_nodes(region_tree)

    print("Spiralling Regions")
    for node in all_nodes:
        if len(node['curves']) > 0:
            if node['type'] == 1:
                start_idx = None
                if 'root' in node['curves']:
                    node['curves'].remove('root')
                if start_pnt:
                    print(node['curves'][0])
                    start_idx = closest_point(start_pnt, rs.DivideCurve(node['curves'][0], precision))
                spiral, indices = spiral_contours(t, node["curves"], precision, start_idx)
                node["fermat_spiral"] = fermat_spiral(t, spiral, indices)
            elif node['type'] == 2:
                node['fermat_spiral'] = rs.DivideCurve(node["curves"][0], precision)

    print("Connecting Spiralled Regions")
    travel_paths = []
    final_spiral = connect_spiralled_nodes(t, region_tree)
    t.pen_up()
    travel_paths.append(rs.AddCurve([t.get_position(), final_spiral[0]]))
    t.set_position(final_spiral[0].X, final_spiral[0].Y, final_spiral[0].Z)
    t.pen_down()
    for p in final_spiral:
        t.set_position(p.X, p.Y, p.Z)

    all_nodes = {node['guid']: node for node in all_nodes}

    return travel_paths


def slice_fermat_fill(t, shape, wall_mode=False, walls=3, fill_bottom=False, bottom_layers=3):
    travel_paths = []
    layers = int(math.floor(get_shape_height(shape) / t.get_layer_height())) + 1
    for l in range(layers):
        if not wall_mode or (wall_mode and fill_bottom and l+1>bottom_layers):
            print("Slicing Layer "+str(l))
            travel_paths = travel_paths + fill_layer_with_fermat_spiral(t, shape, l*t.get_layer_height(), start_pnt=t.get_position(), wall_mode=wall_mode, walls=walls)
        else:
            travel_paths = travel_paths + fill_layer_with_fermat_spiral(t, shape, l*t.get_layer_height(), start_pnt=t.get_position())
 
    return travel_paths


def fill_layer_with_spiral(t, shape, z, start_pnt):
    # get outer curves of shape
    plane = get_plane(z)

    # if there is more than one curve when intersecting
    # with the plane, we will have travel paths between them
    curves = rs.AddSrfContourCrvs(shape, plane)
    curves = get_curve_groupings(curves)

    # slice the shape
    print("Generating Isocontours")
    precision = 300
    isocontours = []
    for curve_group in curves:
        isocontours = isocontours + curve_group
        new_curves = get_isocontours(t, curve_group, 1, precision)
        if new_curves:
            isocontours = isocontours + new_curves

    print("Treeing Isocontours")
    root = tree_isocontours({"guid": "root", "depth": 0, "children":[]}, isocontours)

    region_tree = segment_tree(root)
    set_node_types(region_tree)
    all_nodes = get_all_nodes(region_tree)

    print("Spiralling Regions")
    travel_paths = []
    for node in all_nodes:
        if node['type'] == 1:
            start_idx = None
            if 'root' in node['curves']:
                node['curves'].remove('root')
                if start_pnt:
                    start_idx = closest_point(start_pnt, rs.DivideCurve(node['curves'][0], precision))
            spiral, indices = spiral_contours(t, node["curves"], precision, start_idx)
            t.pen_up()
            travel_paths.append(rs.AddCurve([t.get_position(), spiral[0]]))
            t.set_position(spiral[0].X, spiral[0].Y, spiral[0].Z)
            t.pen_down()
            for p in spiral:
                t.set_position(p.X, p.Y, p.Z)
        elif node['type'] == 2:
            points = rs.DivideCurve(node["curves"][0], precision)
            t.pen_up()
            travel_paths.append(rs.AddCurve([t.get_position(), spiral[0]]))
            t.set_position(points[0].X, points[0].Y, points[0].Z)
            t.pen_down()
            for p in points:
                t.set_position(p.X, p.Y, p.Z)

    return travel_paths


def slice_spiral_fill(t, shape):
    travel_paths = []
    layers = int(math.floor(get_shape_height(shape) / t.get_layer_height())) + 1
    for l in range(layers):
        print("Slicing Layer "+str(l))
        travel_paths = travel_paths + fill_layer_with_spiral(t, shape, l*t.get_layer_height(), t.get_position())

    return travel_paths


def fill_layer_with_contours(t, shape, z):
    # get outer curves of shape
    plane = get_plane(z)
    
    # if there is more than one curve when intersecting
    # with the plane, we will have travel paths between them
    curves = rs.AddSrfContourCrvs(shape, plane)
    curves = get_curve_groupings(curves)

    # slice the shape
    print("Generating Isocontours")
    precision = 300
    isocontours = []
    for curve_group in curves:
        isocontours = isocontours + curve_group
        new_curves = get_isocontours(t, curve_group, 1, precision)
        if new_curves:
            isocontours = isocontours + new_curves

    isocontours = [rs.DivideCurve(i, precision) for i in isocontours]

    travel_paths = []
    start_idx = 0
    for i in range(len(isocontours)):
        start = isocontours[i][start_idx]
        t.pen_up()
        travel_paths.append(rs.AddCurve([t.get_position(), start]))
        t.set_position(start.X, start.Y, start.Z)
        t.pen_down()

        points = [isocontours[i][idx] for idx in range(start_idx, len(isocontours[i])) + range(0, start_idx)]
        for p in points:
                t.set_position(p.X, p.Y, p.Z)

        if i<len(isocontours)-1:
            start_idx = closest_point(start, isocontours[i+1])

    return travel_paths


def slice_contour_fill(t, shape):
    travel_paths = []
    layers = int(math.floor(get_shape_height(shape) / t.get_layer_height())) + 1
    for l in range(layers):
        print("Slicing Layer "+str(l))
        travel_paths = travel_paths + fill_layer_with_contours(t, shape, l*t.get_layer_height())

    return travel_paths
