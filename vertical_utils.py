import math
import Rhino
import rhinoscriptsyntax as rs
import extruder_turtle  
import turtle_utilities as tu
from extruder_turtle import *

import tree_utils
from tree_utils import *

import graph_utils
from graph_utils import *

import geometry_utils
from geometry_utils import *

def test_graph():
    test = Graph()
    a = Graph_Node("A")
    b = Graph_Node("B")
    c = Graph_Node("C")
    d = Graph_Node("D")
    test.add_node(a)
    test.add_node(b)
    test.add_node(c)
    test.add_node(d)
    test.add_edge(Graph_Edge(a, b, 3))
    test.add_edge(Graph_Edge(a, d, 5))
    test.add_edge(Graph_Edge(a, c, 3))
    test.add_edge(Graph_Edge(b, a, 1))
    test.add_edge(Graph_Edge(b, c, 2))
    test.add_edge(Graph_Edge(b, d, 1))
    test.add_edge(Graph_Edge(c, a, 1))
    test.add_edge(Graph_Edge(c, b, 3))
    test.add_edge(Graph_Edge(c, d, 2))
    test.add_edge(Graph_Edge(d, a, 5))
    test.add_edge(Graph_Edge(d, b, 4))
    test.add_edge(Graph_Edge(d, c, 0))

    test.starts.append(a)
    path = test.get_shortest_hamiltonian_path()
    print([node.data for node in path[0]], path[1])

def best_vertical_path(t, shape):
    vert_tree = build_vertical_tree(t, shape)
    all_nodes = vert_tree.get_all_nodes([])

    nozzle_height = 30 #height of nozzle in mm
    height = get_shape_height(shape)
    for h in range(int(math.floor(height / nozzle_height))+1):
        nodes_at_height = [node for node in all_nodes if node.height == h]

        # create a graph for this height chunk
        height_graph = Graph()
        # add nodes to graph for every super node within the height chunk
        for node in nodes_at_height:
            graph_node = Graph_Node(node)
            height_graph.add_node(graph_node)
            # add start nodes to the graph
            if node.sub_nodes[0].height % nozzle_height == 0:
                height_graph.starts.append(graph_node)

        # add edges to graph
        # edges related to height dependency
        for graph_node in height_graph.nodes:
            node1 = graph_node.data
            for child in node1.children:
                if child in nodes_at_height:
                    height_graph.add_edge(Graph_Edge(graph_node, height_graph.get_node(child), 0))

            # edges related to travel between nodes
            siblings_and_counsins = [n for n in nodes_at_height if n not in node1.get_all_descendants([]) + node1.get_all_ancestors([])]
            for node2 in siblings_and_counsins:
                # do not add edge if there is overlap
                if not is_overlapping(node1, node2):
                    # compute travel between start points
                    height_graph.add_edge(Graph_Edge(graph_node, height_graph.get_node(node2), 0))

        num_edges = 0
        for n in height_graph.edges:
            num_edges = num_edges + len(height_graph.edges[n].keys())
        print(len(height_graph.nodes), num_edges)
        #print('shortest path', height_graph.get_shortest_hamiltonian_path())

    return vert_tree

def build_vertical_tree(t, shape):
    layers = int(math.floor(get_shape_height(shape) / t.get_layer_height())) + 1
    root = Node('root')

    previous_nodes = [root]
    for l in range(layers):
        z = l*t.get_layer_height()
        plane = get_plane(z)
        curves = rs.AddSrfContourCrvs(shape, plane)
        new_nodes = []
        for curve in curves:
            node = Node(curve)
            node.depth = l
            node.height = l
            new_nodes.append(node)
            if root in previous_nodes:
                node.parents.append(root)
                root.children.append(node)
            else:
                for prev_n in previous_nodes:
                    if xy_bbox_overlap(prev_n.data, curve):
                        node.parents.append(prev_n)
                        prev_n.children.append(node)
            if len(node.parents) == 0: node.needs_support = True
            else: node.needs_support = False

        previous_nodes = new_nodes

    return segment_tree_by_height(t, root, get_shape_height(shape))

def segment_tree_by_height(t, tree, total_height):
    #nozzle_height = t.get_nozzle_height()
    #nozzle_width = t.get_nozzle_max_width()
    nozzle_height = 30 #height of nozzle in mm
    limit = int(math.floor(nozzle_height / t.get_layer_height()))
    super_root = Node('root')
    super_root.depth = 0
    super_root.height = 0
    idx = 0
    for child in tree.children:
        group_by_height(child, super_root, limit, idx)
        idx = idx + 1

    divide_by_overlap(super_root, total_height)
    return super_root


def group_by_height(node, super_node, height, idx=0):
    s_node = super_node
    if node.depth // height == super_node.height:
        super_node.sub_nodes.append(node)
    elif node.depth // height > super_node.height:
        new_super = Node(str(super_node.data)+'_'+str(idx))
        new_super.parent = super_node
        new_super.depth = super_node.depth + 1
        new_super.height = node.depth // height
        new_super.sub_nodes.append(node)

        super_node.children.append(new_super)

        s_node = new_super
    elif node.depth // height < super_node.height:
        raise ValueError("Error, node should not be below current super_node")

    idx = 0
    if len(node.children) > 1:
        for child in node.children:
            new_new_super = Node(str(s_node.data)+'_'+str(idx))
            new_new_super.parent = s_node
            new_new_super.depth = s_node.depth + 1
            new_new_super.height = node.depth // height
            s_node.children.append(new_new_super)
            group_by_height(child, new_new_super, height, 0)
            idx = idx + 1
    else:
        for child in node.children:
            group_by_height(child, s_node, height, idx)
            idx = idx + 1


def divide_by_overlap(super_root, total_height):
    #nozzle_height = t.get_nozzle_height()
    #nozzle_width = t.get_nozzle_max_width()
    nozzle_width = 8
    nozzle_height = 30 #height of nozzle in mm
    height = int(math.floor(total_height / nozzle_height)) + 1
    nodes = super_root.get_all_nodes([])
    for h in range(height):
        nodes_at_height = [node for node in nodes if node.height == h]
        subdivide_by_overlap(nodes_at_height, nozzle_width)


def subdivide_by_overlap(nodes, width):
    for n1 in range(len(nodes)):
        # check each sub-layer within the node to see if it overlaps with other
        # nodes' layers, if those nodes are siblings or cousins of the node
        for n2 in range(n1+1, len(nodes)):
            node1 = nodes[n1]
            node2 = nodes[n2]
            other_nodes = node1.get_all_ancestors([]) + node1.get_all_descendants([])

            # if node2 within height chunk is a sibling or cousin
            overlap_above = False
            overlap_below = False
            if node2 not in other_nodes:
                for s1 in node1.sub_nodes:
                    for s2 in node2.sub_nodes:
                        if xy_bbox_overlap(s1.data, s2.data, width):
                            if s2 not in s1.overlap:
                                s1.overlap.append(s2)
                            if s1 not in s2.overlap:
                                s2.overlap.append(s1)

    for n1 in range(len(nodes)):
        node1 = nodes[n1]

        splits = []
        for s1 in range(1, len(node1.sub_nodes)):
            sub0 = node1.sub_nodes[s1-1]
            sub1 = node1.sub_nodes[s1]
            if set(sub0.overlap) != set(sub1.overlap):
                splits.append(sub0.height)

        for split in splits:
            split_super_node_at_height(node1, split)


def split_super_node_at_height(node, height):
    split_node = Node(node.data+'_split_'+str(height))
    split_node.depth = node.depth
    split_node.height = node.height
    split_node.children.append(node)
    split_node.parents = [p for p in node.parents]
    for p in split_node.parents:
        p.children.append(split_node)
        p.children.remove(node)
    split_node.sub_nodes = [n for n in node.sub_nodes if n.height <= height]

    node.depth = node.depth + 1
    node.parents = [split_node]
    node.sub_nodes = [n for n in node.sub_nodes if n.height > height]

    descendants = node.get_all_descendants([])
    for d in descendants:
        d.depth = d.depth + 1


def is_overlapping(node1, node2):
    for sub1 in node1.sub_nodes:
        for sub2 in node2.sub_nodes:
            if sub1 in sub2.overlap and sub1.height > sub2.height:
                return True
    return False