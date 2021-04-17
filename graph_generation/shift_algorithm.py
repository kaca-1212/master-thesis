import networkx as nx
from collections import defaultdict
from instances_preliminaries import generate_triangulated_graph, get_canonical_ordering


def combinatorial_embedding_to_pos(g):
    """Assigns every node a (x, y) position based on the given embedding

    The algorithm iteratively inserts nodes of the input graph in a certain
    order and rearranges previously inserted nodes so that the planar drawing
    stays valid. This is done efficiently by only maintaining relative
    positions during the node placements and calculating the absolute positions
    at the end. For more information see [1]_.

    Parameters
    ----------
    embedding : nx.PlanarEmbedding
        This defines the order of the edges

    fully_triangulate : bool
        If set to True the algorithm adds edges to a copy of the input
        embedding and makes it chordal.

    Returns
    -------
    pos : dict
        Maps each node to a tuple that defines the (x, y) position

    References
    ----------
    .. [1] M. Chrobak and T.H. Payne:
        A Linear-time Algorithm for Drawing a Planar Graph on a Grid 1989
        http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.51.6677

    """
    # The following dicts map a node to another node
    # If a node is not in the key set it means that the node is not yet in G_k
    # If a node maps to None then the corresponding subtree does not exist
    left_t_child = {}
    right_t_child = {}

    # The following dicts map a node to an integer
    delta_x = {}
    y_coordinate = {}

    node_list = get_canonical_ordering(g)

    # 1. Phase: Compute relative positions

    # Initialization
    v1, v2, v3 = node_list[0][0], node_list[1][0], node_list[2][0]

    delta_x[v1] = 0
    y_coordinate[v1] = 0
    right_t_child[v1] = v3
    left_t_child[v1] = None

    delta_x[v2] = 1
    y_coordinate[v2] = 0
    right_t_child[v2] = None
    left_t_child[v2] = None

    delta_x[v3] = 1
    y_coordinate[v3] = 1
    right_t_child[v3] = v2
    left_t_child[v3] = None

    for k in range(3, len(node_list)):
        vk, contour_neighbors = node_list[k]
        wp = contour_neighbors[0]
        wp1 = contour_neighbors[1]
        wq = contour_neighbors[-1]
        wq1 = contour_neighbors[-2]
        adds_mult_tri = len(contour_neighbors) > 2

        # Stretch gaps:
        delta_x[wp1] += 1
        delta_x[wq] += 1

        delta_x_wp_wq = sum(delta_x[x] for x in contour_neighbors[1:])

        # Adjust offsets
        delta_x[vk] = (-y_coordinate[wp] + delta_x_wp_wq + y_coordinate[wq]) // 2
        y_coordinate[vk] = (y_coordinate[wp] + delta_x_wp_wq + y_coordinate[wq]) // 2
        delta_x[wq] = delta_x_wp_wq - delta_x[vk]
        if adds_mult_tri:
            delta_x[wp1] -= delta_x[vk]

        # Install v_k:
        right_t_child[wp] = vk
        right_t_child[vk] = wq
        if adds_mult_tri:
            left_t_child[vk] = wp1
            right_t_child[wq1] = None
        else:
            left_t_child[vk] = None

    # 2. Phase: Set absolute positions
    pos = dict()
    pos[v1] = (0, y_coordinate[v1])
    remaining_nodes = [v1]
    while remaining_nodes:
        parent_node = remaining_nodes.pop()

        # Calculate position for left child
        set_position(
            parent_node, left_t_child, remaining_nodes, delta_x, y_coordinate, pos
        )
        # Calculate position for right child
        set_position(
            parent_node, right_t_child, remaining_nodes, delta_x, y_coordinate, pos
        )
    return pos



def set_position(parent, tree, remaining_nodes, delta_x, y_coordinate, pos):
    """Helper method to calculate the absolute position of nodes."""
    child = tree[parent]
    parent_node_x = pos[parent][0]
    if child is not None:
        # Calculate pos of child
        child_x = parent_node_x + delta_x[child]
        pos[child] = (child_x, y_coordinate[child])
        # Remember to calculate pos of its children
        remaining_nodes.append(child)
