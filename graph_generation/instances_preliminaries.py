from collections import defaultdict
import networkx as nx
from random import randint

from errors import NumberOfNodesError

def generate_triangulated_graph(number_of_nodes):
    if number_of_nodes < 3 or number_of_nodes is None:
        raise NumberOfNodesError("number of nodes must be greater than 3")

    G = nx.Graph()
    G.add_nodes_from([1, 2, 3])
    G.add_edge(1, 2)
    G.add_edge(2, 3)
    G.add_edge(3, 1)

    triangles = [(1, 2, 3)]

    for iter in range(4,number_of_nodes):
        chosen_triangle_index = randint(0, len(triangles)-1)
        chosen_triangle = triangles[chosen_triangle_index]
        next_node = iter

        # Add new node and edges
        G.add_edge(chosen_triangle[0], next_node)
        G.add_edge(chosen_triangle[1], next_node)
        G.add_edge(chosen_triangle[2], next_node)

        # Add new triangles
        triangles.append((chosen_triangle[0], chosen_triangle[1], next_node))
        triangles.append((chosen_triangle[1], chosen_triangle[2], next_node))
        triangles.append((chosen_triangle[2], chosen_triangle[0], next_node))

        # Remove old one
        del triangles[chosen_triangle_index]

    G.add_edge(1, number_of_nodes)
    G.add_edge(2, number_of_nodes)
    G.add_edge(3, number_of_nodes)


    return G


def inner_get_canonical_ordering(embedding, outer_face):
    """Returns a canonical ordering of the nodes

    The canonical ordering of nodes (v1, ..., vn) must fulfill the following
    conditions:
    - For the subgraph G_k of the input graph induced by v1, ..., vk it holds:
        - 2-connected
        - internally triangulated
        - the edge (v1, v2) is part of the outer face
    - For a node v(k+1) the following holds:
        - The node v(k+1) is part of the outer face of G_k
        - It has at least two neighbors in G_k
        - All neighbors of v(k+1) in G_k lie consecutively on the outer face of
          G_k (excluding the edge (v1, v2)).

    The algorithm used here starts with G_n (containing all nodes). It first
    selects the nodes v1 and v2. And then tries to find the order of the other
    nodes by checking which node can be removed in order to fulfill the
    conditions mentioned above. This is done by calculating the number of
    chords of nodes on the outer face. For more information see [1]_.

    Parameters
    ----------
    embedding : nx.PlanarEmbedding
        The embedding must be triangulated
    outer_face : list
        The nodes on the outer face of the graph

    Returns
    -------
    ordering : list
        A list of tuples `(vk, wp_wq)`. Here `vk` is the node at this position
        in the canonical ordering. The element `wp_wq` is a list of nodes that
        make up the outer face of G_k.

    """
    v1 = outer_face[0]
    v2 = outer_face[1]
    chords = defaultdict(int)  # Maps nodes to the number of their chords
    marked_nodes = set()
    ready_to_pick = set(outer_face)

    # Initialize outer_face_ccw_nbr (do not include v1 -> v2)
    outer_face_ccw_nbr = {}
    prev_nbr = v2
    for idx in range(2, len(outer_face)):
        outer_face_ccw_nbr[prev_nbr] = outer_face[idx]
        prev_nbr = outer_face[idx]
    outer_face_ccw_nbr[prev_nbr] = v1

    # Initialize outer_face_cw_nbr (do not include v2 -> v1)
    outer_face_cw_nbr = {}
    prev_nbr = v1
    for idx in range(len(outer_face) - 1, 0, -1):
        outer_face_cw_nbr[prev_nbr] = outer_face[idx]
        prev_nbr = outer_face[idx]

    def is_outer_face_nbr(x, y):
        if x not in outer_face_ccw_nbr:
            return outer_face_cw_nbr[x] == y
        if x not in outer_face_cw_nbr:
            return outer_face_ccw_nbr[x] == y
        return outer_face_ccw_nbr[x] == y or outer_face_cw_nbr[x] == y

    def is_on_outer_face(x):
        return x not in marked_nodes and (x in outer_face_ccw_nbr.keys() or x == v1)

    # Initialize number of chords
    for v in outer_face:
        for nbr in embedding.neighbors_cw_order(v):
            if is_on_outer_face(nbr) and not is_outer_face_nbr(v, nbr):
                chords[v] += 1
                ready_to_pick.discard(v)

    # Initialize canonical_ordering
    canonical_ordering = [None] * len(embedding.nodes())  # type: list
    canonical_ordering[0] = (v1, [])
    canonical_ordering[1] = (v2, [])
    ready_to_pick.discard(v1)
    ready_to_pick.discard(v2)

    for k in range(len(embedding.nodes()) - 1, 1, -1):
        # 1. Pick v from ready_to_pick
        v = ready_to_pick.pop()
        marked_nodes.add(v)

        # v has exactly two neighbors on the outer face (wp and wq)
        wp = None
        wq = None
        # Iterate over neighbors of v to find wp and wq
        nbr_iterator = iter(embedding.neighbors_cw_order(v))
        while True:
            nbr = next(nbr_iterator)
            if nbr in marked_nodes:
                # Only consider nodes that are not yet removed
                continue
            if is_on_outer_face(nbr):
                # nbr is either wp or wq
                if nbr == v1:
                    wp = v1
                elif nbr == v2:
                    wq = v2
                else:
                    if outer_face_cw_nbr[nbr] == v:
                        # nbr is wp
                        wp = nbr
                    else:
                        # nbr is wq
                        wq = nbr
            if wp is not None and wq is not None:
                # We don't need to iterate any further
                break

        # Obtain new nodes on outer face (neighbors of v from wp to wq)
        wp_wq = [wp]
        nbr = wp
        while nbr != wq:
            # Get next next neighbor (clockwise on the outer face)
            next_nbr = embedding[v][nbr]["ccw"]
            wp_wq.append(next_nbr)
            # Update outer face
            outer_face_cw_nbr[nbr] = next_nbr
            outer_face_ccw_nbr[next_nbr] = nbr
            # Move to next neighbor of v
            nbr = next_nbr

        if len(wp_wq) == 2:
            # There was a chord between wp and wq, decrease number of chords
            chords[wp] -= 1
            if chords[wp] == 0:
                ready_to_pick.add(wp)
            chords[wq] -= 1
            if chords[wq] == 0:
                ready_to_pick.add(wq)
        else:
            # Update all chords involving w_(p+1) to w_(q-1)
            new_face_nodes = set(wp_wq[1:-1])
            for w in new_face_nodes:
                # If we do not find a chord for w later we can pick it next
                ready_to_pick.add(w)
                for nbr in embedding.neighbors_cw_order(w):
                    if is_on_outer_face(nbr) and not is_outer_face_nbr(w, nbr):
                        # There is a chord involving w
                        chords[w] += 1
                        ready_to_pick.discard(w)
                        if nbr not in new_face_nodes:
                            # Also increase chord for the neighbor
                            # We only iterator over new_face_nodes
                            chords[nbr] += 1
                            ready_to_pick.discard(nbr)
        # Set the canonical ordering node and the list of contour neighbors
        canonical_ordering[k] = (v, wp_wq)

    return canonical_ordering

def get_canonical_ordering(g):
    # get planar embedding
    _, embedding = nx.check_planarity(g)
    max_node = len(g)
    ordering = inner_get_canonical_ordering(embedding, [1, 2, max_node])
    return ordering
