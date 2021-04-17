import networkx as nx
import instances_preliminaries
import io_graph_functions
import matplotlib.pyplot as plt
import os
from shapely.geometry import LineString

UNDEFINED = -1

class GraphDrawingAlgorithm:
    UNDEFINED = -1

    def __init__(self):
        self.g = None
        self.u_set = dict()
        self.DC = dict()
        self.dom = dict()
        self.stable = dict()
        self.pos = dict()
        self.contour = []
        self.debug = False
        self.ordering = None
        self.log_file = None

    # Getting dummy ordering, ordering from original article for debugging purposes.
    def dummy_ordering(self):
        self.ordering = [
            (1, []),
            (2, []),
            (3, [1, 2]),
            (4, [1, 3]),
            (5, [1, 4]),
            (6, [1, 5]),
            (7, [3, 2]),
            (8, [7, 2]),
            (9, [8, 2]),
            (10, [4, 3, 7, 8, 9]),
            (11, [6, 5, 4]),
            (12, [11, 4]),
            (13, [12, 4, 10]),
            (14, [6, 11, 12, 13]),
            (15, [10, 9]),
            (16, [10, 15]),
            (17, [1, 6, 14, 13, 10, 16, 15, 9, 2])
        ]
        return self.ordering


    # Original graph instance used in article. For debugging purposes.
    def test_graph(self):
        g = nx.Graph()
        g.add_edge(1, 2)
        g.add_edge(1, 3)
        g.add_edge(1, 4)
        g.add_edge(1, 5)
        g.add_edge(1, 6)
        g.add_edge(1, 17)

        g.add_edge(2, 3)
        g.add_edge(2, 7)
        g.add_edge(2, 8)
        g.add_edge(2, 9)
        g.add_edge(2, 17)

        g.add_edge(3, 4)
        g.add_edge(3, 7)
        g.add_edge(3, 10)

        g.add_edge(4, 5)
        g.add_edge(4, 10)
        g.add_edge(4, 11)
        g.add_edge(4, 12)
        g.add_edge(4, 13)

        g.add_edge(5, 6)
        g.add_edge(5, 11)

        g.add_edge(6, 17)
        g.add_edge(6, 14)
        g.add_edge(6, 11)

        g.add_edge(7, 10)
        g.add_edge(7, 8)

        g.add_edge(8, 9)
        g.add_edge(8, 10)

        g.add_edge(9, 10)
        g.add_edge(9, 15)
        g.add_edge(9, 17)

        g.add_edge(10, 13)
        g.add_edge(10, 17)
        g.add_edge(10, 16)
        g.add_edge(10, 15)

        g.add_edge(11, 12)
        g.add_edge(11, 14)

        g.add_edge(12, 13)
        g.add_edge(12, 14)

        g.add_edge(13, 14)
        g.add_edge(13, 17)

        g.add_edge(14, 17)

        g.add_edge(15, 16)
        g.add_edge(15, 17)

        g.add_edge(16, 17)

        return g

    # Setting initial position of vertices (first three)
    # for now just for A algorithm.
    def set_inital_pos(self):
        # if there is just 3 nodes (basic case)
        if len(self.ordering) == 3:
            self.pos[self.ordering[0][0]] = (0, 0)
            self.pos[self.ordering[1][0]] = (1, 0)
            self.pos[self.ordering[2][0]] = (0, 1)
            return self.pos

        self.contour = [self.ordering[0][0], self.ordering[2][0],  self.ordering[1][0]]
        # embedding first 3 nodes.
        self.pos[self.ordering[0][0]] = (0, 0)
        self.pos[self.ordering[1][0]] = (2, 0)
        self.pos[self.ordering[2][0]] = (1, 1)

    # Getting z vertex for calculation of Domino chains.
    # k is index of vertex in canonical ordering, v_k is node number, w_p is leftmost neighbour of v_k in G_k-1
    def get_z_for(self, k, v_k, w_p):
        for k in range(k+1, len(self.ordering)):
            if self.g.has_edge(v_k, self.ordering[k][0]) and self.g.has_edge(w_p, self.ordering[k][0]):
                return k
        return self.UNDEFINED

    # Recursive function for calculation domino chain for kth node in ordering
    def find_dc_for_k(self, k):
        v_k, contour_neighbors = self.ordering[k]
        if v_k in self.DC.keys():
            return

        z = self.get_z_for(k, v_k, contour_neighbors[0])
        v_z, z_neighbours = self.ordering[z]
        ind_z_v = z_neighbours.index(v_k)+1

        if ind_z_v == 2:
            self.DC[v_k] = [v_k]
            self.dom[v_k] = v_z
            self.stable[v_k] = False

        if ind_z_v >= 4:
            self.DC[v_k] = [v_k]
            self.dom[v_k] = v_z
            self.stable[v_k] = True

        if ind_z_v == 3:
            if z < len(self.g):
                self.find_dc_for_k(z)
            self.DC[v_k] = self.DC[v_z].copy()
            self.DC[v_k].append(v_k)
            self.dom[v_k] = self.dom[v_z]
            self.stable[v_k] = self.stable[v_z]

    # Parent function for domino chains calculation
    def calculate_DC(self):
        n = len(self.g) - 1
        v_n = self.ordering[n][0]
        v_1 = self.ordering[0][0]
        v_2 = self.ordering[1][0]

        self.DC[v_n] = [v_n]
        self.dom[v_n] = self.UNDEFINED
        self.stable[v_n] = True

        self.stable[v_1] = True
        self.stable[v_2] = True

        for k in range(2, len(self.ordering)):
            self.find_dc_for_k(k)

    # Calculating U-sets for shift operations.
    def calculate_u_sets(self):
        for i in range(3):
            vi = self.ordering[i][0]
            self.u_set[vi] = [vi]

        for k in range(3, len(self.ordering)):
            v_k, contour = self.ordering[k]
            self.u_set[v_k] = [v_k]
            for l in range(1, len(contour)-1):
                for t in self.u_set[contour[l]]:
                    self.u_set[v_k].append(t)

    # Shift vertex logic. wq is rightmost neighbor, or node from which to shift.
    def shift_vertex(self, v_k, wq):
        q = self.contour.index(wq)

        if self.debug:
            print('Shifting {nodes} from node {v_k}'.format(v_k=v_k, nodes=self.contour[q:]))

        for i in range(q, len(self.contour)):
            for t in self.u_set[self.contour[i]]:
                self.pos[t] = (self.pos[t][0] + 1, self.pos[t][1])

    # Updating contour after adding new vertex to Gk
    def update_outer_face(self, k):
        v_k, contour_neighbors = self.ordering[k]
        wp = contour_neighbors[0]
        wq = contour_neighbors[-1]
        p = self.contour.index(wp)
        q = self.contour.index(wq)
        self.contour = self.contour[:p+1] + [v_k] + self.contour[q:]

    # Calculating y_k based on visibility
    # x_k is x coordinate of vertex, begin_y is initial value of y coordinate, k is kth  node in ordering.
    def get_y_prime(self, x_k, begin_y, k):
        y_k = begin_y
        increment = 0
        while increment < 100:
            if self.check_visible(x_k, y_k + increment, k):
                break
            increment = increment + 1

        if increment == 99:
            raise Exception('Problems')

        return y_k + increment

    # Checking if position (x_k, y_k) for kth node is satisfying visibility condition against contour C_k-1
    def check_visible(self, x_k, y_k, k):
        vk, neighbours = self.ordering[k]
        for wi in neighbours:
            line = LineString([(x_k, y_k), self.pos[wi]])
            for i in range(0, len(self.contour)-2):
                if wi == self.contour[i] or wi == self.contour[i+1]:
                    continue
                other = LineString([self.pos[self.contour[i]], self.pos[self.contour[i+1]]])
                if line.intersects(other):
                    return False

        return True

    # Printing calculated values of algorithm
    def print_alg_stats(self, path=None):
        if path is not None:
            os.mkdir('{path}/calculations'.format(path=path))
        print("*********** Canonical Ordering ********************")
        if path is not None:
            f = open('{path}/calculations/ordering.txt'.format(path=path), "w")
        for value in self.ordering:
            print('{0} :  {1}'.format(value[0], value[1]))
            if path is not None:
                f.write('{0} :  {1}\n'.format(value[0], value[1]))

        if path is not None:
            f.close()

        print("*********** Domino chains ********************")
        if path is not None:
            f = open('{path}/calculations/domino_chains.txt'.format(path=path), "w")
        for key, value in self.DC.items():
            print('{0} :  {1}'.format(key, value))
            if path is not None:
                f.write('{0} :  {1}\n'.format(key, value))
        if path is not None:
            f.close()

        print("************* Dominators ******************")
        if path is not None:
            f = open('{path}/calculations/domino_chains.txt'.format(path=path), "w")
        for key, value in self.dom.items():
            print('dom({0}) = {1}'.format(key, value))
            if path is not None:
                f.write('dom({0}) = {1}\n'.format(key, value))
        if path is not None:
            f.close()

        print("************* Stable ******************")
        if path is not None:
            f = open('{path}/calculations/stable.txt'.format(path=path), "w")
        for key, value in self.stable.items():
            print('stable({0}) = {1}'.format(key, value))
            if path is not None:
                f.write('stable({0}) = {1}\n'.format(key, value))
        if path is not None:
            f.close()


        print("*************** U Sets ****************")
        if path is not None:
            f = open('{path}/calculations/u_sets.txt'.format(path=path), "w")
        for key, value in self.u_set.items():
            print('U({0}) = {1}'.format(key, value))
            if path is not None:
                f.write('U({0}) = {1}\n'.format(key, value))
        if path is not None:
            f.close()

        print("*************** Positions ****************")
        if path is not None:
            f = open('{path}/calculations/positions.txt'.format(path=path), "w")
        for key, value in self.pos.items():
            print('{0} = {1}'.format(key, value))
            if path is not None:
                f.write('{0} = {1}\n'.format(key, value))
        if path is not None:
            f.close()

        if self.log_file is not None:
            self.log_file.write('Width: {width}, Height: {height}\n'.format(width=self.pos[2][0], height=self.pos[len(self.pos)][1]))
        print('Width: {width}, Height: {height}'.format(width=self.pos[2][0], height=self.pos[len(self.pos)][1]))

    def init_calculations(self, g, ordering):
        self.__init__()
        self.g = g

        # Calculate embedding for canonical ordering
        _, embedding = nx.check_planarity(self.g)
        max_node = len(g)

        # Calculate canonical ordering
        self.ordering = instances_preliminaries.inner_get_canonical_ordering(embedding, [1, 2,
                                                                                         max_node]) if ordering is None else ordering
        self.set_inital_pos()
        self.calculate_DC()
        self.calculate_u_sets()

    def run(self, g):
        pass


class AAlgorithm(GraphDrawingAlgorithm):
    # Run method of algorithm
    # debug         - indicator for debug values
    # output_path   - path for future saving of snapshots during graph forming
    # ordering      - ordering value for custom orderings (for debugging purposes only).

    def run(self, g, debug=False, output_path=None, ordering=None):
        self.init_calculations(g, ordering)

        if output_path is not None:
            self.log_file = open('{output_path}/runtime.log'.format(output_path=output_path), "w")

        if output_path is not None:
            os.mkdir('{output_path}/steps'.format(output_path=output_path))
            os.mkdir('{output_path}/steps/aalgorithm'.format(output_path=output_path))

        if debug:
            self.debug = debug
            dg = nx.Graph() # debug_graph
            dg.add_edge(self.ordering[0][0], self.ordering[1][0])
            v3, v3_neighbours = self.ordering[2]
            for l in v3_neighbours:
                dg.add_edge(v3, l)

        for k in range(3, len(self.ordering)):
            vk, contour_neighbors = self.ordering[k]
            wp = contour_neighbors[0]
            wp1 = contour_neighbors[1]
            wq = contour_neighbors[-1]
            wq1 = contour_neighbors[-2]

            if debug:
                for l in contour_neighbors:
                    dg.add_edge(vk, l)
                print('Node: {node}, neighbours: {neighbours}, outer_face: {outer_face}'.format(node=vk, neighbours=contour_neighbors, outer_face=self.contour))

            if self.log_file is not None:
                self.log_file.write('Node: {node}, neighbours: {neighbours}, outer_face: {outer_face}\n'.format(node=vk, neighbours=contour_neighbors, outer_face=self.contour))

            y_k = max(self.pos[wp1][1], self.pos[wq1][1])

            deg_v_k = len(contour_neighbors)
            if deg_v_k == 2 :
                # shift
                if not self.stable[vk]:
                    self.shift_vertex(vk, wq)

                # upward
                if self.pos[wp][1] < self.pos[wq][1] and self.pos[wp][0] < self.pos[wq][0]:
                    y_k = self.pos[wq][1]

                # horizontal
                if self.pos[wp][1] == self.pos[wq][1]:
                    y_k = self.pos[wq][1] + 1

                # downward
                if self.pos[wp][1] > self.pos[wq][1] and self.pos[wp][0] < self.pos[wq][0]:
                    if self.stable[vk]:
                        y_k = self.pos[wp][1] + 1
                    else:
                        y_k = self.pos[wp][1]

            if self.stable[vk]:
                x_k = self.pos[wp][0]
            else:
                x_k = self.pos[wp][0] + 1

            if deg_v_k != 2:
                y_k = self.get_y_prime(x_k, y_k, k)
            self.pos[vk] = (x_k, y_k)

            self.update_outer_face(k)

            # if debug:
            #     io_graph_functions.draw_graph(dg, self.pos)

            if debug and output_path:
                plt.clf()
                io_graph_functions.save_graph_pic_to_file(dg, '{path}/steps/aalgorithm/adding_{k}_{vk}.png'.format(vk=vk,k=k, path=output_path),pos=self.pos)

        if debug:
            self.print_alg_stats(path=output_path)

        if self.log_file is not None:
            self.log_file.close()

        return self.pos

# Still in progress...
class BAlgorithm(GraphDrawingAlgorithm):
    # Run method of algorithm
    # debug         - indicator for debug values
    # output_path   - path for future saving of snapshots during graph forming
    # ordering      - ordering value for custom orderings (for debugging purposes only).

    def run(self, g, debug=False, output_path=None, ordering=None):
        self.init_calculations(g, ordering)
        if output_path is not None:
            os.mkdir('{output_path}/steps'.format(output_path=output_path))

        if debug:
            self.debug = debug
            dg = nx.Graph() # debug_graph
            dg.add_edge(self.ordering[0][0], self.ordering[1][0])
            v3, v3_neighbours = self.ordering[2]
            for l in v3_neighbours:
                dg.add_edge(v3, l)

        for k in range(3, len(self.ordering)):
            vk, contour_neighbors = self.ordering[k]
            wp = contour_neighbors[0]
            wp1 = contour_neighbors[1]
            wq = contour_neighbors[-1]
            wq1 = contour_neighbors[-2]

            if self.stable[vk]:
                x_k = self.pos[wp][0]
            else:
                x_k = self.pos[wp][0] + 1

            deg_v_k = len(contour_neighbors)
            if deg_v_k == 2:
                if self.stable[vk]:
                    y_k = max(self.pos[wp][1] + 1, self.pos[wq][1])
                else:
                    self.shift_vertex(vk, wq)
                    if self.pos[wp][1] < self.pos[wq][1] and self.pos[wp][0] < self.pos[wq][0]:
                        y_k = self.pos[wq][1]
                    else:
                        y_k = max(self.pos[wp][1], self.pos[wq][1] + 1)
            else:
                r = self.find_r_for_vk(k)
                wr = contour_neighbors[r]
                wr_1 = contour_neighbors[r-1]
                y_p = self.pos[wr][1] + 4 * (x_k-self.pos[wr][0]) - self.slack(wr_1, wr)

                if vk == 13:
                    d= 'd'
                if r == 1 or (not self.stable[vk] and r == 2):
                    y_p = y_p + 1

                y_k = max(y_p, self.pos[wq][1])

            self.pos[vk] = (x_k, y_k)

            if self.slack(vk, wq) == 0:
                self.shift_vertex(vk, wq)

            if debug:
                for l in contour_neighbors:
                    dg.add_edge(vk, l)
                print('Node: {node}, neighbours: {neighbours}, outer_face: {outer_face}'.format(node=vk, neighbours=contour_neighbors, outer_face=self.contour))

            self.update_outer_face(k)

            if debug:
                io_graph_functions.draw_graph(dg, self.pos)

                if debug and output_path:
                    plt.clf()
                    io_graph_functions.save_graph_pic_to_file(dg, '{path}/steps/aalgorithm/adding_{k}_{vk}.png'.format(
                        vk=vk, k=k, path=output_path), pos=self.pos)

        if debug:
            self.print_alg_stats(path=output_path)

        if self.log_file is not None:
            self.log_file.close()

    def slack(self, u, v):
        return 4 * (self.pos[u][0]-self.pos[v][0]) + (self.pos[u][1]-self.pos[v][1])

    def find_r_for_vk(self, k):
        vk, contour_neighbors = self.ordering[k]
        deg_v_k = len(contour_neighbors)
        if deg_v_k == 2:
            return 1 # Index in contour neighbors vector
        else:
            for r in  range(2, len(contour_neighbors)):
                if self.stable[contour_neighbors[r]]:
                    return r
        return len(contour_neighbors)-1