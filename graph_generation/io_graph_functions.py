import networkx as nx
import matplotlib.pyplot as plt
from shift_algorithm import combinatorial_embedding_to_pos

def save_to_file(graph=None, filepath=''):
    nx.write_gpickle(graph, filepath)

def read_from_file(filepath):
    return nx.read_gpickle(filepath)

def write_graph_format_txt(g, pos, path):
    n = len(g)
    f = open(path, "w")
    f.write('{}\n'.format(n))

    for i in range(1, n+1):
        f.write('{x}, {y}\n'.format(x=pos[i][0], y=pos[i][1]))

    for i in range(1, n+1):
        for k in range(1, n + 1):
            if g.has_edge(i, k):
                f.write('1 ')
            else:
                f.write('0 ')
        f.write('\n')
    f.close()

def load_graph_from_format_txt(path):
    g = nx.Graph()
    pos = dict()

    f = open(path, "r")
    line = f.readline()
    n = int(line)
    for i in range(1, n + 1):
        line = f.readline()
        temp = line.split(',')
        pos[i] = (int(temp[0]), int(temp[1]))

    for i in range(1, n + 1):
        line = f.readline()
        temp = line.split(' ')
        for k in range(1, n + 1):
            if temp[k] == '1' and not g.has_edge(i, k):
                g.add_edge(i, k)

    return g, pos

#==================== Graph plotting ===================
def draw_graph(g, pos = None):
    if pos is None:
        pos = combinatorial_embedding_to_pos(g)

    nx.draw(g, pos)
    nx.draw_networkx_labels(g, pos, labels=dict([(n, str(n)) for n in g.nodes()]))
    plt.show()

def save_graph_pic_to_file(g, output_filename, pos = None, color=None ):
    if pos is None:
        pos = combinatorial_embedding_to_pos(g)
    color_map = []
    for node in g:
        if color is None:
            color_map.append('blue')
        else:
            color_map.append(color)
    nx.draw(g, pos, node_color=color_map)
    nx.draw_networkx_labels(g, pos, labels=dict([(n, str(n)) for n in g.nodes()]))
    plt.savefig(output_filename)