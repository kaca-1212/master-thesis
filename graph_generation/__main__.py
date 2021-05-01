import matplotlib.pyplot as plt
import shutil
import os
from io_graph_functions import save_to_file, read_from_file, draw_graph, write_graph_format_txt, load_graph_from_format_txt,save_graph_pic_to_file
from instances_preliminaries import generate_triangulated_graph, get_canonical_ordering
from shift_algorithm import combinatorial_embedding_to_pos
from graph_algorithms import AAlgorithm,BAlgorithm, GraphDrawingAlgorithm
import time

def make_folder_for_instance(id):
    os.mkdir('instances/{id}/'.format(id=id))
    return 'instances/{id}/'.format(id=id)

def run_instance(instance_id, g, ordering=None):
    path = make_folder_for_instance(instance_id)
    os.mkdir('{path}/steps'.format(path=path))
    algA = AAlgorithm()
    algB = BAlgorithm()

    posAlgA = algA.run(g, debug=True, ordering=ordering, output_path=path)
    posAlgb = algB.run(g, debug=True, ordering=ordering, output_path=path)

    posShift = combinatorial_embedding_to_pos(g)

    plt.clf()
    save_graph_pic_to_file(g, '{path}/aalgorithm.png'.format(path=path),pos=posAlgA)
    write_graph_format_txt(g, posAlgA, '{path}/{instance_id}_aalgorithm.txt'.format(path=path, instance_id=instance_id))

    save_graph_pic_to_file(g, '{path}/comparison_shift.png'.format(path=path), pos=posShift, color='red')
    write_graph_format_txt(g, posShift, '{path}/{instance_id}_shift.txt'.format(path=path, instance_id=instance_id))

    plt.clf()
    save_graph_pic_to_file(g, '{path}/shift.png'.format(path=path), pos=posShift)

    plt.clf()
    save_graph_pic_to_file(g, '{path}/balgorithm.png'.format(path=path), pos=posAlgb)
    write_graph_format_txt(g, posAlgb, '{path}/{instance_id}_balgorithm.txt'.format(path=path, instance_id=instance_id))


def main():
    # clear instances data from before
    try:
        shutil.rmtree('instances/')
    except OSError as e:
        print("Error: %s : %s" % ('instances/', e.strerror))
    os.mkdir('instances/')

    # test on particular graph instance
    alg = AAlgorithm()
    g = alg.test_graph()
    run_instance('test_instance', g)

    # test on random generated graphs
    # for i in range(1, 20):
    #     instance_id = 'graph_{i}'.format(i=i)
    #     g = generate_triangulated_graph(30)
    #     run_instance(instance_id, g)

    instance_id = 'graph_{i}'.format(i=1)
    alg = BAlgorithm()
    g = alg.test_graph()
    run_instance(instance_id, g)

    return

start_time = time.time()
main()
print("--- %s seconds ---" % (time.time() - start_time))