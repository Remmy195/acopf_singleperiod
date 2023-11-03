import igraph as ig
import matplotlib.pyplot as plt

### module for computing the min cut


def mincut(log,all_data):
    
    num_buses = all_data['buses']
    edges = all_data['edges']
    capacities = all_data['capacities']
    
    g = ig.Graph(n_vertices, edges,directed = True)


    flow = g.maxflow(0,3,capacity=g.es["capacity"])

    print("Max flow",flow.value)
    print("Edge assignments:",flow.flow)



    print("--------------")

    cut = g.mincut()

    print("Min cut",cut.value)
    print("vertices in S",cut.partition)
    print("Edge assignments:",cut.cut)


    fig, ax = plt.subplots()
    ig.plot(
        g,
        target=ax,
        layout="circle",
        vertex_label=range(g.vcount()),
        vertex_color="lightblue"
    )
    plt.show()
