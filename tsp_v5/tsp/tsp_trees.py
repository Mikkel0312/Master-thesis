"""
    Computation of 1-trees as relaxation of the symmetric TSP.

    author:
        Andreas Klose
    version:
        1.0 (29/3/2019)
"""
import networkx as nx
   
#---------------------------------------------------------------------

def get1Tree( k, G, elen='weight' ):
    """
    Computes a minimum weight 1-tree of a graph G with special node k.

    Parameters:
        k : int
            The special node of the graph.
        G : object (networkx graph)
            Connected undirected graph.
        elen : str
            The edge attribute to be used as edge lengths when
            computing the tree. Default is the ordinary edge weight.
    
    Returns:
        weight : 
            The weight of the 1-tree.
        tree : list of pairs of int
            The list of edges in the 1-tree. 
    """  
    # Find cheapest edge adjacent to special node k
    kadj = list( G.edges(k))
    wk   = [ G.get_edge_data(*e)[elen] for e in kadj]
    emin = min(kadj, key = lambda e: G.get_edge_data(*e)[elen])
    weight = G.get_edge_data(*emin)[elen] 
    
    # Temporarily set weight of cheapest edge to zero, all other to "big"
    for e in kadj: G[e[0]][e[1]][elen] = 1.0E20
    G[emin[0]][emin[1]][elen] = 0

    # Obtain minimum-weight spanning tree of the graph
    tree = [ (e[0],e[1]) for e in nx.minimum_spanning_edges(G,weight=elen)]
    weight += sum(G.get_edge_data(*e)[elen] for e in tree )

    # Reset the weights of edges incident to k to the old values
    cnt = 0
    for e in kadj:
        G[e[0]][e[1]][elen] = wk[cnt]
        cnt += 1
    
    # Include 2nd cheapest edge incident to special node k
    kadj.remove(emin)
    emin = min(kadj, key = lambda e: G.get_edge_data(*e)[elen])
    weight += G.get_edge_data(*emin)[elen]
    tree.append(emin)

    # Return weight of the 1-tree and the list of edges
    return weight, tree

#---------------------------------------------------------------------

def oneTree( problem, G = None, specNodes=None ):
    """
    Computes the 1-tree lower bound by searching over a set of nodes
    as special nodes.

    Parameters:
        problem : class
            The TSP problem object has returned by function load_problem
            of the package tsplib95.
        G : object (networkx graph)
            if not None, it need to be networkx graph representation of 
            the problem as returned by problem.get_graph(). Default=None.
        specNodes : list of int
            A list of nodes that are used as special nodes; if None
            (the default) all nodes are tried as special nodes.

    Returns:
        weight : 
            The weight of the best 1-tree found.
        tree : list of pairs of int
            The networkx graph object representing this 1-tree. 
    """
    if G is None: 
        G = problem.get_graph()
        G.remove_edges_from(G.selfloop_edges())

    if specNodes is None:
        specNodes = list(G.nodes)

    weight = 0
    for k in specNodes:
        wk, tree_k = get1Tree(k,G)
        if wk > weight:
            weight = wk
            tree = tree_k

    return weight, tree

#---------------------------------------------------------------------

def fast1Tree( problem, G = None ):
    """
    Reinelt's fast method for computing a 1 tree lower bound.

    Parameters:
        problem : class
            The TSP problem object has returned by function load_problem
            of the package tsplib95.
        G : object (networkx graph)
            if not None, it need to be networkx graph representation of 
            the problem as returned by problem.get_graph(). Note that 
            selfloops must be removed from the graph.
            Default=None.
    
    Returns:
        weight :
            The weight of the best 1-tree found.
        tree : list of pairs of int 
            List of the edges in the 1-tree.
    """
    if G is None: 
        G = problem.get_graph()
        G.remove_edges_from(G.selfloop_edges())

    # Compute minimum spanning tree using networkx
    T = nx.minimum_spanning_tree(G)
    tweight = sum(T.get_edge_data(*e)['weight'] for e in T.edges)

    # Initialize weight of best 1-tree
    weight = 0

    # Investigate the leafs of the tree one by one
    leafs = [i for i in T.nodes if T.degree(i)==1]
    for l in leafs:
        # Get the edge in the tree incident to node l
        t_edge = [e for e in T.edges(l)][0] 
        # Get the other edges in G incident to node l
        l_edges = [e for e in G.edges(l) if e != t_edge]
        # Obtain edge of smallest weight from "l_edges"
        emin2 = min(l_edges, key = lambda e: G.get_edge_data(*e)['weight'])
        # Get weight of 1-tree with special node l
        l_weight = tweight + G.get_edge_data(*emin2)['weight']
        # Store the edge if it improves the 1-tree lower bound
        if l_weight > weight:
            weight = l_weight
            ebest  = emin2
    # Collect the edges of the 1-tree in a list
    tree = [e for e in T.edges]
    tree.append( ebest)

    return weight, tree
