"""
    Some simple heuristics for the symmetric TSP.

    author :
        Andreas Klose
    version :
        1.0 (29/3/2019)
"""
from numpy.random import choice
from numpy import Inf
import networkx as nx
from . import tsp_plot
    
#---------------------------------------------------------------------

def nearestNeighbour( problem ):
    """
    Computes solution using nearest neighbour with random start point.

    Parameters:
        problem : class
            TSP problem object has returned by function load_problem
            of the package tsplib95.

    Returns:
        routeLen : int
            Length of the route.
        route : list of int    
            List of points on the route.
    """
    unRouted = [ i for i in problem.get_nodes()]
    routeLen = 0
    # Choose random start point
    node = choice(unRouted)
    route = [ node ]
    unRouted.remove(node)
    while len(unRouted) > 0:
        nextNode = min( unRouted, key = lambda j : problem.wfunc(node,j) )
        unRouted.remove(nextNode)
        route.append(nextNode)
        routeLen += problem.wfunc(node,nextNode)
        node = nextNode
    routeLen += problem.wfunc(node,route[0])
    return routeLen, route

#---------------------------------------------------------------------

def doInsertion( problem, nearest=True ):
    """
    Computes solution using nearest or farthest insertion with random 
    start point.

    Parameters:
        problem : class
            TSP problem object has returned by function load_problem
            of the package tsplib95.
        nearest : bool, optional
            If false, farthest insertion is applied. Otherwise
            nearest insertion is applied (default)

    Returns:
        routeLen : int
            Length of the route.
        route : list of int    
            List of points on the route.
    """
    def insCost( node, i, route ):
        """
        Return cost to insert node at position between i and i+1 in 
        the route.
        """
        j1 = route[i]
        j2 = route[(i+1) % len(route)]
        delta = problem.wfunc(node,j1)+problem.wfunc(node,j2)-problem.wfunc(j1,j2)
        return delta
    
    # Initialize the the route with some randomly selected edge
    unRouted = [ node for node in problem.get_nodes()]
    n1 = choice(unRouted)
    unRouted.remove(n1)
    n2 = choice(unRouted)
    unRouted.remove(n2)
    route = [n1,n2,n1]
    routeLen = problem.wfunc(n1,n2)*2
    sign = 1
    if not nearest : sign = -1
    while len(unRouted) > 0:
        # Determine unrouted node closest (or farthest) to a routed node
        nextNode = None
        minDist  = Inf
        for j in unRouted:
            dist = min( [sign*problem.wfunc(j,node) for node in route] )
            if dist < minDist:
                dist = minDist
                nextNode = j
        # Insert nextNode at cheapest insertion place
        ins_at = min([i for i in range(len(route))],\
                      key = lambda i : insCost( nextNode, i, route ) )
        routeLen += insCost(nextNode,ins_at,route)
        route.insert(ins_at+1,nextNode)
        unRouted.remove( nextNode )

    return routeLen, route

#---------------------------------------------------------------------

def minSpanTree( problem, display=True, animate=True ):
    """
    Apply the minimum spanning tree heuristic.

    Parameters:
        problem : class
            TSP problem object has returned by function load_problem
            of the package tsplib95.
        display : bool, optional
            I true (default) a figure of the minimum spanning tree, 
            the Eulerian cycle and the resulting route is displayed.
        animate : bool, optional
            If true (default) the Eulerian cycle and the route are 
            drawn in an animated way.
              
    Returns:
        routeLen : int
            Length of the route.
        route : list of int    
            List of points on the route.
    """
    # Represent problem data using a networkx graph object
    G = problem.get_graph()

    # Get min. spanning tree and obtain Eulerian graph 
    E = nx.MultiGraph()
    E.add_nodes_from( G.nodes )
    cnt = 0
    tLen = 0
    for e in nx.minimum_spanning_edges(G):
        len = e[2]['weight']
        tLen += len
        E.add_edge(e[0],e[1],key=cnt,weight=len)
        E.add_edge(e[0],e[1],key=cnt+problem.dimension,weight=len)
        cnt+=1

    # Display the minimum spanning tree 
    if display: 
        tsp_plot.plotEdgeList( problem, E.edges, title='Minimum spanning '+\
                              'tree/Eulerian graph. Length= '+str(tLen))
  
    circuit = [u for u,v in nx.eulerian_circuit(E)]
    circuit.append(circuit[0])
    if display:
        tsp_plot.displayRoute(problem.node_coords, circuit, 2*tLen,\
                              animate=animate, title='Eulerian circuit')

    # Build the TSP route by taking shortcuts in the Eulerian cycle
    route = [circuit[0]]
    routeLen = 0
    for node in circuit:
      if node not in route: 
          routeLen += problem.wfunc(route[-1],node)
          route.append(node)
    routeLen += problem.wfunc(route[-1],route[0])
    route.append(route[0])
    if display:
        tsp_plot.displayRoute(problem.node_coords, route, routeLen,\
                              animate=animate,title='Found tour')
    else:
        print("Length of route : ",routeLen)
        print("Route found: ",route )
    
    return routeLen, route    
        
#---------------------------------------------------------------------

def christofides( problem, display=True, animate=True ):
    """
    Apply the Christofides' heuristic.

    Parameters:
        problem : class
            TSP problem object has returned by function load_problem
            of the package tsplib95.
        display : bool, optional
            I true (default) a figure of the minimum spanning tree, 
            the Eulerian cycle and the resulting route is displayed.
        animate : bool, optional
            If true (default) the Eulerian cycle and the route are 
            drawn in an animated way.
                            
    Returns:
        routeLen : int
            Length of the route.
        route : list of int    
            List of points on the route.
    """
    # Represent problem data using a networkx graph object
    G = problem.get_graph()

    # The graph may contain edges connecting a node to itself
    # We first remove all these simple loops
    G.remove_edges_from(G.selfloop_edges())

    # Initialize the graph (that later will be Eulerian)
    # as a multigraph that first contains the edges from
    # a minimim spanning tree of G
    E = nx.MultiGraph()
    E.add_nodes_from( G.nodes )
    cnt = 0
    tLen = 0
    for e in nx.minimum_spanning_edges(G):
        len = e[2]['weight']
        tLen += len
        E.add_edge(e[0],e[1],key=cnt)

    # Find the nodes of odd degree in the minimum spanning tree
    odd = []
    for node in E.nodes:
        if E.degree(node) % 2 > 0: odd.append(node)

    # Display the minimum spanning tree and mark the odd nodes
    if display: 
        tsp_plot.plotEdgeList( problem, E.edges, specialNodes=odd, title=\
                              'Minimum spanning tree (odd nodes are red)')

    # Obtain the subgraph of G induced by the "odd" nodes
    oddG = G.subgraph( odd )

    # Compute a perfect minimum-cost matching in the graph O of 
    # odd nodes. We do that by transforming the problem into the 
    # problem of finding a maximum weight matching of maximum 
    # cardinality
    wmax = max([e[2]['weight'] for e in oddG.edges(data=True)])
    for u, v, w in oddG.edges(data=True):
        ww = w['weight']
        w['weight'] = wmax - ww
    M = nx.max_weight_matching(oddG, maxcardinality=True)

    # Add the matching edges to the minimum spanning tree
    # which thereby is transformed to an Eulerian graph
    cnt = problem.dimension
    for e in M: 
        tLen += problem.wfunc( e[0], e[1] )
        E.add_edge(e[0],e[1],key=cnt)
        cnt+=1
    if display:
        tsp_plot.plotEdgeList( problem, E.edges, specialEdges=M, specialNodes=odd,\
                               title="Spanning tree + matching of odd nodes (red edges)")

    # Obtain an Eulerian circuit in the graph oddG+M
    circuit = [u for u,v in nx.eulerian_circuit(E)]
    circuit.append(circuit[0])
    if display:
        tsp_plot.displayRoute(problem.node_coords, circuit, 2*tLen, title=\
                              'Eulerian circuit', animate=animate)

    # Build the TSP route by taking shortcuts in the Eulerian cycle
    route = [circuit[0]]
    routeLen = 0
    for node in circuit:
      if node not in route: 
          routeLen += problem.wfunc(route[-1],node)
          route.append(node)
    routeLen += problem.wfunc(route[-1],route[0])
    route.append(route[0])
    if display:
        tsp_plot.displayRoute(problem.node_coords, route, routeLen,\
                              title='Found tour',animate=animate)
    else:
        print("Length of route : ",routeLen)
        print("Route found: ",route )

    return routeLen, route
