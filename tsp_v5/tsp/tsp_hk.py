"""
    Held and Karp (Lagrangian) lower bound for the symmetric TSP
    based on 1-trees.

    author:
        Andreas Klose
    version:
        0.0 (29/3/2019). Untested version. Under development!!! 
"""
import networkx as nx
from tsp_v5.tsp.tsp_trees import get1Tree
      
#---------------------------------------------------------------------

def __getGradient( G, tree ):
    """
    Returns the subgradient, i.e. the vector with coordinates
    2 - degree_of_i_in_tree for all nodes i.    
    """
    # degree of node i in the 1 tree
    dg = lambda i : len([ e for e in tree if i in e])
    # subradient
    sg = [ 2 - dg(i) for i in G.nodes ]
    nrm = sum( g*g for g in sg)
    
    return nrm, sg  
    
#---------------------------------------------------------------------

def LR_1tree ( problem, G=None, silent=True ):
    """
    Implements the Lagrangian 1-tree lower bound obtained by relaxing
    degree constraints in Lagrangian manner and (approximately) solving
    the Lagrangian dual by means of subgradient optimization.
    
    Parameters:
        problem : class
            TSP problem object as returned by function load_problem of
            package tsplib95
        G : object, optional
            If not None, the networkx graph of the problem. Default=None.
        silent : bool, optional
            if false, information about the computations is printed.
            Default=True.

    Returns:
        lowBnd : float
            The Lagrangian bound. 
        T_best : list of pairs of int
            The 1-tree (edge list) that obtained when solving the Lagrangian
            subproblem at the best multipliers obtained.
        best_w : list of float
            The list of Lagrangian multiplier values.
    """
    k = min( i for i in G.nodes) # the special node
    
    lowBnd = 0.0
    G_local = G is None 
    if G_local:
        G = problem.get_graph()
        G.remove_edges_from(G.selfloop_edges())

    # Initialize current and best Lagrangian multiplier values
    best_w = [0.0 for i in G.nodes]
    cur_w = [ 0.0 for i in G.nodes]
    best_T = []
        
    # Introduce the Lagrangian multiplier as additional node attribute
    nx.set_node_attributes(G,[],'weight')
    cnt = 0
    for i in G.nodes:
        G.nodes[i]['weight'] = cur_w[cnt]
        cnt += 1
         
    # Introduce the modified edge lengths as additional edge attribute
    eweight = [ G.get_edge_data(*e)['weight'] for e in G.edges]
    nx.set_edge_attributes(G,[],'eweight')
    cnt = 0
    for i in G.edges:
        G[e[0]][e[1]]['eweight'] = eweight[cnt]
        cnt += 1        
        
    iter_max = 10*len(G)
    lam_para = 0.95
    stop = False
    step = 2.0
    iter = 0
    
    # subgradient in previous iteration
    sg_prev = [0.0 for i in G.nodes]
    
    if not silent:
        print("----------------------------------------")
        print("Iter  Lower_Bound  Best_Bound  Grad.norm")
        print("----------------------------------------")
    
    while not stop:
       
        iter += 1
        
        # Compute the 1-tree for the current multiplier values
        cur_bnd, tree = __get1Tree(k, G, elen='eweight' )
        cur_bnd -= 2*sum( cur_w )
        
        # Obtain the subgradient 
        nrm, sg  = __getGradient( G, tree )

        # Check for bound improvement
        if cur_bnd > lowBnd:
            lowBnd = cur_Bnd
            best_w = [ w for w in cur_w]
            T_best = [ e for e in tree]
            
        if nrm < 1.0E-4: break             
          
        # Apply subgradient step
        alpha = 0.7 + 0.3*(iter < 2 )
        for i in range(len(G)): cur_w[i] += step*(alpha*sg[i]+(1.0-alpha)*sg_prev[i])
        sg_prev = sg
        step *= lam_para
        if step < 1.0E-6: break
        if iter >= iter_max: break;
        
        # Display info on current iteration
        if not silent:
            print('{0:4d}  {1:11.2f}  {2:10.2f}  {3:9.2f}\n'.format(iter,cur_bnd,lowBnd,nrm))
        
        # Adjust modified edge length
        for e in G.edges:
            i, j = e[0], e[1]
            G[i][j]['eweight'] = G[i][j]['weight'] \
                                 + G.nodes[i]['weight'] + G.nodes[j]['weight']
                                 
    # Subgradient steps finished
    if not G_local:
        for e in G.edges: del G[e[0]][e[1]]['eweight']
        for i in G.nodes: del G.nodes[i]['weight']                                                           
    
    return lowBnd, T_best, best_w
