"""
Solving the symmetric TSP by means of CPLEX's MIP solver.

author :  
    Andreas Klose
version : 
    1.0 (29/3/2019)
"""
import networkx as nx
from tsp_v5.tsp.tsp_match import twoMatchModel
from tsp_v5.tsp.tsp_plot import plotEdgeList
import cplex
from cplex.callbacks import UserCutCallback
from cplex.callbacks import LazyConstraintCallback
from itertools import combinations as combinat

#---------------------------------------------------------------------

class __SCEcuts( UserCutCallback ):
    """
    Finds sub-cycle elimination (SCE) constraints violated by a 
    fractional solution to the current node LP. Limit maybe to first or 
    first three nodes in search tree!!!
    """
    
    def __call__(self):
        
        G = self.G
        n = len(G)
        
        # Get list of edges with positive value in the current LP solution
        elist = [e for e in G.edges if self.get_values(G[e[0]][e[1]]['var'])>1.0E-5]
        
        # If with animation, show the (fractional) solution
        if self.anim:
            objv = self.get_objective_value()
            xval = [ self.get_values( G[e[0]][e[1]]['var'] ) for e in elist ]
            efrac = [ e for e in elist if self.get_values(G[e[0]][e[1]]['var'])< 0.99 ]
            plotEdgeList(self.problem, elist, specialEdges=efrac,\
                  title="Fractional LP solution of objval="+str(objv) )
        
        # Build the solution's support graph, i.e. the sub-graph induced by
        # the edges in the list "elist" found above.
        supG = nx.edge_subgraph(G, elist)
        
        # If the support graph G is not connected, build the SCE constraint
        # based on the graph's smallest component.
        S = min( nx.connected_components(supG), key=len )
        weight = 1.0
        
        # If the graph is connected, obtain the minimum-weight cut set
        # using the algorithm of Stoer and Wagner
        if len(S)==n:
            # Set the edge weights equal to the variable's solution values
            for e in elist:
                i, j =  e[0], e[1]
                supG[i][j]['weight'] = self.get_values( G[i][j]['var'] )
            weight, part = nx.stoer_wagner(supG)    
            S = part[0] if len(part[0]) < len(part[1]) else part[1]            
        
        # If the cutset constraint is violated, include the cut in 
        # form of the equivalent subcycle elimination constraint
        if ( len(S) > 2 ) and ( weight < 1.98):
            E = list(combinat(S,2))
            cycVars = [ G[e[0]][e[1]]['var'] for e in E ]
            self.add( cplex.SparsePair(cycVars, [1.0]*len(E) ), sense="L", \
                      rhs=len(S)-1)           
        
#---------------------------------------------------------------------

class __LazySCElim( LazyConstraintCallback ):
    """
    Finds sub-cycles in an integer solution and adds the constraint
    to the problem.
    """
              
    def __call__(self):     
        """ 
        Solution is integer. Check for sub-cycles and if there is one
        store the smallest.
        """
    
        G = self.G
        n = len(G)
        
        # List of edges in the solution
        elist = [e for e in G.edges if self.get_values(G[e[0]][e[1]]['var'])>0.99999]
        
        # If with animation, show the (fractional) solution
        if self.anim:
            objv = self.get_objective_value()
            __title = "Checking lazy constraints. Integer solution objval="+str(objv)
            plotEdgeList(self.problem, elist, title=__title )
        
        # Shortest subcycle encountered in the solution
        cyc = None
        cycLen = n
        
        # Search for shortest subcycle
        while len(elist) > 0:
            e = elist.pop(0) 
            cur_cyc = [e[0],e[1]]    
            am_back = False
            while not am_back:
                last = cur_cyc[-1]
                edge = next( e for e in elist if last in e )
                succ = edge[1] if last==edge[0] else edge[0] 
                am_back = succ == cur_cyc[0]
                if not am_back: cur_cyc.append(succ)
                elist.remove(edge)    
            cur_len = len( cur_cyc )
            if cur_len < cycLen:
                cycLen = cur_len
                cyc = cur_cyc
            if cycLen==3: break          
        
        # If subcycle found, include the subcycle elimination constraint
        if cycLen < n:
            edges = list( combinat(cyc,2) )
            cycVars = [ G[e[0]][e[1]]['var'] for e in edges ]
            self.add( cplex.SparsePair(cycVars, [1.0]*len(edges)), sense="L", \
                      rhs=cycLen-1)
                 
#---------------------------------------------------------------------

def __makeRoute( G, edges ):
    """
    Transforms solution given as edge list to a solution in permutation 
    form.
    
    Parameters:
        G : object
            Networkx graph as returned by method get_graph() of the
            problem class in package "tsplib95".
        edges : list of int
            The list of edges selected in the solution.
            
    Returns:
        route : list of int
            The route as a list of nodes.             
    """
    e = edges.pop(0)
    route = [ e[0], e[1] ]
    while len(edges) > 0:
        last = route[-1]
        edge = next( e for e in edges if last in e )
        succ = edge[1] if last==edge[0] else edge[0]
        route.append(succ)
        edges.remove(edge)        
    return( route )
#---------------------------------------------------------------------

def CPXtspSolve( problem, G = None, rootOnly=True, startRoute=None,\
                 nodeLim=None, anim=False ):
    """
    Uses CPLEX's MIP solver for solving the symmetric TSP. Subcycle 
    elimination constraints are included via cut-callback functions.    

    Parameters:
        problem : class
            TSP problem object has returned by function load_problem
            of package tsplib95.
        G : object
            If not None, the networkx graph of the problem.
            Default = None.
        rootOnly : bool, optional  
            If true, only the root node of the branch-and-cut tree
            is processed. Default=True.
        startRoute : list of int, optional
            If not None it specifies a feasible route given as the 
            sequence in which the nodes are visited. Default=None
        nodeLim : int, optional
            Limit on the nodes to be enumerated. Default=None               
        anim : bool, optional
            If true, the (fractional) LP solution obtained before 
            calling a cut callback is displayed (attention: can be 
            time consuming). Default=False.                        

    Returns:
        lowbnd : float
            A lower bound or the optimal solution value.
        routeLen : int
            Length of the computed route.
        route : list of int        
            The route computed.
    """

    G_local = ( G is None )
    if G_local:
        G = problem.get_graph()
        G.remove_edges_from(G.selfloop_edges())
                        
    # Create an instance of the Cplex class from package cplex
    c = cplex.Cplex()
    
    #Initial route provided?
    routeLen = None
    route = startRoute
    if (not route is None):
        routeLen = sum([G[route[i-1]][route[i]]['weight'] \
                        for i in range(1,len(route))] ) 
        routeLen += G[route[0]][route[-1]]['weight']     
        #Pass solution value as upper cutoff value to Cplex
        c.parameters.mip.tolerances.uppercutoff.set(routeLen)    

    #if rootOnly: c.set_results_stream(None) 

    # Build the 2-matching integer program
    twoMatchModel( c, G )
    
    # For using callbacks, we need traditional B&C search
    c.parameters.mip.strategy.search.set( \
        c.parameters.mip.strategy.search.values.traditional)
    
    # Set the node limit if given
    if not nodeLim is None:
        c.parameters.mip.limits.nodes.set( nodeLim )

    # Switch off "advanced" start
    c.parameters.advance.set(0)

    # Install the lazy constraint callback
    c.register_callback( __LazySCElim )
    __LazySCElim.G = G
    __LazySCElim.anim = anim
    if anim: __LazySCElim.problem = problem
    
    # Install the user cut callback
    c.register_callback( __SCEcuts )
    __SCEcuts.G = G
    __SCEcuts.anim = anim
    if anim: __SCEcuts.problem = problem

    # Solve problem with lazy constraints
    c.solve()
    s = c.solution
    
    # Obtain status of the solution
    best_bnd = s.MIP.get_best_objective()
    if ( s.get_status() == s.status.MIP_feasible ) or\
       ( s.get_status() == s.status.MIP_optimal ):
        objval = s.get_objective_value()
        if (route is None) or (objval < routeLen ):
            elist = [e for e in G.edges if s.get_values(G[e[0]][e[1]]['var'])>1.0E-5]
            route = __makeRoute( G, elist )
            routeLen = objval
    c.end()

    # If graph is not a local object, remove the edge attribute "var" again. 
    # Restore also the original weights
    if not G_local:
        for e in G.edges: del G[e[0]][e[1]]['var']
        for e in G.edges: G[e[0]][e[1]]['weight'] = problem.wfunc(e[0],e[1])

    return best_bnd, routeLen, route
