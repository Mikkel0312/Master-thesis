"""
    2-matching and assignment lower bound for the symmetric TSP.

    author: 
        Andreas Klose
    version: 
        1.0 (29/3/2019)
"""
import networkx as nx
import cplex
  
#---------------------------------------------------------------------

def assignLB ( problem ):
    """
    Uses networkx min_cost_flow algorithm for computes the assignment
    lower bound. Alternatively, also networkx network simplex algorithm
    could be used instead.

    Parameters:
        problem : class
            The TSP problem object has returned by function load_problem
            of package tsplib95.

    Returns:
        cost : int
            The lower bound (objective value of the assignment problem).
        elist : list of int
            The list of edges in the assignment.
    """
    # Solve assignment problem as a minimum-cost network flow problem
    n = problem.dimension
    G = nx.DiGraph()
    for i in problem.get_nodes(): G.add_node(i, demand = -1 )
    for j in problem.get_nodes(): G.add_node(j+n, demand = 1 )
    for i in problem.get_nodes():
        for j in problem.get_nodes():
             if i != j:
                 G.add_edge(i,j+n, weight=problem.wfunc(i,j), capacity = n+1 )
    flow = nx.min_cost_flow( G )
    cost = nx.cost_of_flow(G,flow)

    # Retrieve list of edges in the assignment
    elist = [ (e[0],e[1]-n) for e in G.edges if flow[e[0]][e[1]] > 0 ]
    return cost, elist    

#---------------------------------------------------------------------        

def twoMatchModel( c, G, fractional=False ):
    """
    Sets up the (integer or fractional) 2-matching linear program
    using Cplex.

    Parameters:
    G : object (networkx graph)
        The graph holding the TSP data. The edge variables belonging to 
        the model are inserted as additional edge attribute "var" in the 
        graph G.
    fractional : bool, optional
        If false (the default), the integer program is set up; otherwise 
        the linear relaxation.
    """
    # First, set the objective sense to be minimization
    c.objective.set_sense(c.objective.sense.minimize)

    # Add variables, one for each edge
    obj = [ G.get_edge_data(*e)['weight'] for e in G.edges ]
    if fractional:
        vars = c.variables.add ( ub=[1]*len(obj), obj=obj )
    else:
        vars = c.variables.add ( types=["B"]*len(obj), obj=obj )

    # Store which variable belongs to each edge
    nx.set_edge_attributes(G, [], 'var')
    cnt = 0
    for e in G.edges:
        G[e[0]][e[1]]['var'] = vars[cnt]
        cnt += 1
        
    # Add the degree constraints for each node
    for i in G.nodes:
        vars_i = [ G[e[0]][e[1]]['var'] for e in G.edges(i) ]
        lhs_i = [ cplex.SparsePair( vars_i, [1]*len(vars_i) ) ]
        c.linear_constraints.add( lin_expr=lhs_i, senses=['E'], rhs=[2] )

#---------------------------------------------------------------------        

def twoMatch( problem, G = None, fractional=False ):
    """
    Computes the integer or fractional 2-matching lower bound. The 
    (integer) linear program of the 2-matching problem is solved using 
    the Cplex optimizer. 

    Parameters:
        problem : class
            The TSP problem object has returned by function load_problem
            of package tsplib95.
        G : object (networkx graph), optional
            The graph holding the TSP data. Default=None
        fractional : bool, optional
            If true, only the linear relaxation of the 2-matching 
            problem is solved. Default=False.

    Returns:
        weight : int or float
            The (fractional) 2-matching lower bound.
         elist : list of pairs of int
             The list of edges in the 2-matching and a corresponding.
         xvals : list of floats     
            List of (fractional) solution values. xvals is only
            returned if fractional is true.
    """

    G_local = ( G is None )
    if G_local:
        G = problem.get_graph()
        G.remove_edges_from(G.selfloop_edges())

    # Create an instance of the Cplex class from package cplex
    c = cplex.Cplex()
    c.set_results_stream(None) # suppresses output from Cplex

    # Build the 2-matching model
    twoMatchModel( c, G, fractional=fractional )

    # Solve the 2-matching problem
    c.solve()

    # Retrieve objective function value and the edges e with x(e) > 0
    s = c.solution
    weight = s.get_objective_value()
    elist = [ e for e in G.edges if s.get_values( G[e[0]][e[1]]['var'] ) > 1.0E-5 ]
    if fractional:
        xvals = [ s.get_values(G[e[0]][e[1]]['var']) for e in elist ]
    
    # Release the cplex object
    c.end()

    # If graph is not a local object, remove the edge attribute "var" again
    if not G_local:
        for e in G.edges: del G[e[0]][e[1]]['var']

    if fractional:
        return weight, elist, xvals
    else:
        return weight, elist
