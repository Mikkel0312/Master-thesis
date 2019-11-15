#!/usr/bin/python3
"""
    Some simple heuristics for the symmetric TSP

    @author  Andreas Klose
    @version 15/3/2019
"""

from tsplib95 import load_problem

from tsp.tsp_heu import nearestNeighbour
from tsp.tsp_heu import doInsertion
from tsp.tsp_heu import minSpanTree
from tsp.tsp_heu import christofides

from tsp.tsp_plot import displayPoints
from tsp.tsp_plot import displayRoute
from tsp.tsp_plot import plotEdgeList

from tsp.tsp_match import assignLB
from tsp.tsp_match import twoMatch

from tsp.tsp_trees import oneTree
from tsp.tsp_trees import fast1Tree

from tsp.concorde import solveTSP

from tsp.tsp_cpx import CPXtspSolve

#import tsp.tsp_hk

#---------------------------------------------------------------------

def __showSol( problem, route, routeLen, info ):
    """
    Display (if possible graphically) the obtained solution  
    """
    if problem.is_depictable():
        displayRoute(problem.node_coords,route,routeLen,title=info)
    else:
        print("Solution obtained by",info)
        print("Tour length:",routeLen)
        print("Tour       :",route)
    
#---------------------------------------------------------------------

def __heurSolve( problem, method='LK', showResult=False ):
    """
    Obtain a heuristic solution using method "method". Default is
    Concorde's Lin-Kernighan heuristic. If showResult is True
    the solution is shown in a (animated) figure.
    """
    showFig = False
    if method == 'NN':
        header = 'Nearest neighbour route'
        routeLen, route = nearestNeighbour( problem )
    elif method == 'NI':
        header = 'Nearest insertion route'
        routeLen, route = doInsertion( problem )
    elif method == 'FI':
        header = 'Farthest insertion route'
        routeLen, route = doInsertion( problem, nearest=False )
    elif method == 'MST':
        header = 'MST heuristic solution'
        showFig = display=problem.is_depictable() and showResult
        routeLen, route = minSpanTree( problem, display=showFig )
    elif method == 'CHRISTO':
        header = 'Christofides solution'
        showFig = display=problem.is_depictable() and showResult
        routeLen, route = christofides( problem, display=showFig )
    elif method == 'LK':
        header = 'Lin-Kernighan (Concorde) solution'
        routeLen, route = solveTSP( problem, exact=False, logFile="CC_log.log" )
    else:
        routeLen, route = None, None
        
    if showResult and (not showFig) and ( not route is None ):
        __showSol( problem, route, routeLen, header )
            
    return routeLen, route
             
#---------------------------------------------------------------------

if __name__ == "__main__":

    from sys import argv
    # Check if input file given on command line
    if len( argv ) < 3:
        print("Usage: ",argv[0]+" TSPlib data file method")
        print("       where method is on the following ")
        print("       NN      for nearest neighbour")
        print("       NI      for nearest insertion")
        print("       FI      for farthest insertion")
        print("       MST     for min. spanning tree heuristic")
        print("       CHRISTO for Christofides' heuristic")
        print("       OTREE   for one tree lower bound")
        print("       FOTREE  for the fast one treee lower bound")
        print("       ASS     for the assignment lower bound")
        print("       2M      for the 2-matching lower bound")
        print("       F2M     for the fractional 2-matching lower bound")
        print("       LK      for Concorde's Lin-Kernighan heuristic")
        print("       CC      for Concorde's exact branch-and-cut method")
        print("       CPX     for Cplex MIP solver (without start solution")
        print("       CPX-XY  for Cplex MIP solver with start solution XY")
        print("               where XY is NN, NI, FI, MST or LK")
        exit()

    method = argv[2].upper()
    CPXanim = False
    nodeLim = None
    if 'CPX' in method:
        if len(argv) > 3: CPXanim = ('Y' in argv[3].upper() )
        if CPXanim: print("Ok. Trying with animation")
        __answer = input("Put limit on the number of nodes in search tree (y/n)? ").lower()
        if 'y' in __answer:
            __answer = input("Enter node limit:").lower()
            nodeLim = max(0,int(__answer))
            print("Using node limit:",nodeLim)


    # Obtain the problem instance in graph form
    problem = load_problem(argv[1])

    # Plot the points
    if problem.is_depictable(): displayPoints(problem)

    # Apply one of the selected heuristic solution methods
    if method in ['NN','NI','FI','MST','CHRISTO','LK']:
        routeLen, route = __heurSolve(problem, method=method, showResult=True)

    if method=='OTREE':
        # Obtain 1-tree lower bound
        weight, tree = oneTree( problem )
        plotEdgeList( problem, tree, title="1-tree bound = "+str(weight) )

    if method=='FOTREE':
        # Obtain fast 1-tree lower bound
        weight, tree = fast1Tree( problem )
        plotEdgeList( problem, tree,
                      title="Reinelt's fast 1-tree bound = "+str(weight) )

    if method=='ASS':
        # Obtain assignment lower bound
        cost, elist = assignLB( problem )
        plotEdgeList( problem, elist, title="Assignment lower bound = "+str(cost) )

    if method=='2M':
        # Obtain 2-matching lower bound
        weight, elist = twoMatch( problem )
        plotEdgeList( problem, elist, title="2-matching lower bound = "+str(weight) )

    if method=='F2M':
        # Obtain fractional 2-matching lower bound
        weight, elist, xvals = twoMatch( problem, fractional=True )
        efract = [ elist[i] for i in range(len(elist)) if xvals[i] < 0.999999 ]
        plotEdgeList( problem, elist, specialEdges=efract,\
          title="Fract. 2-matching LB (fract.edges are red) = "+str(weight) )

    if method =='CC':
        # Apply Concorde's exact TSP solver
        routeLen, route = solveTSP( problem, logFile="CC_log.log" )
        __showSol( problem, route, routeLen, "Exact solution (Concorde)")

    if 'CPX' in method:
        # Apply Cplex to obtain an (optimal) solution
        route = None
        if '-' in method:
            # Obtain a heuristic start solution
            routeLen, route = __heurSolve( problem, method.split('-')[1], showResult=CPXanim )
        lowbnd, routeLen, route = CPXtspSolve( problem, startRoute=route,\
                                               nodeLim=nodeLim, anim=CPXanim )
        if not route is None:
            header = "Cplex solution. Lower bound= "+str(lowbnd)
            __showSol(problem,route,routeLen,header)
        else:
            print("No feasible solution found. Lower bound is ",lowbnd)
