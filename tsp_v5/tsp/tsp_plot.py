"""
    TSP - plotting routines.

    author:
        Andreas Klose
    version:
        1.0 (29/3/2019)
"""
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

figNum = 0

#--------------------------------------------------------------------- 

def plotEdge( A, B, special=False ):
    """
    Plot line connecting coordinate points A and B.
    
    Parameters:
        A : list of float
            Coordinates of point A
        B : list of float
            Coordinates of point B
        special : bool, optional
            If true, a red line otherwise a blue line is drawn.
            Default=False.
    """
    
    if not special:
        plt.plot( [A[0],B[0]], [A[1],B[1]], 'b-')
    else:
        plt.plot( [A[0],B[0]], [A[1],B[1]], 'r-')

#---------------------------------------------------------------------

def setFigLims( ax, X, Y ):
    """
    Set bounding box for a figure.
    
    Parameters:
        ax : object (matplotlib)
            axis object.
        X  : list of X-coordinates.
        Y  : list of Y-coordinates.   
    """
    xmin, ymin = min(X), min(Y)
    xmax, ymax = max(X), max(Y)
    xsep, ysep = 0.1*(xmax-xmin), 0.1*(ymax-ymin)
    ax.set_xlim(xmin-xsep, xmax+xsep)
    ax.set_ylim(ymin-ysep, ymax+ysep)

#---------------------------------------------------------------------

def displayPoints( problem ):
    """
    Display nodes of the TSP problem.
    
    Parameters:
        problem : object
            TSP problem object as returned by function load_problem
            of the package tsplib95.
    """
    global figNum
    figNum += 1
    fig = plt.figure(figNum)
    ax  = fig.subplots()
    X = [ problem.node_coords[i][0] for i in problem.get_nodes()]
    Y = [ problem.node_coords[i][1] for i in problem.get_nodes()]
    setFigLims( ax, X, Y )
    plt.plot(X,Y,'bo',markersize=2)
    for i in problem.get_nodes():
        plt.text( problem.node_coords[i][0],problem.node_coords[i][1],\
                  str(i),fontsize=6)
    plt.title("Node locations")
    plt.show()

#---------------------------------------------------------------------

def displayRoute( coords, route, routeLen, title=None, animate=True ):
    """ 
    Display (animated) route over nodes.
    
    Parameters:
        coords : list of coordinates
            Node coordinates (attribute problem.node_coords of class
            problem in the package tsplib95).
        route : list of int
            The route to be displayed.
        routeLen : int
            The length of the route.
        title : str
            String used as a header for the figure.
        animate : bool, optional
            If true, the default, the route is displayed in animated way.                
    """
    X = [ coords[i][0] for i in route]
    Y = [ coords[i][1] for i in route]
    X.append( X[0] )
    Y.append( Y[0] )

    global figNum
    figNum += 1
    fig = plt.figure(figNum)
    ax  = fig.subplots()
    ln, = plt.plot([], [], 'b-', markersize=2,animated=True)
    plt.plot(X,Y,'bo',markersize=2)
    for i in range(len(route)):
        plt.text(X[i],Y[i],str(route[i]),fontsize=6)

    def init():
        setFigLims( ax, X, Y )
        return ln,

    def update(j):
        ln.set_data(X[0:j], Y[0:j])
        return ln,

    if title != None:
       plt.title(title+"(Length="+str(routeLen)+")")
    if animate:
        anim = FuncAnimation(fig, update, frames=range(1,len(route)+2),\
                   init_func=init, blit=True, repeat=False)
    else:
        plt.plot( X, Y, 'b-')
    plt.show()

#---------------------------------------------------------------------

def plotEdgeList( problem, edges, specialNodes=None, specialEdges=None,\
                  title=None ):
    """"
    Displays a graph of a list of edges.
    
    Parameters:
        problem : class
            The problem object as returned by function load_problem of
            package tsplib95.
        edges : list of pairs of int
            The edges to be displayed.
        specialNodes : list of int
            A list of nodes seen as special and drawn in red color.
            Default=None.
        specialEdges : list of pairs of int
            List of special edges to be drawn in red color. Default=None.
        title : str
            String used as header for the figure. Default=None.                
    """
    global figNum
    figNum += 1
    fig = plt.figure(figNum)
    ax = fig.subplots()
    X = [ problem.node_coords[i][0] for i in problem.get_nodes()]
    Y = [ problem.node_coords[i][1] for i in problem.get_nodes()]
    setFigLims( ax, X, Y )
    plt.plot(X,Y,'bo',markersize=2)
    for e in edges:        
        plotEdge( problem.node_coords[e[0]], problem.node_coords[e[1]] )
    if specialEdges != None:
        for e in specialEdges:
            plotEdge( problem.node_coords[e[0]], problem.node_coords[e[1]], special=True )
    if specialNodes != None:
        XS = [problem.node_coords[i][0] for i in specialNodes]
        YS = [problem.node_coords[i][1] for i in specialNodes]
        plt.plot(XS,YS,'ro',markersize=4)
    for i in problem.get_nodes():
        plt.text( problem.node_coords[i][0],problem.node_coords[i][1],\
                  str(i),fontsize=6)

    if title != None: plt.title(title)
    plt.show()
