"""
    Python interface to Concorde - A C-library for solving TSPs.

    author:  
        Andreas Klose 
    version: 
        1.0 (29/3/2019)
    
    Remark: 
        Concorde is a TSP solver written in C and created by D. Applegate, 
        R. Bixby, V. Chvatal and W. Cook.  The code is available from
        http://www.math.uwaterloo.ca/tsp/concorde and permission to use it 
        is granted for academic purposes only.
"""
import numpy as np
import ctypes, pathlib, platform, site, time
from tsplib95 import load_problem
import os, sys

__CC_Lib = None
__CPX_Lib = None

#---------------------------------------------------------------------

def __loadConcorde():
    """ 
    Find and load cplex and concorde library.
    """    

    cplexLib = None
    concordeLib = None

    onWindows = 'Windows' in platform.system()
    onMac = 'darwin' in platform.system().lower()
    if onWindows:
        cplexDLL = 'cplex12??.dll'
        concordeDLL = 'libconcorde.dll'

    else:
        cplexDLL = '**/py3?_cplex????.so'
        concordeDLL = 'libconcorde.so'
   
    # Check if the CPLEX callable library is present. On Linux and Mac, we
    # will also have to explictly load this library before Concorde is
    # is loaded, as Concorde depends on Cplex. On Windows, it suffices
    # to load Concorde as the the Concorde Library has been linked with 
    # the Cplex static library. On Linux and Mac, we have to search for the 
    # file "py_3?_cplex????.so"; on Windows, we search cplex????.dll

    if onWindows:
        try:
            cpxEnv = next( p for p in os.environ if 'CPLEX_STUDIO_BINARIES' in p )
        except:
            cpxEnv = 'PATH'
        cpxBins = os.getenv(cpxEnv)
        dirLst = [ p for p in cpxBins.split(';') if 'cplex' in p.lower()]            
    else:    
        usrSitePack = next( p for p in sys.path if 'site-packages' in p)
        dirLst = site.getsitepackages()
        dirLst.append( usrSitePack )
    print(dirLst)
    for direc in dirLst:
        for pathName in pathlib.Path( direc ).rglob( cplexDLL ):
            cplexLib = str(pathlib.Path.resolve(pathName))
            if not cplexLib is None: break
        if not cplexLib is None: break    

    if cplexLib is None:
        print("Could not find CPLEX dynamic link library.")
        return( False )   

    # The dynamic link library libconcorde.so (libconcorde.dll on Windows)
    # is expected to reside in the same directory as this module "concorde.py"
    if onWindows:
        prefix = 'win//'
    elif onMac:
        prefix = 'mac/'
    else:
        prefix = 'linux/'
    concordeHome=pathlib.Path(__file__).resolve().parent
    concordeFile = pathlib.PurePath.joinpath(concordeHome,prefix+concordeDLL)
    if concordeFile.is_file():
        concordeLib = str( concordeFile )
    else:
        print("Could not find Concorde's dynamic link library.")
        return( False )

    # Libraries found. Now try to load them. Cplex first, as Concorde
    # depends on Cplex
    global __CPX_Lib
    global __CC_Lib
    if not onWindows:
        __CPX_Lib = ctypes.CDLL(cplexLib, mode=ctypes.RTLD_GLOBAL)
        if __CPX_Lib is None:
            return( False )
    
    __CC_Lib = ctypes.CDLL(concordeLib, mode=ctypes.RTLD_GLOBAL)
    if __CC_Lib is None:
        print("Error occured when loading library ",concordeLib )
        return( False )

    return( True )
    
#---------------------------------------------------------------------

def solveTSP( problem, route=None, exact=True, logFile=None ):
    """
    Invokes Concorde's TSP solver on the TSP instance "problem".

    Parameters:
        problem : class
            TSP problem object has returned by TSPlib95.
        route : list of int
            If not none it should contain a permutation list of cities
            to be used as initial solution.
        exact : bool
            If true (default) it is tried to solve the instance exactly 
            by means of Concorde's branch-and-cut solver; otherwise 
            Concorde only applies the Lin-Kernighan heuristic.
        logFile : str
            If not None, it need to be the name of file, where to 
            redirect output from Concorde

    Returns:
        routeLen : int 
            Length of the route.
        route : list of int
            Computed route (route[0] is the depot and route[-1] is the
            last customer).
    """
    if __CC_Lib is None:
        print("Concorde Library not loaded!")
        return 0, []

    nn = problem.dimension
    nodeLst = [node for node in problem.get_nodes()]

    n = ctypes.c_int( nn )
    seed = ctypes.c_int( int( time.time() ) )
    tiLim = ctypes.c_double(0.0)
    LKonly = ctypes.c_char(1-int(exact))
    
    # Compute the distance matrix
    dim  = nn*(nn-1)//2
    dist = np.zeros( dim, dtype=ctypes.c_int )
    cnt  = 0
    for i in range(1,nn):
        for j in range(i):
            dist[cnt] = problem.wfunc(nodeLst[i],nodeLst[j])
            cnt += 1
    pdist = dist.ctypes.data_as( ctypes.POINTER( ctypes.c_int ) )        

    # Number the nodes from 0 to n-1
    nodeIdx = {}
    cnt = 0
    for i in nodeLst: 
        nodeIdx[i] = cnt
        cnt += 1    

    # Redirect output from Concorde?
    if logFile is None:
        logPtr = ctypes.c_char_p(0)
    else:
        logPtr = ctypes.c_char_p(logFile.encode('utf-8'))
        old_out = sys.stdout # saver when on windows
    
    # Create integer array representing the tour
    tour = np.zeros(nn, dtype=ctypes.c_int)
    ptour= tour.ctypes.data_as( ctypes.POINTER( ctypes.c_int ) )
    if route is None:
        tLen = ctypes.c_double( 0.0 )
    else:
        tLen = ctypes.c_double( 1.0 )
        cnt = 0
        for i in route:
            if ( cnt < nn ): tour[cnt] = nodeIdx[i]
            cnt += 1
 
    # Call concorde for computing TSP tour
    __CC_Lib.solve_STSP.restype = ctypes.c_int
    status = __CC_Lib.solve_STSP( LKonly, n, seed, tiLim, pdist,\
                 logPtr, ctypes.byref(tLen), ptour )
    
    # Following is safer when on Windows
    if not logFile is None: sys.stdout = old_out
            
    if status < 2:
        routeLen = int(tLen.value)
        route = [ nodeLst[i] for i in tour ]
        return routeLen, route
    else:
        return np.Inf, []


#---------------------------------------------------------------------

def solveTSPdat( n, GRAPH=None, dist=None, symmetric=True, ROUTE=None,\
                 exact=True, logFile = None ):
    """
    Invokes Concorde's TSP solver.

    The problem instance can be specified in three different ways.
    
    (1)
    As a graph with n nodes. In this case, GRAPH must not be None
    and should be a tuple (m, elist, elen), where m is the number
    of edges, elist the list of edges and elen the list of edge
    lengths.

    The list "elist" of edges need to be a numpy array of length 
    2*m containing C integers (ctypes.c_int). For each edge e, 
    (e=1,...,m), the end points of the edge need to be given by 
    elist[2(e-1)] and elist[2e-1]. Moreover, elen need to be a 
    numpy array of length m containing C integers (ctypes.c_int), 
    such that elen[e-1] gives the length of edge e=1,...,m.

    (2)
    Using a full nxn matrix of distances. In this case, GRAPH need
    to be None and symmetric be false. Moreover, dist need to
    be a numpy array of length nx(n-1) containing C integers
    (ctypes.c_int) that stores the distance matrix (excluding the 
    diagonal) row by row.

    (3)
    Using the lower triangle part of a symmetric distance matrix.
    In this case GRAPH need to be None and symmtric be true, and 
    dist an numpy array containing C integers (ctypes.c_int) such
    that dist[(i-1)*i/2 + j] gives the (integer-valued) distance 
    from node i to node j (i,j = 0,...,n-1).

    Parameters:
        n : int
            Number of nodes
        GRAPH : triple
            Tuple (m, elist, elen) specifying number of edges, 
            edge list and edge lengths
        dist : numpy array of c_int
            Full distance matrix or lower triangle part of a 
            symmetric matrix.
        symmetric : bool
            True if distance matrix is symmetric (default).
        ROUTE : tuple of two
            A tuple (r, rlen). If not None, r is an initial route 
            (given by the sequence/permutation of nodes and rlen 
            is the route length.
        exact : bool
            If true (default) it is tried to solve the instance
            exactly by means of Concorde's branch-and-cut solver;
            otherwise Concorde only applies the Lin-Kernighan
            heuristic.
        logFile : str
            If not None, it need to be the name of file, where to 
            redirect output from Concorde

    Returns:
        routeLen : int
            Length of the route (routeLen).
        route : list of int
            Computed route (route[0] is the depot and route[-1] is the
            last customer).
    """
    if __CC_Lib is None:
        print("Concorde Library not loaded!")
        return 0, []
        
    if ( GRAPH is None) and (dist is None): return 0, []

    if GRAPH is None:
        elen = dist
        if symmetric:
            m = n*(n-1)//2
            elist = np.zeros( 2*m, dtype=ctypes.c_int )
            cnt = 0
            for i in range(1,n):
                for j in range(i):
                    elist[cnt], elist[cnt+1] = j, i
                    cnt += 2
        else:
            m = n*(n-1)
            elist = np.zeros( 2*m, dtype=ctypes.c_int )
            cnt = 0
            for i in range(n):
                for j in range(n):
                    if i != j: 
                        elist[cnt], elist[cnt+1] = i, j
                        cnt += 2
    else:                          
        m = GRAPH[0];
        elist = GRAPH[1]
        elen = GRAPH[2]

    p_elen = elen.ctypes.data_as( ctypes.POINTER( ctypes.c_int ) )        
    p_elist = elist.ctypes.data_as( ctypes.POINTER( ctypes.c_int ) )        
    tour = np.zeros(n, dtype=ctypes.c_int) 
    ptour= tour.ctypes.data_as( ctypes.POINTER( ctypes.c_int ) )
    nn = ctypes.c_int( n )
    mm = ctypes.c_int( m )
    seed = ctypes.c_int( int( time.time() ) )
    tiLim = ctypes.c_double(0.0)
    LKonly = ctypes.c_char(1-int(exact))

    minIndx = 0    
    if ROUTE is None:
        tLen = ctypes.c_double( 0.0 )
    else:
        tLen = ctypes.c_double( ROUTE[1] )
        cnt = 0
        minIndx = min( ROUTE[0] )
        for i in range(n):
            tour[cnt] = ROUTE[0][i]-minIndx
            cnt += 1

    # Redirect output from Concorde?
    if logFile is None:
        logPtr = ctypes.c_char_p(0)
    else:
        logPtr = ctypes.c_char_p(logFile.encode('utf-8'))
        old_out = sys.stdout 
    
    # Call concorde for computing TSP tour
    __CC_Lib.solve_STSP.restype = ctypes.c_int
    status = __CC_Lib.solve_sparseTSP( LKonly, nn, mm, seed, tiLim,\
                 p_elist, p_elen, logPtr, ctypes.byref(tLen), ptour )
        
    # Following is safer when on Windows
    if not logFile is None: sys.stdout = old_out
            
    if status < 2:
        routeLen = int(tLen.value)
        route = [ i+minIndx for i in tour ]
        return routeLen, route
    else:
        return np.Inf, []
       
#---------------------------------------------------------------------

def solveTSPLib( fname, exact=True, logFile=None ):
    """
    Solve a TSPLIB instance stored in file "fname" by means of
    Concorde's TSP solver.

    Parameters:
        fname : str
            Name (including) the path to the data file (format TSPLIB).
        exact : bool
            If true (default) it is tried to solve the instance exactly
            by means of Concorde's branch-and-cut solver; otherwise 
            Concorde only applies the Lin-Kernighan heuristic.
        logFile : str
            If not None, it need to be the name of file, where to 
            redirect output from Concorde.

    Returns:
        routeLen : int
            Length of the route (routeLen).
        route : list of int
            Computed route (route[0] is the depot and route[-1] is the
            last customer).
    """

    if __CC_Lib is None:
        print("Concorde Library not loaded!")
        return 0, []
    else:
        # We first create a pointer to an integer array
        problem = load_problem(fname)
        tour = np.zeros(problem.dimension, dtype=ctypes.c_int)
        iptr = ctypes.POINTER( ctypes.c_int )
        ptour= tour.ctypes.data_as( iptr )
        # Redirect output from Concorde?
        if logFile is None:
            logPtr = ctypes.c_char_p(0)
        else:
            logPtr = ctypes.c_char_p(logFile.encode('utf-8'))
            old_out = sys.stdout
        
        # Initialize other parameters of c-function solve_TSLPlib
        seed    = ctypes.c_int( int( time.time() ) )
        status  = ctypes.c_int(0)
        tiLim   = ctypes.c_double(0.0)
        n       = ctypes.c_int(0)
        fnmeptr = ctypes.c_char_p(fname.encode('utf-8'))
        LKonly  = ctypes.c_char(1-int(exact))
        __CC_Lib.solve_TSPlib.restype = ctypes.c_double
        tLen    = __CC_Lib.solve_TSPlib( LKonly, fnmeptr, seed, tiLim,\
                      logPtr, ptour, ctypes.byref(status) );
        routeLen = tLen
        nodeLst = [node for node in problem.get_nodes()]
        route = [ nodeLst[i] for i in tour ]
        
        # Following is safer when on Windows
        if not logFile is None: sys.stdout = old_out
            
        return routeLen, route

#---------------------------------------------------------------------

is_ready = __loadConcorde()

