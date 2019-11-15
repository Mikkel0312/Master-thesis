from ConstructionCTP import *
import math
import pandas as pd
import random
import time



def update_NV_W(distM, dmax, NV, df_W, new_point):

    covering_df = (distM.loc[:,df_W.index]<dmax).loc[new_point]
    NV = NV - set(df_W[covering_df == True].index)
    return NV


def calculate_Pij(distM, dmax, NV_VT, NV_W, last_node_visited, df_T):
    '''
    Calculate the index Pij for each node

    parameters:
        distM: DataFrame
            the distance matrix for the problem.
        dmax: int
            the maximal covering range
        NV_VT: set
            a set with the unvisited vertices from V\T
        NV_W: set
            a set with the uncovered vertices from W
        last_node_visited: int
            the last node visited on the path
        df_T: DataFrame
            a DataFrame with the nodes that must be visited on the route.


    returns:

        df_Pij : Dataframe
            DataFrame with the Pij values for every node in NV_VT

    '''
    listPij = []
    for idx_i in NV_VT:
        d_ij = min(distM.loc[idx_i, last_node_visited], distM.loc[df_T.index, idx_i].min())
        sum_akj = sum((distM[NV_W]<=dmax).astype(int).loc[idx_i])
        if sum_akj != 0:
            listPij.append(d_ij/(math.sqrt(sum_akj)))
        else:
            listPij.append(math.inf)
    df_Pij = pd.DataFrame(listPij, index = NV_VT, columns = ["Pij"])
    return df_Pij

def chose_node(Pij, NV):
    """
    chose 1 node at random from the sqrt(len(NV_W)) nodes with lowest Pij index

    parameters:

        Pij: DataFrame
            DataFrame with a column with the name "Pij" with entries being the Pij values of the index node.
        param NV: int
            integer with the number of nodes a random node should be chosen from.

    returns:

        node: int
            the node that should be added to the route

    """
    number_of_nodes = math.ceil(math.sqrt(NV))
    possible_nodes = Pij.sort_values("Pij").head(number_of_nodes).index.values
    node = random.choice(possible_nodes)
    return node


def MDMC(df, T, V, reduce = True, plot = False):
    """
       Use the Minimal distance-maximal cover heuristic to solve the CTP

        Parameters:
            df : dataframe
               A dataframe with the coordinates for nodes of the problem
            T : int
               Number of nodes that the hamiltonian tour has to visit
            V : int
               number of nodes in the set V, the set V is the nodes that can be visited. W is implicitely defined as the len(df)-T-V point.
            reduce : bool
                if True(default) then the reduction rules is applied to the sets V and W.
            plot : bool
                if True(not default) then the solution is plotted.

        Returns:
            routeLen : int
               Length of the route.
            route : list of int
               Computed route (route[0] is the depot, route[-2] is the
               last customer) and route[-1] is the depot which ends the Hamiltonian cycle.
       """
    cpu_start = time.time()
    df_T = df[:T]
    df_VT = df[T:V]
    df_W = df[V:]
    distM = get_distM(df)
    dmax = max(max(distM.iloc[df_VT.index, df_W.index].min(axis=1)),
               distM.iloc[df_W.index, df_VT.index].apply(lambda x: x.nsmallest(2), axis=1).max(axis=1).max()) #calculate d_max
    A = get_A(distM.iloc[df_W.index, df_VT.index], dmax) #get covering matrix A
    #apply reduction rules to the sets W and V\T
    if reduce == True:
        A = reduce_W_VT(A, df_W, df_T, dmax)
        df_W = df_W.loc[A.index, :]
        df_VT = df_VT.loc[A.columns, :]

    #initialize the MDMC algorithm
    last_visited = 0
    route = []
    NV_W = set(df_W.index.values)
    NV_VT = set(df_VT.index.values)
    #remove sets from NV_W already covered by vertices from T.
    for idx, i in df_T.iterrows():
        NV_W = update_NV_W(distM, dmax, NV_W, df_W, idx) # this is redundant code when the reduction rules is used, since they already remove from W the verties covered by the vertices from T.
        route.append(idx)
    while(len(NV_W)>0): #while all the covering nodes is not covered
        Pij = calculate_Pij(distM, dmax, NV_VT, NV_W, last_visited, df_T) #calculate the Pij
        node = chose_node(Pij, len(NV_W)) # chosen 1 node at random from the sqrt(len(NV_W)) nodes with lowest Pij index
        route.append(node)
        last_visited = node
        #update the sets NV_W and NV_VT
        NV_W = update_NV_W(distM, dmax, NV_W, df_W, last_visited)
        NV_VT = NV_VT - {last_visited}

    #solve a tsp on the nodes to reoptimize the tour and to connect the ends of the tour.
    routeLen, route = solveTSP(distM.loc[route,route], df_T)

    if plot:
        plot_solution(df, df_VT, df_T, df_W, route, dmax) #plots the solution
    cpu_time = time.time()-cpu_start
    return routeLen, route, cpu_time



