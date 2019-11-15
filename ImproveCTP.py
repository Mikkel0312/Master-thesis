import time
from ConstructionCTP import *

def improvement_heuristic(A, distM, df_T, df_VT, route, b, n=1, use_method=1, l=1):

    '''This function removes a node from the route, calculates the insertion costs and solves a SCP and a TSP to find a new CTP solution.

    parameters:

        A: DataFrame
            the covering matrixdist
        M: DataFrame
            the distance matrix for the problem
        df_T: DataFrame
            the nodes that must be visited
        df_VT: DataFrame
            the nodes that can be visited
        route: list of ints
            the route that is currently improved
        n: int (default 1)
            number of nodes to remove from the route (denoted alpha in the thesis)
        b: int (default 1)
            the b’th savings cost to use
        use_method: int (default 1)
            an integer to determine the method to determine the n nodes that should be remove.
                1 : savings cost is used
                2 : a sorted method is used, so it removes the l’th node from the route
                3 : a random node is removed from the route
        l: int (default 1)
            a number to determine which position a node should be removed from in the "sorted" method.

    returns:
        routeLen: int
            route length of the new (improved) solution
        route: list of ints
            the nodes that is in the new route

    ''' #route = best_route
    if use_method == 1: #using savings_cost to remove point(s)
        saving_costs = get_savings_costs(route[:-1], distM, df_T)
        reduced_route = reduce_route_savings(route,saving_costs, n, b)
    if use_method == 2:  #remove l'th (+(l+n)'th point(s) from solution
        reduced_route = reduce_route_l(route,n,l)
    if use_method == 3: #remove random point(s) from solution
        reduced_route = reduce_route_random(route, n)

    weights = get_insertion_costs(distM, df_VT, route, reduced_route, n)

    solution = solve_setcover(A, weights['cost'], None) #route_NoDepot)

    distM_tsp = get_tsp_distM(distM, df_T, solution)


    routeLen, route = solveTSP(distM_tsp, df_T, exact = True)
    #if plot:
    #plot_solution(df, df_VT, df_T, df_W, route, dmax)
    return routeLen, route


def get_savings_costs(route, distM,df_T):
    '''
    gets the saving costs of the nodes on the route. If the nodes is in T, then the saving cost is fixed to 0

    parameters:

        route: list of ints
            the route with a list of nodes the savings cost should be calculated for
        distM: DataFrame
            the distance matrix for the problem.
        df_T:the nodes that must be visited

    returns:

        savings_cost : DataFrame
            DataFrame of len(route) with the saving costs

    '''

    list1 = []
    for i in range(len(route)):
        if i in list(df_T.index.values):
            c_i = 0
        else:
            c_i = distM.iloc[int(route[i-1])-1 , int(route[(i + 1) % len(route)])-1] \
                  - distM.iloc[int(route[i - 1])-1 , int(route[i])-1] \
                  - distM.iloc[int(route[i])-1, int(route[(i + 1) % len(route)])-1]
        list1.append(c_i)

    savings_cost = pd.DataFrame(list1, columns= ["saving cost"])
    savings_cost = savings_cost.set_index(pd.Index(route))
    return savings_cost

def get_insertion_costs(distM, df_VT, route, reduced_route, n):
    '''
    gets the insertion costs of the nodes from V\T.

    parameters:
        distM: DataFrame
            the distance matrix for the problem
        df_VT: DataFrame
            the nodes that can be visited
        route: list of ints
            the route that is currently improved
        reduced_route: list
            the partial route with n nodes removed
        n: int
            number of nodes which is removed from the tour(denoted alpha in the thesis)

    returns:

        insertion_cost : DataFrame
            DataFrame with len(df_VT) with the insertion costs.
    '''
    list1 = []
    here_index = reduced_route.index("here")
    #n_closest = get_n2_closest_tour_points(distM, df_VT, route, n)
    for idx, i in df_VT.iterrows():
        if int(idx) in reduced_route:
            c_i = 0.001
            #distM.iloc[idx,n_closest.loc[idx,0]] \
            #      + distM.iloc[idx,n_closest.loc[idx,1]] \
            #      - distM.iloc[n_closest.loc[idx,0],n_closest.loc[idx,1]]
        elif idx in route and idx not in reduced_route:
            c_i = 99999999
        else:
            c_i = distM.iloc[int(idx)-1, int(reduced_route[here_index - 1])-1]\
                  + distM.iloc[int(idx)-1, int(reduced_route[(here_index + 1) % len(reduced_route)])-1] \
                  - distM.iloc[int(reduced_route[here_index - 1])-1, int(reduced_route[(here_index + 1) % len(reduced_route)])-1]

        list1.append(c_i)

    insertion_cost = pd.DataFrame(list1, columns= ["cost"], index = df_VT.index)

    return insertion_cost

def get_n2_closest_tour_points(distM, df_VT, route, n):
    '''
    returns the m+2 closest nodes for each node of V\T to the nodes on the route.

    parameters:

        distM: DataFrame
            the distance matrix for the problem
        df_VT: DataFrame
            the nodes that can be visited
        route: list of ints
            the route from which the closest points should be calculated, route[0] is the depot and route[-1] is the last node on the tour
        m: int
            number of closests nodes that should be calculated

    returns:

        df : DataFrame
            a Dataframe with m+2 columns and len(df_VT) rows. in the first column is the closests node for the given node. In the second column is the second closest node for the given node. etc.
    '''
    #returns the n+2 closest nodes for each node.
    list1 = []
    for idx, i in df_VT.iterrows():
        list1.append(distM[route].loc[idx].nsmallest(n+2).index.values)
    df = pd.DataFrame(list1)
    df = df.set_index(df_VT.index)
    return df


def reduce_route_savings(route, saving_costs, n=1, b=1):
    '''
    remove n succesive nodes from the route at the node with the b’th highest savings cost

    parameters:

        route : list of ints
            route which a node should be removed from. route[0] = route[-1] is the depot and route[-2] is the last node on the tour.
        savings_cost : DataFrame
            DataFrame with len(df_VT) with the savings cost of all the nodes that can be visited. The column name of the column should be "saving cost".
        n : int (default = 1)
            number of nodes to remove from the tour
        b : int (default = 1)
            the node with the b’th highest saving cost is chosen for removal.

    returns:

        route1: list of ints
            a route with n nodes with the highest saving cost
    '''
    list_s = []
    for i in range(len(route)):
        saving_i=0
        for j in range(n):
            saving_i += saving_costs.loc[route[(i + j) % len(route)]]["saving cost"]
        list_s.append(saving_i)

    max_saving_node = pd.DataFrame(list_s).nsmallest(b,0).tail(1)
    route1 = route.copy()
    if n == 1:
        route1.remove(route[max_saving_node.index.values[0]])
    if n == 2:
        route1.remove(route[max_saving_node.index.values[0]])
        route1.remove(route[(max_saving_node.index.values[0] + 1)% len(route)] )
    if n == 3:
        route1.remove(route[max_saving_node.index.values[0]])
        route1.remove(route[(max_saving_node.index.values[0] + 1)% len(route)])
        route1.remove(route[(max_saving_node.index.values[0] + 2)% len(route)])

    route1.insert(max_saving_node.index.values[0], "here")
    return route1


def reduce_route_random(route, n):
    '''
    This method was used to see if chosen a random node instead of the node with b'th highest saving cost had an impact


    '''
    randomint = random.randint(1,len(route)-2) #exclude the depot
    route1 = route.copy()
    if n == 1:
        route1.remove(route[randomint])
    if n == 2:
        route1.remove(route[randomint])
        route1.remove(route[(randomint + 1)% len(route)])
    if n == 3:
        route1.remove(route[randomint])
        route1.remove(route[(randomint + 1) % len(route)])
        route1.remove(route[(randomint + 2) % len(route)])
    route1.insert(randomint, "here")
    return route1

def reduce_route_l(route,n,l):
    '''
    Method used to check all the positions in the route for local improvement.
    '''
    route1 = route.copy()
    if n == 1:
        route1.remove(route[l])
    if n == 2:
        route1.remove(route[l])
        route1.remove(route[(l + 1) % len(route)])
    if n == 3:
        route1.remove(route[l])
        route1.remove(route[(l + 1) % len(route)])
        route1.remove(route[(l + 2) % len(route)])
    route1.insert(l, "here")
    return route1



#current_obj, route, best_obj, best_route, b, n, method, maxiters, plot1 = val, sol,val, sol, 1, 1, "savings", 10, False
def ImproveCTP(df, T, V, best_route, best_obj, n = 1, method = "savings", maxiters=10, plot = False):
    '''
    This is the function that invokes the ImproveCTP optimization heuristic.

    parameters:

        df: DataFrame
            DataFrame with the coordinates of the points in the problem
        T: int
            number of nodes that must be visited
        V: int
            number of nodes that can be visited
        best_route: list of ints
            the initial route, route[0] = route[-1] is the depot and route[-2] is the last node on the tour.
        best_obj: int
            objective of the initial solution
        n: int (default 1)
            this is number alpha of nodes to remove from the tour
        method: str (default "savings")
            this specifies how the node to remove should be selected, the options are:
                "savings" : calculating the savings cost of each node and select the n successive nodes with the highest savings cost.
                "sorted" : selects the n first nodes on the tour, then the n second↪→first nodes and so on.
                "random" : selects n random nodes on the route to remove.
        maxiters: int
            maximum number of iterations without improvement
        plot: bool (default False)
            if True, then the route obtained is plotted

    returns:
        best_obj : int
            length of the route obtained
        best_route : list of ints
            list of the nodes to visit on the best route obtained.best_route[0] = best_route[-1] is the depot and↪→best_route[-2] is the last node on the tour.
        cpu_time : floatnumber of seconds for the algorithm to run

    '''
    cpu_start = time.time()
    df_T = df[:T]
    df_V = df[:V]
    df_VT = df[T:V]
    df_W = df[V:]
    distM = get_distM(df)
    dmax = max(max(distM.loc[df_V.index.values, df_W.index.values].min(axis=1)),
               distM.loc[df_W.index.values, df_V.index.values].apply(lambda x: x.nsmallest(2), axis=1).max(
                   axis=1).max())
    A = get_A(distM.loc[df_W.index.values, df_VT.index.values], dmax)
    #A = reduce_W_VT(get_A(distM.loc[df_W.index.values, df_VT.index.values], dmax), df_W, df_T, dmax)
    df_W_r = df_W.loc[A.index, :]
    df_VT_r = df_VT.loc[A.columns, :]

    if method == "savings":
        u=0
        if n != None:
            while (u<maxiters and u<len(best_route)-n):
                u+=1
                current_obj , route = improvement_heuristic(A, distM, df_T, df_VT_r, best_route, b=u, n=n, use_method=1)
                print(current_obj, best_obj)
                if current_obj < best_obj:
                    best_obj = current_obj
                    best_route = route
                    u = 0

        else:
            for n in [1,2,3]:
                u = 0
                while (u<maxiters and u<len(best_route)-n):
                    u+=1
                    current_obj , route = improvement_heuristic(A, distM, df_T, df_VT_r, best_route, b=u, n=n, use_method=1)
                    if current_obj < best_obj:
                        best_obj = current_obj
                        best_route = route[:-1]
                        u = 0

    # if method == "sorted":
    #     start = time.time()
    #     while (start < start + 120):
    #         j=1
    #         for i in range(1, len(best_route)-2):
    #             current_obj, route = improvement_heuristic(A, distM, df_T, df_VT_r, df_W_r, dmax, best_route, n, use_method=2, l=i, exact= False)
    #             if current_obj < best_obj:
    #                 best_obj = current_obj
    #                 best_route = route[:-1]
    #                 j=0
    #                 break
    #         if j==1:
    #             break
    #
    # if method == "random":
    #     k = 0
    #     while (k<randomPoints):
    #         current_obj, route = improvement_heuristic(df, nearestN, best_route, n, 0)
    #         k+=1
    #         if current_obj < best_obj:
    #             best_obj = current_obj
    #             best_route = route
    #             k=0
    #             continue

    if plot:
        plot_solution(df, df_VT_r, df_T, df_W_r, route, dmax)
    print("Objective value of the Covering Tour is " + str(best_obj))

    return best_obj, best_route, time.time()-cpu_start
