import pandas as pd
import random
from scipy.spatial import distance_matrix
from docplex.mp.model import Model
import numpy as np
from tsp_v5.tsp.concorde import solveTSPdat
import matplotlib.pyplot as plt
import math
from sklearn import linear_model
import time
import ctypes


figNum = 0


def generate_points(n, seed = None): #generate n points in [0,100]^2 square and a depot in [25,75]^2 square.
    '''
    generate n points in the [0,100]^2 square
    parameters:

        n: int
            number of points to generate
        seed: int (default None)
            to specify a seed that should be used to generate the point
    returns:

        df : DataFrame
            a DataFrame with coordinates of the n points
    '''

    if seed != None:
        random.seed(seed)
    x = []
    y = []
    for i in range(n):
        x.append(random.randint(0,100))
        y.append(random.randint(0,100))
    df = pd.DataFrame(
        {'x' : x,
         'y' : y}
    )


    return df


def plot_points(df, T, VT, W, dmax):
    '''
    plot the points in the problem

    parameters:

        df: DataFrame
            DataFrame with all the points
        T: int
            number of points that must be visited
        VT: int
            number of points that can be visited and no in T
        W: int
            number of covering nodes
        dmax: int
            the covering distance

    returns:

    '''
    #initially plots the points
    global figNum

    figNum += 1
    fig = plt.figure(figNum)
    plt.scatter(VT['x'],VT['y'], c = 'b')
    plt.scatter(T['x'], T['y'], c = 'r', marker = 's')
    plt.scatter(W['x'], W['y'], c='g', marker = 'p')
    ax = fig.gca()
    for i in range(len(VT)+len(T)):
        circle = plt.Circle((df['x'].iloc[i],df['y'].iloc[i]), radius = dmax, alpha = 0.1)
        ax.add_artist(circle)

    plt.show()

def get_distM(df):
    '''
    get the distance matrix of the points in the problem

    parameters:
        df: DataFrame
            the points in the problem

    return:
        distM : DataFrame of with len(df) columns and len(df) rows
            the distance matrix.
    '''
    distM = pd.DataFrame(distance_matrix(df[['x', 'y']], df[['x', 'y']]).round())
    distM.astype(int)
    distM = distM.set_index(df.index.values) #this and the below line is done to set the correct index, if the input is not a dataframe with an ordered index.
    distM.columns = df.index.values
    return distM

def get_A(distM, dmax):
    '''
    Construct the covering (0-1) matrix as a A = [a_ij], where a_ij = 1 if and only if d_ij<dmax

    parameters:

        distM: DataFrame
            distance matrix of the problem
        dmax: int
            the covering distance

    returns:
        A : DataFrame
            A 0-1 DataFrame which specifies what nodes from W is covered by nodes from V\T
    '''
    A = pd.DataFrame((distM<=dmax).astype(int))
    return A


def reduction_rule1(A, df_W, df_T, dmax):
    '''
    Applies reduction rule 1 to the set W, reduction rule 1: Remove i from W if there exists j \in T with A_ij = 1

    parameters:
        A: DataFrame
            the covering matrix
        df_W: DataFrame
            the vertices that must be covered
        df_T: DataFrame
            the vertices that must be visited
        dmax: int
            the covering distance

    returns:
        A: DataFrame
            A DataFrame where the rows has been reduced according to reduction rule 1.
    '''
    #Remove i from W if there exists j \in T with A_ij = 1
    dist = pd.DataFrame(distance_matrix(df_W, df_T))
    dist = dist.set_index(df_W.index)
    dist.columns = df_T.index
    for index, i in A.iterrows():
        for index1, j in df_T.iterrows():
            if ((dist.loc[index] < dmax).iloc[0]):
                A = A.drop(index)
                break
    return A

def reduction_rule2(A):
    '''
    Applies reduction rule 2 to the set W, reduction rule 2: remove from W any vertex i for which A_ij = 1 for all j in V\T

    parameters:
        A: DataFrame
            the covering matrix

    returns:
        A: DataFrame
            A DataFrame where the rows has been reduced according to reduction rule 2.
    '''
    for index, i in A.iterrows():
        if sum(i) == len(A.columns):
            A = A.drop(index)
    return A

def reduction_rule3(A):
    '''
    Applies reduction rule 3 to the set W, reduction rule 3: If several rows were identical, then only one of them was kept.

    parameters:
        A: DataFrame
            the covering matrix

    returns:
        A: DataFrame
            A DataFrame where the rows has been reduced according to reduction rule 3.
    '''

    #
    index_list = []
    for index, i in A.iterrows():
        for index1, j in A.iterrows():
            if (all(i == j) and index != index1):
                A = A.drop(index)
                index_list.append(index)
                break
    return A


def reduction_rule4(A):
    '''
    Applies reduction rule 4 to the set W, reduction rule 4: Remove dominated rows. Remove from W any vertex i for which there exists a vertex j != i, j in W s.t A_ik <= A_jk forall k \in V\T

    parameters:
        A: DataFrame
            the covering matrix

    returns:
        A: DataFrame
            A DataFrame where the rows has been reduced according to reduction rule 4.
    '''

    #
    index_list = []
    for index, i in A.iterrows():
        for index1, j in A.iterrows():
            if (all(i >= j) and index != index1):
                A = A.drop(index)
                index_list.append(index)
                break
    return A



def reduction_rule5(A):
    #
    '''
    Applies reduction rule 5 to the set V\T, reduction rule 5: remove i from V\T if a_ji = 0 forall j \in W

    parameters:
        A: DataFrame
            the covering matrix

    returns:
        A: DataFrame
            A DataFrame where the columns has been reduced according to reduction rule 5.
    '''
    for index, j in A.iteritems():
        if sum(j) == 0:
            A = A.drop(index, axis=1)

    return A


def reduce_W_VT(A, df_W, df_T, dmax):
    '''
    Applies reduction rule 1,2,3,4 and 5 to the sets W and V\T.

    parameters:
        A: DataFrame
            the covering matrix
        df_W: DataFrame
            the nodes that must be covered
        df_T: DataFrame
            the nodes that must be visited
        dmax: int
            the covering distance

    returns:
        A: DataFrame
            A DataFrame where the rows and columns has been reduced according to reduction rule 1,2,3,4 and 5.
    '''
    #iteratively applies all the reduction rules
    df_W_reduced = df_W.loc[A.index, :]
    A = reduction_rule1(A, df_W_reduced, df_T, dmax)
    A = reduction_rule2(A)
    A = reduction_rule3(A)
    A = reduction_rule4(A)
    A = reduction_rule5(A)
    return A


def solve_unweigthed_setcover(A):

    '''
    Solves an unweighted set covering problem based on the covering matrix A.

    parameters:
        A: DataFrame
            the covering matrix

    returns:
        obj: int
            the objective of the set cover, which is the number of nodes included in the solution
        df_solution: DataFrame
            a DataFrame of length len(A.columns) with the solution values of the vertices in V\T.
    '''
    #
    #initialize
    mdl = Model(name = 'Setcover')
    #model paramters
    num_sets = len(A)

    #model variables
    y = {idx_i: mdl.binary_var(name="y" + str(idx_i)) for idx_i,i in A.iteritems()}
    #constraints
    for idx_j, j in A.iterrows():
        mdl.add_constraint(1<=mdl.sum(A.loc[idx_j,idx_i] * y[idx_i] for idx_i, i in A.iteritems()))

    #objective
    mdl.minimize(mdl.sum(y[idx_i] for idx_i, i in A.iteritems()))

    #print solution
    if not mdl.solve():
        print("*** Problem has no solution")
    #else:
       #mdl.print_solution()

    obj = mdl.objective_value
    list1 = []
    for v in mdl.iter_binary_vars():
        list1.append(v.solution_value)
    df_solution = pd.DataFrame(list1, columns=['y_i'])
    df_solution.index = A.columns


    return obj, df_solution


def get_nearest_points(A, distM, dmax):
    '''
    Get the two nearest points j1 and j2 for each node in A.columns

    parameters:

        A: DataFrame
            the covering matrix
        distM: DataFrame
            the distance matrix of the problem.
        dmax: int
            the covering distance

    returns: DataFrame
        A DataFrame with columns "j1" and "j2" with the index of the closest node for each node in A.columns

    '''
    list1 = []
    for idx, i in A.iteritems():
        list1.append(distM[dmax < distM].loc[idx].nsmallest(2).index.values)

    if len(list1[0]) <2:
        print("The problem is too small to get nearest points, use CA-method or RLA-method")
        return None
    df = pd.DataFrame(list1, index=A.columns, columns=['j1', 'j2'])
    return df

def estimate_costs(distM, df_VT, df_T, p, dmax, cust, stemline):
    '''
    Estimate the cost of the nodes in V\T using the "CA-method"

    parameters:

        distM: DataFrame
            distance matrix of the problem
        df_VT: DataFrame
            the nodes of the set V\T
        df_T: DataFrame
            the nodes of the set T
        p: int
            the number of nodes estimated to be in the solution of the CTP
        dmax: int
            the maximal covering range
        cust: int
            integer to specify how the linehaul distance should be estimated:
            0 : d(j) = 0
            1 : d(j) = (djj1 + djj2)/2
            2 : d(j) = (djj1 + djj2+ djh1 + djh2)/4
        stemline: int
            integer to specify how the stemline distance should be estimated:
            0 : D(j) = 0
            1 : D(j) = d(0,j)
            2 : D(j) = (d(0,j1) + d(0,j) + d(0,j2))/3

    returns:

        df_cost : DataFrame
            DataFrame of len(df_VT) with the estimated cost for taking each node on the tour.
    '''

    list1 = []
    for idx, i in df_VT.iterrows():

        if len(df_VT.columns) <4 or (np.isnan(df_VT.loc[idx]['j1']) or np.isnan(df_VT.loc[idx]['j2'])):
            cj = distM.loc[df_T.index, idx].min()/p
        else:
            if stemline == 0:
                Dj = 0
            if stemline == 1:
                d0j = distM.loc[df_T.index, idx].min()
                Dj = d0j
            if stemline == 2:
                d0j = distM.loc[df_T.index, idx].min()
                d0j1 = distM.loc[distM.loc[df_T.index, idx].idxmin(), df_VT['j1'].loc[idx]]
                d0j2 = distM.loc[distM.loc[df_T.index, idx].idxmin(), df_VT['j2'].loc[idx]]
                Dj = (d0j + d0j1 + d0j2) / 3


            if cust == 0:
                dj = 0
            if cust == 1:
                djj1 = distM.loc[idx, df_VT['j1'].loc[idx]]
                if df_VT['j2'].loc[idx]>0:
                    djj2 = distM.loc[idx, df_VT['j2'].loc[idx]]
                    dj = (djj1 + djj2) / 2
                else:
                    dj = djj1
            if cust == 2:
                djj1 = distM.loc[idx, df_VT['j1'].loc[idx]]
                djj2 = distM.loc[idx, df_VT['j2'].loc[idx]]
                h1 = ((distM[dmax > distM].loc[idx].nlargest(1).index.values))[0]
                djh1 = distM.loc[idx, h1]
                h2 = ((distM[dmax > distM].loc[idx].nlargest(2).index.values)[1])
                djh2 = distM.loc[idx, h2]
                dj = (djj1 + djj2 + djh1 + djh2)/4
            if stemline != 0 and cust !=0:
                cj = (1 * Dj - dj) / p + dj
            elif stemline != 0:
                cj = Dj/p
            else:
                cj = dj
        list1.append(cj)
    df_cost = pd.DataFrame(list1, index=df_VT.index, columns=['cj'])

    return df_cost

def solve_setcover(A, weights, prevSolution=None):
    '''
    Solves a weighted setcover

    paramters:

        A: DataFrame
            the covering matrix
        weights: DataFrame
            A DataFrame with the costs of each node
        prevSolution: DataFrame (default None)
            if this is not None, then it should be a 0-1 DataFrame with a solution to the CTP.

    returns:

        df_solution : DataFrame
            A solution to the SCP. A 0-1 DataFrame with len(weights).
    '''
    # this is made, so the improvement heuristic does not find the same solution twice in a row.
    if prevSolution != None:
        J_1 = prevSolution
        str_index = A.columns.values.tolist()
        int_index = [int(i) for i in str_index]
        J_0 = set(int_index) - set(prevSolution)
        J_0 = list(J_0)

    # initialize
    mdl = Model(name='Setcover')

    # model variables
    y = {idx_j: mdl.binary_var(name="y" + str(idx_j)) for idx_j, j in A.iteritems()}
    # constraints
    for idx_j, j in A.iterrows():
        mdl.add_constraint(1 <= mdl.sum(A.loc[idx_j, idx_i] * y[idx_i] for idx_i, i in A.iteritems()))
    if prevSolution != None:
        mdl.add_constraint(1 <= mdl.sum(y[i] for i in J_0) + mdl.sum((y[i]-1) for i in J_1))
    # objective
    mdl.minimize(mdl.sum(y[idx_i] * weights[idx_i] for idx_i, i in A.iteritems()))

    # print solution
    if not mdl.solve():
        print("*** Problem has no solution")
    # else:
    # mdl.print_solution()

    list1 = []
    for v in mdl.iter_binary_vars():
        list1.append(v.solution_value)
    df_solution = pd.DataFrame(list1, columns=['y_i'])
    df_solution.index = A.columns
    return df_solution

def get_tsp_distM(distM, df_T, solution):
    '''
    This method returns a distance matrix of the vertices that should be on the tour

    parameters:

        distM: DataFrame
            the distance matrix
        df_T: DataFrame
            the nodes that must be on the tour
        solution: DataFrame
            a 0-1 DataFrame which specifies which nodes should be on the tour.

    returns:
        distM_tsp : DataFrame
            a distance matrix with only the nodes that should be on the tour.
    '''
    #
    list1 = []
    for i, value in df_T.iterrows():
        list1.append(i)
    for i, value in solution.iterrows():
        if value[0] == 1:
            list1.append(i)
    distM_tsp = distM.loc[list1, list1]
    return distM_tsp

def get_lower_triang_array(distM):
    '''
    produce a lower triangular numpy array from a distance matrix

    parameters:
        distM: DataFrame
            a distance matrix

    returns:
        lower_triang : np.array
            numpy array of C type integers that stores the distance matrix row by row.
    '''
    distance = []
    for i in range(len(distM)):
        for j in range(0, i):
            distance.append(int(distM.iloc[i, j]))
    lower_triang = np.asarray(distance, dtype = ctypes.c_int)
    return lower_triang

def plot_solution(df, VT , T ,W , route_1, dmax):
    '''
    plots the solution using matplotlib.

    parameters:

        df: DataFrame
            the nodes of the problem
        VT: DataFrame
            the nodes of the set V\T
        T: DataFrame
            the nodes of the set T
        W: DataFrame
            the nodes of the set W
        route_1: list of int
            the route that should be plotted, route_1[0] = route_1[-1] is the depot and route_1[-2] is the last node on the tour.
        dmax: int
            the covering distance


    returns:

    '''


    X = [df.iloc[int(i) , 0] for i in route_1]
    Y = [df.iloc[int(i) , 1] for i in route_1]

    # global figNum
    global figNum

    figNum += 1
    fig = plt.figure(figNum)
    ax = fig.subplots()
    ln, = plt.plot([], [], 'b-', markersize=2, animated=True)
    plt.scatter(VT['x'], VT['y'], c='b')
    plt.scatter(T['x'], T['y'], c='r', marker='s')
    plt.scatter(W['x'], W['y'], c='g', marker='p')

    for idx, i in T.iterrows():
        idx = int(idx)-1
        plt.annotate(str(idx), [df.iloc[int(idx), 0], df.iloc[int(idx), 1]])
    for idx, i in VT.iterrows():
        idx = int(idx) - 1
        plt.annotate(str(idx), [df.iloc[int(idx), 0], df.iloc[int(idx), 1]])
    for idx, i in W.iterrows():
        idx = int(idx) - 1
        plt.annotate(str(idx), [df.iloc[int(idx), 0], df.iloc[int(idx), 1]])

    plt.plot(X, Y, 'b-')

    ax = fig.gca()
    for i in range(len(X)):
        circle = plt.Circle((X[i],Y[i]), radius = dmax, alpha = 0.1)
        ax.add_artist(circle)

    plt.show()

def offline_regr(df, T, V, dmax, iterations=10, printit = 0, plot = 0, cust=1):
    '''
    Uses offline regression to determine the regression parameters for the linehaul and the stemline distance.

    parameters:

        df: DataFrame
            the nodes of the problem at hand
        T: int
            the number nodes that must be visited
        V: int
            the number nodes that can be visited
        W:
        dmax: int
            the covering distance
        iterations: int (default 10)
            number of data sets to generate
        plot: bool (default False)
            if True it makes a fitted regression plot.
        cust: int (default 1)
            if 0 the linehaul distance is not used, otherwise the linehaul is used.

    returns:
        a: regression parameter for the linehaul distance (fixed to 0 if cust = 0)
        b: regression parameter for the stemline distance
        c: intercept parameter
    '''
    dS_list = []
    DS_list = []
    LS_list = []
    # calculate dS(linehaul), DS(stemline) and LS(route length) for |iterations| random problems.
    for k in range(iterations):
        df_k = generate_points(len(df))  # generate n random points
        distM_k = get_distM(df_k)  # get the distance matrix of the points
        df_k_T =df[:T]
        df_k_VT =df_k[T:V]
        df_k_W = df_k[V:]
        dmax_k = max(max(distM_k.iloc[df_k_VT.index, df_k_W.index].min(axis=1)),distM_k.iloc[df_k_W.index, df_k_VT.index].apply(lambda x: x.nsmallest(2), axis=1).max(axis=1).max())
        A = get_A(distM_k.iloc[df_k_W.index, df_k_VT.index], dmax_k) # get the A matrix: A = ( a(i,j) ) og a(i,j)=1, naar d(i,j) <= dmax og 0 ellers
        p, y = solve_unweigthed_setcover(A)  # solve the unweighted set covering, with the sets in A and return the number of sets
        nearest_points_k = get_nearest_points(A, distM_k.loc[df_k_VT.index, df_k_VT.index], dmax)
        if (nearest_points_k.isnull().values.sum()) >0:
            continue
        dS, DS, LS = get_dS_DS_LS_offline(y, distM_k, df_k_T, nearest_points_k)
        dS_list.append(dS)
        DS_list.append(DS)
        LS_list.append(LS)
    df_regr = pd.DataFrame()
    df_regr['dS'] = dS_list
    df_regr['DS'] = DS_list
    df_regr['LS'] = LS_list

    if plot == 1:
        global figNum

        figNum += 1
        fig = plt.figure(figNum)
        X = df_regr[['DS']]
        Z = df_regr[['dS']]
        Y = df_regr['LS']
        regr = linear_model.LinearRegression()
        regr.fit(X, Y)

        plt.scatter(X, Y, color='g')
        plt.plot(X, regr.predict(X), color='k')

        regr.fit(Z, Y)
        plt.scatter(Z, Y, color='b')
        plt.plot(Z, regr.predict(Z), color='k')


    if cust == 0: #not including the Linehaul
        X = df_regr[['DS']] #This is changed from  X = df_regr[['ds', 'DS'] to not include the Line-haul distance
        a = 0 #this is added to not include the linehaul distance
        Y = df_regr['LS']

        regr = linear_model.LinearRegression()
        regr.fit(X, Y)

        c = regr.intercept_
        b = regr.coef_
        b = float(b)#change from a,b = regr.coef_ to not include the Line-haul distance
    if cust == 1: #including the Linehaul
        X = df_regr[['dS', 'DS']]
        Y = df_regr['LS']

        regr = linear_model.LinearRegression()
        regr.fit(X, Y)

        c = regr.intercept_
        a,b = regr.coef_  # change from a,b = regr.coef_ to not include the Line-haul distance
    return a,b,c


def get_dS_DS_LS_offline(y,distM_k, df_k_T, nearest_points_k):
    '''
    calculate ds, DS and LS for problem k.

    parameters:

        y: DataFrame
            DataFrame with the solution of the set cover for problem k
        distM_k: DataFrame
            the distance matrix for problem k
        df_k_T: DataFrame
            the nodes that must be visited from problem k
        nearest_points_k: DataFrame
            the nearest points of all the node from the set V\T of problem k.

    returns:
        ds: int
            the sum of linehaul distances for all nodes in V\T of problem k
        DS: int
            the sum of stemline distances for all nodes in V\T of problem k
        LS: int
            the route length of problem k
    '''
    list1 = []
    for i, lol in df_k_T.iterrows():
        list1.append(i)
    for i in range(1, len(y) + 1):
        if y.iloc[i - 1, 0] == 1:
            list1.append(i)

    routeLen, route = solveTSP(distM_k.loc[list1,list1], df_k_T)



    LS = routeLen


    #get d(S) and D(S)
    DS = 0
    dS = 0
    for i in list1:
        if i > len(df_k_T):
            DS += distM_k.loc[distM_k.loc[df_k_T.index, i].idxmin(),i]
            djj1 = distM_k.iloc[i,nearest_points_k['j1'].loc[i]]
            djj2 = distM_k.iloc[i,nearest_points_k['j2'].loc[i]]
            dS += (djj1 + djj2)/2
    #time.sleep(1)
    return dS, DS, LS


def online_regression(distM, df, df_T, df_VT, p, nearest_points, k=1000, cust = 1, plot = 0):
    '''
    Uses online regression to determine parameters a (line-haul) and b(stemline).

    parameters:

        distM: DataFrame
            distance matrix
        df: DataFrame
            coordinates x and y for all points in the problem
        df_T: DataFrame
            coordinates x and y for the points in T
        df_VT: DataFrame
            coordinates x and y for the points in T
        p: int
            estimated number of points that should be on the CTP tour.
        nearest_points: DataFrame
            DataFrame with the nearest points for all node in V\T
        k: number of data points to generate
        cust: int (default 1)
            Determines if the line-haul distance should be used
        plot: int default (0)
            Determines if there should be printed a regression plot.


    returns:
        a : int
            regression parameter for the line-haul distance
        b : int
            regression parameter for the stem-line distance
        c : int
            regression intercept
    '''

    dS_list = []
    DS_list = []
    LS_list = []
    # appends dS, DS, LS of a random sample of size p to the list
    for i in range(k):
        sample_k = list(df_VT.sample(p).index) #appends a sample of the set V\T of size p to the list
        df_sample_k = pd.concat([df_T,df.loc[sample_k]], axis = 0) #appends the depot(s) to the sample
        distM_sample_k = distM.loc[sample_k,sample_k]
        #get the distance matrix of the points
        dS = 0
        DS = 0
        LS = 0
        dS, DS, LS = get_dS_DS_LS_online(distM, sample_k,distM_sample_k, df_sample_k, nearest_points, df_T)
        if dS == None:
            continue
        dS_list.append(dS)
        DS_list.append(DS)
        LS_list.append(LS)


    df_regr = pd.DataFrame()
    df_regr['dS'] = dS_list
    df_regr['DS'] = DS_list
    df_regr['LS'] = LS_list

    if cust == 1:
        X = df_regr[['DS', 'dS']]
        Y = df_regr['LS']
        regr = linear_model.LinearRegression()
        regr.fit(X, Y)
        c = regr.intercept_
        a, b = regr.coef_

    if cust == 0:
        X = df_regr[['DS']]
        Y = df_regr['LS']
        regr = linear_model.LinearRegression()
        regr.fit(X, Y)

        c = regr.intercept_
        b = regr.coef_
        b = float(b)
        a = 0

    if plot:
        X = df_regr[['DS']]
        Z = df_regr[['dS']]
        Y = df_regr['LS']
        regr = linear_model.LinearRegression()
        regr.fit(X, Y)

        plt.scatter(X, Y, color='g') #green is stem-line
        plt.plot(X, regr.predict(X), color='k')

        plt.scatter(Z, Y, color='b') #blue is linehaul
        plt.plot(Z, regr.predict(Z), color='k')

    return a,b,c


def get_dS_DS_LS_online(distM, sample_k,distM_sample_k, df_sample_k, nearest_points, df_T):
    '''
    gets the values dS, dS and lS of the tour S_k

    parameters:
        distM: DataFrame
            distance matrix
        sample_k: list of int
            list of a k points in V\T
        distM_sample_k: DataFrame
            distance matrix for the sample k
        df_sample_k: DataFrame
            coordinates for the sample k
        nearest_points: DataFrame
            DataFrame containing the two nearest points not in range of dmax for all vertices of V\T
        df_T: DataFrame
            DataFrame that contains the coordinates for the points in T

    return:
        dS : int
            sum of the line-haul distances for the sample k
        DS : int
            sum of the stem-line distances for the sample k
        LS : int
            the route length for the sample k
    '''
    #get the objective value of L(S)

    routeLen, route = solveTSP(distM_sample_k, df_T)



    LS = routeLen

    #get the nearest points for each of the customers
    nearest_points_sample = nearest_points.loc[sample_k] #gets the two nearest points, which is not the closest nNearests points, for each customer


    df_sample_k = pd.concat([df_sample_k, nearest_points_sample], axis = 1) #adds the columns j1 and j2, here j1 and j2 columns becomes float, because df is 11 indexes and nearest points is 10.
    #get d(S) and D(S)
    DS = 0
    dS = 0
    for idx_i in sample_k:
        DS += distM.loc[distM.loc[df_T.index, idx_i].idxmin(),idx_i]
        djj1 = 0
        #if not (pd.isna(df_sample_k['j1'].loc[sample_k[i]])):
        djj1 = distM.loc[idx_i,df_sample_k['j1'].loc[idx_i]]
        #if (pd.isna(df_sample_k['j2'].loc[sample_k[idx_i]])):
        #    djj2 = djj1
        #else:
        djj2 = distM.loc[idx_i,df_sample_k['j2'].loc[idx_i]]
        dS += (djj1 + djj2)/2
    #time.sleep(1)
    return dS, DS, LS



def estimate_costs_regression(a,b,c, df_VT, df_T, distM, nearest_points):
    '''
    Uses the regression parameters a and b to determine the cost of taking a node V\T on the tour

    parameters:


        a: int
            regression parameter for the line-haul distance
        b: int
            regression parameter for the stem-line distance
        c: int
            intercept for the regression
        df_VT: DataFrame
            contains the coordinates for the nodes in V\T
        df_T: DataFrame
            contains the coordinates for the nodes in T
        distM: DataFrame
            contains the distance matrix
        nearest_points: DataFrame
            contains the two nearest points not in range of dmax for all vertices of V\T


    returns:

    '''
    list_costs =[]
    for idx_i, i in df_VT.iterrows():
        d0j = distM.loc[distM.loc[df_T.index, idx_i].idxmin(),idx_i]
        #d0j1 = distM.loc[distM.loc[df_T.index, idx_i].idxmin(),nearest_points['j1'].loc[idx_i]]
        #d0j2 = distM.loc[distM.loc[df_T.index, idx_i].idxmin(),nearest_points['j2'].loc[idx_i]]
        djj1 = distM.loc[idx_i,nearest_points['j1'].loc[idx_i]]
        djj2 = distM.loc[idx_i,nearest_points['j2'].loc[idx_i]]
        dj = (djj1 + djj2)/2   #one way to estimate the node-node costs
        #dj = 1.132*dmax[df_A.index[i]-1] #this is another way based on (The Logistics Systems' approximation, reference)
        #Dj = (d0j + d0j1 + d0j2)/3
        Dj = d0j
        cj = a * dj + b*Dj #max (a*df+b*d0j + c/(p+len(df_T)),1)
        list_costs.append(cj)
    df_cost = pd.DataFrame(list_costs, index = df_VT.index, columns=['cj'])
    return df_cost


def Beardwood_set_cover(p, df_T, df_V, A):
    '''
    Uses Bearwoods's RLA to determine which nodes should be in the CTP tour.

    parameters:

        p: int
            number of points estimated to be in the tour
        df_T: DataFrame
            contains the coordinates of the points in T
        df_V: DataFrame
            contains the coordinates of the points in V
        A: DataFrame
            The covering matrix of size W x V\T

    returns:
        df_solution: DataFrame
            The solution that determines what points of V\T that should be on the tour.
    '''

    # initialize

    mdl = Model(name='Bearwood')
    # model paramters
    num_sets = len(df_V)

    # model variables
    y = {idx_i: mdl.binary_var(name="y" + str(idx_i)) for idx_i, i in df_V.iterrows()}
    a_bar = mdl.integer_var(name="a_bar")
    a_min = mdl.integer_var(name="a_min")
    b_min = mdl.integer_var(name="b_min")
    D = mdl.integer_var(name="StemDistance")
    # b_bar = mdl.integer_var(name="b_bar")
    # constraints

    for idx, i in df_T.iterrows():
        mdl.add_constraint(y[idx] == 1)
    for idx_j, j in A.iterrows():
        mdl.add_constraint(1 <= mdl.sum(A.loc[idx_j, idx_i] * y[idx_i] for idx_i, i in A.iteritems()))

    mdl.add_constraint(p + len(df_T) >= mdl.sum(y[idx_i] for idx_i, i in A.iteritems()))

    mdl.add_constraint(a_min <= mdl.min(
        [mdl.max(df_V.loc[idx_i, "x"] * y[idx_i], ((y[idx_i] - 1) * -100)) for idx_i, i in df_V.iterrows()]))

    mdl.add_constraint(a_bar >= mdl.max([(df_V.loc[idx_i, "x"] * y[idx_i] - a_min) for idx_i, i in df_V.iterrows()]))

    mdl.add_constraint(b_min <= mdl.min(
        [mdl.max(df_V.loc[idx_i, "y"] * y[idx_i], ((y[idx_i] - 1) * -100)) for idx_i, i in df_V.iterrows()]))

    mdl.add_constraint(a_bar >= mdl.max([(df_V.loc[idx_i, "y"] * y[idx_i] - b_min) for idx_i, i in df_V.iterrows()]))


    #objective
    mdl.minimize(0.75 * math.sqrt(p) * a_bar)
    #print solution
    if not mdl.solve():
        print("*** Problem has no solution")
    #else:
       #mdl.print_solution()

    obj = mdl.objective_value
    list1 = []
    i=1
    for v in mdl.iter_binary_vars():
        if i <= num_sets:
            list1.append(v.solution_value)

        i+=1
    df_solution = pd.DataFrame(list1, columns=['y_i'])
    df_solution = df_solution.drop(df_T.index)


    return df_solution


def Daganzo_set_cover(distM, p, df_T, df_V, A):
    '''
        Uses Daganzo's RLA to determine which nodes should be in the CTP tour.

        parameters:

            distM : DataFrame
                contains the distance matrix for all the points in the graph.
            p: int
                number of points estimated to be in the tour
            df_T: DataFrame
                contains the coordinates of the points in T
            df_V: DataFrame
                contains the coordinates of the points in V
            A: DataFrame
                The covering matrix of size W x V\T

        returns:
            df_solution: DataFrame
                The solution that determines what points of V\T that should be on the tour.
        '''
    #initialize

    mdl = Model(name = 'Daganzo')
    #model paramters
    num_sets = len(df_V)

    #model variables
    y = {idx_i: mdl.binary_var(name="y" + str(idx_i)) for idx_i,i in df_V.iterrows()}
    a_bar = mdl.integer_var(name = "a_bar")
    a_min = mdl.integer_var(name = "a_min")
    b_min = mdl.integer_var(name = "b_min")
    D = mdl.integer_var(name = "StemDistance")
    #b_bar = mdl.integer_var(name="b_bar")
    #constraints

    for idx, i in df_T.iterrows():
        mdl.add_constraint(y[idx] == 1)
    for idx_j, j in A.iterrows():
        mdl.add_constraint(1<=mdl.sum(A.loc[idx_j,idx_i] * y[idx_i] for idx_i, i in A.iteritems()))

    mdl.add_constraint(p+len(df_T)>=mdl.sum(y[idx_i] for idx_i, i in A.iteritems()))

    mdl.add_constraint(a_min <= mdl.min([mdl.max(df_V.loc[idx_i, "x"] * y[idx_i], ((y[idx_i] - 1) * -100)) for idx_i, i in df_V.iterrows()]))

    mdl.add_constraint(a_bar >= mdl.max([(df_V.loc[idx_i, "x"] * y[idx_i] - a_min) for idx_i, i in df_V.iterrows()]))

    mdl.add_constraint(b_min <= mdl.min([mdl.max(df_V.loc[idx_i, "y"] * y[idx_i], ((y[idx_i] - 1) * -100)) for idx_i, i in df_V.iterrows()]))

    mdl.add_constraint(a_bar >= mdl.max([(df_V.loc[idx_i, "y"] * y[idx_i] - b_min) for idx_i, i in df_V.iterrows()]))

    mdl.add_constraint(D >= mdl.sum(distM.loc[distM.loc[df_T.index, idx_i].idxmin(), idx_i]*y[idx_i] for idx_i, i in A.iteritems())/max(p,0000.1))


    #objective
    mdl.minimize(3.3837*D + 0.5917*math.sqrt(p)*a_bar)
    #print solution
    if not mdl.solve():
        print("*** Problem has no solution")
    #else:
       #mdl.print_solution()

    obj = mdl.objective_value
    list1 = []
    i=1
    for v in mdl.iter_binary_vars():
        if i <= num_sets:
            list1.append(v.solution_value)

        i+=1
    df_solution = pd.DataFrame(list1, columns=['y_i'])
    #df_solution.index = A.columns
    df_solution = df_solution.drop(df_T.index)

    return df_solution


def Chien_set_cover(distM, p, df_T, df_VT, A):
    '''
            Uses Chien's RLA to determine which nodes should be in the CTP tour.

            parameters:

                distM : DataFrame
                    contains the distance matrix for all the points in the graph.
                p: int
                    number of points estimated to be in the tour
                df_T: DataFrame
                    contains the coordinates of the points in T
                df_VT: DataFrame
                    contains the coordinates of the points in V\T
                A: DataFrame
                    The covering matrix of size W x V\T

            returns:
                df_solution: DataFrame
                    The solution that determines what points of V\T that should be on the tour.
            '''
    #initialize
    #mdl.clear()
    mdl = Model(name = 'Chien')
    #model paramters
    num_sets = len(A.columns)
    #model variables
    y = {idx_i: mdl.binary_var(name="y" + str(idx_i)) for idx_i,i in A.iteritems()}
    a_bar = mdl.integer_var(name = "a_bar")
    a_min = mdl.integer_var(name = "a_min")
    b_min = mdl.integer_var(name = "b_min")
    D = mdl.integer_var(name="StemDistance")
    #b_bar = mdl.integer_var(name="b_bar")
    #constraints
    for idx_j, j in A.iterrows():
        mdl.add_constraint(1<=mdl.sum(A.loc[idx_j,idx_i] * y[idx_i] for idx_i, i in A.iteritems()))

    mdl.add_constraint(p>=mdl.sum(y[idx_i] for idx_i, i in A.iteritems()))

    mdl.add_constraint(a_min == mdl.min([mdl.max(df_VT.loc[idx_i, "x"] * y[idx_i],((y[idx_i]-1)*-100)) for idx_i, i in A.iteritems()]))
    mdl.add_constraint(a_bar >= mdl.max([(df_VT.loc[idx_i, "x"] * y[idx_i]- a_min) for idx_i, i in A.iteritems()]))
    mdl.add_constraint(b_min == mdl.min([mdl.max(df_VT.loc[idx_i, "y"] * y[idx_i],((y[idx_i]-1)*-100)) for idx_i, i in A.iteritems()]))
    mdl.add_constraint(a_bar >= mdl.max([(df_VT.loc[idx_i, "y"] * y[idx_i] - b_min) for idx_i, i in A.iteritems()]))


    mdl.add_constraint(
        D >= mdl.sum(distM.loc[distM.loc[df_T.index, idx_i].idxmin(), idx_i] * y[idx_i] for idx_i, i in A.iteritems())/p)

    #for i in mdl.iter_constraints():
    #    print(i)
    #objective
    mdl.minimize(3.9771*D + 0.5775*math.sqrt(p-1)*a_bar)
    #print solution
    if not mdl.solve():
        print("*** Problem has no solution")
    #else:
       #mdl.print_solution()

    obj = mdl.objective_value
    list1 = []
    i=1
    for v in mdl.iter_binary_vars():
        if i <= num_sets:
            list1.append(v.solution_value)

        i+=1
    df_solution = pd.DataFrame(list1, columns=['y_i'])
    df_solution.index = A.columns


    return df_solution



def Kwon1_set_cover(p, df_T, df_V, A, S=1):
    '''
            Uses Kwon's first RLA to determine which nodes should be in the CTP tour.

            parameters:

                p: int
                    number of points estimated to be in the tour
                df_T: DataFrame
                    contains the coordinates of the points in T
                df_V: DataFrame
                    contains the coordinates of the points in V
                A: DataFrame
                    The covering matrix of size W x V\T
                S: int (default 1)
                    The ratio of the width/length ratio of the rectangle that the nodes are within.

            returns:
                df_solution: DataFrame
                    The solution that determines what points of V\T that should be on the tour.
            '''

    mdl = Model(name='Kwon1')
    # model paramters
    num_sets = len(df_V)

    # model variables
    y = {idx_i: mdl.binary_var(name="y" + str(idx_i)) for idx_i, i in df_V.iterrows()}
    a_bar = mdl.integer_var(name="a_bar")
    a_min = mdl.integer_var(name="a_min")
    b_min = mdl.integer_var(name="b_min")

    # constraints
    for idx, i in df_T.iterrows():
        mdl.add_constraint(y[idx] == 1)
    for idx_j, j in A.iterrows():
        mdl.add_constraint(1 <= mdl.sum(A.loc[idx_j, idx_i] * y[idx_i] for idx_i, i in A.iteritems()))

    mdl.add_constraint(p + len(df_T) >= mdl.sum(y[idx_i] for idx_i, i in A.iteritems()))

    mdl.add_constraint(a_min <= mdl.min(
        [mdl.max(df_V.loc[idx_i, "x"] * y[idx_i], ((y[idx_i] - 1) * -100)) for idx_i, i in df_V.iterrows()]))

    mdl.add_constraint(a_bar >= mdl.max([(df_V.loc[idx_i, "x"] * y[idx_i] - a_min) for idx_i, i in df_V.iterrows()]))

    mdl.add_constraint(b_min <= mdl.min(
        [mdl.max(df_V.loc[idx_i, "y"] * y[idx_i], ((y[idx_i] - 1) * -100)) for idx_i, i in df_V.iterrows()]))

    mdl.add_constraint(a_bar >= mdl.max([(df_V.loc[idx_i, "y"] * y[idx_i] - b_min) for idx_i, i in df_V.iterrows()]))


    #objective
    mdl.minimize((0.8740 - 0.0016*p + 1.3403* S/p)*math.sqrt(p)*a_bar)
    #print solution
    if not mdl.solve():
        print("*** Problem has no solution")
    #else:
       #mdl.print_solution()

    list1 = []
    i=1
    for v in mdl.iter_binary_vars():
        if i <= num_sets:
            list1.append(v.solution_value)

        i+=1
    df_solution = pd.DataFrame(list1, columns=['y_i'])

    df_solution = df_solution.drop(df_T.index)


    return df_solution


def Kwon2_set_cover(distM, p, df_T, df_V, A, S=1):
    '''
            Uses Kwon's second RLA to determine which nodes should be in the CTP tour.

            parameters:

                distM : DataFrame
                    contains the distance matrix for all the points in the graph.
                p: int
                    number of points estimated to be in the tour
                df_T: DataFrame
                    contains the coordinates of the points in T
                df_V: DataFrame
                    contains the coordinates of the points in V
                A: DataFrame
                    The covering matrix of size W x V\T
                S: int (default 1)
                    The ratio of the width/length ratio of the rectangle that the nodes are within.

            returns:
                df_solution: DataFrame
                    The solution that determines what points of V\T that should be on the tour.
            '''
    mdl = Model(name='Kwon2')

    # model paramters
    num_sets = len(df_V)

    # model variables
    y = {idx_i: mdl.binary_var(name="y" + str(idx_i)) for idx_i, i in df_V.iterrows()}
    a_bar = mdl.integer_var(name="a_bar")
    a_min = mdl.integer_var(name="a_min")
    b_min = mdl.integer_var(name="b_min")
    D = mdl.integer_var(name="StemDistance")

    # constraints
    for idx, i in df_T.iterrows():
        mdl.add_constraint(y[idx] == 1)
    for idx_j, j in A.iterrows():
        mdl.add_constraint(1 <= mdl.sum(A.loc[idx_j, idx_i] * y[idx_i] for idx_i, i in A.iteritems()))

    mdl.add_constraint(p + len(df_T) >= mdl.sum(y[idx_i] for idx_i, i in A.iteritems()))

    mdl.add_constraint(a_min <= mdl.min(
        [mdl.max(df_V.loc[idx_i, "x"] * y[idx_i], ((y[idx_i] - 1) * -100)) for idx_i, i in df_V.iterrows()]))

    mdl.add_constraint(a_bar >= mdl.max([(df_V.loc[idx_i, "x"] * y[idx_i] - a_min) for idx_i, i in df_V.iterrows()]))

    mdl.add_constraint(b_min <= mdl.min(
        [mdl.max(df_V.loc[idx_i, "y"] * y[idx_i], ((y[idx_i] - 1) * -100)) for idx_i, i in df_V.iterrows()]))

    mdl.add_constraint(a_bar >= mdl.max([(df_V.loc[idx_i, "y"] * y[idx_i] - b_min) for idx_i, i in df_V.iterrows()]))


    mdl.add_constraint(D >= mdl.sum(distM.loc[distM.loc[df_T.index, idx_i].idxmin(), idx_i]*y[idx_i] for idx_i, i in A.iteritems())/max(p,0000.1))


    #objective
    mdl.minimize( (0.7896 - 0.0012 *p + 0.9746 * S / p) * math.sqrt(p) * a_bar + 1.1470 * D)
    #print solution
    if not mdl.solve():
        print("*** Problem has no solution")
    #else:
       #mdl.print_solution()

    list1 = []
    i=1
    for v in mdl.iter_binary_vars():
        if i <= num_sets:
            list1.append(v.solution_value)

        i+=1
    df_solution = pd.DataFrame(list1, columns=['y_i'])

    df_solution = df_solution.drop(df_T.index)

    return df_solution


def solveTSP(distM_tsp, df_T, exact=True):
    '''
    Solves a symmetric Traveling Salesmans Problem with the Concorde solver or with brute force if the problem is small.

    parameters:

        distM_tsp: DataFrame
            Distance matrix of the problem that should be solved
        df_T: DatFrame
            Contains the coordinates of the nodes in T
        exact: bool ( default = True)
            Determines if the Concorde solver should solve the problem exact or not.

    returns:
        routeLen: int
            the length of the route
        route : list of ints
            list that contains the nodes that should be on the tour.
    '''
    lower_triang = get_lower_triang_array(distM_tsp)
    routeLen, route = solveTSPdat(len(distM_tsp), GRAPH=None, dist=lower_triang, symmetric=True, ROUTE=None,
                                  exact=exact, logFile=None)
    if routeLen > 0:

        route_1 = []
        for i in range(len(route)):
            route_1.append(distM_tsp.index[route[i]])
        route_1 = [int(i) for i in route_1]
        route_1.append(route_1[0])

        return routeLen, route_1
    if routeLen == 0:
        route = []
        for idx, i in df_T.iterrows():
            route.append(idx)
        routeLen = 0
        if len(distM_tsp) == 4:
            routeLen = distM_tsp.iloc[0,1] + distM_tsp.iloc[1,2] + distM_tsp.iloc[2,3] + distM_tsp.iloc[3,0]
            route = [0,distM_tsp.iloc[1].name,distM_tsp.iloc[2].name,distM_tsp.iloc[3].name,0]
            if (distM_tsp.iloc[0,2] + distM_tsp.iloc[2,1] + distM_tsp.iloc[1,3] + distM_tsp.iloc[3,0]) < routeLen:
                route = [0,distM_tsp.iloc[2].name,distM_tsp.iloc[1].name,distM_tsp.iloc[3].name,0]
                routeLen = distM_tsp.iloc[0,2] + distM_tsp.iloc[2,1] + distM_tsp.iloc[1,3] + distM_tsp.iloc[3,0]
            if (distM_tsp.iloc[0,1] + distM_tsp.iloc[1, 3] + distM_tsp.iloc[3, 2] + distM_tsp.iloc[2, 0]) < routeLen:
                route = [0,distM_tsp.iloc[1].name,distM_tsp.iloc[3].name,distM_tsp.iloc[2].name,0]
                routeLen = distM_tsp.iloc[0, 1] + distM_tsp.iloc[1, 3] + distM_tsp.iloc[3, 2] + distM_tsp.iloc[2, 0]
        if len(distM_tsp) == 3:
            route = [0,distM_tsp.iloc[1].name,distM_tsp.iloc[2].name, 0]
            routeLen = distM_tsp.iloc[0, 1] + distM_tsp.iloc[1, 2] + distM_tsp.iloc[2, 0]
        if len(distM_tsp) == 2:
            route = [0, distM_tsp.iloc[1].name, 0]
            routeLen = distM_tsp.iloc[0, 1] + distM_tsp.iloc[1, 0]
        return round(routeLen), route


def solve_CTP(df, T, V, W, reduce = False, method = "CA", k=50, cust=0, stemline = 1, plot = False, regplot = False, returnparameter = False, a = None, b= None, c = None):
    '''
    Invoke the CTP on a problem with V+W nodes.

    parameters:

        df: DataFrame
            contains all the points in the problem that should be solved
        T: int
            Number of nodes in T
        V: int
            Number of nodes in V
        W: int
            Number of nodes in W, V+W=len(df)
        reduce: bool (defualt True)
            determines if the reduction rules should be applied to W and V\T
        method: string (default "CA-method")
            This is a string that determines which method that should be used: The options is as follows:
            "CA-method" : uses the cost-allocation method
            "ONR" : uses the regression method with data generated online
            "OFFR" : uses the regression method with data generated online
            "RLA:1" : uses the RLA method with the TSP RLA proposed by Bearwood
            "RLA:2" : uses the RLA method with the TSP RLA proposed by Daganzo
            "RLA:3" : uses the RLA method with the TSP RLA proposed by Chien
            "RLA:4" : uses the RLA method with the first TSP RLA proposed by Kwon et. al
            "RLA:5" : uses the RLA method with the second TSP RLA propose by Kwon et. al
        k: int (default 50)
            number of iterations to use in the regression methods
        cust: int (default 0)
            determines the estimator for the line-haul distance in the "CA-method" and the regression methods
        stemline: int (default 1)
            determines the estimator for the stem-line distance in the "CA-method" and the regression methods
        plot: bool (default False)
            determines if the solution should be plotted
        regplot: bool (default False
            determines if a regression plot should be plotted for the regression methods.
        returnparameter: bool
            determines if the regression parameters a and b should b returned as well

    returns:
        routeLen: int
            length of the route returned
        route : list of int
            list of the nodes that is visited on the tour
        cpu : float
            number of seconds that is used to solve the problem
        a : float ( returned if returnparameter = True)
            regression parameter for the line-haul distance
        b : float ( returned if returnparameter = True)
            regression parameter for the stem-line distance

    '''
    df_T = df[:T]
    df_V = df[:V]
    df_VT = df[T:V]
    df_W = df[V:]
    distM = get_distM(df)
    dmax = max(max(distM.loc[df_V.index.values, df_W.index.values].min(axis=1)),
               distM.loc[df_W.index.values, df_V.index.values].apply(lambda x: x.nsmallest(2), axis=1).max(
                   axis=1).max())
    A = get_A(distM.loc[df_W.index.values, df_VT.index.values], dmax)

    if reduce:
        A = reduce_W_VT(A, df_W, df_T, dmax)
        df_W = df_W.loc[A.index, :]
        df_VT = df_VT.loc[A.columns, :]
        if A.empty:
            if plot:
                plot_points(df, df_T, df_VT, df_W, dmax)
            print(V, W, "A is empty, the depot covers all the nodes")
            return 0, "the depot covers all the nodes", 0
    cpu_start = time.time()
    p, solution = solve_unweigthed_setcover(A)

    if method == "ONR" or method == "OFFR" or cust>0:
        nearest_points = get_nearest_points(A, distM.loc[df_VT.index, df_VT.index], dmax)

        # This should check if the problem is too small for the methods "ONR" and "OFFR" and "CA" with cust>0.
        if type(nearest_points) == type(None) or nearest_points.isnull().values.any():

            if returnparameter:
                return 999, "problem too small to use nearest points, use CA-method or RLA-method", 999, 0, 0
            return 999, "problem too small to use nearest points, use CA-method or RLA-method", 999
        df_VT = pd.concat([df_VT, nearest_points], axis = 1)
    if method == "Setcover":
        df_solution = solution
    if method == "CA":
        df_cost = estimate_costs(distM, df_VT, df_T, p, dmax, cust, stemline)
        df_VT = pd.concat([df_VT, df_cost], axis=1)
        df_solution = solve_setcover(A, df_VT['cj'])
    if method == "ONR":
        a,b,c = online_regression(distM, df, df_T, df_VT, p, nearest_points, k=k, cust= cust, plot = regplot)
        df_cost = estimate_costs_regression(p, a,b,c, df_VT, df_T, distM, nearest_points)
        df_VT = pd.concat([df_VT, df_cost], axis=1)
        df_solution = solve_setcover(A, df_VT['cj'])
    if method == "OFFR":
        if (a == None):
            a,b,c = offline_regr(df, T, V, W, dmax, k, plot = regplot, cust = cust)
        df_cost = estimate_costs_regression(a, b, c, df_VT, df_T, distM, nearest_points)
        df_VT = pd.concat([df_VT, df_cost], axis=1)
        df_solution = solve_setcover(A, df_VT['cj'])
    if method == "RLA:1":
        df_solution = Beardwood_set_cover(p, df_T, df_V, A)
    if method == "RLA:2":
        df_solution = Daganzo_set_cover(distM, p, df_T, df_V, A)
    if method == "RLA:3":
        df_solution = Chien_set_cover(distM, p, df_T, df_VT, A)
    if method == "RLA:4":
        df_solution = Kwon1_set_cover(p, df_T, df_V, A)
    if method == "RLA:5":
        df_solution = Kwon2_set_cover(distM, p, df_T, df_V, A)

    distM_tsp = get_tsp_distM(distM,  df_T, df_solution)

    routeLen, route = solveTSP(distM_tsp, df_T)
    cpu = time.time() - cpu_start
    if plot:
        plot_solution(df, df_VT, df_T, df_W, route, dmax)  # plots the solution

    if returnparameter:
        return routeLen, route, cpu, a, b
    else:
        return routeLen, route, cpu
