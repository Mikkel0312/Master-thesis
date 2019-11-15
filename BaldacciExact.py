import pandas as pd
import numpy as np
from scipy.spatial import distance_matrix
from docplex.mp.model import Model
from ConstructionCTP import *
from ImproveCTP import ImproveCTP


def tspfile_to_df(filename):

    with open(filename) as file:
        dimension = None
        start = None
        end = None
        lines = file.readlines()

        i = 0
        while not end:
            line = lines[i]
            if line.startswith('DIMENSION :'):
                dimension = int(line.split()[2])
            if line.startswith('NODE_COORD_SECTION'):
                start = i
            if line.startswith('EOF'):
                end = i
            i = i+1


        file.seek(0)
        df = pd.read_csv(
            file,
            skiprows=start+1,
            sep=' ',
            names=['number', 'x', 'y'],
            dtype={'number': str, 'x': np.float64, 'y': np.float64},
            nrows=dimension
        )


        df.set_index('number', inplace=True)
        if (df.index[-1] == 'EOF'):
            df = df[:-1]

        return df



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


def create_OPL_file(df, T, V, W, oplfilename):
    df_T = df[:T]
    df_V = df[:V]
    df_W = df[V:V + W]

    distM = get_distM(df_V)

    distM_A = get_distM(df)
    dmax = max(max(distM_A.loc[df_V.index.values, df_W.index.values].min(axis=1)),
               distM_A.loc[df_W.index.values, df_V.index.values].apply(lambda x: x.nsmallest(2), axis=1).max(
                   axis=1).max())

    A = get_A(distM_A.loc[df_W.index.values, df_V.index.values], dmax)
    A = reduce_W_VT(A, df_W, df_T, dmax)
    f = open(oplfilename, "w")
    f.write("NumV = " + str(V - 1) + ";\n\n")
    f.write("Distance = [\n")
    for i in range(len(distM)):
        f.write("[")
        for j in range(len(distM)):
            if j < len(distM):
                f.write(str(distM.iloc[i, j].astype(int)) + ",")
        f.write(str(distM.iloc[i, 0].astype(int)) + "],\n")
    f.write("[")
    for j in range(len(distM)):
        f.write(str(distM.iloc[j, 0].astype(int)) + ",")
    f.write(str(distM.iloc[0, 0].astype(int)))
    f.write("]\n];\n\n")
    f.write("NumW = " + str(len(A.index) - 1) + ";\n\n")
    f.write("Clusters = [\n")

    for i in range(len(A.index)):
        f.write("[")
        for j in range(len(A.columns)):
            if j < len(distM) - 1:
                f.write(str(A.iloc[i, j]) + ",")
            else:
                f.write(str(A.iloc[i, j]) + "]\n")
    f.write("];")

    f.close()


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
    #A = reduction_rule5(A)
    return A


def plot_opl_sol(y):
    df_T = df[:T]
    df_V = df[:V]
    df_VT = df[T:V]
    df_W = df[V:V + W]

    distM = get_distM(df)

    dmax = max(max(distM.loc[df_V.index.values, df_W.index.values].min(axis=1)),
               distM.loc[df_W.index.values, df_V.index.values].apply(lambda x: x.nsmallest(2), axis=1).max(
                   axis=1).max())


    list1 = []
    for i in range(len(y)):
        if y[i]==1:
            list1.append(str(i+1))

    distM_tsp = distM.loc[list1, list1]


    routeLen, route = solveTSP(distM_tsp, df_T)
    plot_solution(df, df_VT, df_T, df_W, route, dmax)


def create_exact_table(file_list, T, V, W):

    results_list = []
    for file in file_list:

        df = tspfile_to_df(file)

        df = df[:V+W]



        val, sol, cpu = solve_CTP(df, T, V, W, method="CA", reduce = True, plot = True)

        Ival, Isol, Icpu = ImproveCTP(df, T, V, sol, val)

        results_list.append([file, T, V, W, val, cpu, Ival, cpu + Icpu])

    df_results = pd.DataFrame(results_list, columns = ["file", "T", "V", "W", "obj", "RT", "obj", "RT"])
    print(df_results.to_latex())

#
# file_list =  ["kroA100.tsp", "kroB100.tsp", "kroC100.tsp"]
#
# T = 1
# V = 10
# W = 90
# create_exact_table(file_list, T, V, W)
#
# filename = "KroB100.tsp"
# df = tspfile_to_df(filename)
# create_OPL_file(df, T, V, W, "CTP.txt")
#
# list_results = []
# for V,W in [[25,75], [50,50], [75,25]]:
#     for file in file_list:
#         df = tspfile_to_df(file)
#         val_1, sol, cpu_l = solve_CTP(df, T, V, W, method="OFFR", a= 0.375, b=0.114)
#         val, sol, cpu = solve_CTP(df, T, V, W, method="CA")
#         list_results.append([file, T, V, W, val, cpu, val_1,cpu_l])
#
# df_results = pd.DataFrame(list_results, columns = ["file", "T", "V", "W", "obj1", "RT1", "obj2", "RT2"])
# list_results.append(["Average", "", "", "", df_results["obj1"].mean(),df_results["RT1"].mean(),df_results["obj2"].mean(),df_results["RT2"].mean()])
# df_results = pd.DataFrame(list_results, columns = ["file", "T", "V", "W", "obj1", "RT1", "obj2", "RT2"])
#
# print(df_results.to_latex())