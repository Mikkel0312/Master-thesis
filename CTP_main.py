from ConstructionCTP import *
from ImproveCTP import ImproveCTP
from MDMC import MDMC




# This main method solves the CTP with one of the solution methods on a random set of points. Prints the solution value and plots the solution.

if __name__ == "__main__":
    from sys import argv

    if len(argv)<6:
        print("Usage: ", argv[0] + " T V W bool method seed")
        print("       where T, V, W is integers with the number of vertices in T, V, W respectively")
        print("       and where bool is a boolean that determines if the ImproveCTP heuristic should be used")
        print("       and where method is string with one of the follwing")
        print("       CA               for cost-allocation method")
        print("       OFFR             for regression method with data generated offline")
        print("       ONR              for regression method with data generated online")
        print("       RLA1             for RLA-method that uses Bearwoods route length approximation")
        print("       RLA2             for RLA-method that uses Daganzos route length approximation")
        print("       RLA3             for RLA-method that uses Chiens route length approximation")
        print("       RLA4             for RLA-method that uses Kwon et. al 's first route length approximation")
        print("       RLA5             for RLA-method that uses Kwon et. al 's second route length approximation")
        print("       MDMC             for the Minimum-distance Maximum Covering method")
        exit()

    method = argv[5].upper()
    T = int(argv[1])
    V = int(argv[2])
    W = int(argv[3])

    df = generate_points(V+W, 10)
    if method == "CA":
        routeLen, route, cpu = solve_CTP(df, T,V,W, method = "CA", stemline = 1, cust = 0,plot = True)
    if method == "OFFR":
        routeLen, route, cpu = solve_CTP(df, T,V,W, method = "OFFR", stemline = 1, cust = 1, k=25 ,plot = True)
    if method == "ONR":
        routeLen, route, cpu = solve_CTP(df, T, V, W, method="ONR", stemline=1, cust=1, k =200, plot=True)
    if method == "RLA1":
        routeLen, route, cpu = solve_CTP(df, T, V, W, method="RLA:1", plot=True)
    if method == "RLA2":
        routeLen, route, cpu = solve_CTP(df, T, V, W, method="RLA:2", plot=True)
    if method == "RLA3":
        routeLen, route, cpu = solve_CTP(df, T, V, W, method="RLA:3", plot=True)
    if method == "RLA4":
        routeLen, route, cpu = solve_CTP(df, T, V, W, method="RLA:4", plot=True)
    if method == "RLA5":
        routeLen, route, cpu = solve_CTP(df, T, V, W, method="RLA:5", plot=True)
    if method == "MDMC":
        routeLen, route, cpu = MDMC(df,T,V,W, plot = True)

    if bool(argv[4]) == True:
        routeLen_i, route, cpu= ImproveCTP(df, T, V, route, routeLen, plot = True)
    print("The length of the Covering Tour solution obatined with the construction heuristic " + str(method) + " is " + str(routeLen))
    if bool(argv[4]) == True:
        print("The length of the Covering Tour solution obatined with the construction heuristic " + str(method) + " and the ImproveCTP method is " + str(routeLen_i))




